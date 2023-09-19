#!/usr/bin/env Rscript
# Copyright 2023 Martin Sikora <martin.sikora@sund.ku.dk>
#
#  This file is free software: you may copy, redistribute and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation, either version 2 of the License, or (at your
#  option) any later version.
#
#  This file is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(lsei))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))


## --------------------------------------------------
## functions

get_palette_matrix <- function(d, groups, g_x, g_r) {
    ## set up vars
    g_x1 <- enquo(g_x)
    g_r1 <- enquo(g_r)

    ## get grouped summary and reshape
    r1 <- d %>%
        group_by(!!g_x1, !!g_r1) %>%
        summarise(
            value = sum(ibd),
            .groups = "drop_last"
        ) %>%
        mutate(p = value / sum(value)) %>%
        ungroup() %>%
        select(!!g_x1, !!g_r1, p) %>%
        pivot_wider(
            values_from = "p",
            names_from = !!g_x1,
            values_fill = 0
        )

    m1 <- r1 %>%
        select(-!!g_r1) %>%
        as.matrix()

    rownames(m1) <- r1 %>%
        select(!!g_r1) %>%
        pull(!!g_r1)

    ## set up result matrix and return results
    m2 <- matrix(0,
        ncol = ncol(m1),
        nrow = length(groups)
    )
    rownames(m2) <- groups
    colnames(m2) <- colnames(m1)
    m2[rownames(m1), ] <- m1
    m2
}


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("files",
    nargs = "+",
    help = "Files with ibd sharing data for each chromosome"
)

parser$add_argument("-s", "--sample_file",
    action = "store",
    dest = "sample_file",
    help = "File with sample to group mapping for model"
)

parser$add_argument("-g", "--group_file",
    action = "store",
    dest = "group_file",
    help = "File with sample to group mapping"
)

parser$add_argument("-o", "--out",
    action = "store",
    dest = "out_file",
    help = "Output filename"
)

parser$add_argument("-t", "--threads",
    action = "store",
    dest = "threads",
    default = 1L,
    help = "Number of threads [default %(default)s]"
)

parser$add_argument("-f", "--fixed",
    action = "store",
    dest = "fixed",
    default = NULL,
    help = "Source groups to be fixed in every model"
)

args <- parser$parse_args()

registerDoParallel(as.integer(args$threads))


## --------------------------------------------------
## read input data

cat("__ reading IBD data __\n")
ibd_pop <- map_dfr(args$files, ~ {
    r <- read_tsv(.x,
        col_types = "icccdi"
    )
    r
})

cat("__ reading metadata __\n")
sample_map <- read_tsv(args$sample_file,
    col_types = "cc"
)

group_map <- read_tsv(args$group_file,
    col_types = "cc"
)


## --------------------------------------------------
## set up helpers for models

cat("__ setting up models__\n")

sample_info <- sample_map %>%
    left_join(group_map,
        by = c("sample_id" = "sample_id")
    )

## target samples to model
target_samples <- sample_info %>%
    filter(
        group == "target",
        !pop_id == "exclude"
    ) %>%
    pull(sample_id)

## source populations, combinations and individuals
source_pops <- sample_info %>%
    filter(group == "source") %>%
    distinct(pop_id) %>%
    pull(pop_id)

source_samples <- sample_info %>%
    filter(group == "source") %>%
    distinct(sample_id) %>%
    pull(sample_id)

if (!is.null(args$fixed)) {
    source_pops_fixed <- args$fixed %>%
        strsplit(",") %>%
        unlist()
} else {
    source_pops_fixed <- character()
}

source_pops_rotate <- source_pops[!source_pops %in% source_pops_fixed]

source_pop_combn <- map_dfr(seq_along(length(source_pops_rotate)), ~ {
    r1 <- combn(source_pops_rotate, .x)
    r2 <- matrix(rep(source_pops_fixed, ncol(r1)), ncol = ncol(r1))
    r3 <- rbind(r2, r1)
    r4 <- tibble(
        n_pop = nrow(r3),
        comb_idx = paste(n_pop,
            rep(seq_len(ncol(r3)), each = n_pop),
            sep = "_"
        ),
        pop_id = as.vector(r3)
    )
})

## full set of populations
all_pops <- sample_info %>%
    filter(pop_id != "exclude") %>%
    distinct(pop_id) %>%
    pull(pop_id)


## --------------------------------------------------
## genome-wide estimates

cat("__ estimating genome-wide coefficients __\n")

## target sharing matrix (populations x target samples)
ibd_pop_target <- ibd_pop %>%
    filter(sample1 %in% target_samples) %>%
    get_palette_matrix(., all_pops, sample1, pop_id2)

## source sharing matrix (populations x source populations)
ibd_pop_source <- ibd_pop %>%
    filter(sample1 %in% source_samples) %>%
    get_palette_matrix(., all_pops, pop_id1, pop_id2)

## NNLS for all combinatios
idx <- unique(source_pop_combn$comb_idx)

p_genome <- foreach(i = seq_along(idx)) %dopar% {
    cat(length(idx) - i, "\r")

    s1 <- source_pop_combn %>%
        filter(comb_idx == idx[i])
    n_param <- nrow(s1)

    map_dfr(target_samples, function(x) {
        r1 <- pnnls(
            ibd_pop_source[, s1$pop_id, drop = FALSE],
            ibd_pop_target[, x],
            sum = 1
        )

        r2 <- tibble(
            sample_id = x,
            source_pop = s1$pop_id,
            p = r1$x,
            res_norm = r1$rnorm,
            comb_idx = idx[i],
            n_pop = n_param,
        )
        r2
    })
}


## --------------------------------------------------
## combine and estimate standard errors, write output

p_full <- p_genome %>%
    bind_rows() %>%
    left_join(sample_map,
        by = c("sample_id" = "sample_id")
    ) %>%
    select(sample_id, pop_id, source_pop, p, res_norm, n_pop, comb_idx)

cat("\n__ writing output __\n")

write_tsv(p_full,
    file = args$out_file
)

cat("__ done! __\n")
