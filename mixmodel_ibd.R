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

se_jk <- function(theta_hat, theta_j, mi) {
    ## jackknife SE given resampled estimates and block sizes (Busing 1999)
    g <- length(mi) ## number of blocks
    n <- sum(mi) ## total number of observations
    hi <- n / mi

    t2 <- sum((1 - mi / n) * theta_j) ## constant sum term 2 within chrom sums

    t1 <- (hi * theta_hat - (hi - 1) * theta_j - g * theta_hat + t2)^2 / (hi - 1) # nolint: line_length_linter.

    se <- sqrt(sum(t1) / g)
    se
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

parser$add_argument("-l", "--length_file",
    action = "store",
    dest = "length_file",
    help = "File with number of markers for each chromosome"
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
)

n_markers <- read_tsv(args$length_file,
    col_types = "ii"
)


## --------------------------------------------------
## set up helpers for models

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

## source individuals
source_samples <- sample_info %>%
    filter(group == "source") %>%
    distinct(sample_id) %>%
    pull(sample_id)

## full set of populations
all_pops <- sample_info %>%
    filter(pop_id != "exclude") %>%
    distinct(pop_id) %>%
    pull(pop_id)

## chromosomes
chroms <- n_markers$chrom


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

## NNLS
n_param <- ncol(ibd_pop_source)

p_genome <- map_dfr(target_samples, function(x) {
    r <- pnnls(
        ibd_pop_source,
        ibd_pop_target[, x],
        sum = 1
    )

    r1 <- tibble(
        sample_id = x,
        source_pop = colnames(ibd_pop_source),
        p = r$x,
        res_norm = r$rnorm,
        chrom = 0L
    ) ## full genome gets 0

    r1
}, .progress = TRUE)


## --------------------------------------------------
## per chromosome JK estimates

cat("__ estimating JK coefficients __\n")

p_jk <- foreach(i = chroms) %dopar% {
    cat(i, "\r")

    ## target sharing matrix excluding chrom i (populations x target samples)
    ibd_pop_target <- ibd_pop %>%
        filter(
            sample1 %in% target_samples,
            chrom != i
        ) %>%
        get_palette_matrix(., all_pops, sample1, pop_id2)

    ## source sharing matrix excluding chrom i (populations x source pops)
    ibd_pop_source <- ibd_pop %>%
        filter(
            sample1 %in% source_samples,
            chrom != i
        ) %>%
        get_palette_matrix(., all_pops, pop_id1, pop_id2)

    p_chrom <- map_dfr(target_samples, function(x) {
        r <- pnnls(
            ibd_pop_source,
            ibd_pop_target[, x],
            sum = 1
        )

        r1 <- tibble(
            sample_id = x,
            source_pop = colnames(ibd_pop_source),
            p = r$x,
            res_norm = r$rnorm,
            chrom = i
        )

        r1
    })
    p_chrom
}

## --------------------------------------------------
## combine and estimate standard errors, add sources, write output

p_full <- bind_rows(p_genome, p_jk) %>%
    arrange(chrom) %>%
    left_join(n_markers,
        by = c("chrom" = "chrom")
    ) %>%
    group_by(sample_id, source_pop) %>%
    summarise(
        se = se_jk(
            p[chrom == 0],
            p[chrom > 0],
            n[chrom > 0]
        ),
        p = p[chrom == 0],
        res_norm = res_norm[chrom == 0],
        .groups = "drop"
    ) %>%
    left_join(sample_map,
        by = c("sample_id" = "sample_id")
    ) %>%
    mutate(group = "target") %>%
    select(sample_id, pop_id, group, source_pop, p, se, res_norm)

p_source <- group_map %>%
    filter(
        group == "source",
    ) %>%
    select(sample_id, group) %>%
    left_join(sample_map) %>%
    mutate(
        source_pop = pop_id,
        p = 1,
        se = 0)
 
o <- bind_rows(p_full, p_source)

cat("__ writing output __\n")

write_tsv(o,
    file = args$out_file
)

cat("__ done! __\n")
