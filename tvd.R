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
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phytools))


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("files",
    nargs = "+",
    help = "Files with ibd sharing data for each chromosome"
)

parser$add_argument("-c", "--color_file",
    action = "store",
    dest = "color_file",
    help = "File with color mapping"
)

parser$add_argument("-s", "--sample_file",
    action = "store",
    dest = "sample_file",
    help = "File with sample to population mapping"
)

parser$add_argument("-o", "--out",
    action = "store",
    dest = "out_file",
    help = "output file"
)

parser$add_argument("-p", "--plot_file",
    action = "store",
    dest = "plot_file",
    help = "Plot file"
)

args <- parser$parse_args()


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

ibd_pop <- map_dfr(args$files, ~ {
    r <- read_tsv(.x,
        col_types = "icccdi"
    )
    r
})

sample_map <- read_tsv(args$sample_file,
    col_types = "ccc"
)

color_map <- read_tsv(args$color_file,
    col_types = "ccci"
)


## --------------------------------------------------
##  TVD between groups

cat("__ calculating TVD __\n")

## aggregate IBD sharing across chromosomes and groups
d <- ibd_pop %>%
    group_by(pop_id1, pop_id2) %>%
    summarise(ibd = sum(ibd)) %>%
    mutate(p_ibd = ibd / sum(ibd)) %>%
    ungroup()

## reshape into matrix
m <- d %>%
    select(-ibd) %>%
    pivot_wider(
        names_from = pop_id2,
        values_from = p_ibd,
        values_fill = 0
    ) %>%
    column_to_rownames("pop_id1") %>%
    as.matrix() %>%
    t()


## TVD for each group
tvd <- map_dfr(colnames(m), ~ {
    r <- colSums(abs(m[, .x] - m) / 2)

    r1 <- tibble(
        pop_id1 = .x,
        pop_id2 = names(r),
        tvd = r
    )
    r1
})


## --------------------------------------------------
## plot tree

cat("__ generating plot __\n")

## get NJ tree
m1 <- tvd %>%
    pivot_wider(
        names_from = pop_id2,
        values_from = tvd
    ) %>%
    column_to_rownames("pop_id1") %>%
    as.matrix()

tr <- njs(m1) %>%
    midpoint.root()

## set up and plot
pal_c <- color_map$color
names(pal_c) <- color_map$pop_id

pal_f <- color_map$fill
names(pal_f) <- color_map$pop_id

pal_s <- color_map$shape
names(pal_s) <- color_map$pop_id

w <- tvd %>%
    pull(pop_id1) %>%
    unique() %>%
    length() %/% 20 + 3

pdf(args$plot_file,
    width = 5,
    height = w
)
p <- ggtree(tr,
    aes(
        x = x,
        y = y
    ),
    size = 0.1,
)
p %<+% color_map +
    geom_tippoint(
        aes(
            colour = label,
            fill = label,
            shape = label
        ),
        size = 1,
        alpha = 1
    ) +
    geom_tiplab(
        size = 1.5
    ) +
    scale_color_manual(values = pal_c) +
        scale_fill_manual(values = pal_f) +
        scale_shape_manual(values = pal_s) +
        xlim(c(0, 1)) + 
    theme_tree() +
    theme(legend.position = "none")
dev.off()


## --------------------------------------------------
## write table

cat("__ writing output __\n")

write_tsv(tvd,
    file = args$out_file
)

cat("__ done! __\n")
