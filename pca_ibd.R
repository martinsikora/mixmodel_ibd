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

parser$add_argument("-g", "--group_file",
    action = "store",
    dest = "group_file",
    help = "File with sample to group mapping"
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

parser$add_argument("-n", "--number_pcs",
    action = "store",
    dest = "n_pcs",
    default = 20,
    help = "Number of pcs to plot [default %(default)s]"
)

parser$add_argument("--project",
    action = "store_true",
    dest = "project",
    default = FALSE,
    help = "Project target samples on PCA from source groups"
)

arg <- parser$parse_args()


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

ibd_pop <- map_dfr(arg$files, ~ {
    r <- read_tsv(.x,
        col_types = "icccdi"
    )
    r
})

sample_map <- read_tsv(arg$sample_file,
    col_types = "cc"
)

color_map <- read_tsv(arg$color_file,
    col_types = "ccci"
)

group_map <- read_tsv(arg$group_file,
    col_types = "cc"
)


## --------------------------------------------------
## calculate PCs

cat("__ calculating PCA __\n")

## aggregate IBD sharing across chromosomes
d <- ibd_pop %>%
    group_by(sample1, pop_id2) %>%
    summarise(ibd = sum(ibd)) %>%
    mutate(p_ibd = ibd / sum(ibd)) %>%
    ungroup()


## reshape into matrix
m <- d %>%
    filter(sample1 %in% group_map$sample_id[group_map$group != "exclude"]) %>%
    select(-ibd) %>%
    pivot_wider(
        names_from = pop_id2,
        values_from = p_ibd,
        values_fill = 0
    ) %>%
    column_to_rownames("sample1") %>%
    as.matrix() %>%
    t()

if (arg$project) {
    ## split matrix in ref and project
    samples_ref <- group_map %>%
        filter(group == "source") %>%
        pull(sample_id)

    samples_proj <- group_map %>%
        filter(group == "target") %>%
        pull(sample_id)

    ## SVD
    sv <- svd(m[, samples_ref])

    proj <- t(m[, samples_proj]) %*% sv$u
    proj <- t(t(proj) / sv$d)

    ## results table
    pca_ref <- sv$v %>%
        as_tibble(.name_repair = ~ paste("PC",
            seq_len(ncol(sv$v)),
            sep = ""
        )) %>%
        mutate(
            sample_id = samples_ref,
            group = "source"
        ) %>%
        left_join(sample_map) %>%
        select(sample_id, pop_id, group, everything())

    pca_proj <- proj %>%
        as_tibble(.name_repair = ~ paste("PC",
            seq_len(ncol(proj)),
            sep = ""
        )) %>%
        mutate(
            sample_id = samples_proj,
            group = "target"
        ) %>%
        left_join(sample_map) %>%
        select(sample_id, pop_id, group, everything())

    pca_res <- bind_rows(pca_ref, pca_proj)
} else {
    ## SVD
    sv <- svd(m)

    ## results table
    pca_res <- sv$v %>%
        as_tibble(.name_repair = ~ paste("PC",
            seq_len(ncol(sv$v)),
            sep = ""
        )) %>%
        mutate(sample_id = colnames(m)) %>%
        left_join(sample_map) %>%
        left_join(group_map) %>%
        select(sample_id, pop_id, group, everything())
}

write_tsv(pca_res,
    file = arg$out_file,
    col_names = TRUE
)


## --------------------------------------------------
## set up plot

cat("__ preparing plot data __\n")

## helpers
pal_c <- color_map$color
names(pal_c) <- color_map$pop_id

pal_f <- color_map$fill
names(pal_f) <- color_map$pop_id

pal_s <- color_map$shape
names(pal_s) <- color_map$pop_id

th <- theme_bw() +
    theme(
        panel.grid.major = element_line(
            linewidth = 0.25,
            linetype = "dotted"
        ),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(0.0015, "npc"),
        legend.text = element_text(size = 6)
    )

var_explained <- sv$d^2 / sum(sv$d^2) * 100

labs <- matrix(
    paste("PC", 1:arg$n_pcs,
        " (",
        formatC(var_explained[1:arg$n_pcs], format = "f", digits = 1),
        "%)",
        sep = ""
    ),
    nrow = 2
)
pcs <- matrix(paste("PC", 1:arg$n_pcs, sep = ""),
    nrow = 2
)

d <- pca_res %>%
    mutate(pop_id = factor(pop_id, levels = color_map$pop_id))

d_source <- d %>%
    filter(group == "source")

d_target <- d %>%
    filter(group == "target")

## plot
cat("__ generating plot __\n")

pdf(arg$plot_file,
    width = 12,
    height = 6
)

walk(aeq_len(ncol(pcs)), ~ {
    p <- ggplot(d, aes(
        x = !!sym(pcs[1, .x]),
        y = !!sym(pcs[2, .x]),
        color = pop_id,
        fill = pop_id,
        shape = pop_id
    ))
    print(p +
        geom_point(
            size = 2,
            alpha = 1,
            stroke = 0.75,
            data = d_source,
            show.legend = FALSE
        ) +
        geom_point(
            size = 1,
            alpha = 0.8,
            data = d_target,
            show.legend = FALSE
        ) +
        scale_color_manual(values = pal_c) +
        scale_fill_manual(values = pal_f) +
        scale_shape_manual(values = pal_s) +
        xlab(labs[1, .x]) +
        ylab(labs[2, .x]) +
        guides(color = guide_legend(
            ncol = 4,
            override.aes = list(size = 2)
        )) +
        th)
})
dev.off()

cat("__ done! __\n")
