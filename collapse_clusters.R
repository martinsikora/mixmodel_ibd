#!/usr/bin/env Rscript
# Copyright 2024 Martin Sikora <martin.sikora@sund.ku.dk>
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

suppressPackageStartupMessages({
    library(argparse)
    library(readr)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(igraph)
    library(data.tree)
    library(ggraph)
    library(tidygraph)
})


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--input_file",
    action = "store",
    dest = "input_file",
    help = "Input file with full clustering result"
)

parser$add_argument("-s", "--sample_file",
    action = "store",
    dest = "sample_file",
    help = "File with sample to group mapping for model"
)

parser$add_argument("--out_tsv",
    action = "store",
    dest = "out_file_tsv",
    help = "Output table filename"
)

parser$add_argument("--out_pdf",
    action = "store",
    dest = "out_file_pdf",
    help = "Output plot filename"
)

parser$add_argument("--p_val",
    action = "store",
    dest = "p_val",
    type = "double",
    default = 0.05,
    help = "P-value cut off for terminal cluster assignment [default %(default)s]" # nolint
)

parser$add_argument("--modularity",
    action = "store",
    dest = "modularity",
    type = "double",
    default = 0,
    help = "Minimum modularity for terminal cluster assignment [default %(default)s]" # nolint
)

args <- parser$parse_args()


## --------------------------------------------------
## read input data

cat("__ reading clustering data __\n")
cl_res <- read_tsv(
    file = args$input_file,
    col_types = "cciccdddd"
)

cat("__ reading metadata __\n")
sample_label <- read_tsv(args$sample_file,
    col_types = "cc"
)


## --------------------------------------------------
## collapse into final terminal clusters based on p-value and modularity

cl_level_terminal <- cl_res %>%
    arrange(cluster_level) %>%
    group_by(sample_id) %>%
    summarise(
        cluster_level_terminal = max(
            which(
                cluster_p_val <= args$p_val & cluster_modularity > args$modularity
            )
        )
    ) %>%
    mutate(cluster_level_terminal = case_when(
        is.na(cluster_level_terminal) ~ 1,
        TRUE ~ cluster_level_terminal
    ))

cl_final <- cl_res %>%
    left_join(cl_level_terminal) %>%
    mutate(cluster_type = case_when(
        cluster_level == cluster_level_terminal ~ "terminal",
        cluster_level < cluster_level_terminal ~ "ancestor",
        TRUE ~ "descendant_nonsig"
    ))


## --------------------------------------------------
## plot cluster hierarchies

cat("__ plotting __\n")


## convert to graph
g_plots_terminal <- cl_final %>%
    filter(cluster_type != "descendant_nonsig") %>%
    group_by(sample_id) %>%
    summarise(pathString = paste(
        "0",
        paste(cluster_id, collapse = "/"),
        sample_id[1],
        sep = "/"
    ), .groups = "drop") %>%
    as.Node() %>%
    as.igraph(directed = TRUE) %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(
        cluster_modularity = cl_final$cluster_modularity[match(name, cl_final$cluster_id)],
        ibd_avg = cl_final$ibd_avg[match(name, cl_final$cluster_id)],
        label = paste(
            name, sample_label$label[match(name, sample_label$sample_id)],
            sep = " / "
        )
    )

## plot clustering
l_terminal <- create_layout(g_plots_terminal,
    "partition",
    circular = TRUE
)

h <- vcount(g_plots_terminal) %/% 50

pdf(args$out_file_pdf,
    width = h + 9,
    height = h + 7
)

p <- ggraph(l_terminal)
p +
    geom_edge_diagonal0(
        color = "grey40",
        lineend = "round",
        edge_width = 0.25
    ) +
    geom_node_point(
        aes(
            filter = !leaf,
            size = ibd_avg,
            fill = cluster_modularity
        ),
        shape = 21,
    ) +
    geom_node_text(
        aes(
            filter = leaf,
            angle = node_angle(x, y),
            label = label
        ),
        hjust = "outward",
        size = 1.5
    ) +
    scale_edge_width(range = c(0.25, 0.5)) +
    scale_fill_gradient2(
        low = "royalblue3",
        high = "firebrick3"
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
        legend.position = "right",
        plot.margin = unit(rep(50, 4), "points")
    )
dev.off()


## --------------------------------------------------
## write results table

cat("__ writing output __\n")

write_tsv(cl_final,
    file = args$out_file_tsv
)

cat("__ done! __\n")
