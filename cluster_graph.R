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
    library(multidplyr)
})


## --------------------------------------------------
## functions

get_pval_perm <- function(g_full, nodes, cl_membership, n_perm, res_p) {
    ## extract sub-graph
    g_sub <- induced_subgraph(g_full, nodes)
    cl_idx <- cl_membership %>%
        as.factor() %>%
        as.integer()

    ## if sub-graph has 3 or less nodes, return 1 for p-value
    if (vcount(g_sub) <= 3) {
        p_val <- 1
        p_val
    } else {
        ## modularity for sub-graph clustering result
        m_obs <- modularity(g_sub, cl_idx, E(g_sub)$weight)

        ## modularities for clustering of rewired graphs
        m_perm <- map_dbl(seq_len(n_perm), ~ {
            ## rewire network and redistribute edge weights
            g_perm <- rewire(g_sub, each_edge(p = 1))
            E(g_perm)$weight <- sample(E(g_sub)$weight)
            cl_perm <- cluster_leiden(
                graph = g_perm,
                objective_function = "modularity",
                resolution_parameter = res_p
            ) 
            ## modularity of permuted clustering
            r <- modularity(g_perm, membership(cl_perm), E(g_perm)$weight)
            r
        })
        p_val <- 1 - ecdf(m_perm)(m_obs)
        p_val
    }
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

parser$add_argument("--ibd_min",
    action = "store",
    dest = "ibd_min",
    type = "double",
    default = 5,
    help = "Minimum total IBD length for edge to keep [default %(default)s]"
)

parser$add_argument("--res_min",
    action = "store",
    dest = "res_min",
    type = "double",
    default = 1,
    help = "Minimum resolution parameter for leiden algorithm [default %(default)s]" # nolint
)

parser$add_argument("--res_max",
    action = "store",
    dest = "res_max",
    type = "double",
    default = 2,
    help = "Maximum resolution parameter for leiden algorithm [default %(default)s]" # nolint
)

parser$add_argument("--res_step",
    action = "store",
    dest = "res_step",
    type = "double",
    default = 0.5,
    help = "Resolution parameter increment size [default %(default)s]"
)

parser$add_argument("--level_max",
    action = "store",
    dest = "level_max",
    type = "integer",
    default = 10L,
    help = "Maximum clustering levels [default %(default)s]"
)

parser$add_argument("--permutations",
    action = "store",
    dest = "n_perm",
    type = "integer",
    default = 100L,
    help = "Number of network permutations for p-value [default %(default)s]"
)

parser$add_argument("-t", "--threads",
    action = "store",
    dest = "threads",
    type = "integer",
    default = 1L,
    help = "Number of threads [default %(default)s]"
)

args <- parser$parse_args()

cluster_worker <- new_cluster(args$threads)

## --------------------------------------------------
## read input data

cat("__ reading IBD data __\n")
d <- map_dfr(args$files, ~ {
    r <- read_tsv(.x,
        col_names = c("sample1", "sample2", "chrom", "ibd"),
        col_types = "cccd"
    )
    r
})

ibd <- d %>%
    group_by(sample1, sample2) %>%
    summarise(ibd = sum(ibd), .groups = "drop")


cat("__ reading metadata __\n")
sample_label <- read_tsv(args$sample_file,
    col_types = "cc"
)


## --------------------------------------------------
## main clustering

## set up graph
g <- ibd %>%
    filter(
        ibd >= args$ibd_min,
        !sample1 %in% sample_label$sample_id[sample_label$label == "exclude"],
        !sample2 %in% sample_label$sample_id[sample_label$label == "exclude"]
    ) %>%
    rename("weight" = "ibd") %>%
    graph_from_data_frame(directed = FALSE)

cat("__ clustering __\n")

## set up helpers
n_nodes <- vcount(g)
node_ids <- V(g)$name

m_cl <- matrix(character(),
    ncol = args$level_max + 1,
    nrow = n_nodes
) ## matrix for cluster membership indicator at each level
rownames(m_cl) <- node_ids

## set up loop variables and result list
i <- 0
is_clustered <- FALSE
m_cl[, 1] <- "0"
l_cl <- vector("list", length = args$level_max)

## loop over clustering until either all samples are clustered
## or maximum number of levels reached

while (i < args$level_max && !is_clustered) {
    cat("clustering level", i + 1, "\r")

    ## split node IDs into current clusters
    idx_cl <- split(node_ids, m_cl[, i + 1])

    ## loop over all sub-clusters
    r <- map2_dfr(idx_cl, names(idx_cl), ~ {
        ## if sub-cluster has only one sample, assign final cluster ID
        if (length(.x) == 1) {
            tibble(
                cluster_index = as.character(1),
                cluster_id_parent = .y,
                cluster_id = paste(cluster_id_parent, cluster_index, sep = "_"),
                cluster_level = i + 1
            )
        } else {
            ## extract sub-graph of current cluster nodes
            g_cl <- induced_subgraph(g, .x)

            ## set up loop variables
            n_clusters <- 1
            res <- args$res_min

            ## loop until either number of new clusters is > 1
            ## or maximum resolution parameter has been reached
            while (n_clusters == 1 & res < args$res_max) {
                ## cluster graph using leiden algorithim
                cl_membership <- cluster_leiden(
                    graph = g_cl,
                    objective_function = "modularity",
                    resolution_parameter = res
                ) %>%
                    membership()

                n_clusters <- cl_membership %>%
                    unique() %>%
                    length()

                ## increase resolution if no clustering found
                if (n_clusters == 1) {
                    res <- res + args$res_step
                }
            }

            ## final clustering for subgraph and resolution parameter
            cl_membership <- cluster_leiden(
                graph = g_cl,
                objective_function = "modularity",
                resolution_parameter = res
            ) %>%
                membership()

            ## average IBD sharing within subgraphs
            cl_ibd <- map_dfr(1:max(cl_membership), function(i){
                g1 <- induced_subgraph(
                    graph = g_cl, 
                    vids = names(cl_membership[cl_membership == i])
                )
                tibble(
                    cluster_index = as.character(i),
                    ibd_avg = mean(E(g1)$weight)
                )
            })

            ## return tibble with results
            tibble(
                cluster_index = as.character(cl_membership),
                cluster_id_parent = .y,
                cluster_id = paste(cluster_id_parent, cluster_index, sep = "_"),
                cluster_level = i + 1,
                resolution_parameter = res,
                cluster_modularity = modularity(g_cl, cl_membership, E(g_cl)$weight), # nolint
            ) %>%
                left_join(cl_ibd)
        }
    })

    ## add sample id and number of new sub-clusters
    r1 <- r %>%
        mutate(sample_id = unlist(idx_cl)) %>%
        group_by(cluster_id_parent) %>%
        mutate(
            n_clusters = length(unique(cluster_id)),
        ) %>%
        ungroup()

    ## add new clustering level to result indicator matrix
    ## terminal sub-clusters (N=1) set to "NA" to stop further sub-clustering
    ## split function at start of the loop won't return those

    idx <- r1$cluster_id
    idx[r1$n_clusters == 1] <- NA
    m_cl[r1$sample_id, i + 2] <- idx
    l_cl[[i + 1]] <- r1

    ## check if there is no more new sub-clusters
    ## TRUE if all parent clusters only 1 new subclusters
    is_clustered <- all(r1$n_clusters == 1)
    i <- i + 1
}

## assemble final clustering table, get p-value and add sample labels
cat("__ calculating p-values __\n")

cluster_copy(cluster_worker, c("get_pval_perm", "g", "args"))
cluster_library(cluster_worker, c("igraph", "purrr"))

cl_final <- bind_rows(l_cl) %>%
    group_by(cluster_id_parent) %>%
    partition(cluster_worker) %>%
    mutate(
        cluster_p_val = get_pval_perm(g, sample_id, cluster_id, args$n_perm, resolution_parameter[1]) # nolint
    ) %>%
    collect() %>%
    arrange(cluster_id_parent, cluster_id) %>%
    left_join(sample_label) %>%
    select(
        sample_id, label, cluster_level, cluster_id, cluster_id_parent,
        resolution_parameter, cluster_p_val, cluster_modularity, ibd_avg
    )


## --------------------------------------------------
## plot cluster hierarchies

cat("__ plotting __\n")

## find highest clustering level with all p-val < cutoff
## and modularity > 0 for each sample

## convert to graph
g_plots_all <- cl_final %>%
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
        cluster_modularity = cl_final$cluster_modularity[match(name, cl_final$cluster_id_parent)],
        ibd_avg = cl_final$ibd_avg[match(name, cl_final$cluster_id)],
        label = paste(
            name, sample_label$label[match(name, sample_label$sample_id)],
            sep = " / "
        )
    )

## plot clustering
l_all <- create_layout(
    g_plots_all,
    "partition",
    circular = TRUE
)
h <- vcount(g_plots_all) %/% 50

pdf(args$out_file_pdf,
    width = h + 9,
    height = h + 7
)

p <- ggraph(l_all)
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
        color = "grey",
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
    geom_node_text(
        aes(
            filter = !is.na(cluster_modularity),
            label = formatC(cluster_modularity, 
                format = "f", 
                digits = 4)
        ),
        hjust = "center",
        size = 1
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
 
cat("__ writing output __\n")

write_tsv(cl_final,
    file = args$out_file_tsv
)

cat("__ done! __\n")
