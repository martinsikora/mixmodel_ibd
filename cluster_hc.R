#!/usr/bin/env Rscript
# Copyright 2025 Martin Sikora <martin.sikora@sund.ku.dk>
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
    library(tibble)
    library(furrr)
    library(data.tree)
    library(parallelDist)
    library(dynamicTreeCut)
    library(igraph)
    library(ggraph)
    library(tidygraph)
    library(scico)
    library(heatmap3)
})

options(future.globals.maxSize = 1e9)

## --------------------------------------------------
## functions

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

parser$add_argument("--out_cl",
    action = "store",
    dest = "out_file_cl",
    help = "Output cluster plot filename"
)

parser$add_argument("--out_hm",
    action = "store",
    dest = "out_file_hm",
    help = "Output heatmap filename"
)

parser$add_argument("--clust_method",
    action = "store",
    dest = "clust_method",
    type = "character",
    default = "ward.D2",
    help = "Clustering method for hclust [default %(default)s]"
)

parser$add_argument("--heights",
    action = "store",
    dest = "heights",
    type = "character",
    default = "1000,10000",
    help = "Height values for adaptive tree cutting [default %(default)s]" # nolint
)

parser$add_argument("--cl_size",
    action = "store",
    dest = "cl_size",
    type = "integer",
    default = 2,
    help = "Minimum cluster size for adaptive tree cutting [default %(default)s]" # nolint
)

parser$add_argument("--deep_split",
    action = "store",
    dest = "deep_split",
    type = "integer",
    default = 3,
    help = "Deep split parameter for adaptive tree cutting [default %(default)s]" # nolint
)

parser$add_argument("-t", "--threads",
    action = "store",
    dest = "threads",
    type = "integer",
    default = 1L,
    help = "Number of threads [default %(default)s]"
)

args <- parser$parse_args()

plan(multisession, workers = args$threads)


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

cat("__ reading metadata __\n")
sample_label <- read_tsv(args$sample_file,
    col_types = "ccc"
)


## --------------------------------------------------
## main clustering

cat("__ clustering __\n")

## helpers
inds <- sample_label %>%
    filter(group != "exclude") %>%
    pull(sample_id)

heights <- args$heights %>%
    strsplit(",") %>%
    unlist() %>%
    as.numeric()

## set up distance matrix
d1 <- d %>%
    filter(
        sample1 %in% inds,
        sample2 %in% inds
    ) %>%
    group_by(sample1, sample2) %>%
    summarise(ibd = sum(ibd), .groups = "drop")

d2 <- d1 %>%
    select(sample2, sample1, ibd)
colnames(d2) <- colnames(d1)
d3 <- tibble(sample1 = inds, sample2 = inds, ibd = 0)

m <- bind_rows(d1, d2, d3) %>%
    arrange(sample1, sample2) %>%
    pivot_wider(
        names_from = "sample2",
        values_from = "ibd",
        values_fill = 0
    ) %>%
    column_to_rownames("sample1") %>%
    as.matrix()
m <- m[inds, inds]
m_d <- parDist(m, threads = args$threads) %>%
    as.matrix()

## clustering of full clustering indivduals
inds_cl_full <- sample_label %>%
    filter(group == "cluster_full") %>%
    pull(sample_id)

res_hc <- hclust(
    d = as.dist(m_d[inds_cl_full, inds_cl_full]),
    method = args$clust_method
)

res_hc_cut <- future_map_dfr(heights, ~ {
    r1 <- cutreeDynamic(res_hc,
        distM = m_d[inds_cl_full, inds_cl_full],
        method = "hybrid",
        deepSplit = args$deep_split,
        minClusterSize = args$cl_size,
        cutHeight = .x
    )
    r2 <- tibble(
        sample_id = inds_cl_full,
        cluster_terminal = r1,
        cut_height = .x
    ) %>%
        arrange(cluster_terminal, sample_id)
    r2
})


## --------------------------------------------------
## reformat

## generate cluster hierarchy strings
res_hc_el <- res_hc %>%
    as.dendrogram() %>%
    as.Node(name = "0") %>%
    ToDataFrameNetwork() %>%
    as_tibble()

cl_hier <- res_hc_el %>%
    mutate(sample_id = gsub(".*/", "", to)) %>%
    filter(sample_id %in% inds) %>%
    mutate(
        cluster_id = gsub("/", "_", from),
    ) %>%
    select(sample_id, cluster_id) %>%
    arrange(cluster_id, sample_id)

## expand cluster hierarchy strings with ancestral strings
cl_ids <- cl_hier %>%
    distinct(cluster_id) %>%
    pull(cluster_id)

cl_ids_expand <- map_dfr(cl_ids, function(k) {
    r <- strsplit(k, "_") %>%
        unlist()

    r1 <- map_chr(seq_along(r), ~ {
        paste(r[1:.x], collapse = "_")
    })
    tibble(
        cluster_id = k,
        cluster_id_anc = r1,
        cluster_level = seq_along(r1) - 1L
    )
})

## add cluster cut results
cl_expand <- cl_hier %>%
    left_join(res_hc_cut) %>%
    left_join(cl_ids_expand, relationship = "many-to-many") %>%
    arrange(cut_height, cluster_id, cluster_level)

cl_term <- cl_expand %>%
    count(cut_height, cluster_terminal, cluster_level, cluster_id_anc) %>%
    group_by(cut_height, cluster_terminal) %>%
    slice_max(n) %>%
    slice_max(cluster_level)

cl_1 <- cl_expand %>%
    semi_join(cl_term) %>%
    select(sample_id, cut_height, cluster_id_anc, cluster_level) %>%
    rename("cluster_id" = "cluster_id_anc")

cl_2 <- cl_hier %>%
    left_join(cl_ids_expand, relationship = "many-to-many") %>%
    filter(cluster_id == cluster_id_anc) %>%
    mutate(cut_height = -1, cluster_id) %>%
    select(sample_id, cut_height, cluster_id, cluster_level)

cl_final <- bind_rows(cl_2, cl_1) %>%
    left_join(sample_label) %>%
    arrange(sample_id)


## add min dist clustering individuals to cluster with lowest distance
inds_cl_min <- sample_label %>%
    filter(group == "cluster_min_dist") %>%
    pull(sample_id)

if (length(inds_cl_min) > 0) {
    cl_final_inds_cl_min <- map_dfr(inds_cl_min, ~ {
        idx_min_dist <- which(m_d[.x, inds_cl_full] == min(m_d[.x, inds_cl_full]))
        if (length(idx_min_dist) > 1) {
            idx_min_dist <- sample(idx_min_dist, 1)
        }
        cl_final %>%
            filter(sample_id == inds_cl_full[idx_min_dist]) %>%
            select(sample_id:cluster_level) %>%
            mutate(sample_id = .x) %>%
            left_join(sample_label)
    })
    cl_final <- bind_rows(cl_final, cl_final_inds_cl_min)
}

## --------------------------------------------------
## plot cluster hierarchies

cat("__ plotting clusters __\n")

## convert to graph
g_plots_all <- res_hc %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(
        label1 = paste(
            label, sample_label$label[match(label, sample_label$sample_id)],
            sep = " / "
        )
    )

## plot clustering
w <- length(inds) %/% 50

pdf(args$out_file_cl,
    width = w + 7,
    height = w + 7
)

walk(unique(heights), ~ {
    cl <- cl_final %>%
        filter(cut_height == .x)

    cl_ids <- cl %>%
        count(cluster_id, sort = T) %>%
        pull(cluster_id)

    pal_c1 <- scico(length(cl_ids), palette = "batlow")
    names(pal_c1) <- cl_ids

    d_title <- tibble(
        x = 0,
        y = 0,
        label = paste("cut height", .x)
    )

    g1 <- g_plots_all %>%
        left_join(cl, by = c("label" = "sample_id"))

    l1 <- create_layout(
        g1,
        "dendrogram",
        height = height,
        circular = TRUE
    )

    p <- ggraph(l1)
    print(p +
        geom_edge_elbow(
            color = "grey40",
            lineend = "round",
            edge_width = 0.25
        ) +
        geom_node_point(aes(
            filter = leaf,
            color = cluster_id,
        )) +
        geom_node_text(
            aes(
                filter = leaf,
                angle = node_angle(x, y),
                label = label1
            ),
            hjust = "outward",
            size = 1.5
        ) +
        geom_node_text(
            aes(label = label),
            data = d_title
        ) +
        coord_cartesian(clip = "off") +
        scale_color_manual(values = pal_c1) +
        theme_void() +
        theme(
            plot.margin = unit(rep(100, 4), "points"),
            legend.position = "none",
        ))
})
dev.off()


## --------------------------------------------------
## plot heatmap

cat("__ plotting heatmap __\n")

w <- length(inds) %/% 30 + 6
pal_c <- colorRampPalette(c(rev(scico(20, palette = "davos")), "black"))(1000)
cl_dendro <- res_hc %>%
    as.dendrogram()

m1 <- m[inds_cl_full, inds_cl_full]
lab <- paste(inds_cl_full,
    sample_label$label[match(inds_cl_full, sample_label$sample_id)],
    sep = " / "
)
colnames(m1) <- lab
rownames(m1) <- lab

m_l <- log10(m1)
m_l[m_l == -Inf] <- 0

pdf(
    file = args$out_file_hm,
    width = w,
    height = w
)
walk(unique(heights), ~ {
    cl <- cl_final %>%
        filter(
            cut_height == .x,
            group == "cluster_full"
        )

    cl_ids <- cl %>%
        count(cluster_id, sort = T) %>%
        pull(cluster_id)

    pal_c1 <- scico(length(cl_ids), palette = "batlow")
    names(pal_c1) <- cl_ids

    pal_c2 <- pal_c1[cl$cluster_id]
    names(pal_c2) <- cl$sample_id

    heatmap3(m1,
        scale = "none",
        col = pal_c,
        bg = pal_c,
        symm = TRUE,
        ColSideLabs = paste("cluster height ", .x),
        ColSideColors = pal_c2[inds_cl_full],
        RowSideLabs = "",
        useRaster = FALSE,
        cexRow = 0.25,
        cexCol = 0.25,
        Rowv = cl_dendro,
        Colv = cl_dendro,
        margins = c(10, 10)
    )
    heatmap3(m_l,
        scale = "none",
        col = pal_c,
        bg = pal_c,
        ColSideLabs = paste("cluster height ", .x),
        ColSideColors = pal_c2[inds_cl_full],
        RowSideLabs = "",
        useRaster = FALSE,
        cexRow = 0.25,
        cexCol = 0.25,
        Rowv = cl_dendro,
        Colv = cl_dendro,
        margins = c(10, 10)
    )
})
dev.off()

cat("__ writing output __\n")

write_tsv(cl_final,
    file = args$out_file_tsv
)

cat("__ done! __\n")
