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
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-i", "--in",
    action = "store",
    dest = "in_file",
    help = "File with mixmodel results"
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

parser$add_argument("-m", "--min_p",
    action = "store",
    default = 0.01,
    help = "Minimum ancestry proportion to plot SE  [default %(default)s]"
)

parser$add_argument("--source_grid",
    action = "store_true",
    default = FALSE,
    help = "Plot as grid with sources in separate rows"
)

args <- parser$parse_args()


## --------------------------------------------------
## read input data

cat("__ reading data __\n")

d_mixmodel <- read_tsv(args$in_file,
    col_types = "ccccddd"
)

sample_map <- read_tsv(args$sample_file,
    col_types = "cc"
)

color_map <- read_tsv(args$color_file,
    col_types = "cc"
)


## --------------------------------------------------
## set up helpers and plot data

## plotting helpers

cat("__ preparing plot data __\n")

pal_c <- color_map$color
names(pal_c) <- color_map$pop_id

pal_f <- color_map$fill
names(pal_f) <- color_map$pop_id

th <- theme_bw() +
    theme(
        panel.border = element_rect(linewidth = 0.1),
        axis.text.x = element_text(
            angle = 90,
            size = 5,
            hjust = 1,
            vjust = 0.5
        ),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(
            angle = 90,
            size = 7,
            hjust = 0,
            vjust = 0.5
        ),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")
    )


## plot data
d1 <- d_mixmodel %>%
    mutate(
        sample_id = factor(sample_id, levels = sample_map$sample_id),
        pop_id = factor(pop_id, levels = color_map$pop_id),
        source_pop = factor(source_pop, levels = color_map$pop_id)
    )

d_s1 <- d_mixmodel %>%
    filter(group == "source") %>%
    mutate(
        sample_id = factor(sample_id, levels = sample_map$sample_id),
        pop_id = factor(pop_id, levels = color_map$pop_id),
        source_pop = factor(source_pop, levels = color_map$pop_id)
    )

w <- d1 %>%
    pull(sample_id) %>%
    unique() %>%
    length() %/% 15 + 3

h <- d1 %>%
    pull(sample_id) %>%
    str_length() %>%
    max() %/% 20


## --------------------------------------------------
## plot barplot or grid

cat("__ generating plot __\n")

if (args$source_grid) {
    h1 <- d_mixmodel %>%
        pull(source_pop) %>%
        unique() %>%
        length() %/% 1.5 + 3

    pdf(args$out_file,
        width = w,
        height = h + h1
    )

    p <- ggplot(d1, aes(
        x = sample_id,
        y = p
    ))
    print(p +
        geom_col(
            aes(
                color = source_pop,
                fill = source_pop
            ),
            linewidth = 0.25
        ) +
        geom_errorbar(
            aes(
                ymax = p - se,
                ymin = p + se,
                group = source_pop
            ),
            linewidth = 0.25,
            width = 0.25,
        ) +
        geom_text(
            y = 0.5,
            label = "S",
            size = 2,
            data = d_s1,
        ) +
        geom_hline(
            yintercept = c(0, 1),
            linewidth = 0.25
        ) +
        facet_grid(source_pop ~ pop_id,
            space = "free_x",
            scales = "free_x"
        ) +
        scale_fill_manual(
            name = "Source population",
            values = pal_f
        ) +
        scale_color_manual(
            name = "Source population",
            values = pal_c
        ) +
        scale_size_manual(values = c(0.2, 0)) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("") +
        ylab("") +
        guides(fill = guide_legend(nrow = 2)) +
        th +
        theme(strip.text.y = element_text(
            angle = 0,
            size = 7,
            hjust = 0,
            vjust = 0.5
        ), ))

    dev.off()
} else {
    d1_m <- d1 %>%
        filter(p >= args$min_p) %>%
        select(sample_id, source_pop)

    d2 <- d1 %>%
        group_by(sample_id) %>%
        arrange(desc(source_pop)) %>%
        mutate(
            p_sum = cumsum(p),
            p_low = p_sum - se
        ) %>%
        ungroup() %>%
        semi_join(d1_m)

    pdf(args$out_file,
        width = w,
        height = h + 4
    )

    p <- ggplot(d1, aes(
        x = sample_id,
        y = p
    ))
    print(p +
        geom_col(
            aes(
                color = source_pop,
                fill = source_pop
            ),
            linewidth = 0.25
        ) +
        geom_errorbar(
            aes(
                ymax = p_sum,
                ymin = p_low,
                group = source_pop
            ),
            linewidth = 0.25,
            width = 0.25,
            data = d2
        ) +
        geom_text(
            y = 0.5,
            label = "S",
            size = 2,
            data = d_s1,
        ) +
        geom_point(
            aes(
                size = res_norm,
                alpha = res_norm
            ),
            fill = "grey",
            color = "grey40",
            y = 1.05,
            shape = 22
        ) +
        geom_hline(
            yintercept = c(0, 1),
            linewidth = 0.25
        ) +
        facet_grid(. ~ pop_id,
            space = "free_x",
            scales = "free_x",
        ) +
        scale_fill_manual(
            name = "Source population",
            values = pal_f
        ) +
        scale_color_manual(
            name = "Source population",
            values = pal_c
        ) +
        scale_size_continuous(range = c(0.1, 2)) +
        scale_alpha_continuous(range = c(0.1, 1)) +
        coord_cartesian(ylim = c(0, 1.07)) +
        xlab("") +
        ylab("") +
        guides(fill = guide_legend(nrow = 1)) +
        th)

    dev.off()
}

cat("__ done! __\n")
