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


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-f", "--first",
    action = "store_true",
    default = FALSE,
    help = "Remove first sample from between population sharing calculation instead of randomly selected sample [default]"
)

parser$add_argument("-i", "--in",
    action = "store",
    dest = "ibd_file",
    help = "File with pairwise IBD sharing values"
)

parser$add_argument("-o", "--out",
    action = "store",
    dest = "out_file",
    help = "Output filename"
)

parser$add_argument("-s", "--sample_file",
    action = "store",
    dest = "sample_file",
    help = "File with sample to population mapping"
)

args <- parser$parse_args()


## --------------------------------------------------
## read input data

cat("__ reading IBD data __\n")
ibd <- read_tsv(args$ibd_file,
    col_names = c("sample1", "sample2", "chrom", "ibd"),
    col_types = "ccid"
)


cat("__ reading metadata __\n")
sample_map <- read_tsv(args$sample_file,
    col_types = "cc"
)


## --------------------------------------------------
## set up ibd data for aggregate calculations

cat("__ aggregating data __\n")

## add reverse order combination for each pair
## remove exclude flagged samples
ibd1 <- ibd %>%
    select(sample2, sample1, chrom, ibd)
colnames(ibd1) <- colnames(ibd)

ibd <- bind_rows(ibd, ibd1) %>%
    filter(
        !sample1 %in% sample_map$sample_id[sample_map$pop_id == "exclude"],
        !sample2 %in% sample_map$sample_id[sample_map$pop_id == "exclude"]
    ) %>%
    left_join(sample_map,
        by = c("sample1" = "sample_id")
    ) %>%
    left_join(sample_map,
        by = c("sample2" = "sample_id"),
        suffix = c("1", "2")
    )

## count number of samples per population
pop_size <- sample_map %>%
    filter(pop_id != "exclude") %>%
    count(pop_id)

## select samples to use for between population sharing
## has to be n-1 to match sample sizes for within population sharing
## defaults to random sample per population
## if args$first, remove first sample for population in sample_map table instead
if (args$first) {
    samples_donor_between <- sample_map %>%
        filter(pop_id %in% pop_size$pop_id[pop_size$n > 1]) %>%
        group_by(pop_id) %>%
        slice(2:n())
} else {
    samples_donor_between <- sample_map %>%
        filter(pop_id %in% pop_size$pop_id[pop_size$n > 1]) %>%
        group_by(pop_id) %>%
        sample_n(size = n() - 1)
}


## --------------------------------------------------
## aggregate IBD sharing

## total IBD shared with n-1 individuals of same population
ibd_within <- ibd %>%
    filter(pop_id1 == pop_id2) %>%
    group_by(chrom, sample1, pop_id1, pop_id2) %>%
    summarise(
        ibd = sum(ibd),
        .groups = "drop"
    ) %>%
    left_join(pop_size,
        by = c("pop_id2" = "pop_id")
    ) %>%
    mutate(n_inds = n - 1) %>%
    select(-n)


## total IBD shared with n-1 individuals of different populations
ibd_between <- ibd %>%
    filter(
        pop_id1 != pop_id2,
        sample2 %in% samples_donor_between$sample_id
    ) %>%
    group_by(chrom, sample1, pop_id1, pop_id2) %>%
    summarise(
        ibd = sum(ibd),
        .groups = "drop"
    ) %>%
    left_join(pop_size,
        by = c("pop_id2" = "pop_id")
    ) %>%
    mutate(n_inds = n - 1) %>%
    select(-n)


## write final table
cat("__ writing output __\n")

ibd_pop <- ibd_within %>%
    bind_rows(ibd_between) %>%
    arrange(sample1, pop_id2)

write_tsv(ibd_pop,
    file = args$out_file,
    col_names = TRUE
)

cat("__ done! __\n")
