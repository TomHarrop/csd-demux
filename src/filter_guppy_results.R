#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

guppy_summary_file <- snakemake@input[["guppy_results"]]
guppy_filtered_file <- snakemake@output[["guppy_results"]]

# dev
# guppy_summary_file <- "output/030_guppy-barcoder/barcoding_summary.txt"

guppy_summary <- fread(guppy_summary_file)

guppy_summary[, barcode_rear_bc := unlist(
    strsplit(barcode_rear_id, "_"))[[1]],
    by = barcode_rear_id]
guppy_summary[, barcode_front_bc := unlist(
    strsplit(barcode_front_id, "_"))[[1]],
    by = barcode_front_id]
guppy_summary[, barcode_full_bc := unlist(
    strsplit(barcode_full_arrangement, "_"))[[1]],
    by = barcode_front_id]

guppy_filtered <- guppy_summary[barcode_arrangement != "unclassified"]

# guppy_filtered <- guppy_summary[
#     barcode_front_score >= median(barcode_front_score) |
#         barcode_rear_score >= median(barcode_rear_score)][
#         barcode_front_bc == barcode_rear_bc]

fwrite(guppy_filtered, guppy_filtered_file, sep = '\t')

# Log
sessionInfo()


