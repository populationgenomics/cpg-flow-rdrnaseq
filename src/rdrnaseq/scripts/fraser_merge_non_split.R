#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Merge Non-Split Counts and Calculate PSI")
parser$add_argument("--fds_path", required = TRUE, help = "Path to FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--filtered_ranges_path", required = TRUE, help = "Path to g_ranges_non_split_counts.RDS")
parser$add_argument("--working_dir", default = "output", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

fds <- loadFraserDataSet(dir = args$working_dir, name = args$cohort_id)
register(MulticoreParam(workers = args$nthreads))

options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

#1. Merge Split

fds <- getSplitReadCountsForAllSamples(fds, recount = FALSE)

#2. Merge Non-Split

split_count_ranges <- readRDS(args$filtered_ranges_path)
fds <- getNonSplitReadCountsForAllSamples(
  fds = fds,
  splitCountRanges = split_count_ranges,
  minAnchor = 5,
  recount = FALSE
)

#3. FRASER 2.0: Calculate Jaccard & PSI

fds <- calculatePSIValues(fds)
fds <- saveFraserDataSet(fds)
