#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Count Non-Split Reads for a Single Sample")
parser$add_argument("--fds_path", required = TRUE, help = "Path to FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--sample_id", required = TRUE, help = "Sample ID")
parser$add_argument("--coords_path", required = TRUE, help = "Path to splice_site_coords.RDS")
parser$add_argument("--working_dir", default = "output", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

fds <- loadFraserDataSet(dir = args$working_dir, name = args$cohort_id)
splice_site_coords <- readRDS(args$coords_path)

options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

countNonSplicedReads(
  sampleID = args$sample_id,
  fds = fds,
  spliceSiteCoords = splice_site_coords,
  splitCountRanges = NULL,
  minAnchor = 5,
  NcpuPerSample = args$nthreads
)
