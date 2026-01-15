#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Count Split Reads for a Single Sample")
parser$add_argument("--fds_path", required = TRUE, help = "Path to FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--sample_id", required = TRUE, help = "Sample ID to count")
parser$add_argument("--working_dir", default = "output", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

#Setup FRASER directory structure

dir.create(file.path(args$working_dir, "savedObjects", paste0("FRASER_", args$cohort_id)), recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(args$working_dir, "savedObjects", paste0("FRASER_", args$cohort_id), "fds-object.RDS"))

#Force HDF5 to save RAM

options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

fds <- loadFraserDataSet(dir = args$working_dir, name = args$cohort_id)

fds <- countSplitReads(
  sampleID = args$sample_id,
  fds = fds,
  NcpuPerSample = args$nthreads,
  recount = TRUE
)
