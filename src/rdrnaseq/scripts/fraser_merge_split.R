#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Merge Split Read Counts")
parser$add_argument("--fds_path", required = TRUE, help = "Path to FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--working_dir", default = "output", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

dir.create(file.path(args$working_dir, "savedObjects", paste0("FRASER_", args$cohort_id)), recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(args$working_dir, "savedObjects", paste0("FRASER_", args$cohort_id), "fds-object.RDS"))

options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

fds <- loadFraserDataSet(dir = args$working_dir, name = args$cohort_id)

bp <- MulticoreParam(workers = args$nthreads)
register(bp)


fds <- countRNAData(fds, BPPARAM = bp)

#Merge split counts

fds <- getSplitReadCountsForAllSamples(fds = fds, recount = FALSE)

#Prepare ranges for non-split counting

split_counts <- asSE(fds, type = "j")
split_count_ranges <- rowRanges(split_counts)
split_count_ranges <- FRASER:::annotateSpliceSite(split_count_ranges)
saveRDS(split_count_ranges, "g_ranges_split_counts.RDS")

#Filter for expression to reduce non-split counting overhead

minExpressionInOneSample <- 20
max_count <- rowMaxs(assay(split_counts, "rawCountsJ"))
passed <- max_count >= minExpressionInOneSample
filtered_ranges <- split_count_ranges[passed, ]
saveRDS(filtered_ranges, "g_ranges_non_split_counts.RDS")

splice_site_coords <- FRASER:::extractSpliceSiteCoordinates(filtered_ranges, fds)
saveRDS(splice_site_coords, "splice_site_coords.RDS")

saveFraserDataSet(fds)
