#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)

parser <- ArgumentParser(description = "FRASER 2.0 Statistical Analysis")
parser$add_argument("--fds_dir", required = TRUE, help = "Directory containing savedObjects")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--pval_cutoff", type = "double", default = 0.05)
parser$add_argument("--z_cutoff", type = "double")
parser$add_argument("--delta_psi_cutoff", type = "double", default = 0.3)
parser$add_argument("--min_count", type = "integer", default = 5)
parser$add_argument("--nthreads", type = "integer", default = 1)
args <- parser$parse_args()

#Force HDF5

options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

z_cutoff <- if (is.null(args$z_cutoff)) NA else args$z_cutoff
bp <- MulticoreParam(workers = args$nthreads)
register(bp)

#Load

fds <- loadFraserDataSet(dir = args$fds_dir, name = args$cohort_id)

#Filter

fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0.0, filter = FALSE)
dir.create("plots/misc", recursive = TRUE, showWarnings = FALSE)

png("plots/misc/filter_expression.png", width = 2000, height = 2000, res = 300)
plotFilterExpression(fds, bins = 100)
dev.off()

fds_filtered <- fds[mcols(fds, type = "j")[, "passed"], ]
psi_types <- c("psi5", "psi3", "jaccard")

optimal_qs <- c(
  psi5 = bestQ(fds_filtered, type = "psi5"),
  psi3 = bestQ(fds_filtered, type = "psi3"),
  jaccard = bestQ(fds_filtered, type = "jaccard")
)

fds_fit <- FRASER(fds_filtered, q = optimal_qs, BPPARAM = bp)

#Annotation

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds_fit <- annotateRangesWithTxDb(fds_fit, txdb = txdb, orgDb = orgDb)

#Results

res <- results(fds_fit, padjCutoff = args$pval_cutoff, deltaPsiCutoff = args$delta_psi_cutoff,
               zScoreCutoff = z_cutoff, minCount = args$min_count)
res_all <- results(fds_fit, padjCutoff = 1, deltaPsiCutoff = 0, minCount = 0)

write_csv(as.data.frame(res), "results.significant.csv")
write_csv(as.data.frame(res_all), "results.all.csv")

#Save final state

saveFraserDataSet(fds_fit, dir = getwd(), name = "FraserDataSet")
