library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Initialize FRASER Data Set")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--sample_ids", required = TRUE, help = "Comma-separated sample IDs")
parser$add_argument("--bam_files", required = TRUE, help = "Comma-separated BAM file paths")
parser$add_argument("--working_dir", default = "output", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

sample_ids <- unlist(strsplit(args$sample_ids, ","))
bam_files <- unlist(strsplit(args$bam_files, ","))

sample_table <- DataFrame(
  sampleID = sample_ids,
  bamFile = bam_files,
  group = seq_len(length(bam_files)),
  pairedEnd = TRUE
)
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

fds <- FraserDataSet(
  colData = sample_table,
  workingDir = args$working_dir,
  name = args$cohort_id
)

#Setup parallel execution
bp <- MulticoreParam(workers = args$nthreads)
register(bp)


#Initialize counts and metadata
fds <- countRNAData(fds, BPPARAM = bp)
fds <- calculatePSIValues(fds,BPPARAM = bp)
fds <- saveFraserDataSet(fds)

#Print location for Python to capture if needed
fds_save_path <- file.path(
  args$working_dir,
  "savedObjects",
  paste0("FRASER_", args$cohort_id),
  "fds-object.RDS"
)

cat(fds_rds_path, "\n")
