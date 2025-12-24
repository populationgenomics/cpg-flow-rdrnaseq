## Pipeline Overview

This pipeline implements a comprehensive RNA-seq workflow with four main stages:

**Per-sample stages:**
- **TrimAlignRNA**: Trims paired FASTQ files using fastp and aligns reads to the genome with STAR. Outputs a CRAM file with index (`<dataset_prefix>/cram/<sequencing_group_id>.cram` and `.crai`) and temporary BAM file with index (`<tmp_prefix>/bam/<sequencing_group_id>.bam` and `.bai`).

- **Count**: Quantifies gene/transcript read counts from aligned reads using featureCounts. Produces count files (`<dataset_prefix>/count/<sequencing_group_id>.count`) and summary statistics (`<dataset_prefix>/count/<sequencing_group_id>.count.summary`).

**Cohort-level stages:**
- **Fraser**: Performs aberrant splicing analysis across samples in a cohort. Consumes BAM files (preferred) or CRAM files and generates an FDS archive (`<dataset_prefix>/fraser/<cohort_id>.fds.tar.gz`).

- **Outrider**: Conducts outlier gene expression analysis using count data from all samples in a cohort. Outputs results as R data files (`<dataset_prefix>/outrider/<cohort_id>.outrider.RData`).
## Planned Future Improvements
- Integration of additional QC metrics and visualization tools. (Integrating PICARD.)
- Updating Fraser to its latest version for improved splicing analysis.
- Optimization to improve scalability and efficiency.
## Usage

```bash
analysis-runner \
    --dataset seqr \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images-dev/rdrnaseq:0.2.1-1 \
    --skip-repo-checkout \
    --description "RNA-seq analysis" \
    -o "output-description" \
    --access-level full \
    --config src/rdrnaseq/config_template.toml \
    python3 src/rdrnaseq/run_workflow.py
