## Pipeline Overview

This pipeline implements a comprehensive RNA-seq workflow with six stages:

**Per-sample stages:**
- **TrimAlignRNA**: Trims paired FASTQ files using fastp and aligns reads to the genome with STAR. Outputs a CRAM file with index (`<dataset_prefix>/cram/<sequencing_group_id>.cram` and `.crai`) and temporary BAM file with index (`<tmp_prefix>/bam/<sequencing_group_id>.bam` and `.bai`).

- **Count**: Quantifies gene/transcript read counts from aligned reads using featureCounts. Produces count files (`<dataset_prefix>/count/<sequencing_group_id>.count`) and summary statistics (`<dataset_prefix>/count/<sequencing_group_id>.count.summary`).

**Cohort-level stages:**
- **Fraser**: Performs aberrant splicing analysis across samples in a cohort. Consumes BAM files (preferred) or CRAM files and generates an FDS archive (`<dataset_prefix>/fraser/<cohort_id>.fds.tar.gz`).

- **Outrider**: Conducts outlier gene expression analysis using count data from all samples in a cohort. Outputs results as R data files (`<dataset_prefix>/outrider/<cohort_id>.outrider.RData`).

- **VariantSpliceMatch**: Annotates FRASER significant splicing regions with rare variants from a seqr-loader MatrixTable. Produces per-dataset BED and TSV files mapping variant positions to aberrant splicing events (`<dataset_prefix>/variant_splice_match/<cohort_id>_<cell_type>_<library_type>/variants_of_interest.{bed,tsv}`).

- **Dashboard**: Generates per-dataset interactive HTML dashboards combining FRASER, OUTRIDER, and variant annotation results. Includes volcano plots, searchable tables, IGV links, and seqr links. Outputs are published to the web bucket (`<web_prefix>/dashboard/<cohort_id>_<cell_type>_<library_type>/rna_dashboard.html`).

## Planned Future Improvements
- Integration of additional QC metrics and visualization tools. (Integrating PICARD.)
- Optimization to improve scalability and efficiency.
## Usage

```bash
analysis-runner \
    --dataset seqr \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/rdrnaseq:2.2.10-1 \
    --skip-repo-checkout \
    --description "RNA-seq analysis" \
    -o "output-description" \
    --access-level full \
    --config src/rdrnaseq/config_template.toml \
    python3 src/rdrnaseq/run_workflow.py

analysis-runner \
    --dataset seqr \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/rdrnaseq:2.2.10-1 \
    --skip-repo-checkout \
    --description "RNA-seq analysis" \
    -o "output-description" \
    --access-level full \
    --config src/rdrnaseq/config_template.toml \
    python3 src/rdrnaseq/run_workflow_alignment_only.py
