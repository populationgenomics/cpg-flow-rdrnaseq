"""
Re-implementation of a production-pipelines RNAseq pipeline, using CPG-Flow
"""

from collections import defaultdict

from hailtop.batch.job import Job
from loguru import logger

from cpg_flow import stage, targets, utils
from cpg_flow.filetypes import (
    BamPath,
    CramPath,
    FastqPair,
    FastqPairs,
)
from cpg_utils import Path, config

from rdrnaseq.jobs import align_rna, bam_to_cram, count, fraser, outrider, rna_dashboard, trim, variant_splice_match


def get_trim_inputs(sequencing_group: targets.SequencingGroup) -> FastqPairs | None:
    """
    Get the input FASTQ file pairs for trimming
    """
    alignment_input = sequencing_group.alignment_input
    if not alignment_input or not alignment_input.exists() or not isinstance(alignment_input, FastqPair | FastqPairs):
        return None
    if isinstance(alignment_input, FastqPair):
        alignment_input = FastqPairs([alignment_input])
    return alignment_input


# an object to collect all samples where a CRAM->BAM job was already scheduled
# this allows cross-stage job dependencies to exist, not very CPG-Flow
samples_needing_bams: dict[str, Job] = {}


@stage.stage()
class TrimAlignRNA(stage.SequencingGroupStage):
    """
    Trim and align RNA-seq FASTQ reads with fastp and STAR
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path | str]:
        """
        Expect a pair of CRAM and CRAI files, one per set of input FASTQ files

        Creates an output BAM in temporary storage, location is a String not a path, to prevent existence checking from
        always running this stage if BAM isn't required
        """
        return {
            'cram': sequencing_group.dataset.prefix() / 'cram' / f'{sequencing_group.id}.cram',
            'bam': str(sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam'),
        }

    def queue_jobs(
        self,
        sequencing_group: targets.SequencingGroup,
        inputs: stage.StageInput,
    ) -> stage.StageOutput | None:
        """
        Queue a job to align the input FASTQ files to the genome using STAR
        """

        outputs = self.expected_outputs(sequencing_group)
        attributes = self.get_job_attrs(sequencing_group)

        jobs = []

        # Run trim
        input_fq_pairs = get_trim_inputs(sequencing_group)
        if not input_fq_pairs:
            return self.make_outputs(target=sequencing_group, error_msg='No FASTQ input found')
        if not isinstance(input_fq_pairs, FastqPairs):
            raise Exception(f'Invalid FASTQ input for {sequencing_group}')
        trimmed_fastq_pairs = []
        for fq_pair in input_fq_pairs:
            j, out_fqs = trim.trim(
                input_fq_pair=fq_pair,
                job_attrs=attributes,
            )
            if j:
                if not isinstance(j, Job):
                    raise TypeError(f"Expected 'j' to be a Job, got {type(j).__name__}")
                jobs.append(j)
            if not out_fqs or not isinstance(out_fqs, FastqPair):
                raise Exception(f'Error trimming FASTQs for {sequencing_group}')
            trimmed_fastq_pairs.append(out_fqs)

        # Run alignment
        trimmed_fastq_pairs = FastqPairs(trimmed_fastq_pairs)

        aligned_bam = BamPath(
            path=outputs['bam'],
            index_path=f'{outputs["bam"]}.bai',
        )
        aligned_cram = CramPath(
            path=outputs['cram'],
            index_path=f'{outputs["cram"]!s}.crai',
        )
        try:
            align_jobs = align_rna.align(
                fastq_pairs=trimmed_fastq_pairs,
                sample_name=sequencing_group.id,
                output_bam=aligned_bam,
                output_cram=aligned_cram,
                job_attrs=attributes,
            )
            logger.debug(f'Generating BAM for {sequencing_group.id} (Align stage)')

            # during this run, this SG will have a BAM created
            # If this trim align job is running then a bam is already being created for this sample,
            # so we can add it to the dict of samples needing bams
            samples_needing_bams[sequencing_group.id] = align_jobs[-1]

            if align_jobs:
                jobs.extend(align_jobs)
        except Exception as e:
            logger.error(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')
            raise RuntimeError(f'Error aligning RNA-seq reads for {sequencing_group}') from e

        # Create outputs and return jobs
        return self.make_outputs(sequencing_group, data=outputs, jobs=jobs)


@stage.stage(required_stages=TrimAlignRNA)
class Count(stage.SequencingGroupStage):
    """
    Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path]:
        """
        Generate a text file output containing read counts.
        """
        return {
            'count': sequencing_group.dataset.prefix() / 'count' / f'{sequencing_group.id}.count',
            'summary': sequencing_group.dataset.prefix() / 'count' / f'{sequencing_group.id}.count.summary',
        }

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Queue a job to count the reads with featureCounts.
        """
        outputs = self.expected_outputs(sequencing_group)

        jobs = []

        cram_and_bam_paths = inputs.as_dict(sequencing_group, TrimAlignRNA)

        # if this stage is running, this sample needs to have a BAM
        if not (utils.exists(cram_and_bam_paths['bam']) or (sequencing_group.id in samples_needing_bams)):
            bam_job = bam_to_cram.cram_to_bam(
                input_cram_path=cram_and_bam_paths['cram'],
                output_bam=cram_and_bam_paths['bam'],
                job_attrs=self.get_job_attrs(target=sequencing_group),
            )
            logger.info(f'Generating BAM for {sequencing_group.id} (Count stage)')
            samples_needing_bams[sequencing_group.id] = bam_job
            jobs.append(bam_job)

        count_job = count.count(
            input_bam=cram_and_bam_paths['bam'],
            output_path=outputs['count'],
            summary_path=outputs['summary'],
            sg_id=sequencing_group.id,
            job_attrs=self.get_job_attrs(sequencing_group),
        )

        # if there was a non-alignment BAM creation job, this job must wait for that to conclude
        if sequencing_group.id in samples_needing_bams:
            count_job.depends_on(samples_needing_bams[sequencing_group.id])

        return self.make_outputs(sequencing_group, data=outputs, jobs=count_job)


@stage.stage(
    required_stages=TrimAlignRNA, analysis_type='fraser', analysis_keys=['Rds_data', 'seqr_data', 'sig_results']
)
class Fraser(stage.CohortStage):
    """
    Perform aberrant splicing analysis with FRASER.
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        """
        Generate FRASER outputs.
        """
        return {
            'Rds_data': str(cohort.dataset.tmp_prefix() / 'fraser' / f'{cohort.id}.fds.tar.gz'),
            'seqr_data': cohort.dataset.prefix() / 'fraser' / f'{cohort.id}.results.all.csv.gz',
            'sig_results': cohort.dataset.prefix() / 'fraser' / f'{cohort.id}.results.significant.csv',
            'temp_data': str(cohort.dataset.tmp_prefix() / 'fraser' / f'{cohort.id}.fraser_temp_data'),
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Queue a job to run the refactored FRASER analysis.
        """
        output = self.expected_outputs(cohort)

        bam_inputs: list[tuple[str, str]] = []

        for sequencing_group in cohort.get_sequencing_groups():
            cram_and_bam_paths = inputs.as_dict(sequencing_group, TrimAlignRNA)
            # if this stage is running, this sample needs to have a BAM
            if not (utils.exists(cram_and_bam_paths['bam']) or (sequencing_group.id in samples_needing_bams)):
                bam_job = bam_to_cram.cram_to_bam(
                    input_cram_path=cram_and_bam_paths['cram'],
                    output_bam=cram_and_bam_paths['bam'],
                    job_attrs=self.get_job_attrs(target=sequencing_group),
                )
                logger.info(f'Generating BAM for {sequencing_group.id} (FRASER stage)')
                samples_needing_bams[sequencing_group.id] = bam_job

            bam_inputs.append((sequencing_group.id, cram_and_bam_paths['bam']))

        jobs = fraser.fraser_pipeline(
            input_bams=bam_inputs,
            output_fds_path=output,
            cohort_id=cohort.id,
            job_attrs=self.get_job_attrs(),
            output_prefix=output['temp_data'],
        )

        jobs = [x for x in jobs if x is not None]
        # if there was a non-alignment BAM creation job, this job must wait for that to conclude
        for sequencing_group in cohort.get_sequencing_groups():
            parent_job = samples_needing_bams.get(sequencing_group.id)
            if parent_job and jobs:
                for j in jobs:
                    # We already cleaned 'jobs', but a final check doesn't hurt
                    if j is not None:
                        j.depends_on(parent_job)

        return self.make_outputs(cohort, data=output, jobs=jobs)


@stage.stage(required_stages=Count, analysis_type='outrider', analysis_keys=['RData', 'seqr_out'])
class Outrider(stage.CohortStage):
    """
    Perform outlier gene expression analysis with Outrider.
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        """
        Generate outrider outputs.
        """
        return {
            'RData': cohort.dataset.prefix() / 'outrider' / f'{cohort.id}.outrider.RData',
            'seqr_out': cohort.dataset.prefix() / 'outrider' / f'{cohort.id}.outrider.aberrant_genes_per_sample.csv',
            'outrider_csv': cohort.dataset.prefix() / 'outrider' / f'{cohort.id}.outrider.results.all.csv',
            'outrider_sig_csv': cohort.dataset.prefix() / 'outrider' / f'{cohort.id}.outrider.results.csv',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Queue a job to run outrider.
        """
        requested_nthreads: int = config.config_retrieve(['cohort_job_resources', 'ncpu'])
        output = self.expected_outputs(cohort)
        count_inputs = [
            inputs.as_path(sequencing_group, Count, 'count') for sequencing_group in cohort.get_sequencing_groups()
        ]
        j = outrider.outrider(
            input_counts=count_inputs,
            output_rdata_path=output['RData'],
            cohort_id=cohort.id,
            job_attrs=self.get_job_attrs(),
            requested_nthreads=requested_nthreads,
        )
        return self.make_outputs(cohort, data=output, jobs=j)


@stage.stage(required_stages=Fraser, analysis_type='custom', analysis_keys=['bed', 'tsv'])
class VariantSpliceMatch(stage.CohortStage):
    """
    Annotate FRASER significant regions (determined using rna data)
     with rare variants from a seqr-loader MatrixTable (by finding their corresponding genome data).
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        return {
            'bed': cohort.dataset.prefix() / 'variant_splice_match' / f'{cohort.id}_variants_of_interest.bed',
            'tsv': cohort.dataset.prefix() / 'variant_splice_match' / f'{cohort.id}_variants_of_interest.tsv',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        output = self.expected_outputs(cohort)

        fraser_csv = inputs.as_path(cohort, Fraser, 'sig_results')


        sg_ids_by_dataset: dict[str, list[str]] = defaultdict(list)
        # todo this doesn't correctly identify `dataset-test`
        for sg in cohort.get_sequencing_groups():
            sg_ids_by_dataset[sg.dataset.name].append(sg.id)
            # TODO: if this is more than one dataset per cohort, this will need to be refactored

        jobs = variant_splice_match.match_variants_and_splicing(
            fraser_csv=fraser_csv,
            sg_ids_by_dataset=sg_ids_by_dataset,
            output=str(output['bed']).removesuffix('.bed'),
            cohort_id=cohort.id,
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=output, jobs=jobs)


@stage.stage(
    required_stages=[Fraser, Outrider, VariantSpliceMatch], analysis_type='web', analysis_keys=['dashboard_html']
)
class Dashboard(stage.CohortStage):
    """
    Create an interactive HTML dashboard from FRASER and OUTRIDER results.
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        prefix = cohort.dataset.web_prefix() / 'rna_dashboard'
        base = f'{cohort.id}.rna_dashboard'
        return {
            'dashboard_html': prefix / f'{base}.html',
            'fraser_csv': prefix / f'{base}.fraser.csv',
            'outrider_csv': prefix / f'{base}.outrider.csv',
            'variant_bed': prefix / f'{base}.variants.bed',
            'variant_tsv': prefix / f'{base}.variants.tsv',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        output = self.expected_outputs(cohort)

        fraser_csv = inputs.as_path(cohort, Fraser, 'sig_results')
        outrider_csv = inputs.as_path(cohort, Outrider, 'outrider_sig_csv')
        variant_bed = inputs.as_path(cohort, VariantSpliceMatch, 'bed')
        variant_tsv = inputs.as_path(cohort, VariantSpliceMatch, 'tsv')

        sg_ids_by_dataset: dict[str, list[str]] = defaultdict(list)
        for sg in cohort.get_sequencing_groups():
            sg_ids_by_dataset[sg.dataset.name].append(sg.id)
            # TODO: if this is more than one dataset per cohort, this will need to be refactored

        logger.info(f'Sequence group IDs by dataset for dashboard: {dict(sg_ids_by_dataset)}')

        jobs = rna_dashboard.make_dashboards(
            fraser_csv=fraser_csv,
            outrider_csv=outrider_csv,
            variant_bed=variant_bed,
            variant_tsv=variant_tsv,
            output_paths=output,
            sg_ids_by_dataset=sg_ids_by_dataset,
            cohort_id=cohort.id,
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=output, jobs=jobs)
