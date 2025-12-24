"""
Re-implementation of a production-pipelines RNAseq pipeline, using CPG-Flow
"""

import logging
import re
from dataclasses import dataclass
from os.path import basename

from cpg_flow import stage, targets, utils
from cpg_flow.filetypes import (
    BamPath,
    CramPath,
    FastqPair,
    FastqPairs,
)
from cpg_utils import Path
from hailtop.batch.job import Job

from rdrnaseq.jobs import align_rna, count, fraser, outrider, trim


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


@dataclass
class InOutFastqPair:
    """
    Represents a single set of input and output paired FASTQ files
    """

    id: str
    input_pair: FastqPair
    output_pair: FastqPair


def get_input_output_pairs(sequencing_group: targets.SequencingGroup) -> list[InOutFastqPair]:
    """
    Get the input FASTQ pairs, determine the trimmed output FASTQ file pair paths and
    output a list of InOutFastqPair objects
    """
    inputs = get_trim_inputs(sequencing_group)
    if not inputs or not isinstance(inputs, FastqPairs):
        return []
    prefix = sequencing_group.dataset.tmp_prefix() / 'trim'
    trim_suffix = '.trimmed.fastq.gz'
    input_output_pairs = []
    for i, pair in enumerate(inputs, 1):
        if not isinstance(pair, FastqPair):
            raise Exception(f'Invalid FASTQ pair in sequencing group {sequencing_group.id}')
        input_r1_bn = re.sub('.f(ast)?q.gz', '', basename(str(pair.r1)))
        input_r2_bn = re.sub('.f(ast)?q.gz', '', basename(str(pair.r2)))
        output_r1 = prefix / f'{input_r1_bn}{trim_suffix}'
        output_r2 = prefix / f'{input_r2_bn}{trim_suffix}'
        input_output_pairs.append(
            InOutFastqPair(
                str(i),  # ID
                pair,  # input FASTQ pair
                FastqPair(output_r1, output_r2),  # output FASTQ pair
            ),
        )
    return input_output_pairs


@stage.stage()
class TrimAlignRNA(stage.SequencingGroupStage):
    """
    Trim and align RNA-seq FASTQ reads with fastp and STAR
    """

    def expected_outputs(self, sequencing_group: targets.SequencingGroup) -> dict[str, Path | str]:
        """
        Expect a pair of CRAM and CRAI files, one per set of input FASTQ files
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
            if align_jobs:
                jobs.extend(align_jobs)
        except Exception as e:
            logging.error(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')
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

    def queue_jobs(
        self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput
    ) -> stage.StageOutput | None:
        """
        Queue a job to count the reads with featureCounts.
        """
        outputs = self.expected_outputs(sequencing_group)

        cram_path = inputs.as_str(sequencing_group, TrimAlignRNA, 'cram')
        input_cram_or_bam: BamPath | CramPath = CramPath(cram_path, f'{cram_path!s}.crai')

        bam_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'bam')
        if utils.exists(bam_path):
            input_cram_or_bam = BamPath(bam_path, index_path=f'{bam_path!s}.bai')

        jobs = count.count(
            input_cram_or_bam=input_cram_or_bam,
            cram_to_bam_path=bam_path,
            output_path=outputs['count'],
            summary_path=outputs['summary'],
            sg_id=sequencing_group.id,
            job_attrs=self.get_job_attrs(sequencing_group),
        )

        return self.make_outputs(sequencing_group, data=outputs, jobs=jobs)


@stage.stage(required_stages=TrimAlignRNA)
class Fraser(stage.CohortStage):
    """
    Perform aberrant splicing analysis with FRASER.
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        """
        Generate FRASER outputs.
        """
        return {
            'Rds_data': cohort.dataset.prefix() / 'fraser' / f'{cohort.id}.fds.tar.gz',
            'seqr_data': cohort.dataset.prefix() / 'fraser' / f'{cohort.id}.results.all.csv',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Queue a job to run FRASER.
        """

        output = self.expected_outputs(cohort)

        sequencing_groups = cohort.get_sequencing_groups()

        bam_or_cram_inputs: list[tuple[str, BamPath, None] | tuple[str, CramPath, Path]] = []
        for sequencing_group in sequencing_groups:
            cram_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'cram')
            bam_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'bam')
            if utils.exists(bam_path):
                bam_or_cram_inputs.append((sequencing_group.id, BamPath(bam_path, f'{bam_path}.bai'), None))
            else:
                bam_or_cram_inputs.append((sequencing_group.id, CramPath(cram_path, f'{cram_path!s}.crai'), bam_path))

        j = fraser.fraser(
            input_bams_or_crams=bam_or_cram_inputs,
            output_fds_path=output,
            cohort_id=cohort.id,
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=output, jobs=j)


@stage.stage(required_stages=Count)
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
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput | None:
        """
        Queue a job to run outrider.
        """
        output = self.expected_outputs(cohort)
        count_inputs = [
            inputs.as_path(sequencing_group, Count, 'count') for sequencing_group in cohort.get_sequencing_groups()
        ]
        j = outrider.outrider(
            input_counts=count_inputs,
            output_rdata_path=output['RData'],
            cohort_id=cohort.id,
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=output, jobs=j)
