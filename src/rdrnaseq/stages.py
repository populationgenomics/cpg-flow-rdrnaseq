"""
This file exists to define all the Stages for the workflow.
The logic for each stage can be contained here (if it is not too complex),
or can be delegated to a separate file in jobs.

Naming conventions for Stages are not enforced, but a series of recommendations have been made here:

https://cpg-populationanalysis.atlassian.net/wiki/spaces/ST/pages/185597962/Pipeline+Naming+Convention+Specification

A suggested naming convention for a stages is:
  - PascalCase (each word capitalized, no hyphens or underscores)
  - If the phrase contains an initialism (e.g. VCF), only the first character should be capitalised
  - Verb + Subject (noun) + Preposition + Direct Object (noun)  TODO(anyone): please correct my grammar is this is false
  e.g. AlignShortReadsWithBowtie2, or MakeSitesOnlyVcfWithBcftools
  - This becomes self-explanatory when reading the code and output folders

Each Stage should be a Class, and should inherit from one of
  - SequencingGroupStage
  - DatasetStage
  - CohortStage
  - MultiCohortStage
"""

from typing import TYPE_CHECKING

from workflow_name.jobs.DoSomethingGenericWithBash import echo_statement_to_file
from workflow_name.jobs.PrintPreviousJobOutputInAPythonJob import print_file_contents

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from cpg_flow.stage import MultiCohortStage, stage

if TYPE_CHECKING:
    # Path is a classic return type for a Stage, and is a shortcut for [CloudPath | pathlib.Path]
    from cpg_utils import Path
    from cpg_flow.targets import MultiCohort, StageInput, StageOutput
from cpg_utils import Path
from cpg_workflows import get_batch
from cpg_workflows.jobs import outrider
from cpg_workflows.stages.count import Count
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    StageInput,
    StageOutput,
    stage,
)



"""
Perform outlier gene expression analysis with Outrider.
"""

"""
Perform aberrant splicing analysis with FRASER.
"""

from cpg_utils import Path
from cpg_workflows import get_batch
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
)
from cpg_workflows.jobs import fraser
from cpg_workflows.stages.trim_align import TrimAlignRNA
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    StageInput,
    StageOutput,
    stage,
)

"""
Align RNA-seq reads to the genome using STAR.
"""

import logging
import re
from dataclasses import dataclass
from os.path import basename

from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import BamPath, CramPath, FastqPair, FastqPairs
from cpg_workflows.jobs import align_rna, trim
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


import logging

from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
)
from cpg_workflows.jobs import bam_to_cram, count
from cpg_workflows.stages.trim_align import TrimAlignRNA
from cpg_workflows.utils import can_reuse
from cpg_workflows.workflow import (
    SequencingGroup,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


def get_trim_inputs(sequencing_group: SequencingGroup) -> FastqPairs | None:
    """
    Get the input FASTQ file pairs for trimming
    """
    alignment_input = sequencing_group.alignment_input
    if (
        not alignment_input
        or (get_config()['workflow'].get('check_inputs', True) and not alignment_input.exists())
        or not isinstance(alignment_input, (FastqPair, FastqPairs))
    ):
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


def get_input_output_pairs(sequencing_group: SequencingGroup) -> list[InOutFastqPair]:
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
        assert isinstance(pair, FastqPair)
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


@stage
class TrimAlignRNA(SequencingGroupStage):
    """
    Trim and align RNA-seq FASTQ reads with fastp and STAR
    """

    def expected_tmp_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Expect a pair of BAM and BAI files, one per set of input FASTQ files
        """
        return {
            suffix: sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.{extension}'
            for suffix, extension in [
                ('bam', 'bam'),
                ('bai', 'bam.bai'),
            ]
        }

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Expect a pair of CRAM and CRAI files, one per set of input FASTQ files
        """
        expected_outs = {
            suffix: sequencing_group.dataset.prefix() / 'cram' / f'{sequencing_group.id}.{extension}'
            for suffix, extension in [
                ('cram', 'cram'),
                ('crai', 'cram.crai'),
            ]
        }
        # Also include the temporary BAM and BAI files, but only if the CRAM and CRAI files don't exist
        if not (expected_outs['cram'].exists() and expected_outs['crai'].exists()):
            expected_outs.update(self.expected_tmp_outputs(sequencing_group))
        return expected_outs

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to align the input FASTQ files to the genome using STAR
        """
        jobs = []

        # Run trim
        input_fq_pairs = get_trim_inputs(sequencing_group)
        if not input_fq_pairs:
            return self.make_outputs(target=sequencing_group, error_msg='No FASTQ input found')
        assert isinstance(input_fq_pairs, FastqPairs)
        trimmed_fastq_pairs = []
        for fq_pair in input_fq_pairs:
            j, out_fqs = trim.trim(
                b=get_batch(),
                sequencing_group=sequencing_group,
                input_fq_pair=fq_pair,
                job_attrs=self.get_job_attrs(sequencing_group),
                overwrite=sequencing_group.forced,
            )
            if j:
                assert isinstance(j, Job)
                jobs.append(j)
            if not out_fqs or not isinstance(out_fqs, FastqPair):
                raise Exception(f'Error trimming FASTQs for {sequencing_group}')
            trimmed_fastq_pairs.append(out_fqs)

        # Run alignment
        trimmed_fastq_pairs = FastqPairs(trimmed_fastq_pairs)
        aligned_bam_dict = self.expected_tmp_outputs(sequencing_group)
        aligned_bam = BamPath(
            path=aligned_bam_dict['bam'],
            index_path=aligned_bam_dict['bai'],
        )
        aligned_cram_dict = self.expected_outputs(sequencing_group)
        aligned_cram = CramPath(
            path=aligned_cram_dict['cram'],
            index_path=aligned_cram_dict['crai'],
        )
        try:
            align_jobs = align_rna.align(
                b=get_batch(),
                fastq_pairs=trimmed_fastq_pairs,
                sample_name=sequencing_group.id,
                genome_prefix=get_config()['references']['star'].get('ref_dir'),
                mark_duplicates=True,
                output_bam=aligned_bam,
                output_cram=aligned_cram,
                job_attrs=self.get_job_attrs(sequencing_group),
                overwrite=sequencing_group.forced,
            )
            if align_jobs:
                assert isinstance(align_jobs, list)
                assert all([isinstance(j, Job) for j in align_jobs])
                jobs.extend(align_jobs)
        except Exception as e:
            logging.error(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')
            raise Exception(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')

        # Create outputs and return jobs
        return self.make_outputs(sequencing_group, data=aligned_cram_dict, jobs=jobs)

"""
Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
"""


@stage(
    required_stages=TrimAlignRNA,
)
class Count(SequencingGroupStage):
    """
    Count reads with featureCounts.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Generate a text file output containing read counts.
        """
        return {
            'count': sequencing_group.dataset.prefix() / 'count' / f'{sequencing_group.id}.count',
            'summary': sequencing_group.dataset.prefix() / 'count' / f'{sequencing_group.id}.count.summary',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to count the reads with featureCounts.
        """
        cram_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'cram')
        crai_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'crai')
        potential_bam_path = sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam'
        potential_bai_path = sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam.bai'
        input_cram_or_bam: BamPath | CramPath | None = None
        try:
            bam_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'bam')
            bai_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'bai')
            input_cram_or_bam = BamPath(bam_path, bai_path)
        except KeyError:
            if potential_bam_path.exists() and potential_bai_path.exists():
                input_cram_or_bam = BamPath(potential_bam_path, potential_bai_path)
            else:
                input_cram_or_bam = CramPath(cram_path, crai_path)

        output_path = self.expected_outputs(sequencing_group)['count']
        summary_path = self.expected_outputs(sequencing_group)['summary']

        jobs = count.count(
            b=get_batch(),
            input_cram_or_bam=input_cram_or_bam,
            cram_to_bam_path=potential_bam_path,
            output_path=output_path,
            summary_path=summary_path,
            sample_name=sequencing_group.id,
            job_attrs=self.get_job_attrs(sequencing_group),
            overwrite=sequencing_group.forced,
        )

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)


@stage(
    required_stages=TrimAlignRNA,
)
class Fraser(CohortStage):
    """
    Perform aberrant splicing analysis with FRASER.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Generate FRASER outputs.
        """
        dataset_prefix = cohort.get_sequencing_groups()[0].dataset.prefix()
        return {cohort.id: dataset_prefix / 'fraser' / f'{cohort.id}.fds.tar.gz'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to run FRASER.
        """
        sequencing_groups = cohort.get_sequencing_groups()

        bam_or_cram_inputs: list[tuple[BamPath, None] | tuple[CramPath, Path]] = []
        for sequencing_group in sequencing_groups:
            cram_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'cram')
            crai_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'crai')
            input_bam_or_cram: BamPath | CramPath | None = None
            try:
                bam_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'bam')
                bai_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'bai')
                input_bam_or_cram = BamPath(bam_path, bai_path)
            except KeyError:
                potential_bam_path = sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam'
                potential_bai_path = sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam.bai'
                if potential_bam_path.exists() and potential_bai_path.exists():
                    input_bam_or_cram = BamPath(potential_bam_path, potential_bai_path)
                else:
                    input_bam_or_cram = CramPath(cram_path, crai_path)
            if isinstance(input_bam_or_cram, (BamPath)):
                bam_or_cram_inputs.append((input_bam_or_cram, None))
            elif isinstance(input_bam_or_cram, (CramPath)):
                bam_or_cram_inputs.append((input_bam_or_cram, potential_bam_path))

        j = fraser.fraser(
            b=get_batch(),
            input_bams_or_crams=bam_or_cram_inputs,
            output_fds_path=list(self.expected_outputs(cohort).values())[0],
            cohort_id=cohort.id,
            job_attrs=self.get_job_attrs(),
            overwrite=cohort.forced,
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=j)




@stage(required_stages=Count)
class Outrider(CohortStage):
    """
    Perform outlier gene expression analysis with Outrider.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Generate outrider outputs.
        """
        dataset_prefix = cohort.get_sequencing_groups()[0].dataset.prefix()
        return {cohort.id: dataset_prefix / 'outrider' / f'{cohort.id}.outrider.RData'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to run outrider.
        """
        count_inputs = [
            inputs.as_path(sequencing_group, Count, 'count') for sequencing_group in cohort.get_sequencing_groups()
        ]
        j = outrider.outrider(
            b=get_batch(),
            input_counts=count_inputs,
            output_rdata_path=list(self.expected_outputs(cohort).values())[0],
            cohort_id=cohort.id,
            job_attrs=self.get_job_attrs(),
            overwrite=cohort.forced,
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=j)


