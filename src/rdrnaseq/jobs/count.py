"""
Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
"""

import hailtop.batch as hb
from cpg_flow.filetypes import BamPath
from cpg_flow.resources import STANDARD
from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import Job


class FeatureCounts:
    """
    Construct a featureCounts command for counting reads.
    """

    def __init__(
        self,
        input_bam: BamPath | str | Path,
        gtf_file: hb.ResourceFile,
        output_path: str | Path,
        summary_path: str | Path,
        paired_end: bool = True,
        feature_type: str = 'exon',
        attribute: str = 'gene_id',
        strandness: str = 'reverse',  # one of: 'none', 'forward', 'reverse'
        multi_mapping: bool = False,
        min_quality: int | None = None,
        primary_only: bool = False,
        ignore_duplicates: bool = False,
        count_pairs: bool = True,
        both_ends_mapped: bool = True,
        both_ends_same_chr: bool = True,
        threads: int = 1,
    ) -> None:
        self.command = [
            'featureCounts',
            *('-t', feature_type),
            *('-g', attribute),
            *('-s', {'none': '0', 'forward': '1', 'reverse': '2'}[strandness]),
            *('-a', str(gtf_file)),
            *('-T', str(threads)),
        ]

        if paired_end:
            self.command.append('-p')

        if multi_mapping:
            self.command.append('-M')

        if min_quality is not None:
            self.command.extend(['-Q', str(min_quality)])

        if primary_only:
            self.command.append('--primary')

        if ignore_duplicates:
            self.command.append('--ignoreDup')

        if paired_end and count_pairs:
            self.command.append('--countReadPairs')

        if both_ends_mapped:
            self.command.append('-B')

        if not both_ends_same_chr:
            self.command.append('-C')

        self.tmp_output = '$BATCH_TMPDIR/count_out/count'
        self.tmp_output_summary = f'{self.tmp_output}.summary'

        self.command.extend(['-o', self.tmp_output, str(input_bam)])

        self.make_tmpdir_command = 'mkdir -p $BATCH_TMPDIR/count_out'

        self.finalise_outputs_command = (
            f'ln {self.tmp_output} {output_path} && ln {self.tmp_output_summary} {summary_path}'
        )

    def __str__(self):
        return ' && '.join(
            [
                self.make_tmpdir_command,
                ' '.join(self.command),
                self.finalise_outputs_command,
            ],
        )

    def __repr__(self):
        return self.__str__()


def count(
    input_bam: BamPath,
    output_path: Path,
    summary_path: Path,
    job_attrs: dict[str, str],
    sg_id: str,
) -> Job:
    """
    Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
    """

    b = hail_batch.get_batch()

    # Localise input
    input_bam_reads = b.read_input_group(
        **{
            'bam': input_bam,
            'bam.bai': f'{input_bam}.bai',
        }
    ).bam

    counting_reference = b.read_input(config.config_retrieve(['references', 'star', 'gtf']))

    # Create job
    j = b.new_bash_job(f'count_{sg_id}', attributes=job_attrs | {'tool': 'featureCounts'})
    j.image(config.config_retrieve(['images', 'subread']))

    # Set resource requirements
    res = STANDARD.set_resources(
        j=j,
        ncpu=config.config_retrieve(['workflow', 'count_cpu'], 8),
        storage_gb=50,  # TODO: make configurable
    )

    # Declare output resource group
    j.declare_resource_group(
        count_output={
            'count': '{root}.count',
            'count.summary': '{root}.count.summary',
        },
    )

    # Create counting command
    fc = FeatureCounts(
        input_bam=input_bam_reads,
        gtf_file=counting_reference,
        output_path=j.count_output['count'],
        summary_path=j.count_output['count.summary'],
        paired_end=True,
        feature_type='exon',
        attribute='gene_id',
        strandness='reverse',  # TODO: confirm this
        multi_mapping=False,  # TODO: determine default value
        min_quality=None,  # TODO: determine default value
        primary_only=True,  # TODO: determine default value
        ignore_duplicates=False,  # TODO: determine default value
        count_pairs=True,  # TODO: determine default value
        both_ends_mapped=True,  # TODO: determine default value
        both_ends_same_chr=True,  # TODO: determine default value
        threads=res.get_nthreads(),
    )

    # Add command to job
    j.command(hail_batch.command(str(fc), monitor_space=True))

    # Write output to file
    b.write_output(j.count_output['count'], str(output_path))
    b.write_output(j.count_output['count.summary'], str(summary_path))

    return j
