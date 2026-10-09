"""
Run samtools stats on a CRAM file to produce post-alignment QC metrics.
"""

from hailtop.batch.job import Job

from cpg_flow.resources import STANDARD
from cpg_utils import Path, config
from cpg_utils.hail_batch import command, get_batch


def samtools_stats(
    input_cram: str | Path,
    output_stats: Path,
    job_attrs: dict[str, str],
) -> Job:
    """Run samtools stats on a CRAM and write the output."""
    b = get_batch()
    ref_fasta = config.config_retrieve(['references', 'ref_fasta'])

    j = b.new_bash_job('SamtoolsStats', job_attrs | {'tool': 'samtools'})
    j.image(config.config_retrieve(['images', 'samtools']))
    STANDARD.set_resources(j=j, ncpu=4, storage_gb=30)

    cram_input = b.read_input_group(**{'cram': str(input_cram), 'cram.crai': f'{input_cram}.crai'})
    ref_input = b.read_input_group(**{'fasta': ref_fasta, 'fasta.fai': f'{ref_fasta}.fai'})

    j.command(
        command(
            f'samtools stats --reference {ref_input.fasta} --threads 3 {cram_input.cram} > {j.stats}',
            monitor_space=True,
        ),
    )
    b.write_output(j.stats, str(output_stats))

    return j
