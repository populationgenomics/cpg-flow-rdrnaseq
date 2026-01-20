"""
Convert BAM to CRAM.
"""

# ruff: noqa: E501
from cpg_flow.resources import STANDARD
from cpg_utils import Path, config
from cpg_utils.hail_batch import command, get_batch
from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job


def bam_to_cram(
    input_bam: ResourceGroup,
    job_attrs: dict,
    requested_nthreads: int | None = None,
    reference_fasta_path: str | None = None,
) -> tuple[Job, ResourceGroup]:
    """
    Convert a BAM file to a CRAM file.
    """

    b = get_batch()

    j = b.new_bash_job(name='bam_to_cram', attributes=job_attrs | {'tool': 'samtools'})
    j.image(config.config_retrieve(['images', 'samtools']))

    # Get fasta file
    fasta = b.read_input_group(
        fasta=reference_fasta_path,
        fasta_fai=f'{reference_fasta_path}.fai',
    )

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j=j,
        ncpu=nthreads,
        storage_gb=config.config_retrieve(['resource_overrides', 'bam_to_cram', 'storage_gb'], 50),
    )

    j.declare_resource_group(
        sorted_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        },
    )

    cmd = f'samtools view -@ {res.get_nthreads() - 1} -T {fasta.fasta} -C {input_bam.bam} \
        | tee {j.sorted_cram["cram"]} \
        | samtools index -@ {res.get_nthreads() - 1} - {j.sorted_cram["cram.crai"]}'
    j.command(command(cmd, monitor_space=True))

    return j, j.sorted_cram


def cram_to_bam(
    input_cram_path: Path,
    output_bam: str,
    job_attrs: dict[str, str],
) -> Job:
    """
    Convert a CRAM file to a BAM file.
    """
    b = get_batch()

    reference_fasta_path = config.config_retrieve(['references', 'ref_fasta'])

    input_cram = b.read_input_group(
        **{
            'cram': str(input_cram_path),
            'cram.crai': f'{input_cram_path}.crai',
        },
    )

    # Get fasta file
    fasta = b.read_input_group(
        fasta=reference_fasta_path,
        fasta_fai=f'{reference_fasta_path}.fai',
    )

    j = b.new_bash_job(name='cram_to_bam', attributes=job_attrs | {'tool': 'samtools'})
    j.image(config.config_retrieve(['images', 'samtools']))

    # Set resource requirements
    res = STANDARD.set_resources(
        j=j,
        ncpu=config.config_retrieve(['workflow', 'cram_to_bam_cpu'], 8),
        storage_gb=50,
    )

    j.declare_resource_group(
        sorted_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai',
        },
    )

    cmd = f"""samtools view -@ {res.get_nthreads() - 1} -T {fasta.fasta} -b {input_cram.cram} > {j.sorted_bam['bam']} && \
    samtools index -@ {res.get_nthreads() - 1} {j.sorted_bam['bam']} {j.sorted_bam['bam.bai']}"""
    j.command(command(cmd, monitor_space=True))

    # Write BAM if requested
    b.write_output(j.sorted_bam, output_bam.removesuffix('.bam'))

    return j
