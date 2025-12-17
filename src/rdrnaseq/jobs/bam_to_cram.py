"""
Convert BAM to CRAM.
"""

# ruff: noqa: E501
from cpg_flow.resources import STANDARD
from cpg_utils import Path
from cpg_utils.config import config_retrieve, image_path
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

    if not isinstance(input_bam, ResourceGroup):
        raise TypeError(f'Expected input_bam to be a ResourceGroup, but got {type(input_bam).__name__}')

    j = b.new_job(name='bam_to_cram', attributes=job_attrs | {'tool': 'samtools'})
    j.image(image_path('samtools'))

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
        storage_gb=config_retrieve(['resource_overrides', 'bam_to_cram', 'storage_gb'], 50),
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
    input_cram: ResourceGroup,
    output_bam: Path | None = None,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
    reference_fasta_path: str | None = None,
) -> tuple[Job, ResourceGroup]:
    """
    Convert a CRAM file to a BAM file.
    """
    b = get_batch()

    # Get fasta file
    fasta = b.read_input_group(
        fasta=reference_fasta_path,
        fasta_fai=f'{reference_fasta_path}.fai',
    )
    if not isinstance(input_cram, ResourceGroup):
        raise TypeError(f'Expected input_cram to be a ResourceGroup, but got {type(input_cram).__name__}')

    job_name = 'cram_to_bam'
    if extra_label:
        job_name += f' {extra_label}'

    convert_tool = 'samtools_view_cram_to_bam'
    j_attrs = (job_attrs or {}) | {'label': job_name, 'tool': convert_tool}
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j=j,
        ncpu=nthreads,
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
    if output_bam:
        b.write_output(j.sorted_bam, str(output_bam.with_suffix('')))

    return j, j.sorted_bam
