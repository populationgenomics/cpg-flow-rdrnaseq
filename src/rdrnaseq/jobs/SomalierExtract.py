"""
Job for CRAM fingerprinting using Somalier.
"""

from cpg_utils import Path, config, hail_batch
from hailtop.batch.job import Job
from pygments.lexer import default


def extract_somalier(
    cram_path: str,
    output: Path,
    job_attrs: dict[str, str],
) -> Job:
    """Run `somalier extract` to generate a fingerprint (i.e. a `*.somalier` file)."""

    batch_instance = hail_batch.get_batch()
    default_job_storage= 50

    job = batch_instance.new_job('Somalier extract', job_attrs | {'tool': 'somalier'})

    job.image(config.config_retrieve(['images', 'somalier']))
    storage_gb = default_job_storage
    job.storage(f'{storage_gb}GB')

    ref = hail_batch.fasta_res_group(batch_instance)

    # read in the somalier sites VCF
    sites = batch_instance.read_input(config.config_retrieve(['references', 'somalier_sites']))

    # read in the CRAM and index
    cram_localised = batch_instance.read_input_group(
        cram=cram_path,
        crai=f'{cram_path}.crai',
    ).cram

    job.command(f"""
    somalier extract -d extracted/ --sites {sites} -f {ref.base} {cram_localised}
    mv extracted/*.somalier {job.output_file}
    """)
    batch_instance.write_output(job.output_file, output)
    return job
