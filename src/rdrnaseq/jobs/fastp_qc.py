"""
Run fastp in QC-only mode and check the output against thresholds.
"""

from hailtop.batch.job import Job

from cpg_flow.filetypes import FastqPairs
from cpg_flow.resources import STANDARD
from cpg_utils import Path, config
from cpg_utils.config import image_path
from cpg_utils.hail_batch import command, get_batch

from rdrnaseq.jobs.trim import AdapterPairs


def fastp_qc(
    input_fq_pairs: FastqPairs,
    sg_id: str,
    qc_json_path: Path,
    qc_html_path: Path,
    status_path: Path,
    job_attrs: dict[str, str],
) -> list[Job]:
    """Run fastp in QC-only mode and check thresholds.

    Returns a list of jobs. The final job writes the status file.
    """
    b = get_batch()
    trim_config = config.config_retrieve('trim', {})
    adapter_type = trim_config.get('adapter_type', '')
    nthreads = trim_config.get('nthreads', 8)

    try:
        adapters = AdapterPairs[adapter_type].value
    except KeyError as e:
        raise ValueError(f'Invalid adapter type: {adapter_type}') from e

    jobs: list[Job] = []
    fastp_json_resources = []

    for idx, fq_pair in enumerate(input_fq_pairs):
        fq_resources = fq_pair.as_resources(b)

        suffix = f'_{idx}' if len(input_fq_pairs) > 1 else ''
        fastp_j = b.new_bash_job(f'FastpQC{suffix}', job_attrs | {'tool': 'fastp'})
        fastp_j.image(image_path('fastp'))
        res = STANDARD.set_resources(j=fastp_j, ncpu=nthreads, storage_gb=30)

        fastp_cmd = (
            f'fastp'
            f' --in1 {fq_resources.r1}'
            f' --in2 {fq_resources.r2}'
            f' --json {fastp_j.qc_json}'
            f' --html {fastp_j.qc_html}'
            f' --adapter_sequence {adapters.r1.sequence}'
            f' --adapter_sequence_r2 {adapters.r2.sequence}'
            f' --thread {res.get_nthreads()}'
        )
        fastp_j.command(command(fastp_cmd, monitor_space=True))
        jobs.append(fastp_j)
        fastp_json_resources.append(fastp_j.qc_json)

        if idx == 0:
            b.write_output(fastp_j.qc_json, str(qc_json_path))
            b.write_output(fastp_j.qc_html, str(qc_html_path))

    check_j = b.new_bash_job('FastpQC check', job_attrs | {'tool': 'python'})
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))
    STANDARD.set_resources(j=check_j, ncpu=2)
    for j in jobs:
        check_j.depends_on(j)

    json_args = ' '.join(f'--fastp-json {rsc}' for rsc in fastp_json_resources)
    check_j.command(
        f'python3 -m rdrnaseq.scripts.check_fastp_qc {json_args} --sample-id {sg_id} --output-status {check_j.status}'
    )
    b.write_output(check_j.status, str(status_path))
    jobs.append(check_j)

    return jobs
