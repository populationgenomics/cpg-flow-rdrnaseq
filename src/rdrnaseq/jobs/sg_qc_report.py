"""
Batch job to generate the SequencingGroup QC HTML report from Metamist QC flags.
"""

from hailtop.batch.job import Job

from cpg_utils import Path, config, hail_batch


def sg_qc_report_job(
    dataset: str,
    outputs: dict[str, Path],
    out_html_url: str,
    multiqc_url: str,
    job_attrs: dict,
) -> Job:
    """Create a Hail Batch job that queries Metamist for all RNA QC flags and renders a report."""
    batch = hail_batch.get_batch()

    j = batch.new_bash_job(f'SG QC Report: {dataset}', job_attrs | {'tool': 'python'})
    j.image(config.config_retrieve(['workflow', 'driver_image'])).memory('standard').cpu(2)

    j.command(
        f"""\
    python3 -m rdrnaseq.scripts.sg_qc_report \
        --dataset {dataset} \
        --fixed-output {outputs['html']} \
        --timestamped-output {outputs['timestamped']} \
        --html-url {out_html_url} \
        --multiqc-url {multiqc_url}
    """
    )

    return j
