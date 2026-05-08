"""
Annotate FRASER significant regions with rare variants from a seqr-loader MatrixTable.

Creates one Hail Batch job per dataset in the cohort. Each job runs the
variant_of_interest_subset script, which queries Metamist to map RNA SG IDs
to genome SG IDs and subsets the MT to matching rare variants.
"""

from hailtop.batch.job import Job
from loguru import logger

from cpg_utils import Path, config
from cpg_utils.hail_batch import command, get_batch


def annotate_variants(
    fraser_csv: str | Path,
    mt_path: str,
    sg_ids_by_dataset: dict[str, list[str]],
    output: str,
    cohort_id: str,
    job_attrs: dict,
) -> list[Job]:
    """
    Create one Hail Batch job per dataset that runs variant_of_interest_subset.

    Each job localises the Fraser significant CSV, passes the MT path (read
    directly from GCS by Hail), and writes a BED + TSV output.
    """
    b = get_batch()
    fraser_input = b.read_input(fraser_csv)

    jobs: list[Job] = []
    for dataset_name, sg_ids in sg_ids_by_dataset.items():
        logger.info(f'Variant annotation for dataset {dataset_name} with {len(sg_ids)} RNA SG IDs')

        j = b.new_job(
            f'variant_annotation_{dataset_name}_{cohort_id}',
            attributes=job_attrs | {'tool': 'variant_annotation'},
        )
        j.image(config.config_retrieve('workflow')['driver_image'])

        rna_ids_str = ' '.join(sg_ids)

        j.command(
            command(f"""\
python3 -m rdrnaseq.scripts.variant_of_interest_subset \
    --mt {mt_path} \
    --csv {fraser_input} \
    --rna_ids {rna_ids_str} \
    --query_dataset {dataset_name} \
    --output {output}
"""),
        )
        jobs.append(j)

    return jobs
