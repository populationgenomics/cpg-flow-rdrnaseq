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
from metamist import graphql

LONG_READ_STRING = 'LongRead'
METAMIST_ANALYSIS_QUERY = graphql.gql(
    """
    query MyQuery($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
                meta
            }
        }
    }
""",
)


def query_for_latest_analysis(
    dataset: str,
    analysis_type: str,
    sequencing_type: str = 'all',
    long_read: bool = False,
    stage_name: str | None = None,
) -> str | None:
    """
    Query for the latest analysis object of a given type in the requested project.

    Analysis entries for Talos all have unique types, so we can use this generic query method

    Args:
        dataset (str):         project to query for
        analysis_type (str):   analysis type to query for - rd_combiner writes MTs to metamist as 'matrixtable',
                               seqr_loader used 'custom': using a config entry we can decide which type to use
        sequencing_type (str): optional, if set, only return entries with meta.sequencing_type == this
        long_read (bool):      if True, will skip over any entries that are not LongRead (SNPsIndels/SV)
        stage_name (str):      optional, if set, will only return entries with meta.stage == this
    Returns:
        str, the path to the latest object for the given type, or log a warning and return None
    """

    # swapping to a string we can freely modify
    query_dataset = dataset
    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    logger.info(f'Querying for {analysis_type} in {query_dataset}')

    result = graphql.query(METAMIST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': analysis_type})

    # get all the relevant entries, and bin by date
    analysis_by_date = {}
    for analysis in result['project']['analyses']:
        if analysis['output'] and (sequencing_type in {'all', analysis['meta'].get('sequencing_type')}):
            # skip over the partial-cohort AnnotateDataset objects
            if '_families-' in analysis['output']:
                logger.debug(
                    f'Skipping analysis {analysis["output"]} for dataset {query_dataset}. '
                    f'It is a partial-cohort AnnotateDataset object',
                )
                continue

            # manually implementing an XOR check - long read (bool) and LongRead in output must match
            if long_read != (LONG_READ_STRING in analysis['output']):
                logger.debug(
                    f'Skipping analysis {analysis["output"]} for dataset {query_dataset}. '
                    f'It does not match query parameter long_read={long_read}',
                )
                continue

            if stage_name is not None and analysis['meta'].get('stage') != stage_name:
                continue

            analysis_by_date[analysis['timestampCompleted']] = analysis['output']

    if not analysis_by_date:
        logger.warning(f'No Analysis Entries found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return analysis_by_date[sorted(analysis_by_date)[-1]]


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
        mt_path = config.config_retrieve(['variant_annotation', 'mt_path', str(dataset_name)], default=None)
        logger.info(f'Config lookup for variant_annotation.mt_path.{dataset_name}: {mt_path!r}')
        if mt_path is None:
            analysis_type = config.config_retrieve(
                ['workflow', 'variant_annotation_mt_analysis_type'], default='matrixtable'
            )
            seq_type = config.config_retrieve(['workflow', 'sequencing_type'], default='genome')
            long_read = config.config_retrieve(['workflow', 'long_read'], default=False)
            stage_name = config.config_retrieve(
                ['workflow', 'variant_annotation_stage_name'], default='AnnotateDataset'
            )
            logger.info(
                f'Config lookup returned None, querying metamist for latest analysis: '
                f'analysis_type={analysis_type!r}, sequencing_type={seq_type!r}, '
                f'long_read={long_read!r}, stage_name={stage_name!r}',
            )
            mt_path = query_for_latest_analysis(
                dataset=dataset_name,
                analysis_type=analysis_type,
                sequencing_type=seq_type,
                long_read=long_read,
                stage_name=stage_name,
            )
            logger.info(f'query_for_latest_analysis returned: {mt_path!r}')
        if mt_path is None:
            logger.error(f'No MT path found for dataset {dataset_name} — skipping job')
            continue

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
