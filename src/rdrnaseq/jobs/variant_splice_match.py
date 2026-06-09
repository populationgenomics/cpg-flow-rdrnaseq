"""
Annotate FRASER significant regions with rare variants from a seqr-loader MatrixTable.

Creates one Hail Batch job per dataset in the cohort. Each job runs the
variant_of_interest_subset script, which queries Metamist to map RNA SG IDs
to genome SG IDs and subsets the MT to matching rare variants.
"""

from hailtop.batch.job import Job
from loguru import logger

from cpg_flow.status import complete_analysis_job
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


    Args:
        dataset (str):         project to query for
        analysis_type (str):   analysis type to query for - rd_combiner writes MTs to metamist as 'matrixtable',
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


def register_multiple_analyses(outputs, analysis_type, cohort_ids, sg_ids, project_name, meta):
    for output in outputs:
        complete_analysis_job(output, analysis_type, cohort_ids, sg_ids, project_name, meta)


def match_variants_and_splicing(
    fraser_csv: str | Path,
    sg_ids_by_dataset: dict[str, list[str]],
    output_by_dataset: dict[str, str],
    cohort_id: str,
    job_attrs: dict,
) -> list[Job]:
    """
    Create one Hail Batch job per dataset that runs variant_of_interest_subset.

    Each job localises the Fraser significant CSV, passes the MT path (read
    directly from GCS by Hail), and writes a BED + TSV output per dataset.
    """
    b = get_batch()
    fraser_input = b.read_input(fraser_csv)

    jobs: list[Job] = []
    for dataset_name, sg_ids in sg_ids_by_dataset.items():
        logger.info(f'Variant annotation for dataset {dataset_name} with {len(sg_ids)} RNA SG IDs')
        mt_path = config.config_retrieve(['variant_splice_match', 'mt_path', str(dataset_name)], default=None)
        logger.info(f'Config lookup for variant_splice_match.mt_path.{dataset_name}: {mt_path!r}')
        if mt_path is None:
            analysis_type = 'matrixtable'

            seq_type = config.config_retrieve(['workflow', 'variant_splice_match_sequencing_type'], default='genome')
            long_read = config.config_retrieve(['workflow', 'long_read'], default=False)
            stage_name = config.config_retrieve(
                ['workflow', 'variant_splice_match_stage_name'], default='AnnotateDataset'
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

        rna_ids_str = ' '.join(sg_ids)
        dataset_output = output_by_dataset[dataset_name]
        coarse_output = output_by_dataset[f'coarse_{dataset_name}']

        j = b.new_job(
            f'variant_splice_match_{dataset_name}_{cohort_id}',
            attributes=job_attrs | {'tool': 'variant_splice_match'},
        )
        j.image(config.config_retrieve('workflow')['driver_image'])
        j.command(
            command(f"""\
python3 -m rdrnaseq.scripts.variant_of_interest_subset \
    --mt {mt_path} \
    --csv {fraser_input} \
    --rna_ids {rna_ids_str} \
    --query_dataset {dataset_name} \
    --output {coarse_output} \
    --output-tsv {j.output_tsv} \
    --output-bed {j.output_bed}
"""),
        )
        b.write_output(j.output_tsv, f'{dataset_output}.tsv')
        b.write_output(j.output_bed, f'{dataset_output}.bed')
        jobs.append(j)

        # --- Registration job (now depends on verify) ---
        registration_job = b.new_python_job(
            f'register_variant_splice_match_{dataset_name}_{cohort_id}',
            attributes=job_attrs | {'tool': 'metamist'},
        )
        registration_job.image(config.config_retrieve('workflow')['driver_image'])
        registration_job.call(
            register_multiple_analyses,
            outputs=[str(dataset_output) + '.bed', str(dataset_output) + '.tsv', str(coarse_output) + '.tsv'],
            analysis_type='variantsplicematch',
            cohort_ids=[cohort_id],
            sg_ids=sg_ids,
            project_name=dataset_name,
            meta={'stage': 'VariantSpliceMatch', 'dataset': dataset_name, 'cohort_id': cohort_id},
        )
        registration_job.depends_on(verify_j)
        jobs.append(registration_job)

    return jobs
