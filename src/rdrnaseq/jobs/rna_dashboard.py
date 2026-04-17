"""
Create an interactive HTML dashboard for RNA-seq FRASER/OUTRIDER results.

This job module creates one Hail Batch job per dataset in the cohort. Each job
localises the cohort-level Fraser/Outrider CSVs, subsets to the dataset's SG IDs
via family-mapping metadata, and runs the dashboard CLI script to produce a
self-contained HTML file.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import command, get_batch
from loguru import logger
from metamist.graphql import gql, query

if TYPE_CHECKING:
    from cpg_utils import Path
    from hailtop.batch.job import Job


METADATA_QUERY = gql(
    """
    query Pedigree($project: String!, $sgIds: [String!]!) {
        project(name: $project) {
            sequencingGroups(id: {in_: $sgIds}) {
                id
                sample {
                    participant {
                        externalId
                        families {
                            externalId
                        }
                        familyParticipants {
                            affected
                        }
                    }
                }
            }
        }
    }
    """
)


def get_cpg_metadata(dataset_name: str, sg_ids: list[str]) -> dict[str, dict[str, str | int]]:
    """
    Query metamist for SG ID -> family/participant metadata.

    Returns a dict mapping cpg_id to {family_id, external_id, affected}.
    Runs at stage construction time (not inside the batch job).
    """
    query_dataset = dataset_name

    logger.info(f'Querying metamist project={query_dataset!r} for {len(sg_ids)} SG IDs: {sg_ids}')
    result = query(METADATA_QUERY, variables={'project': query_dataset, 'sgIds': sg_ids})

    returned_sgs = result.get('project', {}).get('sequencingGroups', [])
    logger.info(f'Metamist returned {len(returned_sgs)} sequencing groups')

    cpg_metadata: dict[str, dict[str, str | int]] = {}
    for group in returned_sgs:
        cpg_id = group.get('id')
        try:
            participant = group['sample']['participant']
            cpg_metadata[cpg_id] = {
                'family_id': participant['families'][0]['externalId'],
                'external_id': participant['externalId'],
                'affected': participant['familyParticipants'][0]['affected'],
            }
        except (KeyError, IndexError, TypeError) as e:
            logger.warning(f'Missing metadata for {cpg_id}: {type(e).__name__}: {e}')
            continue

    logger.info(f'Built metadata for {len(cpg_metadata)}/{len(sg_ids)} SG IDs')
    return cpg_metadata


def make_dashboards(
    fraser_csv: str | Path,
    outrider_csv: str | Path,
    output_paths: dict[str, str | Path],
    sg_ids_by_dataset: dict[str, list[str]],
    cohort_id: str,
    job_attrs: dict,
) -> list[Job]:
    """
    Create one Hail Batch job per dataset that renders an interactive dashboard.

    Each job localises the full cohort-level Fraser and Outrider CSVs, writes a
    family-mapping CSV scoped to that dataset's SG IDs, and runs the dashboard
    CLI script. The script writes an HTML file plus companion FRASER/OUTRIDER
    CSVs in the same directory — the HTML loads them at runtime via relative
    fetch() + PapaParse.

    output_paths keys: 'dashboard_html', 'fraser_csv', 'outrider_csv'
    """
    b = get_batch()
    access_level = config_retrieve(['workflow', 'access_level'], 'main')

    if access_level == 'test':
        fraser_csv = str(fraser_csv).replace('main', 'test')
        outrider_csv = str(outrider_csv).replace('main', 'test')

    fraser_input = b.read_input(str(fraser_csv))
    outrider_input = b.read_input(str(outrider_csv))

    ensg_to_symbol_path = config_retrieve(['references', 'ensg_to_symbol'])
    ensg_input = b.read_input(ensg_to_symbol_path)

    jobs: list[Job] = []
    for dataset_name, sg_ids in sg_ids_by_dataset.items():
        logger.info(f'Processing dataset {dataset_name} with SG IDs: {sg_ids}')
        cpg_metadata = get_cpg_metadata(dataset_name, sg_ids)

        j = b.new_job(f'rna_dashboard_{dataset_name}_{cohort_id}', attributes=job_attrs | {'tool': 'rna_dashboard'})
        j.image(config_retrieve(['images', 'rdrnaseq']))
        j.declare_resource_group(
            out={
                'dashboard_html': f'{cohort_id}.rna_dashboard.html',
                'fraser_csv': f'{cohort_id}.rna_dashboard.fraser.csv',
                'outrider_csv': f'{cohort_id}.rna_dashboard.outrider.csv',
            },
        )

        # Build a family-mapping CSV scoped to this dataset's SG IDs
        csv_lines = ['sequencing_group.id,family.external_ids']
        for cpg_id, meta in cpg_metadata.items():
            csv_lines.append(f'{cpg_id},{meta["family_id"]}')
        family_csv_content = '\n'.join(csv_lines)

        j.command(
            command(f"""\
cat > /tmp/family_mapping.csv << 'FAMILY_EOF'
{family_csv_content}
FAMILY_EOF

python3 -m rdrnaseq.scripts.create_interactive_dashboard \
    --fraser {fraser_input} \
    --outrider {outrider_input} \
    --family-mapping /tmp/family_mapping.csv \
    --ensg-to-symbol {ensg_input}\
    --output {j.out.dashboard_html} \
    --output-fraser-csv {j.out.fraser_csv} \
    --output-outrider-csv {j.out.outrider_csv}
"""),
        )

        for k, p in output_paths.items():
            b.write_output(j.out[k], str(p))

        web_path = f'https://{access_level}-web.populationgenomics.org.au/{dataset_name}/transcriptome/rna_dashboard/'
        logger.info(f'Dashboard job created for dataset {dataset_name}: {web_path}')
        jobs.append(j)

    return jobs
