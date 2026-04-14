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
    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(METADATA_QUERY, variables={'project': query_dataset, 'sgIds': sg_ids})

    cpg_metadata: dict[str, dict[str, str | int]] = {}
    for group in result.get('project', {}).get('sequencingGroups', []):
        cpg_id = group.get('id')
        try:
            participant = group['sample']['participant']
            cpg_metadata[cpg_id] = {
                'family_id': participant['families'][0]['externalId'],
                'external_id': participant['externalId'],
                'affected': participant['familyParticipants'][0]['affected'],
            }
        except (KeyError, IndexError, TypeError):
            if cpg_id in sg_ids:
                logger.warning(f'Missing metadata for requested ID {cpg_id}')
            continue

    return cpg_metadata


def make_dashboards(
    fraser_csv: str | Path,
    outrider_csv: str | Path,
    output_html: str | Path,
    sg_ids_by_dataset: dict[str, list[str]],
    cohort_id: str,
    job_attrs: dict,
) -> list[Job]:
    """
    Create one Hail Batch job per dataset that renders an interactive dashboard.

    Each job localises the full cohort-level Fraser and Outrider CSVs, writes a
    family-mapping CSV scoped to that dataset's SG IDs, and runs the dashboard
    CLI script. The script itself filters results to matching sample IDs via
    the family mapping.
    """
    b = get_batch()
    access_level = config_retrieve(['workflow', 'access_level'], 'main')

    if access_level == 'test':
        fraser_csv = str(fraser_csv).replace('main', 'test')
        outrider_csv = str(outrider_csv).replace('main', 'test')

    fraser_input = b.read_input(str(fraser_csv))
    outrider_input = b.read_input(str(outrider_csv))

    jobs: list[Job] = []
    for dataset_name, sg_ids in sg_ids_by_dataset.items():
        cpg_metadata = get_cpg_metadata(dataset_name, sg_ids)

        j = b.new_job(f'rna_dashboard_{dataset_name}_{cohort_id}', attributes=job_attrs | {'tool': 'rna_dashboard'})

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
    --output {j.output}
"""),
        )

        b.write_output(j.output, str(output_html))

        web_path = (
            f'https://{access_level}-web.populationgenomics.org.au'
            f'/{dataset_name}/rna_dashboard/{cohort_id}.rna_dashboard.html'
        )
        logger.info(f'Dashboard job created for dataset {dataset_name}: {web_path}')
        jobs.append(j)

    return jobs
