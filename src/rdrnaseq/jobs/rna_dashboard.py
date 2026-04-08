
from typing import TYPE_CHECKING
from collections import defaultdict
import loguru
from cpg_flow import targets
from cpg_utils import Path, config, hail_batch, to_path
from metamist.graphql import gql, query

from pyarrow.dataset import dataset
from cpg-flow-rdrnaseq.scripts import create_interactive_dashboard, dashboard_utilities.py

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob

"A query that checks which project each of the sequencing groups in my cohort belong to"

project_query = gql(
"""Pedigree($cohort: String!) { 
  cohorts(id: {eq: "$cohort"}) {
    sequencingGroups {
      id
      sample {
        project {
          name
        }
      }
    }
  }
}
"""
)




"""
Create Hail Batch jobs to run STRipy
"""


# not used, needs correction
REPORT_TEMPLATE_PATH = './cpg_flow_stripy/stripy_report_template.html'

COMBINED_QUERY = gql(
    """
    query Pedigree($project: String!, $sgIds: [String!]!) {
        project(name: $project) {
            sequencingGroups(id: {in_: $sgIds}}) {
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


def get_cpg_metadata(dataset: str, relevant_ids: list[str]) -> dict[str, dict[str, str | int]]:
    """
    Returns a dictionary mapping cpgID to metadata:
    {cpgID: {"family_id": str, "external_id": str, "affected": int or str}}
    """

    # Handle test environment naming conventions
    query_dataset = dataset
    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    variables = {'project': query_dataset, 'sgIds': relevant_ids}
    result = query(COMBINED_QUERY, variables=variables)

    cpg_metadata = {}

    sequencing_groups = result.get('project', {}).get('sequencingGroups', [])

    for group in sequencing_groups:
        cpg_id = group.get('id')
        try:
            participant = group['sample']['participant']
            ext_id = participant['externalId']
            family_id = participant['families'][0]['externalId']
            affected = participant['familyParticipants'][0]['affected']
            cpg_metadata[cpg_id] = {
                'family_id': family_id,
                'external_id': ext_id,
                'affected': affected,
            }
        except (KeyError, IndexError, TypeError):
            if cpg_id in relevant_ids:
                print(f'Warning: Missing metadata for requested ID {cpg_id}')
            continue

    return cpg_metadata

def make_dashboard(
    fraser_results: Path,
    outrider_results: Path,
    output_dashboard_path: Path,
    cohort_id: str,
    output_latest: Path,
    job_attrs: dict,
) -> 'BashJob':
    """Makes an Dashboard for all rnaseq results."""
    batch_instance = hail_batch.get_batch()

    variables = {'cohort': cohort_id}
    result = query(COMBINED_QUERY, variables=variables)

    dataset_dict = defaultdict(list)
    for sg in result['data']['cohorts'][0]['sequencingGroups']:
        project = sg['sample']['project']['name']
        dataset_dict[project].append(sg['id'])
    list_of_jobs = []

    for dataset_name, sg_ids in dataset_dict:
        j = batch_instance.new_bash_job(name=f'Make RNA-seq Results Dashboard for {dataset_name}', attributes=job_attrs)
        j.image(config.config_retrieve(['workflow', 'driver_image']))


        cpg_metadata = get_cpg_metadata(dataset_name, sg_ids)

        file_prefix = config.config_retrieve(['storage', dataset_name, 'web'])
        html_prefix = config.config_retrieve(['storage', dataset_name, 'web_url'])

        # an object to store all the content we need to write
        collected_lines: list[str] = []
        for cpg_id, output_dict in inputs.items():
            fam_id = cpg_metadata[cpg_id]['family_id']
            external_id = cpg_metadata[cpg_id]['external_id']
            affected = cpg_metadata[cpg_id]['affected']
            # possible values for affected:0(unknown), 1(unaffected), 2(affected) -9(unknown)
            match affected:
                case 1:
                    affected_status = 'Unaffected'
                case 2:
                    affected_status = 'Affected'
                case 0 | -9 | 'Unknown':
                    affected_status = 'Unknown'
                case _:
                    affected_status = 'Unknown'
                    loguru.logger.warning(f'Unexpected status {affected} for ID {cpg_id}.')

            for report_type, report_path in output_dict.items():
                # substitute the report HTML path for a proxy-rendered path
                corrected_path = str(report_path).replace(file_prefix, html_prefix)
                collected_lines.append(
                    f'{cpg_id}\t{fam_id}\t{external_id}\t{report_type}\t{corrected_path}\t{affected_status}'
                )

        # write all reports to a single temp file, instead of passing an arbitrary number of CLI/script arguments
        with to_path(all_reports).open('w') as f:
            f.write('\n'.join(collected_lines))

        # localise that file
        mega_input_file = hail_batch.get_batch().read_input(all_reports)
        # --- Job Command (SINGLE STEP) ---
        # Runs your script, telling it to write to the local VM path
        j.command(f"""
            python3 -m cpg_flow_stripy.scripts.make_stripy_index \\
            --manifest {mega_input_file} \\
            --dataset {dataset_name} \\
            --output {j.output} \\
            --logfile {j.biglog}
        """)
        batch_instance.write_output(j.output, output_archive)
        batch_instance.write_output(j.output, output_latest)

        corrected_path_index = str(output_archive).replace(file_prefix, html_prefix)
        corrected_path_latest = str(output_latest).replace(file_prefix, html_prefix)

        loguru.logger.info(f'Index page job created for dataset {dataset_name} at {corrected_path_index}')
        loguru.logger.info(f'latest page job created for dataset {dataset_name} at {corrected_path_latest}')

        list_of_jobs.append(j)

    return j
