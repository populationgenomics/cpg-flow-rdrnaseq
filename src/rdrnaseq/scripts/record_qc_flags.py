import json
from dataclasses import asdict
from datetime import datetime

import click
from loguru import logger

from metamist.graphql import gql, query

from rdrnaseq.utils import QcFlag

DATASET_SG_META_QUERY = gql(
    """
    query datasetSgMeta($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups {
                id
                meta
            }
        }
    }
    """
)

SG_META_MUTATION = gql(
    """
    mutation updateSgMeta($dataset: String!, $sgId: String!, $sgMeta: JSON!) {
        sequencingGroup {
            updateSequencingGroup(
                project: $dataset
                sequencingGroup: {id: $sgId, meta: $sgMeta}
            ) {
                id
                meta
            }
        }
    }
    """
)


def compare_qc_flag(current_flag: dict, new_flag: dict) -> bool:
    """Returns True if the two flags refer to the same QC issue and the current flag is still unresolved."""
    return (
        current_flag.get('flag') == new_flag.get('flag')
        and current_flag.get('comparison') == new_flag.get('comparison')
        and current_flag.get('threshold') == new_flag.get('threshold')
        and current_flag.get('section') == new_flag.get('section')
        and not current_flag.get('resolved', False)
    )


def reconcile_sg_qc_flags(
    sg: dict,
    new_flags_by_sg: dict[str, list[dict]],
    dataset: str,
    label: str,
    today: datetime,
) -> None:
    """Reconcile current and new QC flags for a single SG and update Metamist."""
    sg_id = sg['id']
    report = f'{label}MultiQC'

    qc_flags_key = f'{label.lower()}_qc_flags'
    current_qc_flags: list[dict] = (sg['meta'] or {}).get(qc_flags_key, [])
    unresolved_current_flags = [flag for flag in current_qc_flags if not flag.get('resolved', False)]

    new_qc_flags: list[dict] = new_flags_by_sg.get(sg_id, [])

    if not current_qc_flags and not new_qc_flags:
        logger.info(f'{sg_id} :: No existing or new {report} flags for this SG, skipping.')
        return

    if not unresolved_current_flags and not new_qc_flags:
        logger.info(f'{sg_id} :: No unresolved existing or new {report} flags for this SG, skipping.')
        return

    new_qc_flags_by_key = {(flag['section'], flag['flag']): flag for flag in new_qc_flags}
    existing_flag_keys = {(flag['section'], flag['flag']) for flag in current_qc_flags}

    final_flags: list[QcFlag] = []
    stats = {'resolved': 0, 'retained': 0, 'updated': 0, 'added': 0}

    if current_qc_flags:
        logger.info(f'{sg_id} :: Found {len(current_qc_flags)} existing {report} flags. Reconciling.')
        for flag in current_qc_flags:
            key = (flag['section'], flag['flag'])
            present_in_new = key in new_qc_flags_by_key
            if not present_in_new:
                if not flag['resolved']:
                    flag['resolved'] = True
                    flag['resolution_date'] = today.isoformat(timespec='seconds')
                    logger.info(f"{sg_id} :: Marking {report} flag '{flag['flag']}' as resolved.")
                    stats['resolved'] += 1
            elif compare_qc_flag(flag, new_qc_flags_by_key[key]):
                new_flag = new_qc_flags_by_key[key]
                flag['value'] = new_flag['value']
                flag['severity'] = new_flag.get('severity', flag.get('severity', 'fail'))
                logger.info(f"{sg_id} :: {report} flag '{flag['flag']}' remains unresolved (value refreshed).")
                stats['retained'] += 1
            else:
                old_severity = flag.get('severity', 'fail')
                flag.update(new_qc_flags_by_key[key])
                transition = f' ({old_severity} -> {flag["severity"]})' if old_severity != flag.get('severity') else ''
                logger.info(f"{sg_id} :: {report} flag '{flag['flag']}' updated{transition}.")
                stats['updated'] += 1
            final_flags.append(QcFlag(**flag))
    else:
        logger.info(f'{sg_id} :: No existing {report} flags, adding {len(new_qc_flags)} new flags.')

    for flag in new_qc_flags:
        if (flag['section'], flag['flag']) in existing_flag_keys:
            continue
        logger.info(f"{sg_id} :: Adding new {report} flag '{flag['flag']}'.")
        final_flags.append(QcFlag(**flag))
        stats['added'] += 1

    query(
        SG_META_MUTATION,
        variables={
            'dataset': dataset,
            'sgId': sg_id,
            'sgMeta': {qc_flags_key: [asdict(flag) for flag in final_flags]},
        },
    )
    logger.info(
        f'{sg_id} :: Recorded {len(final_flags)} {report} flags. '
        f'Resolved: {stats["resolved"]}, Retained: {stats["retained"]}, '
        f'Updated: {stats["updated"]}, Added: {stats["added"]}'
    )


@click.command()
@click.option('--dataset', required=True, help='Dataset name')
@click.option('--label', required=True, help='Report type label')
@click.option('--qc-flags-json', 'qc_flags_json_path', required=True, help='Path to QC flags JSON')
@click.option('--sequencing-group-ids-map', 'sg_id_mapping_file', required=True, help='Path to SG IDs mapping TSV')
def main(
    dataset: str,
    label: str,
    qc_flags_json_path: str,
    sg_id_mapping_file: str,
):
    """Read QC flags JSON and update SG meta in Metamist."""
    today = datetime.now()  # noqa: DTZ005

    with open(sg_id_mapping_file) as f:
        sg_id_map = dict(line.strip().split('\t') for line in f if line.strip())

    with open(qc_flags_json_path) as f:
        qc_flags_data = json.load(f)

    response = query(DATASET_SG_META_QUERY, variables={'dataset': dataset})
    sequencing_groups = response['project']['sequencingGroups']

    new_flags_by_sg = qc_flags_data['qc_flags']
    for sg in sequencing_groups:
        if sg['id'] in sg_id_map:
            reconcile_sg_qc_flags(sg, new_flags_by_sg, dataset, label, today)


if __name__ == '__main__':
    main()
