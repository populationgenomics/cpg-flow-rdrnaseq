import argparse
from datetime import datetime, timezone

from loguru import logger

from cpg_utils import to_path
from metamist.apis import AnalysisApi
from metamist.graphql import gql, query
from metamist.models import Analysis, AnalysisStatus

RNASEQ_SG_QUERY = gql(
    """
    query rnaSgs($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups(type: {in_: ["polyarna", "totalrna"]}) {
                id
                analyses(type: {eq: "cram"}) {
                    id
                }
            }
        }
    }
    """
)


def main():
    parser = argparse.ArgumentParser(description='Register analyses for CRAM files')
    parser.add_argument('dataset', type=str, help='Dataset name to register analyses for')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    # Query for sequencing groups and their analyses
    result = query(RNASEQ_SG_QUERY, variables={'dataset': args.dataset})
    sequencing_groups = result['project']['sequencingGroups']

    analysis_api = AnalysisApi()

    for sg in sequencing_groups:
        sg_id = sg['id']
        analyses = sg['analyses']
        if not analyses:
            cram_path = to_path(f'gs://cpg-{args.dataset}-main/transcriptome/cram/{sg_id}.cram')
            if not cram_path.exists():
                logger.info(f'CRAM file for sequencing group {sg_id} does not exist at expected path: {cram_path}')
                continue
            timestamp_completed = datetime.fromtimestamp(cram_path.stat().st_mtime, tz=timezone.utc).strftime(
                '%Y-%m-%dT%H:%M:%SZ'
            )
            if not args.dry_run:
                logger.info(f'Registering analysis for sequencing group {sg_id} with CRAM path: {cram_path}')
                analysis_api.create_analysis(
                    project=args.dataset,
                    analysis=Analysis(
                        type='cram',
                        status=AnalysisStatus('completed'),
                        output=str(cram_path),
                        sequencing_group_ids=[sg_id],
                        timestamp_completed=timestamp_completed,
                    ),
                )
            else:
                logger.info(f'Dry run: would register analysis for sequencing group {sg_id}:')
                logger.info(f' {cram_path} with timestamp {timestamp_completed}')
        else:
            logger.info(f'Sequencing group {sg_id} already has an analysis registered')


if __name__ == '__main__':
    main()
