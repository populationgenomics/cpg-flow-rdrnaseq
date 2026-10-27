"""
Batch jobs to run MultiQC, check thresholds, and record QC flags.
"""

from collections import defaultdict
from os.path import basename

from hailtop.batch import Batch, Resource, ResourceFile
from hailtop.batch.job import Job
from loguru import logger

from cpg_flow import resources, targets, utils
from cpg_utils import Path, config, hail_batch
from metamist.graphql import gql, query

QC_ANALYSES_QUERY = gql(
    """
    query QcAnalyses($dataset: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: "qc"}, status: {eq: COMPLETED}) {
                output
                meta
                timestampCompleted
                sequencingGroups {
                    id
                }
            }
        }
    }
""",
)

MULTIQC_QC_STAGES = {'FastpQC', 'FastqScreen', 'SamtoolsStats', 'PicardRnaSeqMetrics'}


def _get_registered_qc_paths(dataset: targets.Dataset) -> dict[str, list[str]]:
    """Query Metamist for all registered QC outputs in the dataset.

    Returns a dict mapping SG ID -> list of output paths (one per QC stage).
    Deduplicates by taking the latest analysis per SG per stage.
    """
    dataset_name = dataset.name
    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in dataset_name:
        dataset_name += '-test'

    result = query(QC_ANALYSES_QUERY, variables={'dataset': dataset_name})

    latest: dict[tuple[str, str], tuple[str, str]] = {}
    for analysis in result.get('project', {}).get('analyses', []):
        meta = analysis.get('meta', {})
        stage_name = meta.get('stage')
        if stage_name not in MULTIQC_QC_STAGES:
            continue
        sg_ids = [sg['id'] for sg in analysis.get('sequencingGroups', [])]
        if len(sg_ids) != 1:
            continue
        sg_id = sg_ids[0]
        ts = analysis.get('timestampCompleted', '')
        key = (sg_id, stage_name)
        if key not in latest or ts > latest[key][0]:
            latest[key] = (ts, analysis['output'])

    sg_paths: dict[str, list[str]] = defaultdict(list)
    for (sg_id, _stage_name), (_ts, output) in latest.items():
        sg_paths[sg_id].append(output)

    return dict(sg_paths)


def _add_fastp_paths(sg_paths: dict[str, list[str]], dataset: targets.Dataset) -> None:
    """Derive fastp JSON paths by convention and add to sg_paths if the file exists."""
    fastp_prefix = dataset.prefix() / 'qc' / 'fastp'
    for sg_id in list(sg_paths):
        fastp_path = str(fastp_prefix / f'{sg_id}.fastp.json')
        if utils.can_reuse(fastp_path):
            sg_paths[sg_id].append(fastp_path)


def multiqc(
    dataset: targets.Dataset,
    outputs: dict[str, Path],
    out_html_url: str,
    out_checks_path: Path,
    label: str,
    ending_to_trim: set[str],
    modules_to_trim_endings: set[str],
    job_attrs: dict,
) -> list[Job]:
    """Run MultiQC, then check thresholds, then record flags in Metamist."""
    batch_instance = hail_batch.get_batch()
    tmp_prefix = dataset.tmp_prefix() / 'multiqc' / label
    sequencing_group_id_map = dataset.rich_id_map()

    sg_paths = _get_registered_qc_paths(dataset)
    _add_fastp_paths(sg_paths, dataset)
    logger.info(f'MultiQC for {dataset.name}: {len(sg_paths)} SGs with registered QC outputs')

    all_paths = {p for sg_files in sg_paths.values() for p in sg_files}
    if not all_paths:
        logger.warning(f'No QC outputs found for dataset {dataset.name}, skipping MultiQC')
        return []

    title = f'MultiQC [{label}]'

    mqc_j = batch_instance.new_bash_job(title, job_attrs | {'tool': 'MultiQC'})
    mqc_j.image(config.config_retrieve(['images', 'multiqc']))
    mqc_j.cpu(8)
    mqc_j.storage('20Gi')

    sg_id_mapping_file_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_rename-sg-map.tsv'

    dry_run = config.config_retrieve(['workflow', 'dry_run'], False)
    if not dry_run:
        with sg_id_mapping_file_path.open('w') as fh:
            for sgid, new_sgid in sequencing_group_id_map.items():
                fh.write('\t'.join([sgid, new_sgid]) + '\n')

    sg_id_mapping_file = batch_instance.read_input(sg_id_mapping_file_path)

    qc_inputs = batch_instance.read_input_group(**{basename(p): p for p in sorted(all_paths)})

    joined_endings = ', '.join(ending_to_trim)
    joined_modules = ', '.join(modules_to_trim_endings)

    mqc_j.command(
        f"""
        multiqc -f $(dirname {qc_inputs}) -o output \\
        --replace-names {sg_id_mapping_file} \\
        --title "{title} for dataset <b>{dataset.name}</b>" \\
        --filename report.html \\
        --cl-config "extra_fn_clean_exts: [{joined_endings}]" \\
        --cl-config "max_table_rows: 10000" \\
        --cl-config "use_filename_as_sample_name: [{joined_modules}]" \\

        cp output/report.html {mqc_j.html}
        cp output/report_data/multiqc_data.json {mqc_j.json}
        """
    )
    mqc_j.command(f'echo "HTML URL: {out_html_url}"')

    batch_instance.write_output(mqc_j.html, outputs['html'])
    batch_instance.write_output(mqc_j.html, outputs['latest'])
    batch_instance.write_output(mqc_j.json, outputs['json'])
    all_jobs: list[Job] = [mqc_j]

    if config.config_retrieve(['qc_thresholds'], {}):
        check_j = _check_report_job(
            b=batch_instance,
            multiqc_json_file=mqc_j.json,
            multiqc_html_url=out_html_url,
            rich_id_map=dataset.rich_id_map(),
            dataset_name=dataset.name,
            label=label,
            out_checks_path=out_checks_path,
            job_attrs=job_attrs,
        )
        check_j.depends_on(mqc_j)
        all_jobs.append(check_j)

        record_j = _record_qc_flags_job(
            b=batch_instance,
            dataset_name=dataset.name,
            label=label,
            sg_id_mapping_file=sg_id_mapping_file,
            check_multiqc_json_file=check_j.output,
            job_attrs=job_attrs,
        )
        record_j.depends_on(check_j)
        all_jobs.append(record_j)

    return all_jobs


def _rich_sequencing_group_id_seds(
    rich_id_map: dict[str, str],
    file_names: list[str | Resource | ResourceFile],
) -> str:
    """Extend sequencing group IDs in files with external IDs via sed."""
    cmd = ''
    for sgid, rich_sgid in rich_id_map.items():
        for fname in file_names:
            cmd += f"sed -iBAK 's/{sgid}/{rich_sgid}/g' {fname}\n"
    return cmd


def _check_report_job(
    b: Batch,
    multiqc_json_file: Resource | ResourceFile,
    dataset_name: str,
    multiqc_html_url: str,
    label: str,
    out_checks_path: Path,
    rich_id_map: dict[str, str] | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """Check MultiQC JSON and send Slack notification."""
    title = f'MultiQC [{label}]' if label else 'MultiQC'
    check_j = b.new_job(f'{title} check', (job_attrs or {}) | {'tool': 'python'})
    resources.STANDARD.set_resources(j=check_j, ncpu=2)
    check_j.image(config.config_retrieve(['workflow', 'driver_image']))

    cmd = f"""\
    {_rich_sequencing_group_id_seds(rich_id_map, [multiqc_json_file]) if rich_id_map else ''}

    python3 -m rdrnaseq.scripts.check_multiqc \\
    --multiqc-json {multiqc_json_file} \\
    --html-url {multiqc_html_url} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --output-json {check_j.output}
    """

    check_j.command(cmd)
    b.write_output(check_j.output, str(out_checks_path))
    return check_j


def _record_qc_flags_job(
    b: Batch,
    dataset_name: str,
    label: str,
    sg_id_mapping_file: Resource | ResourceFile,
    check_multiqc_json_file: Resource | ResourceFile,
    job_attrs: dict | None = None,
) -> Job:
    """Record QC flags in Metamist."""
    record_j = b.new_job('Record QC flags', (job_attrs or {}) | {'tool': 'python'})
    resources.STANDARD.set_resources(j=record_j, ncpu=2)
    record_j.image(config.config_retrieve(['workflow', 'driver_image']))

    cmd = f"""\
    python3 -m rdrnaseq.scripts.record_qc_flags \\
    --dataset {dataset_name} \\
    --label {label} \\
    --qc-flags-json {check_multiqc_json_file} \\
    --sequencing-group-ids-map {sg_id_mapping_file}
    """

    record_j.command(cmd)
    return record_j
