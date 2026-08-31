"""
Batch jobs to run MultiQC, check thresholds, and record QC flags.
"""

from hailtop.batch import Batch, ResourceFile
from hailtop.batch.job import Job

from cpg_flow import resources, targets
from cpg_utils import Path, config, hail_batch


def multiqc(
    dataset: targets.Dataset,
    tmp_prefix: Path,
    paths: list[Path | str],
    outputs: dict[str, Path],
    out_html_url: str | None,
    out_checks_path: Path | None,
    label: str,
    ending_to_trim: set[str],
    modules_to_trim_endings: set[str],
    job_attrs: dict,
    sequencing_group_id_map: dict[str, str],
    extra_config: dict | None = None,
    send_to_slack: bool = True,
) -> list[Job]:
    """Run MultiQC, then check thresholds, then record flags in Metamist."""
    batch_instance = hail_batch.get_batch()
    extra_config = extra_config or {}

    title = f'MultiQC [{label}]'

    mqc_j = batch_instance.new_bash_job(title, job_attrs | {'tool': 'MultiQC'})
    mqc_j.image(config.config_retrieve(['images', 'multiqc']))
    mqc_j.cpu(8)
    mqc_j.storage('20Gi')

    file_list_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_multiqc-file-list.txt'
    sg_id_mapping_file_path = tmp_prefix / f'{dataset.get_alignment_inputs_hash()}_rename-sg-map.tsv'

    dry_run = config.config_retrieve(['workflow', 'dry_run'], False)
    if not dry_run:
        with file_list_path.open('w') as f:
            f.writelines([f'{p}\n' for p in paths])
        with sg_id_mapping_file_path.open('w') as fh:
            for sgid, new_sgid in sequencing_group_id_map.items():
                fh.write('\t'.join([sgid, new_sgid]) + '\n')

    file_list = batch_instance.read_input(file_list_path)
    sg_id_mapping_file = batch_instance.read_input(sg_id_mapping_file_path)

    joined_endings = ', '.join(ending_to_trim)
    joined_modules = ', '.join(modules_to_trim_endings)

    extra_config_param = ''
    if extra_config:
        for k, v in extra_config.items():
            serialised = f'{k}: {v}'
            extra_config_param += f'--cl-config "{serialised}" \\\n            '

    mqc_j.command(
        f"""
        mkdir inputs
        cat {file_list} | gcloud storage cp -I inputs/

        multiqc -f inputs -o output \\
        --replace-names {sg_id_mapping_file} \\
        --title "{title} for dataset <b>{dataset.name}</b>" \\
        --filename report.html \\
        --cl-config "extra_fn_clean_exts: [{joined_endings}]" \\
        --cl-config "max_table_rows: 10000" \\
        --cl-config "use_filename_as_sample_name: [{joined_modules}]" \\
        {extra_config_param}

        cp output/report.html {mqc_j.html}
        cp output/report_data/multiqc_data.json {mqc_j.json}
        """
    )
    if out_html_url:
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
            send_to_slack=send_to_slack,
        )
        check_j.depends_on(mqc_j)
        all_jobs.append(check_j)

        dataset_name = config.dataset_for_access_level(dataset.name)
        record_j = _record_qc_flags_job(
            b=batch_instance,
            dataset_name=dataset_name,
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
    file_names: list[str | ResourceFile],
) -> str:
    """Extend sequencing group IDs in files with external IDs via sed."""
    cmd = ''
    for sgid, rich_sgid in rich_id_map.items():
        for fname in file_names:
            cmd += f"sed -iBAK 's/{sgid}/{rich_sgid}/g' {fname}\n"
    return cmd


def _check_report_job(
    b: Batch,
    multiqc_json_file: ResourceFile,
    dataset_name: str,
    multiqc_html_url: str | None = None,
    label: str | None = None,
    rich_id_map: dict[str, str] | None = None,
    out_checks_path: Path | None = None,
    job_attrs: dict | None = None,
    send_to_slack: bool = True,
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
    --output-json {check_j.output} \\
    --{'no-' if not send_to_slack else ''}send-to-slack
    """

    check_j.command(cmd)
    if out_checks_path:
        b.write_output(check_j.output, str(out_checks_path))
    return check_j


def _record_qc_flags_job(
    b: Batch,
    dataset_name: str,
    label: str,
    sg_id_mapping_file: ResourceFile,
    check_multiqc_json_file: ResourceFile,
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
