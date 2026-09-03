"""
Per-sample fastp QC threshold check.

Reads fastp JSON output and compares key metrics against configurable
thresholds. Writes a plain-text status file used by the alignment gate.
"""

import json
import logging
from dataclasses import asdict
from datetime import datetime
from typing import Any

import click

from cpg_utils import config

from rdrnaseq.utils import DIRECTIONS, QcFlag, load_thresholds, worst_breach

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

# Maps config metric keys to (json_path, description) in fastp JSON.
# json_path is a dot-separated path into the fastp JSON structure.
FASTP_METRIC_MAP: dict[str, tuple[str, str]] = {
    'after_filtering_total_reads': ('summary.after_filtering.total_reads', 'Total reads after filtering'),
    'after_filtering_total_bases': ('summary.after_filtering.total_bases', 'Total bases after filtering'),
    'after_filtering_q20_rate': ('summary.after_filtering.q20_rate', 'Q20 rate after filtering'),
    'after_filtering_q30_rate': ('summary.after_filtering.q30_rate', 'Q30 rate after filtering'),
    'after_filtering_gc_content': ('summary.after_filtering.gc_content', 'GC content after filtering'),
    'before_filtering_total_reads': ('summary.before_filtering.total_reads', 'Total reads before filtering'),
    'before_filtering_q30_rate': ('summary.before_filtering.q30_rate', 'Q30 rate before filtering'),
    'duplication_rate': ('duplication.rate', 'Duplication rate'),
    'adapter_trimmed_reads_rate': (
        'adapter_cutting.adapter_trimmed_reads_rate',
        'Adapter trimmed reads rate',
    ),
    'insert_size_peak': ('insert_size.peak', 'Insert size peak'),
}


def resolve_json_path(data: dict, path: str) -> float | None:
    """Walk a dot-separated path into a nested dict, returning None if any key is missing."""
    current: Any = data
    for key in path.split('.'):
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    try:
        return float(current)
    except (TypeError, ValueError):
        return None


def extract_metrics(fastp_data: dict) -> dict[str, float]:
    """Extract configured metrics from a fastp JSON dict."""
    metrics = {}
    for metric_key, (json_path, _desc) in FASTP_METRIC_MAP.items():
        val = resolve_json_path(fastp_data, json_path)
        if val is not None:
            metrics[metric_key] = val
    return metrics


@click.command()
@click.option('--fastp-json', 'fastp_json_paths', multiple=True, required=True, help='Path(s) to fastp JSON output')
@click.option('--output-status', required=True, help='Path to write status file (PASS/FAIL)')
@click.option('--output-json', default=None, help='Optional path to write structured flags JSON')
@click.option('--sample-id', required=True, help='Sample/sequencing group ID')
def main(
    fastp_json_paths: tuple[str, ...],
    output_status: str,
    output_json: str | None,
    sample_id: str,
):
    """Check fastp QC metrics against configured thresholds."""
    today = datetime.now()  # noqa: DTZ005
    seq_type = config.config_retrieve(['workflow', 'sequencing_type'])
    thresholds = load_thresholds(seq_type)

    if not any(thresholds[d] for d in DIRECTIONS):
        logging.warning('No qc_thresholds configured; all samples will pass.')

    all_flags: list[QcFlag] = []

    for json_path in fastp_json_paths:
        with open(json_path) as f:
            fastp_data = json.load(f)

        metrics = extract_metrics(fastp_data)

        for direction, metric_tiers in thresholds.items():
            sign = DIRECTIONS[direction][0]
            for metric, tiers in metric_tiers.items():
                if metric not in metrics:
                    if metric in FASTP_METRIC_MAP:
                        logging.warning(f'{sample_id}: metric {metric!r} not found in fastp JSON')
                    continue
                val = metrics[metric]
                verdict = worst_breach(val, tiers, direction)
                if verdict is None:
                    logging.info(f'{sample_id}: {metric}={val:.4f} within thresholds')
                    continue
                severity, threshold = verdict
                icon = '!' if severity == 'fail' else '?'
                logging.info(f'{icon} {sample_id}: {metric}={val:.4f}{sign}{threshold:.4f} [{severity}]')
                all_flags.append(
                    QcFlag(
                        flag=metric,
                        value=val,
                        comparison=sign,
                        threshold=threshold,
                        section='fastp',
                        date=today.isoformat(timespec='seconds'),
                        ar_guid=config.try_get_ar_guid(),
                        severity=severity,
                    ),
                )

    has_fail = any(f.severity == 'fail' for f in all_flags)
    block_failed: bool = config.config_retrieve(['workflow', 'fastp_qc', 'block_failed_samples'], True)
    skip_gate_for: list[str] = config.config_retrieve(['workflow', 'fastp_qc', 'skip_gate_for'], [])

    if has_fail and not block_failed:
        logging.warning(f'{sample_id}: would FAIL but block_failed_samples is false — forcing PASS')
        status = 'PASS'
    elif has_fail and sample_id in skip_gate_for:
        logging.warning(f'{sample_id}: would FAIL but is in skip_gate_for override list — forcing PASS')
        status = 'PASS'
    else:
        status = 'FAIL' if has_fail else 'PASS'

    with open(output_status, 'w') as f:
        f.write(f'{status}\n')
        for flag in all_flags:
            f.write(f'{flag.flag}={flag.value:.4f} {flag.comparison} {flag.threshold:.4f} [{flag.severity}]\n')

    if output_json:
        result = {
            'sample_id': sample_id,
            'status': status,
            'qc_flags': {sample_id: [asdict(flag) for flag in all_flags]},
        }
        with open(output_json, 'w') as f:
            json.dump(result, f, indent=2)

    n_fail = sum(1 for f in all_flags if f.severity == 'fail')
    logging.info(f'{sample_id}: QC status = {status} ({len(all_flags)} flags, {n_fail} fail)')


if __name__ == '__main__':
    main()
