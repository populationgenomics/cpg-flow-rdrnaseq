"""
Queries Metamist for RNA QC flags across a dataset's sequencing groups
and renders them into an HTML report.
"""

import re
from argparse import ArgumentParser
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import jinja2
from loguru import logger

from cpg_utils import to_path
from metamist.graphql import gql, query

from rdrnaseq.utils import QcFlag

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'

DATASET_SGS_QUERY = gql(
    """
    query datasetSgs($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups(type: {eq: "transcriptome"}) {
                id
                meta
            }
        }
    }
    """
)

METRIC_LABELS: dict[str, tuple[str, str, int]] = {
    'PCT_MRNA_BASES': ('mRNA bases', '%', 1),
    'PCT_INTRONIC_BASES': ('Intronic bases', '%', 1),
    'MEDIAN_5PRIME_TO_3PRIME_BIAS': ("5'-3' bias", '', 1),
    'MEDIAN_CV_COVERAGE': ('Coverage CV', '', 1),
    'reads_mapped_percent': ('Reads mapped', '%', 1),
}

SECTION_LABELS: dict[str, str] = {
    'fastp': 'fastp',
    'fastq_screen': 'FastQ Screen',
    'samtools': 'Samtools',
    'picard': 'Picard',
}

SEVERITY_RANK: dict[str, int] = {'fail': 0, 'warn': 1}


def _metric_label(metric: str) -> tuple[str, str, int]:
    return METRIC_LABELS.get(metric, (metric, '', 1))


def _section_label(section: str) -> str:
    base = re.sub(r'_\d+$', '', section or '')
    return SECTION_LABELS.get(base.lower(), base.replace('_', ' ').title() or '—')


def _fmt_num(n: float | str) -> str:
    if isinstance(n, bool) or not isinstance(n, int | float):
        return str(n)
    f = float(n)
    if f.is_integer():
        return str(int(f))
    s = f'{f:.2f}' if abs(f) >= 1 else f'{f:.2g}'
    if '.' in s and 'e' not in s.lower():
        s = s.rstrip('0').rstrip('.')
    return s


def _value_display(value: float, comparison: str, threshold: float, unit: str) -> str:
    v = f'{_fmt_num(value)}{unit}'
    t = f'{_fmt_num(threshold)}{unit}'
    if comparison == '<':
        return f'{v} (below minimum {t})'
    if comparison == '>':
        return f'{v} (above maximum {t})'
    return f'{v} {comparison} {t}'


@dataclass(frozen=True)
class SGReport:
    sg_id: str
    flags: list[QcFlag]


def _flag_to_dict(flag: QcFlag) -> dict:
    label, unit, multiplier = _metric_label(flag.flag)
    date_full = (flag.resolution_date if flag.resolved else flag.date) or ''
    severity = flag.severity or 'fail'
    return {
        'flag': flag.flag,
        'metric_label': label,
        'section': flag.section,
        'section_label': _section_label(flag.section),
        'resolved': flag.resolved,
        'severity': severity,
        'severity_label': 'Fail' if severity == 'fail' else 'Warn',
        'value_display': _value_display(flag.value * multiplier, flag.comparison, flag.threshold * multiplier, unit),
        'date_full': date_full,
        'date_short': date_full[:10],
    }


def build_sections(reports: list[SGReport]) -> tuple[list[dict], list[dict]]:
    unresolved, resolved = [], []
    for report in reports:
        all_flags = [_flag_to_dict(f) for f in report.flags]
        active = [f for f in all_flags if not f['resolved']]
        past = [f for f in all_flags if f['resolved']]

        if active:
            active.sort(key=lambda f: (SEVERITY_RANK.get(f['severity'], 0), f['date_full']))
            n_fail = sum(1 for f in active if f['severity'] == 'fail')
            n_warn = sum(1 for f in active if f['severity'] == 'warn')
            unresolved.append(
                {
                    'sg_id': report.sg_id,
                    'flags': active,
                    'count_summary': f'{len(active)} active flag' + ('' if len(active) == 1 else 's'),
                    'n_fail': n_fail,
                    'n_warn': n_warn,
                    'row_severity': 'fail' if n_fail else 'warn',
                    'sort_key': (0 if n_fail else 1, report.sg_id),
                }
            )
        if past:
            past.sort(key=lambda f: f['date_full'], reverse=True)
            resolved.append(
                {
                    'sg_id': report.sg_id,
                    'flags': past,
                    'count_summary': f'{len(past)} resolved flag' + ('' if len(past) == 1 else 's'),
                    'n_fail': 0,
                    'n_warn': 0,
                    'row_severity': '',
                    'sort_key': (2, report.sg_id),
                }
            )

    def _sort_key(row: dict) -> tuple[int, str]:
        return row['sort_key']

    unresolved.sort(key=_sort_key)
    resolved.sort(key=_sort_key)
    return unresolved, resolved


def collect_qc_flags(sequencing_groups: list[dict]) -> list[SGReport]:
    results = []
    for sg in sequencing_groups:
        meta = sg.get('meta') or {}
        raw_flags = meta.get('rna_qc_flags', [])
        if not raw_flags:
            continue
        flags = []
        for f in raw_flags:
            known = {k: f[k] for k in QcFlag.__dataclass_fields__ if k in f}
            try:
                flags.append(QcFlag(**known))
            except TypeError:
                logger.warning(f'Skipping malformed flag on {sg["id"]}: {f}')
        if flags:
            results.append(SGReport(sg_id=sg['id'], flags=flags))
    return results


def summarise_flags(reports: list[SGReport], total_sgs: int) -> dict:
    all_flags = [f for r in reports for f in r.flags]
    active = [f for f in all_flags if not f.resolved]
    return {
        'total_sgs': total_sgs,
        'active_flags': len(active),
        'active_fail': sum(1 for f in active if (f.severity or 'fail') == 'fail'),
        'active_warn': sum(1 for f in active if (f.severity or 'fail') == 'warn'),
        'sgs_affected': sum(1 for r in reports if any(not f.resolved for f in r.flags)),
        'resolved_flags': sum(1 for f in all_flags if f.resolved),
    }


def metric_histogram(rows: list[dict]) -> list[dict]:
    counts: dict[str, dict] = {}
    for row in rows:
        for key in {f['flag'] for f in row['flags']}:
            entry = counts.setdefault(key, {'key': key, 'label': _metric_label(key)[0], 'count': 0})
            entry['count'] += 1
    return sorted(counts.values(), key=lambda d: (-d['count'], d['label']))


def severity_histogram(rows: list[dict]) -> list[dict]:
    counts: dict[str, int] = {}
    for row in rows:
        for severity in {f['severity'] for f in row['flags']}:
            counts[severity] = counts.get(severity, 0) + 1
    labels = {'fail': 'Failing', 'warn': 'Warnings'}
    return [{'key': k, 'label': labels[k], 'count': counts[k]} for k in ('fail', 'warn') if k in counts]


def render_report(dataset: str, reports: list[SGReport], *, summary: dict, multiqc_url: str) -> str:
    unresolved, resolved = build_sections(reports)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('sg_qc_overview.html.jinja')
    return template.render(
        dataset=dataset,
        generated_at=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),  # noqa: DTZ005
        summary=summary,
        multiqc_url=multiqc_url,
        unresolved=unresolved,
        resolved=resolved,
        active_metrics=metric_histogram(unresolved),
        active_severities=severity_histogram(unresolved),
    )


def main(dataset: str, output: str, timestamped_output: str, out_html_url: str, multiqc_url: str):
    logger.info(f'{dataset} :: Querying Metamist for RNA QC flags')
    response = query(DATASET_SGS_QUERY, variables={'dataset': dataset})
    sequencing_groups = response['project']['sequencingGroups']
    logger.info(f'{dataset} :: Found {len(sequencing_groups)} sequencing groups')

    reports = collect_qc_flags(sequencing_groups)
    summary = summarise_flags(reports, total_sgs=len(sequencing_groups))
    logger.info(f'{dataset} :: {len(reports)} sequencing groups have QC flags')

    html = render_report(dataset, reports, summary=summary, multiqc_url=multiqc_url)

    with to_path(output).open('w') as f:
        f.write(html)
    logger.info(f'{dataset} :: Wrote SG QC report to {output}')

    with to_path(timestamped_output).open('w') as f:
        f.write(html)
    logger.info(f'{dataset} :: Wrote timestamped report to {timestamped_output}')

    logger.info(f'{dataset} :: Report available at {out_html_url}')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--fixed-output', required=True)
    parser.add_argument('--timestamped-output', required=True)
    parser.add_argument('--html-url', required=True)
    parser.add_argument('--multiqc-url', required=True)
    args = parser.parse_args()
    main(
        dataset=args.dataset,
        output=args.fixed_output,
        timestamped_output=args.timestamped_output,
        out_html_url=args.html_url,
        multiqc_url=args.multiqc_url,
    )
