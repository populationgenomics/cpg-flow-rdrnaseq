"""
Create an interactive dashboard with volcano plot and IGV browser.

Produces an HTML file plus companion FRASER/OUTRIDER CSV files in the same
directory.  The HTML loads the CSVs at runtime via relative fetch() + PapaParse,
keeping the HTML lightweight regardless of data size.

Usage:
    python create_interactive_dashboard.py --fraser fraser_results.csv --output dashboard.html
    python create_interactive_dashboard.py --fraser fraser.csv --outrider outrider.csv --output dashboard.html

Output files (all in the same directory as --output):
    <stem>.html          — the dashboard
    <stem>.fraser.csv    — processed FRASER results
    <stem>.outrider.csv  — processed OUTRIDER results (if --outrider given)
"""

import argparse
import ast
import csv
import logging
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from rdrnaseq.scripts.dashboard_utilities import (
    SeqrVariantLinkEngine,
    add_family_ids,
    build_ensg_to_hgnc_subset,
    enrich_with_gene_mapping,
    load_cpg_to_family_mapping,
    load_ensg_to_symbol,
    load_fraser_data,
    load_outrider_data,
)

# =============================================================================
# HTML Generation
# =============================================================================


def render_dashboard(
    fraser_csv_filename: str,
    outrider_csv_filename: str | None,
    output_path: str,
    family_map: dict,
    ensg_to_hgnc: dict,
    pvalue_threshold: float = 0.05,
    deltapsi_threshold: float = 0.2,
    zscore_threshold: float = 2.0,
    variant_bed_filename: str | None = None,
    variant_tsv_filename: str | None = None,
    seqr_links: dict[str, str] | None = None,
    template_name: str = 'interactive_dashboard.html.j2',
    display_name: str = '',
    cell_library_type: str = '',
    sg_ids: list[str] | None = None,
) -> None:
    """Render the dashboard HTML using Jinja2 template.

    Only small configuration data is embedded in the HTML.  The bulk data
    (FRASER / OUTRIDER rows) is loaded client-side from companion CSV files.
    """
    template_dir = Path(__file__).parent / 'templates'
    env = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
    template = env.get_template(template_name)

    display_cell_library_type = cell_library_type.replace('_', ' ').replace('-', ' ').title() if cell_library_type else ''

    html_content = template.render(
        fraser_csv_filename=fraser_csv_filename,
        outrider_csv_filename=outrider_csv_filename,
        variant_bed_filename=variant_bed_filename,
        variant_tsv_filename=variant_tsv_filename,
        family_map=family_map,
        ensg_to_hgnc=ensg_to_hgnc,
        default_pvalue_threshold=pvalue_threshold,
        default_deltapsi_threshold=deltapsi_threshold,
        default_zscore_threshold=zscore_threshold,
        seqr_links=seqr_links or {},
        genome='hg38',
        display_name=display_name,
        cell_library_type=display_cell_library_type,
        sg_ids=sg_ids or [],
    )

    with open(output_path, 'w') as f:
        f.write(html_content)

    print(f'  Dashboard written to: {output_path}')


# =============================================================================
# CLI Interface
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Create interactive FRASER/OUTRIDER results dashboard',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with FRASER only
    python create_interactive_dashboard.py --fraser results.csv --output dashboard.html

    # With OUTRIDER results
    python create_interactive_dashboard.py --fraser fraser.csv --outrider outrider.csv --output dashboard.html

    # With family mapping
    python create_interactive_dashboard.py --fraser fraser.csv --outrider outrider.csv \\
        --family-mapping families.csv --output dashboard.html

    # Explicit CSV output paths (used by Hail Batch jobs)
    python create_interactive_dashboard.py --fraser fraser.csv --output dashboard.html \\
        --output-fraser-csv /tmp/out.fraser.csv --output-outrider-csv /tmp/out.outrider.csv
        """,
    )

    parser.add_argument('--fraser', '-f', required=True, help='Path to FRASER results file (CSV or gzipped TSV)')
    parser.add_argument('--outrider', '-r', default=None, help='Path to OUTRIDER results file (CSV or gzipped TSV)')
    parser.add_argument('--output', '-o', required=True, help='Output HTML file path')

    parser.add_argument(
        '--output-fraser-csv',
        default=None,
        help='Output path for processed FRASER CSV (default: derived from --output)',
    )
    parser.add_argument(
        '--output-outrider-csv',
        default=None,
        help='Output path for processed OUTRIDER CSV (default: derived from --output)',
    )

    parser.add_argument(
        '--pvalue-threshold', type=float, default=0.05, help='Default p-value threshold (default: 0.05)'
    )
    parser.add_argument(
        '--deltapsi-threshold', type=float, default=0.2, help='Default |Delta PSI| threshold (default: 0.2)'
    )
    parser.add_argument(
        '--zscore-threshold', type=float, default=2.0, help='Default |Z-Score| threshold for OUTRIDER (default: 2.0)'
    )

    parser.add_argument(
        '--family-mapping', default=None, help='Path to CSV with sequencing_group.id and family.external_ids columns'
    )
    parser.add_argument(
        '--ensg-to-symbol', required=True, help='Path to ENSG-to-HGNC-symbol TSV mapping file (two columns, no header)'
    )

    parser.add_argument(
        '--variant-bed-filename',
        default=None,
        help='Filename of variant annotation BED (already in output dir, for IGV.js track)',
    )
    parser.add_argument(
        '--variant-tsv-filename',
        default=None,
        help='Filename of variant annotation TSV (already in output dir, for data table)',
    )
    parser.add_argument(
        '--private-output',
        default=None,
        help='Output path for private dashboard HTML ',
    )
    parser.add_argument('--dataset-name', required=True, help='CPG dataset name (e.g. rdnow)')
    parser.add_argument('--cohort-id', required=True, help='Cohort ID (e.g. COH10509)')
    parser.add_argument('--display-name', default='', help='Project display name from metamist')
    parser.add_argument('--cell-library-type', default='', help='Cell/library type (e.g. fibroblast_polyA)')
    parser.add_argument('--sg-ids', nargs='*', default=[], help='Sequencing group IDs analysed')

    return parser.parse_args()


def apply_family_mapping(mapping_file, fraser_df, outrider_df):
    """Load family mapping, filter DataFrames to mapped SG IDs, and add metadata columns."""
    cpg_to_family = load_cpg_to_family_mapping(mapping_file)
    sg_ids = set(cpg_to_family.keys())
    print(f'\nFiltering to {len(sg_ids)} sample IDs from family mapping...')

    pre = len(fraser_df)
    fraser_df = fraser_df[fraser_df['sampleID'].isin(sg_ids)].copy()
    print(f'  FRASER: {pre} -> {len(fraser_df)} rows')

    if outrider_df is not None:
        pre = len(outrider_df)
        outrider_df = outrider_df[outrider_df['sampleID'].isin(sg_ids)].copy()
        print(f'  OUTRIDER: {pre} -> {len(outrider_df)} rows')

    fraser_df = add_family_ids(fraser_df, cpg_to_family)
    if outrider_df is not None:
        outrider_df = add_family_ids(outrider_df, cpg_to_family)

    return cpg_to_family, fraser_df, outrider_df


def main() -> None:
    """Main entry point."""
    args = parse_args()

    # ── Load data ────────────────────────────────────────────────────────────
    fraser_df = load_fraser_data(args.fraser)
    outrider_df = load_outrider_data(args.outrider) if args.outrider else None

    # ── Load and apply family mapping ────────────────────────────────────────
    cpg_to_family: dict[str, dict] = {}
    if args.family_mapping:
        cpg_to_family, fraser_df, outrider_df = apply_family_mapping(args.family_mapping, fraser_df, outrider_df)

    # ── Load ENSG-to-symbol mapping and enrich DataFrames ──────────────────
    print('\nLoading ENSG-to-symbol mapping...')
    ensg_to_symbol = load_ensg_to_symbol(args.ensg_to_symbol)
    print('\nEnriching DataFrames with gene name/ID columns...')
    fraser_df, outrider_df = enrich_with_gene_mapping(fraser_df, outrider_df, ensg_to_symbol)

    # ── Determine output paths ───────────────────────────────────────────────
    output_html = Path(args.output)
    output_dir = output_html.parent
    base_name = output_html.stem  # e.g. "COH10509.rna_dashboard"

    fraser_csv_path = args.output_fraser_csv or str(output_dir / f'{base_name}.fraser.csv')

    outrider_csv_path: str | None = None
    if outrider_df is not None:
        outrider_csv_path = args.output_outrider_csv or str(output_dir / f'{base_name}.outrider.csv')

    # ── Write processed CSVs ─────────────────────────────────────────────────
    print('\nWriting processed CSVs...')
    fraser_df.to_csv(fraser_csv_path, index=False)
    print(f'  FRASER CSV: {fraser_csv_path} ({len(fraser_df)} rows)')

    if outrider_df is not None and outrider_csv_path:
        outrider_df.to_csv(outrider_csv_path, index=False)
        print(f'  OUTRIDER CSV: {outrider_csv_path} ({len(outrider_df)} rows)')

    # ── Build embedded data for JS ───────────────────────────────────────────
    seqr_links: dict[str, str] = {}
    if args.variant_tsv_filename:
        link_engine = SeqrVariantLinkEngine(args.dataset_name)
        with open(str(output_dir / args.variant_tsv_filename)) as variants:
            reader = csv.DictReader(variants, delimiter='\t')
            for variant in reader:
                # grab locus
                loc_parts = variant['locus'].split(':')
                chrom = loc_parts[0].replace('chr', '')
                pos = loc_parts[1]
                # grab alleles and sanitize for use in URL
                alleles = ast.literal_eval(variant['alleles'])
                variant_id = f'{chrom}-{pos}-{"-".join(alleles)}'
                family_id = variant['family_id']
                link = link_engine.build_link(family_id, variant_id)
                if link is not None:
                    seqr_links[f'{variant_id}|{family_id}'] = link
    else:
        logging.warning('\nNo variant TSV provided, skipping Seqr link generation.')

    ensg_to_hgnc: dict[str, str] = {}
    if outrider_df is not None and ensg_to_symbol:
        ensg_to_hgnc = build_ensg_to_hgnc_subset(outrider_df, ensg_to_symbol)

    # ── Render dashboard HTML ────────────────────────────────────────────────
    render_kwargs = {
        'fraser_csv_filename': Path(fraser_csv_path).name,
        'outrider_csv_filename': Path(outrider_csv_path).name if outrider_csv_path else None,
        'family_map': cpg_to_family,
        'ensg_to_hgnc': ensg_to_hgnc,
        'pvalue_threshold': args.pvalue_threshold,
        'deltapsi_threshold': args.deltapsi_threshold,
        'zscore_threshold': args.zscore_threshold,
        'variant_bed_filename': args.variant_bed_filename,
        'variant_tsv_filename': args.variant_tsv_filename,
        'seqr_links': seqr_links,
        'display_name': args.display_name,
        'cell_library_type': args.cell_library_type,
        'sg_ids': args.sg_ids,
    }

    render_dashboard(output_path=args.output, **render_kwargs)

    render_dashboard(
        output_path=args.private_output,
        template_name='private_interactive_dashboard.html.j2',
        **render_kwargs,
    )
    print(
        f'  HTML: https://main-web.populationgenomics.org.au/{args.dataset_name}/transcriptome/rna_dashboard/{args.cohort_id}.rna_dashboard.html'
    )

    print(f'  FRASER CSV:  {fraser_csv_path}')
    if outrider_csv_path:
        print(f'  OUTRIDER CSV: {outrider_csv_path}')


if __name__ == '__main__':
    main()
