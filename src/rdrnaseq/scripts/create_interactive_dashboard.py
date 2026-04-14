"""
Create an interactive dashboard with volcano plot and IGV browser.

Supports FRASER and OUTRIDER results with optional TABIX indexing for efficient
region-based queries. Optionally annotates FRASER results with proximity to rare
variants from a Hail table.

Usage:
    python create_interactive_dashboard.py --fraser fraser_results.csv --output dashboard.html
    python create_interactive_dashboard.py --fraser fraser.tsv.gz --outrider outrider.tsv.gz --output dashboard.html
Input file formats:
    - CSV files: Will be converted to tabix-indexed format automatically
    - .gz files with .tbi index: Used directly for tabix queries
    - Hail tables (.ht): Used for variant proximity annotation

Required FRASER columns: seqnames, start, end, pValue, padjust, deltaPsi, psiValue, type, hgncSymbol
Required OUTRIDER columns: seqnames, start, end, pValue, padjust, (optional: zScore, log2fc, geneID/hgncSymbol)

Coordinate system: 1-based, inclusive (BED-like but 1-based for consistency with R/bioconductor)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader

from rdrnaseq.scripts.dashboard_utilities import (
    _col,
    _genes_match,
    add_family_ids,
    load_cpg_to_family_mapping,
    load_fraser_data,
    load_outrider_data,
    prepare_tabix_file,
)

# =============================================================================
# Common Results Identification
# =============================================================================


def _build_common_entry(sample: str, fraser_row: pd.Series, outrider_row: pd.Series) -> dict:
    """Build a common-significant result dict from a FRASER and OUTRIDER row pair."""
    entry = {
        'gene': fraser_row.get('hgncSymbol', 'NA'),
        'ensembl_id': outrider_row.get('geneID', outrider_row.get('hgncSymbol', 'NA')),
        'sampleID': sample,
        'chr': fraser_row.get('seqnames', 'NA'),
        'start': fraser_row.get('start', 0),
        'end': fraser_row.get('end', 0),
        'fraser_pvalue': fraser_row['pValue'],
        'fraser_padjust': fraser_row['padjust'],
        'fraser_deltaPsi': fraser_row['deltaPsi'],
        'fraser_psiValue': fraser_row.get('psiValue', 'NA'),
        'fraser_type': fraser_row.get('type', 'NA'),
        'outrider_pvalue': outrider_row['pValue'],
        'outrider_padjust': outrider_row['padjust'],
    }
    if 'zScore' in outrider_row:
        entry['outrider_zScore'] = outrider_row['zScore']
    if 'l2fc' in outrider_row:
        entry['outrider_l2fc'] = outrider_row['l2fc']
    elif 'log2fc' in outrider_row:
        entry['outrider_l2fc'] = outrider_row['log2fc']
    if 'familyID' in fraser_row.index:
        entry['familyID'] = fraser_row['familyID']
    return entry


def find_common_significant(
    df_fraser: pd.DataFrame, df_outrider: pd.DataFrame, pval_threshold: float = 0.05, deltapsi_threshold: float = 0.2
) -> tuple:
    """Find genes that are significant in both FRASER and OUTRIDER"""
    print('\nFinding common significant results...')

    # FRASER significant: meets p-value and deltaPSI thresholds
    fraser_sig = df_fraser[
        (df_fraser['padjust'] <= pval_threshold) & (np.abs(df_fraser['deltaPsi']) >= deltapsi_threshold)
    ].copy()

    # OUTRIDER significant: meets p-value and z-score thresholds (if zScore column exists)
    if 'zScore' in df_outrider.columns:
        outrider_sig = df_outrider[
            (df_outrider['padjust'] <= pval_threshold) & (np.abs(df_outrider['zScore']) >= 2)
        ].copy()
    else:
        outrider_sig = df_outrider[df_outrider['padjust'] <= pval_threshold].copy()

    print(f'  FRASER significant: {len(fraser_sig)}')
    print(f'  OUTRIDER significant: {len(outrider_sig)}')

    common_samples = set(fraser_sig['sampleID']) & set(outrider_sig['sampleID'])
    print(f'  Samples with both FRASER and OUTRIDER hits: {len(common_samples)}')

    common_by_sample = []
    for sample in common_samples:
        fraser_rows = fraser_sig[fraser_sig['sampleID'] == sample]
        outrider_rows = outrider_sig[outrider_sig['sampleID'] == sample]

        # Find gene-matched pairs
        found_gene_match = False
        for _, f_row in fraser_rows.iterrows():
            for _, o_row in outrider_rows.iterrows():
                if _genes_match(f_row, o_row):
                    common_by_sample.append(_build_common_entry(sample, f_row, o_row))
                    found_gene_match = True

        # If no gene match, pair the top hits from each to show the sample has both aberration types
        if not found_gene_match:
            f_top = fraser_rows.nsmallest(1, 'pValue').iloc[0]
            o_top = outrider_rows.nsmallest(1, 'pValue').iloc[0]
            common_by_sample.append(_build_common_entry(sample, f_top, o_top))

    df_common = pd.DataFrame(common_by_sample)
    print(f'  Common significant results: {len(df_common)}')
    print('  (Samples showing aberrations in both splicing and expression)')

    return df_common, fraser_sig, outrider_sig


def prepare_table_data(df: pd.DataFrame, columns_subset: list | None = None) -> pd.DataFrame:
    """Prepare data for searchable tables with formatted numeric values"""
    available_columns = [col for col in columns_subset if col in df.columns] if columns_subset else df.columns.tolist()

    table_df = df[available_columns].copy()

    # Format numeric columns
    for col in table_df.columns:
        if 'pValue' in col or 'padjust' in col or 'pvalue' in col:
            table_df[col] = table_df[col].apply(
                lambda x: f'{x:.2e}' if pd.notna(x) and isinstance(x, int | float) else str(x)
            )
        elif 'psi' in col.lower() or 'Score' in col or 'score' in col.lower() or 'l2fc' in col or 'deltaPsi' in col:
            table_df[col] = table_df[col].apply(
                lambda x: f'{x:.3f}' if pd.notna(x) and isinstance(x, int | float) else str(x)
            )

    # Fill NaN values
    return table_df.fillna('NA')


# =============================================================================
# Data Preparation for Dashboard
# =============================================================================


def get_top_positions(df: pd.DataFrame, n: int = 100) -> list[dict]:
    """Get top N positions by p-value for IGV."""
    top_df = df.nsmallest(n, 'pValue').copy()

    positions = []
    for _, row in top_df.iterrows():
        positions.append(
            {
                'chr': str(row['seqnames']),
                'start': int(row['start']),
                'end': int(row['end']),
                'gene': str(row['hgncSymbol']) if pd.notna(row['hgncSymbol']) else 'NA',
                'pValue': float(row['pValue']),
                'deltaPsi': float(row['deltaPsi']),
                'psiValue': float(row['psiValue']),
            }
        )

    return positions


def prepare_volcano_data(df: pd.DataFrame) -> dict:
    """Prepare data for volcano plot."""
    return {
        'deltaPsi': df['deltaPsi'].tolist(),
        'log10pValue': df['-log10(pValue)'].tolist(),
        'padjust': df['padjust'].tolist(),
        'gene': _col(df, 'hgncSymbol', fill='NA'),
        'chr': df['seqnames'].astype(str).tolist(),
        'start': df['start'].tolist(),
        'end': df['end'].tolist(),
        'psiValue': df['psiValue'].tolist(),
        'type': df['type'].tolist(),
        'sampleID': df['sampleID'].tolist(),
        'familyID': _col(df, 'familyID', default='Unknown', fill='Unknown'),
    }


def prepare_outrider_volcano_data(df: pd.DataFrame) -> dict:
    """Prepare data for OUTRIDER volcano plot (zScore vs -log10(pValue))."""
    # Add -log10(pValue) column if not present
    if '-log10(pValue)' not in df.columns:
        df['-log10(pValue)'] = -np.log10(df['pValue'].replace(0, 1e-300))

    # Resolve log2fc from either 'l2fc' or 'log2fc' column
    log2fc = _col(df, 'l2fc', default=0, fill=0)
    if 'l2fc' not in df.columns:
        log2fc = _col(df, 'log2fc', default=0, fill=0)

    return {
        'zScore': _col(df, 'zScore', default=0, fill=0),
        'log10pValue': df['-log10(pValue)'].tolist(),
        'padjust': df['padjust'].tolist(),
        'gene': _col(df, 'hgncSymbol', fill='NA') if 'hgncSymbol' in df.columns else _col(df, 'geneID', default='NA'),
        'geneID': _col(df, 'geneID', default='NA'),
        'log2fc': log2fc,
        'sampleID': df['sampleID'].tolist(),
        'familyID': _col(df, 'familyID', default='Unknown', fill='Unknown'),
    }


def prepare_outrider_top_genes_data(df: pd.DataFrame, n: int = 50) -> dict:
    """
    Prepare data for OUTRIDER top genes plot.
    Shows the top N genes by significance.
    """
    df = df.copy()
    df['gene_display'] = (
        df['hgncSymbol'].fillna(df.get('geneID', 'NA')) if 'hgncSymbol' in df.columns else df.get('geneID', 'NA')
    )

    top_df = df.nsmallest(n, 'pValue').copy()
    if '-log10(pValue)' not in top_df.columns:
        top_df['-log10(pValue)'] = -np.log10(top_df['pValue'].replace(0, 1e-300))

    log2fc = _col(top_df, 'l2fc', default=0, fill=0)
    if 'l2fc' not in top_df.columns:
        log2fc = _col(top_df, 'log2fc', default=0, fill=0)

    return {
        'genes': top_df['gene_display'].tolist(),
        'log10pValue': top_df['-log10(pValue)'].tolist(),
        'padjust': top_df['padjust'].tolist(),
        'zScore': _col(top_df, 'zScore', default=0, fill=0),
        'log2fc': log2fc,
        'sampleID': top_df['sampleID'].tolist(),
        'familyID': _col(top_df, 'familyID', default='Unknown', fill='Unknown'),
    }


def prepare_top_genes_data(df: pd.DataFrame, n: int = 100) -> dict[str, list]:
    """Prepare data for top genes plot."""
    top_df = df.nsmallest(n, 'pValue').copy()

    # Add -log10(pValue) if not present
    if '-log10(pValue)' not in top_df.columns:
        top_df['-log10(pValue)'] = -np.log10(top_df['pValue'].replace(0, 1e-300))

    return {
        'gene': _col(top_df, 'hgncSymbol', fill='NA'),
        'log10pValue': top_df['-log10(pValue)'].tolist(),
        'padjust': top_df['padjust'].tolist(),
        'deltaPsi': top_df['deltaPsi'].tolist(),
        'psiValue': top_df['psiValue'].tolist(),
        'type': top_df['type'].tolist(),
        'sampleID': top_df['sampleID'].tolist(),
        'familyID': _col(top_df, 'familyID', default='Unknown', fill='Unknown'),
    }


def prepare_tabix_data_for_js(df: pd.DataFrame, data_type: str = 'fraser') -> list[dict]:
    """
    Prepare data for JavaScript-side region queries.

    This creates a lightweight index that can be used for client-side
    region lookups without requiring server-side tabix queries.

    For OUTRIDER data, genomic coordinates (seqnames, start, end) are optional.
    Only records with coordinates will be included for region queries.
    """
    records = []

    for _, row in df.iterrows():
        # Check if genomic coordinates are present (required for region queries)
        if 'seqnames' not in row or pd.isna(row.get('seqnames')):
            continue
        if 'start' not in row or pd.isna(row.get('start')):
            continue
        if 'end' not in row or pd.isna(row.get('end')):
            continue

        record = {
            'chr': str(row['seqnames']),
            'start': int(row['start']),
            'end': int(row['end']),
            'gene': str(row.get('hgncSymbol', 'NA')) if pd.notna(row.get('hgncSymbol')) else 'NA',
            'pValue': float(row['pValue']) if pd.notna(row['pValue']) else None,
            'padjust': float(row['padjust']) if pd.notna(row['padjust']) else None,
        }

        if data_type == 'fraser':
            record['deltaPsi'] = float(row['deltaPsi']) if pd.notna(row.get('deltaPsi')) else None
            record['psiValue'] = float(row['psiValue']) if pd.notna(row.get('psiValue')) else None
            record['type'] = str(row['type']) if pd.notna(row.get('type')) else 'NA'

        elif data_type == 'outrider':
            record['zScore'] = float(row['zScore']) if pd.notna(row.get('zScore')) else None
            record['log2fc'] = float(row['log2fc']) if pd.notna(row.get('log2fc')) else None

        records.append(record)

    return records


def calculate_stats(fraser_df: pd.DataFrame, outrider_df: pd.DataFrame | None, pvalue_threshold: float) -> dict:
    """Calculate summary statistics for the dashboard."""
    stats = {
        'total_variants': len(fraser_df),
        'chromosomes': fraser_df['seqnames'].nunique(),
        'unique_genes': fraser_df['hgncSymbol'].nunique(),
        'significant_count': int((fraser_df['padjust'] < pvalue_threshold).sum()),
    }

    if outrider_df is not None:
        stats['outrider_total'] = len(outrider_df)

    return stats


# =============================================================================
# HTML Generation
# =============================================================================


def render_dashboard(  # noqa: PLR0912,PLR0915
    fraser_df: pd.DataFrame,
    outrider_df: pd.DataFrame | None,
    output_path: str,
    pvalue_threshold: float = 0.05,
    deltapsi_threshold: float = 0.2,
    zscore_threshold: float = 2.0,
) -> None:
    """Render the dashboard HTML using Jinja2 template."""

    print('Preparing data for dashboard...')

    # Prepare all data structures
    volcano_data = prepare_volcano_data(fraser_df)
    top_positions = get_top_positions(fraser_df, n=100)

    # Prepare OUTRIDER plot data if available
    outrider_volcano_data = None
    outrider_top_genes_data = None
    if outrider_df is not None and len(outrider_df) > 0:
        outrider_volcano_data = prepare_outrider_volcano_data(outrider_df)
        outrider_top_genes_data = prepare_outrider_top_genes_data(outrider_df, n=50)

    # Prepare tabix-like data for JS queries
    fraser_tabix_data = prepare_tabix_data_for_js(fraser_df, 'fraser')
    outrider_tabix_data = prepare_tabix_data_for_js(outrider_df, 'outrider') if outrider_df is not None else []

    # Find common significant results if OUTRIDER data is available
    df_common = pd.DataFrame()
    if outrider_df is not None and len(outrider_df) > 0:
        df_common, _, _ = find_common_significant(fraser_df, outrider_df, pvalue_threshold, deltapsi_threshold)

    # Prepare table data
    print('Preparing table data...')

    # FRASER table
    fraser_table_cols = [
        'hgncSymbol',
        'seqnames',
        'start',
        'end',
        'type',
        'pValue',
        'padjust',
        'psiValue',
        'deltaPsi',
        'sampleID',
    ]
    if 'familyID' in fraser_df.columns:
        fraser_table_cols.append('familyID')

    fraser_table = prepare_table_data(fraser_df, fraser_table_cols)
    fraser_table_data = fraser_table.to_dict('records')
    fraser_table_columns = [{'data': col, 'title': col} for col in fraser_table.columns]

    # OUTRIDER table
    outrider_table_data = []
    outrider_table_columns = []
    if outrider_df is not None and len(outrider_df) > 0:
        outrider_table_cols = ['sampleID', 'pValue', 'padjust']
        if 'geneID' in outrider_df.columns:
            outrider_table_cols.insert(0, 'geneID')
        elif 'hgncSymbol' in outrider_df.columns:
            outrider_table_cols.insert(0, 'hgncSymbol')
        if 'zScore' in outrider_df.columns:
            outrider_table_cols.append('zScore')
        if 'l2fc' in outrider_df.columns:
            outrider_table_cols.append('l2fc')
        elif 'log2fc' in outrider_df.columns:
            outrider_table_cols.append('log2fc')
        if 'familyID' in outrider_df.columns:
            outrider_table_cols.append('familyID')

        outrider_table = prepare_table_data(outrider_df, outrider_table_cols)
        outrider_table_data = outrider_table.to_dict('records')
        outrider_table_columns = [{'data': col, 'title': col} for col in outrider_table.columns]

    # Common results table
    common_table_data = []
    common_table_columns = []
    if len(df_common) > 0:
        common_table = prepare_table_data(df_common)
        common_table_data = common_table.to_dict('records')
        common_table_columns = [{'data': col, 'title': col} for col in common_table.columns]

    # Calculate stats
    stats = calculate_stats(fraser_df, outrider_df, pvalue_threshold)
    stats['common_count'] = len(df_common)

    # Get unique families for filtering (combine from both FRASER and OUTRIDER)
    families = set()

    # Safety check for FRASER
    if 'fraser_df' in locals() and fraser_df is not None and 'familyID' in fraser_df.columns:
        # dropna() removes actual nulls; unique() gets the rest
        families.update(fraser_df['familyID'].dropna().unique())

    # Safety check for OUTRIDER
    if 'outrider_df' in locals() and outrider_df is not None and 'familyID' in outrider_df.columns:
        families.update(outrider_df['familyID'].dropna().unique())

    # Filter 'Unknown', ensure everything is a string (to avoid sort errors), and sort
    sorted_families = sorted([str(f) for f in families if str(f) != 'Unknown'])

    # Set up Jinja2 environment
    template_dir = Path(__file__).parent / 'templates'
    env = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
    template = env.get_template('interactive_dashboard.html.j2')

    # Render template
    print('Rendering HTML template...')
    html_content = template.render(
        volcano_data=volcano_data,
        outrider_volcano_data=outrider_volcano_data,
        outrider_top_genes_data=outrider_top_genes_data,
        top_positions=top_positions,
        fraser_tabix_data=fraser_tabix_data,
        outrider_tabix_data=outrider_tabix_data,
        fraser_table_data=fraser_table_data,
        fraser_table_columns=fraser_table_columns,
        outrider_table_data=outrider_table_data,
        outrider_table_columns=outrider_table_columns,
        common_table_data=common_table_data,
        common_table_columns=common_table_columns,
        stats=stats,
        families=sorted_families,
        default_pvalue_threshold=pvalue_threshold,
        default_deltapsi_threshold=deltapsi_threshold,
        default_zscore_threshold=zscore_threshold,
    )

    # Write output
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

    # Complete example with all features
    python create_interactive_dashboard.py --fraser fraser.csv --outrider outrider.csv \\
        --family-mapping families.csv --output dashboard.html

    # With pre-indexed tabix files
    python create_interactive_dashboard.py --fraser fraser.tsv.gz --outrider outrider.tsv.gz --output dashboard.html

    # Custom thresholds
    python create_interactive_dashboard.py --fraser results.csv --output dashboard.html \\
        --pvalue-threshold 0.01 --deltapsi-threshold 0.3
        """,
    )

    parser.add_argument('--fraser', '-f', required=True, help='Path to FRASER results file (CSV or tabix-indexed .gz)')

    parser.add_argument(
        '--outrider',
        '-r',
        required=False,
        default=None,
        help='Path to OUTRIDER results file (CSV or tabix-indexed .gz)',
    )

    parser.add_argument('--output', '-o', required=True, help='Output HTML file path')

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
        '--bam-mapping',
        default=None,
        help='Path to TSV file mapping sample IDs to BAM URLs (columns: sampleID, bam_url, bai_url)',
    )

    parser.add_argument(
        '--prepare-tabix',
        action='store_true',
        help='Convert CSV files to tabix-indexed format (keeps the indexed files)',
    )

    parser.add_argument(
        '--tabix-output-dir', default=None, help='Directory for tabix output files (default: same as input)'
    )

    parser.add_argument(
        '--family-mapping', default=None, help='Path to CPG-to-Family mapping CSV (rdnow-export-project-summary format)'
    )

    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()

    print('=' * 60)
    print('FRASER/OUTRIDER Interactive Dashboard Generator')
    print('=' * 60)

    # Load CPG to Family mapping if provided
    cpg_to_family = {}
    if args.family_mapping:
        cpg_to_family = load_cpg_to_family_mapping(args.family_mapping)

    # Optionally prepare tabix files
    fraser_path = args.fraser
    outrider_path = args.outrider

    if args.prepare_tabix:
        print('\nPreparing tabix-indexed files...')
        if args.fraser.endswith('.csv'):
            fraser_path = prepare_tabix_file(args.fraser, args.tabix_output_dir)
        if args.outrider and args.outrider.endswith('.csv'):
            outrider_path = prepare_tabix_file(args.outrider, args.tabix_output_dir)

    # Load data
    print('\nLoading data...')
    fraser_df = load_fraser_data(args.fraser)  # Use original for full data load
    outrider_df = load_outrider_data(args.outrider) if args.outrider else None

    # If a family mapping was provided, use its SG IDs to filter both DataFrames
    # down to only the samples belonging to this dataset, then annotate with family IDs.
    if cpg_to_family:
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

    # Get top positions info
    if len(fraser_df) > 0:
        top_positions = get_top_positions(fraser_df, n=1)
        if top_positions:
            print(
                f'\nTop FRASER hit: {top_positions[0]["chr"]}:{top_positions[0]["start"]}-{top_positions[0]["end"]} '
                f'(p={top_positions[0]["pValue"]:.2e}, gene={top_positions[0]["gene"]})'
            )

    # Create dashboard
    print('\nCreating interactive dashboard...')
    render_dashboard(
        fraser_df=fraser_df,
        outrider_df=outrider_df,
        output_path=args.output,
        pvalue_threshold=args.pvalue_threshold,
        deltapsi_threshold=args.deltapsi_threshold,
        zscore_threshold=args.zscore_threshold,
    )
    print(' Dashboard created successfully!')
    print(f'\nOutput: {args.output}')
    if args.prepare_tabix:
        print('\nTabix files created:')
        if fraser_path != args.fraser:
            print(f'  - {fraser_path}')
        if outrider_path and outrider_path != args.outrider:
            print(f'  - {outrider_path}')


if __name__ == '__main__':
    main()
