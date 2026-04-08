"""Utility functions for the interactive dashboard."""

import os
from typing import Optional

import numpy as np
import pandas as pd
import pysam

# =============================================================================
# Constants and Configuration
# =============================================================================

FRASER_REQUIRED_COLUMNS = ['seqnames', 'start', 'end', 'pValue', 'padjust',
                           'deltaPsi', 'psiValue', 'type', 'sampleID']
FRASER_OPTIONAL_COLUMNS = ['hgncSymbol', 'counts', 'totalCounts', 'nonsplitCounts']

OUTRIDER_REQUIRED_COLUMNS = ['pValue', 'padjust', 'sampleID']
OUTRIDER_OPTIONAL_COLUMNS = ['zScore', 'l2fc', 'log2fc', 'geneID', 'hgncSymbol', 'seqnames', 'start', 'end']

# Chromosome sort order
CHR_ORDER = {
    **{str(i): i for i in range(1, 23)},
    **{f'chr{i}': i for i in range(1, 23)},
    'X': 23, 'chrX': 23,
    'Y': 24, 'chrY': 24,
    'M': 25, 'MT': 25, 'chrM': 25, 'chrMT': 25
}


# =============================================================================
# Column Helper
# =============================================================================

def _col(df: pd.DataFrame, col: str, default=None, fill=None) -> list:
    """Extract a column as a list, with fallback if the column is missing."""
    if col not in df.columns:
        return [default] * len(df)
    s = df[col]
    if fill is not None:
        s = s.fillna(fill)
    return s.tolist()


# =============================================================================
# Gene Matching
# =============================================================================

def _genes_match(fraser_row: pd.Series, outrider_row: pd.Series) -> bool:
    """Check if a FRASER row's HGNC symbol matches an OUTRIDER row's gene ID."""
    fraser_gene = str(fraser_row.get('hgncSymbol', '')).upper() if pd.notna(fraser_row.get('hgncSymbol')) else ''
    outrider_gene = str(outrider_row.get('geneID', outrider_row.get('hgncSymbol', ''))).upper()
    if not fraser_gene:
        return False
    return fraser_gene in outrider_gene or any(fraser_gene in part for part in outrider_gene.split(';'))


# =============================================================================
# Chromosome Utilities
# =============================================================================

def normalize_chromosome(chrom: str) -> str:
    """Normalize chromosome name to include 'chr' prefix."""
    chrom = str(chrom)
    if not chrom.startswith('chr'):
        return f'chr{chrom}'
    return chrom


def chr_sort_key(chr_name: str) -> int:
    """Get sort key for chromosome ordering."""
    return CHR_ORDER.get(chr_name, CHR_ORDER.get(str(chr_name).replace('chr', ''), 26))


# =============================================================================
# Data Validation
# =============================================================================

def validate_dataframe(df: pd.DataFrame, required_cols: list, optional_cols: list, name: str) -> None:
    """Validate that a DataFrame has required columns."""
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(
            f"{name} data is missing required columns: {missing}\n"
            f"Required: {required_cols}\n"
            f"Found: {list(df.columns)}"
        )

    # Check for optional columns that are present
    present_optional = set(optional_cols) & set(df.columns)
    if present_optional:
        print(f"  Found optional columns: {present_optional}")


# =============================================================================
# CPG to Family ID Mapping
# =============================================================================

def load_cpg_to_family_mapping(mapping_file: str) -> dict:
    """
    Load CPG ID to Family ID mapping from project summary CSV.

    Args:
        mapping_file: Path to rdnow-export-project-summary CSV file

    Returns:
        Dictionary mapping CPG IDs (sequencing_group.id) to family.external_ids
    """
    print(f"Loading CPG to Family mapping from {mapping_file}...")

    df = pd.read_csv(mapping_file)

    # Create mapping from sequencing_group.id to family.external_ids
    mapping = {}
    for _, row in df.iterrows():
        cpg_id = row['sequencing_group.id']
        family_id = row['family.external_ids']
        mapping[cpg_id] = family_id

    unique_families = len(set(mapping.values()))
    print(f"  Loaded {len(mapping)} CPG IDs mapping to {unique_families} families")

    return mapping


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_fraser_data(filepath: str) -> pd.DataFrame:
    """Load and prepare FRASER results data."""
    print(f"Loading FRASER data from {filepath}...")

    if filepath.endswith('.gz'):
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(filepath)

    validate_dataframe(df, FRASER_REQUIRED_COLUMNS, FRASER_OPTIONAL_COLUMNS, "FRASER")

    # Add hgncSymbol if missing
    if 'hgncSymbol' not in df.columns:
        df['hgncSymbol'] = 'NA'

    # Calculate -log10(pValue) for plotting
    df['-log10(pValue)'] = -np.log10(df['pValue'].replace(0, np.finfo(float).tiny))

    print(f"  Loaded {len(df)} FRASER rows")
    return df


def load_outrider_data(filepath: str) -> Optional[pd.DataFrame]:
    """Load and prepare OUTRIDER results data."""
    if filepath is None:
        return None

    print(f"Loading OUTRIDER data from {filepath}...")

    if filepath.endswith('.gz'):
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(filepath)

    validate_dataframe(df, OUTRIDER_REQUIRED_COLUMNS, OUTRIDER_OPTIONAL_COLUMNS, "OUTRIDER")

    # Standardize gene column name
    if 'hgncSymbol' not in df.columns and 'geneID' in df.columns:
        df['hgncSymbol'] = df['geneID']
    elif 'hgncSymbol' not in df.columns:
        df['hgncSymbol'] = 'NA'

    print(f"  Loaded {len(df)} OUTRIDER rows")
    return df


# =============================================================================
# Tabix Preparation Functions
# =============================================================================

def is_tabix_indexed(filepath: str) -> bool:
    """Check if a file is already tabix-indexed (.gz with .tbi)."""
    if not filepath.endswith('.gz'):
        return False
    tbi_path = filepath + '.tbi'
    return os.path.exists(tbi_path)


def csv_to_tabix(csv_path: str, output_dir: Optional[str] = None) -> str:
    """
    Convert a CSV file to a sorted, bgzipped, tabix-indexed TSV.

    Args:
        csv_path: Path to input CSV file
        output_dir: Directory for output files (default: same as input)

    Returns:
        Path to the bgzipped file (.gz)

    Coordinate system: The output uses 1-based coordinates with columns:
        - Column 1: seqnames (chromosome)
        - Column 2: start (1-based)
        - Column 3: end (1-based, inclusive)
    """
    print(f"  Converting {csv_path} to tabix format...")

    # Read CSV
    df = pd.read_csv(csv_path)

    # Ensure numeric coordinates
    df['start'] = pd.to_numeric(df['start'], errors='coerce').astype('Int64')
    df['end'] = pd.to_numeric(df['end'], errors='coerce').astype('Int64')

    # Sort by chromosome and start position
    df['_chr_sort'] = df['seqnames'].apply(chr_sort_key)
    df = df.sort_values(['_chr_sort', 'start']).drop(columns=['_chr_sort'])

    # Determine output path
    if output_dir is None:
        output_dir = os.path.dirname(csv_path) or '.'

    base_name = os.path.splitext(os.path.basename(csv_path))[0]
    tsv_path = os.path.join(output_dir, f"{base_name}.sorted.tsv")
    gz_path = tsv_path + '.gz'

    # Ensure seqnames, start, end are first columns for tabix
    cols = ['seqnames', 'start', 'end'] + [c for c in df.columns if c not in ['seqnames', 'start', 'end']]
    df = df[cols]

    # Write TSV with header
    df.to_csv(tsv_path, sep='\t', index=False)

    # Bgzip compress
    print(f"  Compressing with bgzip...")
    pysam.tabix_compress(tsv_path, gz_path, force=True)
    os.remove(tsv_path)

    # Create tabix index
    print(f"  Creating tabix index...")
    pysam.tabix_index(gz_path, seq_col=0, start_col=1, end_col=2, meta_char='#', line_skip=1, force=True)

    print(f"  Created: {gz_path} and {gz_path}.tbi")
    return gz_path


def prepare_tabix_file(filepath: str, output_dir: Optional[str] = None) -> str:
    """
    Prepare a file for tabix queries. Converts CSV to tabix if needed.

    Args:
        filepath: Path to input file (CSV or .gz)
        output_dir: Directory for output files if conversion needed

    Returns:
        Path to tabix-ready .gz file
    """
    if is_tabix_indexed(filepath):
        print(f"  File already tabix-indexed: {filepath}")
        return filepath

    if filepath.endswith('.csv'):
        return csv_to_tabix(filepath, output_dir)

    if filepath.endswith('.gz'):
        # .gz but no .tbi - need to index it
        print(f"  Creating tabix index for {filepath}...")
        try:
            pysam.tabix_index(filepath, seq_col=0, start_col=1, end_col=2, meta_char='#', line_skip=1, force=True)
        except Exception:
            raise RuntimeError(f"Failed to index {filepath}. Ensure it's a sorted, bgzipped TSV.")
        return filepath

    raise ValueError(f"Unsupported file format: {filepath}. Expected .csv or .gz")