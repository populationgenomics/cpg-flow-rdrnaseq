"""Utility functions for the interactive dashboard."""

import os

import numpy as np
import pandas as pd
import pysam
from loguru import logger

from metamist.apis import WebApi
from metamist.graphql import gql

# =============================================================================
# Constants and Configuration
# =============================================================================

FRASER_REQUIRED_COLUMNS = ['seqnames', 'start', 'end', 'pValue', 'padjust', 'deltaPsi', 'psiValue', 'type', 'sampleID']
FRASER_OPTIONAL_COLUMNS = ['hgncSymbol', 'counts', 'totalCounts', 'nonsplitCounts']

OUTRIDER_REQUIRED_COLUMNS = ['pValue', 'padjust', 'sampleID']
OUTRIDER_OPTIONAL_COLUMNS = ['zScore', 'l2fc', 'log2fc', 'geneID', 'hgncSymbol', 'seqnames', 'start', 'end']

# Chromosome sort order
CHR_ORDER = {
    **{str(i): i for i in range(1, 23)},
    **{f'chr{i}': i for i in range(1, 23)},
    'X': 23,
    'chrX': 23,
    'Y': 24,
    'chrY': 24,
    'M': 25,
    'MT': 25,
    'chrM': 25,
    'chrMT': 25,
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


def load_ensg_to_symbol(filepath: str) -> dict[str, str]:
    """
    Load the ENSG-to-HGNC-symbol mapping from a two-column TSV (no header).

    Returns a dict mapping ENSG IDs to HGNC symbols.  Each entry is stored
    twice — with and without the version suffix — so lookups work regardless
    of whether the caller's gene IDs include versions.

    Example return: {"ENSG00000001461.17": "NIPAL3", "ENSG00000001461": "NIPAL3", ...}
    """
    mapping: dict[str, str] = {}
    with open(filepath) as fh:
        for raw_line in fh:
            stripped = raw_line.strip()
            if not stripped:
                continue
            parts = stripped.split('\t')
            if len(parts) < 2:
                continue
            ensg_id, symbol = parts[0], parts[1]
            mapping[ensg_id] = symbol
            base_id = ensg_id.split('.')[0]
            if base_id != ensg_id:
                mapping[base_id] = symbol
    print(f'  Loaded {len(mapping)} ENSG-to-symbol entries from {filepath}')
    return mapping


def build_ensg_to_hgnc_subset(df_outrider: pd.DataFrame, ensg_to_symbol: dict[str, str]) -> dict[str, str]:
    """
    Build a small ENSG→HGNC mapping limited to genes present in the OUTRIDER data.

    This subset is embedded in the HTML template as a JS object so the browser
    can match OUTRIDER (ENSG) genes with FRASER (HGNC) genes for common
    significant results.
    """
    if 'geneID' not in df_outrider.columns:
        return {}
    unique_gene_ids = df_outrider['geneID'].dropna().unique()
    subset: dict[str, str] = {}
    for gid in unique_gene_ids:
        gid_str = str(gid)
        if gid_str in ensg_to_symbol:
            subset[gid_str] = ensg_to_symbol[gid_str]
    print(f'  Built ENSG-to-HGNC subset: {len(subset)} of {len(unique_gene_ids)} OUTRIDER genes mapped')
    return subset


def enrich_with_gene_mapping(
    fraser_df: pd.DataFrame,
    outrider_df: pd.DataFrame | None,
    ensg_to_symbol: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """
    Add both gene name and gene ID columns to FRASER and OUTRIDER DataFrames.

    - OUTRIDER: adds ``hgncSymbol`` by looking up ``geneID`` in the mapping.
    - FRASER: adds ``geneID`` using the reverse (HGNC→ENSG) mapping.
    """
    # Build reverse mapping: HGNC symbol -> ENSG ID (prefer versioned IDs)
    symbol_to_ensg: dict[str, str] = {}
    for ensg_id, symbol in ensg_to_symbol.items():
        # Only use the versioned (longer) ID so reverse lookups are unambiguous
        if '.' in ensg_id:
            symbol_to_ensg[symbol] = ensg_id

    # OUTRIDER: geneID -> hgncSymbol
    if outrider_df is not None and 'geneID' in outrider_df.columns:
        outrider_df['hgncSymbol'] = outrider_df['geneID'].map(ensg_to_symbol).fillna('')
        mapped = (outrider_df['hgncSymbol'] != '').sum()
        print(f'  OUTRIDER: mapped {mapped}/{len(outrider_df)} rows to HGNC symbols')

    # FRASER: hgncSymbol -> geneID
    if 'hgncSymbol' in fraser_df.columns:
        fraser_df['geneID'] = fraser_df['hgncSymbol'].map(symbol_to_ensg).fillna('')
        mapped = (fraser_df['geneID'] != '').sum()
        print(f'  FRASER: mapped {mapped}/{len(fraser_df)} rows to ENSG gene IDs')
    else:
        fraser_df['geneID'] = ''

    return fraser_df, outrider_df


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
            f'{name} data is missing required columns: {missing}\nRequired: {required_cols}\nFound: {list(df.columns)}'
        )

    # Check for optional columns that are present
    present_optional = set(optional_cols) & set(df.columns)
    if present_optional:
        print(f'  Found optional columns: {present_optional}')


# =============================================================================
# Metamist Pedigree / SG-ID Mapping
# =============================================================================

PEDIGREE_QUERY = gql(
    """
    query Pedigree($project: String!, $RnaSequencingGroupIds: [String!]!) {
      project(name: $project) {
        sequencingGroups(id: {in_: $RnaSequencingGroupIds}) {
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
              samples {
                sequencingGroups(type: {eq: "genome"}, technology: {eq: "short-read"}, activeOnly: {eq: true}) {
                  id
                  technology
                  type
                }
              }
            }
          }
        }
      }
    }
    """
)


def build_rna_to_genome_map(query_result: dict) -> dict[str, set[str]]:
    """Parse metamist query result into a mapping of RNA SG IDs to genome SG IDs."""
    rna_to_genome_ids: dict[str, set[str]] = {}
    for group in query_result['project']['sequencingGroups']:
        rna_id = group['id']
        participant = group['sample']['participant']
        genome_ids = set()
        for sample in participant['samples']:
            for sg in sample['sequencingGroups']:
                genome_ids.add(sg['id'])
        rna_to_genome_ids[rna_id] = genome_ids
    return rna_to_genome_ids


def build_genome_to_rna_map(rna_to_genome_ids: dict[str, set[str]]) -> dict[str, str]:
    """Reverse the RNA-to-genome mapping: genome SG ID -> RNA SG ID."""
    genome_to_rna: dict[str, str] = {}
    for rna_id, genome_ids in rna_to_genome_ids.items():
        for gid in genome_ids:
            if gid in genome_to_rna:
                raise ValueError(f'Genome SG: {gid} linked to multiple RNA SGs: {genome_to_rna[gid]}, {rna_id}')
            genome_to_rna[gid] = rna_id
    return genome_to_rna


def build_rna_to_metadata_map(query_result: dict) -> dict[str, dict[str, str | int]]:
    """Parse metamist query result into a mapping of RNA SG ID to participant metadata."""
    rna_to_metadata: dict[str, dict[str, str | int]] = {}
    for group in query_result['project']['sequencingGroups']:
        rna_id = group['id']
        participant = group['sample']['participant']
        try:
            rna_to_metadata[rna_id] = {
                'participant_external_id': participant['externalId'],
                'family_id': participant['families'][0]['externalId'],
                'affected': AFFECTED_LABELS.get(participant['familyParticipants'][0]['affected'], 'Unknown'),
            }
        except (KeyError, IndexError, TypeError):
            logger.warning(f'Incomplete metadata for RNA SG {rna_id}, skipping')
    return rna_to_metadata


# =============================================================================
# CPG to Family ID Mapping
# =============================================================================
def load_cpg_to_family_mapping(mapping_file: str) -> dict:
    """Load CPG ID to Family ID mapping from a CSV."""
    print(f'Loading CPG to Family mapping from {mapping_file}...')
    col_renames = {'family.external_ids': 'familyID', 'participant.external_id': 'participantExternalId'}
    meta_df = pd.read_csv(mapping_file).rename(columns=col_renames).set_index('sequencing_group.id')
    for col in ['familyID', 'participantExternalId', 'affected']:
        if col not in meta_df.columns:
            meta_df[col] = ''
    mapping = meta_df[['familyID', 'participantExternalId', 'affected']].fillna('').to_dict('index')
    unique_families = len({m['familyID'] for m in mapping.values()})
    print(f'  Loaded {len(mapping)} CPG IDs mapping to {unique_families} families')
    return mapping


AFFECTED_LABELS: dict[int, str] = {0: 'Unknown', 1: 'Unaffected', 2: 'Affected'}


def affected_status_label(value) -> str:
    """Map numeric affected status to a human-readable label."""
    try:
        numeric = int(value)
    except (ValueError, TypeError):
        return 'Unknown'
    return AFFECTED_LABELS.get(numeric, 'Unknown')


def add_family_ids(df: pd.DataFrame, cpg_to_metadata: dict) -> pd.DataFrame:
    """Add familyID, participantExternalId, and affected columns by mapping sampleID."""
    meta_df = pd.DataFrame.from_dict(cpg_to_metadata, orient='index')
    joined = df.join(meta_df, on='sampleID', how='left')
    joined['familyID'] = joined['familyID'].fillna('Unknown')
    joined = joined.fillna({'participantExternalId': '', 'affected': ''})
    joined['affected'] = joined['affected'].apply(affected_status_label)
    return joined


# =============================================================================
# Data Loading Functions
# =============================================================================


def load_fraser_data(filepath: str) -> pd.DataFrame:
    """Load and prepare FRASER results data."""
    print(f'Loading FRASER data from {filepath}...')

    df = pd.read_csv(filepath, sep='\t', compression='gzip') if filepath.endswith('.gz') else pd.read_csv(filepath)

    validate_dataframe(df, FRASER_REQUIRED_COLUMNS, FRASER_OPTIONAL_COLUMNS, 'FRASER')

    # Add hgncSymbol if missing
    if 'hgncSymbol' not in df.columns:
        df['hgncSymbol'] = 'NA'

    # Calculate -log10(pValue) for plotting
    df['-log10(pValue)'] = -np.log10(df['pValue'].replace(0, np.finfo(float).tiny))

    print(f'  Loaded {len(df)} FRASER rows')
    return df


def load_outrider_data(filepath: str) -> pd.DataFrame | None:
    """Load and prepare OUTRIDER results data."""
    if filepath is None:
        return None

    print(f'Loading OUTRIDER data from {filepath}...')

    df = pd.read_csv(filepath, sep='\t', compression='gzip') if filepath.endswith('.gz') else pd.read_csv(filepath)

    validate_dataframe(df, OUTRIDER_REQUIRED_COLUMNS, OUTRIDER_OPTIONAL_COLUMNS, 'OUTRIDER')

    # Ensure hgncSymbol column exists (will be properly populated by
    # enrich_with_gene_mapping() when the ENSG-to-symbol mapping is available)
    if 'hgncSymbol' not in df.columns:
        df['hgncSymbol'] = ''

    print(f'  Loaded {len(df)} OUTRIDER rows')
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


def csv_to_tabix(csv_path: str, output_dir: str | None = None) -> str:
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
    print(f'  Converting {csv_path} to tabix format...')

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
    tsv_path = os.path.join(output_dir, f'{base_name}.sorted.tsv')
    gz_path = tsv_path + '.gz'

    # Ensure seqnames, start, end are first columns for tabix
    cols = ['seqnames', 'start', 'end'] + [c for c in df.columns if c not in ['seqnames', 'start', 'end']]
    df = df[cols]

    # Write TSV with header
    df.to_csv(tsv_path, sep='\t', index=False)

    # Bgzip compress
    print('  Compressing with bgzip...')
    pysam.tabix_compress(tsv_path, gz_path, force=True)
    os.remove(tsv_path)

    # Create tabix index
    print('  Creating tabix index...')
    pysam.tabix_index(gz_path, seq_col=0, start_col=1, end_col=2, meta_char='#', line_skip=1, force=True)

    print(f'  Created: {gz_path} and {gz_path}.tbi')
    return gz_path


def prepare_tabix_file(filepath: str, output_dir: str | None = None) -> str:
    """
    Prepare a file for tabix queries. Converts CSV to tabix if needed.

    Args:
        filepath: Path to input file (CSV or .gz)
        output_dir: Directory for output files if conversion needed

    Returns:
        Path to tabix-ready .gz file
    """
    if is_tabix_indexed(filepath):
        print(f'  File already tabix-indexed: {filepath}')
        return filepath

    if filepath.endswith('.csv'):
        return csv_to_tabix(filepath, output_dir)

    if filepath.endswith('.gz'):
        # .gz but no .tbi - need to index it
        print(f'  Creating tabix index for {filepath}...')
        pysam.tabix_index(filepath, seq_col=0, start_col=1, end_col=2, meta_char='#', line_skip=1, force=True)

        return filepath

    raise ValueError(f'Unsupported file format: {filepath}. Expected .csv or .gz')


# =============================================================================
# Link Preparation Functions
# =============================================================================


SEQR_VARIANT_TEMPLATE = 'https://seqr.populationgenomics.org.au/variant_search/variant/{variant}/family/{sample}'


class SeqrVariantLinkEngine:
    """Build seqr variant-search URLs from a family ID and variant string."""

    def __init__(self, project: str, variant_template: str = SEQR_VARIANT_TEMPLATE):
        """querying metamist for external family id mapping to seqr sample mapping
        the configurabile things are
        - project to query
        - sequencing type which will always be genome here, but is changeable in case we want to use
        this for exome data in the future"""
        api = WebApi()
        self.guid_map: dict[str, str] = api.get_seqr_family_guid_map('genome', project)
        self.variant_template = variant_template

    def build_link(self, family_id: str, variant_string: str) -> str | None:
        """Return a seqr URL, or None if the family ID has no GUID mapping."""
        guid = self.guid_map.get(family_id)
        if not guid:
            return None
        return self.variant_template.format(variant=variant_string, sample=guid)
