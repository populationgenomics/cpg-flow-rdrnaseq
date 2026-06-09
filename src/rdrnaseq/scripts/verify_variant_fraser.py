"""Verify variant-FRASER sample matches and correct false positives from interval merging.

Reads the TSV produced by variant_of_interest_subset, cross-references each
variant-sample match against the original FRASER CSV, removes false positives,
recalculates fraser_distance (indel-span-aware), and attaches deltaPsi.
"""

import argparse
import json

import pandas as pd
from loguru import logger

from metamist.graphql import query

from rdrnaseq.scripts.dashboard_utilities import (
    PEDIGREE_QUERY,
    build_genome_to_rna_map,
    build_rna_to_genome_map,
    build_rna_to_metadata_map,
)

BUFFER_BP = 200


def compute_span_distance(variant_start: int, variant_end: int, region_start: int, region_end: int) -> int:
    "Compute distance between a variant span and a region span, treating any overlap as distance 0."
    if variant_start <= region_end and variant_end >= region_start:
        return 0
    return min(abs(variant_start - region_end), abs(variant_end - region_start))


def build_fraser_lookup(fraser_df: pd.DataFrame) -> dict[str, dict[str, list[tuple[int, int, float, float]]]]:
    """Index FRASER regions as sample -> chrom -> list of (start, end, deltaPsi, psiValue)."""
    lookup: dict[str, dict[str, list[tuple[int, int, float, float]]]] = {}
    for _, row in fraser_df.iterrows():
        sample_regions = lookup.setdefault(row['sampleID'], {})
        delta_psi = float(row['deltaPsi'])
        psi_value = float(row['psiValue'])
        sample_regions.setdefault(row['seqnames'], []).append(
            (int(row['start']), int(row['end']), delta_psi, psi_value)
        )
    return lookup


def verify_variant_fraser_matches(
    tsv_path: str,
    fraser_csv_path: str,
    genome_to_rna: dict[str, str],
    rna_to_metadata: dict[str, dict],
    output_tsv_path: str,
    output_bed_path: str,
    buffer_bp: int = BUFFER_BP,
) -> None:
    """Re-check each variant-sample match against the original FRASER regions.

    Removes false positives from interval merging, recalculates per-sample
    fraser_distance, attaches deltaPsi, and denormalizes to one row per
    verified variant-sample pair with correct participant metadata.
    """
    variant_df = pd.read_csv(tsv_path, sep='\t')
    fraser_df = pd.read_csv(fraser_csv_path)
    original_count = len(variant_df)
    logger.info(f'Verifying {original_count} variant rows against {len(fraser_df)} FRASER regions')

    variant_df['_genome_ids'] = variant_df['matching_genome_sg_ids'].apply(json.loads)
    variant_df['_alleles'] = variant_df['alleles'].apply(json.loads)
    locus_split = variant_df['locus'].str.split(':', n=1)
    variant_df['_chrom'] = locus_split.str[0]
    variant_df['_pos'] = locus_split.str[1].astype(int)

    fraser_lookup = build_fraser_lookup(fraser_df)

    # Columns from the Hail TSV to carry through (drop the old array ID and metadata columns)
    drop_cols = {'matching_genome_sg_ids', 'rna_sg_ids', 'fraser_distance'}
    base_cols = [c for c in variant_df.columns if c not in drop_cols and not c.startswith('_')]

    output_rows = []
    for _, row in variant_df.iterrows():
        chrom = row['_chrom']
        pos = row['_pos']
        ref, alt = row['_alleles'][0], row['_alleles'][1]
        variant_end = pos + max(len(ref), len(alt)) - 1

        for genome_id in row['_genome_ids']:
            rna_id = genome_to_rna.get(genome_id)
            if rna_id is None:
                continue

            regions = fraser_lookup.get(rna_id, {}).get(chrom, [])
            best_dist = None
            best_delta_psi = None
            best_psi_value = None
            for region_start, region_end, delta_psi, psi_value in regions:
                dist = compute_span_distance(pos, variant_end, region_start, region_end)
                if dist <= buffer_bp and (best_dist is None or dist < best_dist):
                    best_dist = dist
                    best_delta_psi = delta_psi
                    best_psi_value = psi_value

            if best_dist is not None:
                sample_row = {col: row[col] for col in base_cols}
                sample_row['matching_genome_sg_id'] = genome_id
                sample_row['rna_sg_id'] = rna_id
                sample_row['fraser_distance'] = best_dist
                sample_row['deltaPsi'] = best_delta_psi
                sample_row['psiValue'] = best_psi_value
                meta = rna_to_metadata.get(rna_id, {})
                sample_row['participant_external_id'] = meta.get('participant_external_id', '')
                sample_row['family_id'] = meta.get('family_id', '')
                sample_row['affected'] = meta.get('affected', 'Unknown')
                output_rows.append(sample_row)

    result_df = pd.DataFrame(output_rows)
    verified_variants = result_df['locus'].nunique() if len(result_df) > 0 else 0
    logger.info(
        f'Verification complete: {len(result_df)} variant-sample rows '
        f'({verified_variants} unique variants) from {original_count} input variants'
    )

    result_df.to_csv(output_tsv_path, sep='\t', index=False)
    export_bed(result_df, output_bed_path)


def export_bed(df: pd.DataFrame, bed_path: str) -> None:
    """Write a minimal IGV-compatible BED from the verified variant DataFrame."""
    bed = df[['bed_chrom', 'bed_start', 'bed_end']].copy()
    names = []
    for _, r in df.iterrows():
        tid = r.get('transcript_id', '')
        symbol = r.get('gene_symbol', '')
        hgvsc = r.get('hgvsc', '')
        if pd.notna(tid) and tid:
            suffix = hgvsc.split(':')[1] if pd.notna(hgvsc) and ':' in str(hgvsc) else hgvsc
            names.append(f'{tid}({symbol}):{suffix}')
        else:
            alleles = json.loads(r['alleles'])
            names.append(f'{r["bed_chrom"]}:{r["bed_start"]}:{alleles[0]}>{alleles[1]}')
    bed['name'] = names
    bed.to_csv(bed_path, sep='\t', header=False, index=False)
    logger.info(f'Exported verified BED ({len(bed)} rows) to {bed_path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tsv', required=True, help='Input TSV from variant_of_interest_subset')
    parser.add_argument('--fraser-csv', required=True, help='FRASER significant results CSV')
    parser.add_argument('--rna-ids', required=True, nargs='+', help='RNA sequencing group IDs')
    parser.add_argument('--query-dataset', required=True, help='Metamist project name')
    parser.add_argument('--output-tsv', required=True, help='Output path for verified TSV')
    parser.add_argument('--output-bed', required=True, help='Output path for verified BED')
    parser.add_argument('--buffer-bp', type=int, default=BUFFER_BP, help='Buffer in bp (default 200)')
    args = parser.parse_args()

    # Query Metamist for genome-to-RNA mapping
    result = query(
        PEDIGREE_QUERY,
        variables={
            'project': args.query_dataset,
            'RnaSequencingGroupIds': list(set(args.rna_ids)),
        },
    )
    rna_to_genome_ids = build_rna_to_genome_map(result)
    genome_to_rna = build_genome_to_rna_map(rna_to_genome_ids)
    rna_to_metadata = build_rna_to_metadata_map(result)

    verify_variant_fraser_matches(
        tsv_path=args.tsv,
        fraser_csv_path=args.fraser_csv,
        genome_to_rna=genome_to_rna,
        rna_to_metadata=rna_to_metadata,
        output_tsv_path=args.output_tsv,
        output_bed_path=args.output_bed,
        buffer_bp=args.buffer_bp,
    )


if __name__ == '__main__':
    main()
