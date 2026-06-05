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
)


BUFFER_BP = 200


def compute_span_distance(variant_start: int, variant_end: int, region_start: int, region_end: int) -> int:
    "Compute distance between a variant span and a region span, treating any overlap as distance 0."
    if variant_start <= region_end and variant_end >= region_start:
        return 0
    return min(abs(variant_start - region_end), abs(variant_end - region_start))


def build_fraser_lookup(fraser_df: pd.DataFrame) -> dict[str, dict[str, list[tuple[int, int, float]]]]:
    """Index FRASER regions as sample -> chrom -> list of (start, end, deltaPsi)."""
    lookup: dict[str, dict[str, list[tuple[int, int, float]]]] = {}
    for _, row in fraser_df.iterrows():
        sample_regions = lookup.setdefault(row['sampleID'], {})
        delta_psi = float(row['deltaPsi'])
        sample_regions.setdefault(row['seqnames'], []).append((int(row['start']), int(row['end']), delta_psi))
    return lookup


def verify_variant_fraser_matches(
    tsv_path: str,
    fraser_csv_path: str,
    genome_to_rna: dict[str, str],
    output_tsv_path: str,
    output_bed_path: str,
    buffer_bp: int = BUFFER_BP,
) -> None:
    """Re-check each variant-sample match against the original FRASER regions.
    The main point is to Remove false positives, recalculate distance, and attach deltaPsi."""
    df = pd.read_csv(tsv_path, sep='\t')
    fraser_df = pd.read_csv(fraser_csv_path)
    original_count = len(df)
    logger.info(f'Verifying {original_count} variant rows against {len(fraser_df)} FRASER regions')
    # Local testing revealed weird hail savin
    # Parse Hail JSON-encoded array columns
    df['_genome_ids'] = df['matching_genome_sg_ids'].apply(json.loads)
    df['_rna_ids'] = df['rna_sg_ids'].apply(json.loads)
    df['_alleles'] = df['alleles'].apply(json.loads)

    # Parse locus -> chrom, pos
    locus_split = df['locus'].str.split(':', n=1)
    df['_chrom'] = locus_split.str[0]
    df['_pos'] = locus_split.str[1].astype(int)

    fraser_lookup = build_fraser_lookup(fraser_df)

    verified_genome_ids_col = []
    verified_rna_ids_col = []
    verified_distance_col = []
    verified_delta_psi_col = []
    keep_mask = []

    for _, row in df.iterrows():
        chrom = row['_chrom']
        pos = row['_pos']
        alleles = row['_alleles']
        ref, alt = alleles[0], alleles[1]
        variant_end = pos + max(len(ref), len(alt)) - 1

        verified_genome_ids = []
        verified_rna_ids = []
        min_distance = None
        closest_delta_psi = None

        for genome_id in row['_genome_ids']:
            rna_id = genome_to_rna.get(genome_id)
            if rna_id is None:
                continue

            regions = fraser_lookup.get(rna_id, {}).get(chrom, [])
            sample_verified = False
            for region_start, region_end, delta_psi in regions:
                dist = compute_span_distance(pos, variant_end, region_start, region_end)
                if dist <= buffer_bp:
                    sample_verified = True
                    if min_distance is None or dist < min_distance:
                        min_distance = dist
                        closest_delta_psi = delta_psi

            if sample_verified:
                verified_genome_ids.append(genome_id)
                verified_rna_ids.append(rna_id)

        if verified_genome_ids:
            keep_mask.append(True)
            verified_genome_ids_col.append(json.dumps(verified_genome_ids))
            verified_rna_ids_col.append(json.dumps(verified_rna_ids))
            verified_distance_col.append(min_distance)
            verified_delta_psi_col.append(closest_delta_psi)
        else:
            keep_mask.append(False)
            verified_genome_ids_col.append(None)
            verified_rna_ids_col.append(None)
            verified_distance_col.append(None)
            verified_delta_psi_col.append(None)

    df['matching_genome_sg_ids'] = verified_genome_ids_col
    df['rna_sg_ids'] = verified_rna_ids_col
    df['fraser_distance'] = verified_distance_col
    df['deltaPsi'] = verified_delta_psi_col

    df = df[keep_mask].drop(columns=[c for c in df.columns if c.startswith('_')])
    verified_count = len(df)
    removed = original_count - verified_count
    logger.info(
        f'Verification complete: {verified_count}/{original_count} variants retained, '
        f'{removed} false positives removed'
    )

    df.to_csv(output_tsv_path, sep='\t', index=False)

    # Export BED from verified results
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
    bed.to_csv(output_bed_path, sep='\t', header=False, index=False)
    logger.info(f'Exported verified BED ({len(bed)} rows) to {output_bed_path}')


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
    result = query(PEDIGREE_QUERY, variables={
        'project': args.query_dataset,
        'RnaSequencingGroupIds': list(set(args.rna_ids)),
    })
    rna_to_genome_ids = build_rna_to_genome_map(result)
    genome_to_rna = build_genome_to_rna_map(rna_to_genome_ids)

    verify_variant_fraser_matches(
        tsv_path=args.tsv,
        fraser_csv_path=args.fraser_csv,
        genome_to_rna=genome_to_rna,
        output_tsv_path=args.output_tsv,
        output_bed_path=args.output_bed,
        buffer_bp=args.buffer_bp,
    )


if __name__ == '__main__':
    main()