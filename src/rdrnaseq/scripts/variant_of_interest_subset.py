"""Subset a seqr-loader MatrixTable to FRASER significant regions and extract variants of interest.

For each FRASER significant region (±200bp buffer), finds variants where at least one
carrier (het or hom-alt) has a genome SG ID matching the RNA SG ID flagged by FRASER (which has a rna sg id and
is mapped to the same participant via Metamist).
Exports as a BED-like TSV with all variant annotations.
"""

import argparse
import json

import pandas as pd
from loguru import logger

import hail as hl

from cpg_utils import hail_batch
from metamist.graphql import query

from rdrnaseq.scripts.dashboard_utilities import (
    PEDIGREE_QUERY,
    build_genome_to_rna_map,
    build_rna_to_genome_map,
    build_rna_to_metadata_map,
)

BUFFER_BP = 200
REFERENCE_GENOME = 'GRCh38'
MAX_POPMAX_AF = 0.01


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
    genome_to_rna: dict[str, set[str]],
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
    with hl.current_backend().fs.open(tsv_path, 'r') as f:
        variant_df = pd.read_csv(f, sep='\t')
    fraser_df = pd.read_csv(fraser_csv_path)
    original_count = len(variant_df)
    logger.info(f'Verifying {original_count} variant rows against {len(fraser_df)} FRASER regions')

    variant_df['_genome_ids'] = variant_df['matching_genome_sg_ids'].apply(json.loads)
    variant_df['_alleles'] = variant_df['alleles'].apply(json.loads)
    locus_split = variant_df['locus'].str.split(':', n=1)
    variant_df['_chrom'] = locus_split.str[0]
    variant_df['_pos'] = locus_split.str[1].astype(int)

    fraser_lookup = build_fraser_lookup(fraser_df)

    drop_cols = {'matching_genome_sg_ids', 'rna_sg_ids', 'fraser_distance'}
    base_cols = [c for c in variant_df.columns if c not in drop_cols and not c.startswith('_')]

    output_rows = []
    for _, row in variant_df.iterrows():
        chrom = row['_chrom']
        pos = row['_pos']
        ref, alt = row['_alleles'][0], row['_alleles'][1]
        variant_end = pos + max(len(ref), len(alt)) - 1

        for genome_id in row['_genome_ids']:
            rna_ids = genome_to_rna.get(genome_id)
            if rna_ids is None:
                continue
            for rna_id in rna_ids:
                regions = fraser_lookup.get(rna_id, {}).get(chrom, [])
                meta = rna_to_metadata.get(rna_id, {})
                for region_start, region_end, delta_psi, psi_value in regions:
                    dist = compute_span_distance(pos, variant_end, region_start, region_end)
                    if dist <= buffer_bp:
                        sample_row = {col: row[col] for col in base_cols}
                        sample_row['matching_genome_sg_id'] = genome_id
                        sample_row['rna_sg_id'] = rna_id
                        sample_row['fraser_distance'] = dist
                        sample_row['deltaPsi'] = delta_psi
                        sample_row['psiValue'] = psi_value
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
    if len(result_df) > 0:
        export_bed(result_df, output_bed_path)
    else:
        open(output_bed_path, 'w').close()
        logger.info(f'No verified variants; wrote empty BED to {output_bed_path}')


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


def merge_overlapping_intervals(df: pd.DataFrame) -> list[dict]:
    """Merge overlapping buffered intervals, unioning their sequencing group ID sets.

    Takes a DataFrame with columns: seqnames, start, end, genome_sg_ids (set per row).
    Returns a list of dicts with keys: chrom, start, end, genome_sg_ids (set of genome IDs).
    """
    df = df.copy()
    df['buf_start'] = df['start'] - BUFFER_BP
    df['buf_end'] = df['end'] + BUFFER_BP
    df = df.sort_values(['seqnames', 'buf_start']).reset_index(drop=True)

    merged = []
    for chrom, group in df.groupby('seqnames'):
        rows = group.to_dict('records')
        current = {
            'chrom': chrom,
            'start': rows[0]['buf_start'],
            'end': rows[0]['buf_end'],
            'genome_sg_ids': set(rows[0]['genome_ids']),
            'fraser_regions': [(rows[0]['start'], rows[0]['end'])],
        }
        for row in rows[1:]:
            if row['buf_start'] <= current['end']:
                current['end'] = max(current['end'], row['buf_end'])
                current['genome_sg_ids'].update(row['genome_ids'])
                current['fraser_regions'].append((row['start'], row['end']))
            else:
                merged.append(current)
                current = {
                    'chrom': chrom,
                    'start': row['buf_start'],
                    'end': row['buf_end'],
                    'genome_sg_ids': set(row['genome_ids']),
                    'fraser_regions': [(row['start'], row['end'])],
                }
        merged.append(current)

    return merged


def build_hail_intervals(merged_intervals: list[dict]) -> tuple[list[hl.Interval], hl.Table]:
    """Build Hail interval list for filtering and a keyed Table for SG ID annotation."""
    hail_intervals = [
        hl.parse_locus_interval(
            f'[{iv["chrom"]}:{iv["start"]}-{iv["end"]}]',
            reference_genome=REFERENCE_GENOME,
        )
        for iv in merged_intervals
    ]
    fraser_region_schema = hl.tstruct(start=hl.tint32, end=hl.tint32)
    rows = [
        hl.Struct(
            chrom=iv['chrom'],
            start=iv['start'],
            end=iv['end'],
            genome_sg_ids=sorted(iv['genome_sg_ids']),
            fraser_regions=[hl.Struct(start=s, end=e) for s, e in iv['fraser_regions']],
        )
        for iv in merged_intervals
    ]
    ht = hl.Table.parallelize(
        rows,
        schema=hl.tstruct(
            chrom=hl.tstr,
            start=hl.tint32,
            end=hl.tint32,
            genome_sg_ids=hl.tarray(hl.tstr),
            fraser_regions=hl.tarray(fraser_region_schema),
        ),
    )
    ht = ht.annotate(
        interval=hl.locus_interval(ht.chrom, ht.start, ht.end, includes_end=True, reference_genome=REFERENCE_GENOME),
        genome_sg_ids=hl.set(ht.genome_sg_ids),
    )
    ht = ht.select('interval', 'genome_sg_ids', 'fraser_regions')
    return hail_intervals, ht.key_by('interval')


def map_to_genome_ids(rna_id: str, rna_to_genome_ids: dict[str, set[str]]) -> set[str]:
    """Look up genome SG IDs for a given RNA SG ID."""
    genome_ids = rna_to_genome_ids.get(rna_id)
    if genome_ids is None:
        return set()
    return genome_ids


def subset_mt_to_variants_of_interest(
    mt_path: str,
    hail_intervals: list[hl.Interval],
    interval_ht: hl.Table,
    genome_to_rna: dict[str, str],
) -> hl.Table:
    """Subset MT to FRASER regions, filter to carrier-matched rare variants, and annotate BED fields."""
    logger.info(f'Reading MatrixTable: {mt_path}')
    mt = hl.read_matrix_table(mt_path)

    logger.info('Filtering to FRASER significant regions')
    mt = hl.filter_intervals(mt, hail_intervals)

    relevant_genome_ids = hl.literal(set(genome_to_rna.keys()))
    mt = mt.filter_cols(relevant_genome_ids.contains(mt.s))

    mt = mt.annotate_rows(carriers=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.collect_as_set(mt.s)))

    ht = mt.rows()

    logger.info('Annotating variants with CSV region SG IDs')
    interval_annot = interval_ht.index(ht.locus)
    ht = ht.annotate(
        genome_sg_ids=interval_annot.genome_sg_ids,
        fraser_regions=interval_annot.fraser_regions,
    )

    ht = ht.annotate(
        matching_samples=ht.carriers.intersection(ht.genome_sg_ids),
    )
    ht = ht.filter(ht.matching_samples.length() > 0)

    gnomad_af = ht.gnomad_joint.AF_POPMAX_OR_GLOBAL
    ht = ht.filter(hl.is_missing(gnomad_af) | (gnomad_af < MAX_POPMAX_AF))

    pos = ht.locus.position
    ht = ht.annotate(
        fraser_distance=hl.min(
            ht.fraser_regions.map(
                lambda r: hl.if_else(
                    (pos >= r.start) & (pos <= r.end),
                    0,
                    hl.min(hl.abs(pos - r.start), hl.abs(pos - r.end)),
                )
            )
        ),
        bed_chrom=ht.locus.contig,
        bed_start=ht.locus.position - 1,
        bed_end=ht.locus.position - 1 + hl.max(ht.alleles[0].length(), ht.alleles[1].length()),
    )

    genome_to_rna_hl = hl.literal(genome_to_rna)
    ht = ht.annotate(
        rna_sg_ids=ht.matching_samples.map(genome_to_rna_hl.get),
    )

    main_tid = ht.mainTranscript.transcript_id
    matching_tc = ht.vep.transcript_consequences.filter(lambda tc: tc.transcript_id == main_tid)
    mane_raw = hl.or_missing(matching_tc.length() > 0, matching_tc[0].mane_select)
    refseq_tid = hl.or_missing(hl.is_defined(mane_raw) & (mane_raw != ''), mane_raw.split(':')[0])

    fields = {
        'gene_symbol': ht.mainTranscript.gene_symbol,
        'transcript_id': hl.or_else(refseq_tid, main_tid),
        'hgvsc': ht.mainTranscript.hgvsc,
        'major_consequence': ht.mainTranscript.major_consequence,
        'splice_ai_delta_score': ht.splice_ai.delta_score,
        'splice_ai_consequence': ht.splice_ai.splice_consequence,
        'matching_genome_sg_ids': ht.matching_samples,
        'rna_sg_ids': ht.rna_sg_ids,
        'fraser_distance': ht.fraser_distance,
        'bed_chrom': ht.bed_chrom,
        'bed_start': ht.bed_start,
        'bed_end': ht.bed_end,
    }
    if 'avis_phred' in ht.row:
        fields['avis_phred'] = ht.avis_phred
    if 'avis' in ht.row:
        fields['avi_score'] = ht.avis

    return ht.select(**fields)


def read_fraser_csv(csv_path: str, rna_to_genome_ids: dict[str, set[str]]) -> pd.DataFrame:
    """Read FRASER CSV, validate columns, and remap sampleIDs (which are rna_sg_ids) to genome SG IDs."""
    logger.info(f'Reading FRASER CSV: {csv_path}')
    df = pd.read_csv(csv_path)
    required_cols = {'seqnames', 'start', 'end', 'sampleID'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f'CSV missing required columns: {missing}')
    df['genome_ids'] = df['sampleID'].apply(map_to_genome_ids, rna_to_genome_ids=rna_to_genome_ids)
    return df[df['genome_ids'].map(len) > 0].reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mt', required=True, help='GCS path to the seqr-loader MatrixTable')
    parser.add_argument('--csv', required=True, help='Path to FRASER significant results CSV')
    parser.add_argument('--rna_ids', required=True, help='RNA sequencing group IDs', nargs='+')
    parser.add_argument('--query_dataset', required=True, help='Metamist project name')
    parser.add_argument('--output', required=True, help='Output root for coarse TSV')
    parser.add_argument('--output-tsv', required=True, help='Output path for verified TSV')
    parser.add_argument('--output-bed', required=True, help='Output path for verified BED')
    args = parser.parse_args()

    hail_batch.init_batch(driver_memory='highmem', driver_cores=2)
    # Step 0: RNA SG IDs to a set of genome SG IDs from the same participants

    relevant_ids = list(set(args.rna_ids))
    logger.info(f'Input RNA Sequencing Group IDs: {relevant_ids}')
    variables = {'project': args.query_dataset, 'RnaSequencingGroupIds': relevant_ids}

    result = query(PEDIGREE_QUERY, variables=variables)
    rna_to_genome_ids = build_rna_to_genome_map(result)
    genome_to_rna: dict[str, str] = build_genome_to_rna_map(rna_to_genome_ids)
    genome_to_rna_full: dict[str, set[str]] = {}
    for rna_id, genome_ids in rna_to_genome_ids.items():
        for gid in genome_ids:
            genome_to_rna_full.setdefault(gid, set()).add(rna_id)
    rna_to_metadata = build_rna_to_metadata_map(result)

    coarse_tsv_path = f'{args.output}.tsv'

    if hl.current_backend().fs.exists(coarse_tsv_path):
        logger.info(f'Coarse TSV already exists at {coarse_tsv_path}, skipping Hail subset')
    else:
        df = read_fraser_csv(args.csv, rna_to_genome_ids)

        cpg_ids = set().union(*df['genome_ids'])
        logger.info(f'Found {len(cpg_ids)} unique genome sequencing group IDs after mapping from RNA SG IDs')
        logger.info(f'Found {len(df)} FRASER significant regions')

        merged = merge_overlapping_intervals(df)
        logger.info(f'Merged into {len(merged)} non-overlapping intervals (±{BUFFER_BP}bp buffer)')

        hail_intervals, interval_ht = build_hail_intervals(merged)

        ht = subset_mt_to_variants_of_interest(args.mt, hail_intervals, interval_ht, genome_to_rna)

        logger.info(f'Exporting coarse TSV to {coarse_tsv_path}')
        ht.export(coarse_tsv_path, delimiter='\t')

    logger.info('Running verification against original FRASER regions')
    verify_variant_fraser_matches(
        tsv_path=coarse_tsv_path,
        fraser_csv_path=args.csv,
        genome_to_rna=genome_to_rna_full,
        rna_to_metadata=rna_to_metadata,
        output_tsv_path=args.output_tsv,
        output_bed_path=args.output_bed,
    )

    logger.info('Done.')


if __name__ == '__main__':
    main()
