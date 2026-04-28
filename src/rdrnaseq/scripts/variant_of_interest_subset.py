"""Subset a seqr-loader MatrixTable to FRASER significant regions and extract variants of interest.

For each FRASER significant region (±200bp buffer), finds variants where at least one
carrier sample (het or hom-alt) matches the FRASER sampleID for that region.
Exports as a BED-like TSV with all variant annotations.
"""

import argparse

import pandas as pd
from cpg_utils import hail_batch
from loguru import logger
from metamist.graphql import gql, query

import hail as hl

BUFFER_BP = 200
REFERENCE_GENOME = 'GRCh38'
MAX_POPMAX_AF = 0.01

query_ids = gql(
    """
        query Pedigree($project: String!, $sample_external_IDs: [String!]!) {
      project(name: $project) {
    sequencingGroups(
      id: {in_: $sample_external_IDs}) {
      id
      meta
      sample {
        meta
        id
        participant {
          externalId
          samples {
            id
            externalId
            sequencingGroups(type: {eq: "genome"}, technology: {eq: "short-read"}) {
              id
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


def merge_overlapping_intervals(df: pd.DataFrame) -> list[dict]:
    """Merge overlapping buffered intervals, unioning their sample ID sets.

    Takes a DataFrame with columns: seqnames, start, end, genome_ids (set per row).
    Returns a list of dicts with keys: chrom, start, end, sample_ids (set of genome IDs).
    """
    df = df.copy()
    df['buf_start'] = df['start'] - BUFFER_BP
    df['buf_end'] = df['end'] + BUFFER_BP
    df = df.sort_values(['seqnames', 'buf_start']).reset_index(drop=True)

    merged = []
    for chrom, group in df.groupby('seqnames'):
        rows = group.sort_values('buf_start').to_dict('records')
        current = {
            'chrom': chrom,
            'start': rows[0]['buf_start'],
            'end': rows[0]['buf_end'],
            'sample_ids': set(rows[0]['genome_ids']),
        }
        for row in rows[1:]:
            if row['buf_start'] <= current['end']:
                current['end'] = max(current['end'], row['buf_end'])
                current['sample_ids'].update(row['genome_ids'])
            else:
                merged.append(current)
                current = {
                    'chrom': chrom,
                    'start': row['buf_start'],
                    'end': row['buf_end'],
                    'sample_ids': set(row['genome_ids']),
                }
        merged.append(current)

    return merged


def build_interval_table(merged_intervals: list[dict]) -> hl.Table:
    """Build a Hail Table keyed by non-overlapping intervals with associated sample ID sets."""
    rows = [
        hl.Struct(chrom=iv['chrom'], start=iv['start'], end=iv['end'], sample_ids=sorted(iv['sample_ids']))
        for iv in merged_intervals
    ]
    ht = hl.Table.parallelize(
        rows,
        schema=hl.tstruct(chrom=hl.tstr, start=hl.tint32, end=hl.tint32, sample_ids=hl.tarray(hl.tstr)),
    )
    ht = ht.annotate(
        interval=hl.locus_interval(ht.chrom, ht.start, ht.end, includes_end=True, reference_genome=REFERENCE_GENOME),
        csv_sample_ids=hl.set(ht.sample_ids),
    )
    ht = ht.select('interval', 'csv_sample_ids')
    return ht.key_by('interval')


def main():  # noqa: PLR0915
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mt', required=True, help='GCS path to the seqr-loader MatrixTable')
    parser.add_argument('--csv', required=True, help='Path to FRASER significant results CSV')
    parser.add_argument('--rna_ids', required=True, help='RNA sample IDs', nargs='+')
    parser.add_argument('--query_dataset', help='Metamist project name', default='rdnow')
    parser.add_argument('--output', required=True, help='Output path for BED-like TSV')
    args = parser.parse_args()

    hail_batch.init_batch()
    # Step 0: RNA sample ids to a set of their Genome sample IDs

    relevant_ids = list(set(args.rna_ids))
    logger.info(f'Input RNA sample IDs: {relevant_ids}')
    variables = {'project': args.query_dataset, 'sample_external_IDs': relevant_ids}

    # build a dictionary mapping from RNA sample ID to set of genome sample IDs for that participant
    result = query(query_ids, variables=variables)
    rna_to_genome_ids = {}
    for group in result['project']['sequencingGroups']:
        rna_id = group['id']
        participant = group['sample']['participant']
        genome_ids = set()
        for sample in participant['samples']:
            for sg in sample['sequencingGroups']:
                if sg['type'] == 'genome':
                    genome_ids.add(sg['id'])
        rna_to_genome_ids[rna_id] = genome_ids

    # --- Step 1: Parse CSV and build merged intervals ---
    logger.info(f'Reading FRASER CSV: {args.csv}')
    df = pd.read_csv(args.csv)
    required_cols = {'seqnames', 'start', 'end', 'sampleID'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f'CSV missing required columns: {missing}')

    # remap the sampleID column to be the set of genome sample IDs for that RNA sample ID
    def map_to_genome_ids(rna_id) -> set[str]:
        genome_ids = rna_to_genome_ids.get(rna_id)
        if genome_ids is None:
            logger.warning(f'RNA sample ID {rna_id} not found in query results, skipping')
            return set()
        return genome_ids

    df['genome_ids'] = df['sampleID'].apply(map_to_genome_ids)
    df = df[df['genome_ids'].map(len) > 0].reset_index(drop=True)

    cpg_ids = set().union(*df['genome_ids'])
    logger.info(f'Found {len(cpg_ids)} unique genome sample IDs after mapping from RNA IDs')
    logger.info(f'Found {len(df)} FRASER significant regions')

    merged = merge_overlapping_intervals(df)
    logger.info(f'Merged into {len(merged)} non-overlapping intervals (±{BUFFER_BP}bp buffer)')

    # Build Hail interval objects for filter_intervals
    hail_intervals = [
        hl.parse_locus_interval(
            f'[{iv["chrom"]}:{iv["start"]}-{iv["end"]}]',
            reference_genome=REFERENCE_GENOME,
        )
        for iv in merged
    ]

    interval_ht = build_interval_table(merged)

    # --- Step 2: Subset MT to regions ---
    logger.info(f'Reading MatrixTable: {args.mt}')
    mt = hl.read_matrix_table(args.mt)

    logger.info('Filtering to FRASER significant regions')
    mt = hl.filter_intervals(mt, hail_intervals)

    # --- Step 3: Drop to rows table ---
    ht = mt.rows()

    # --- Step 4: Identify variants of interest ---
    logger.info('Annotating variants with CSV region sample IDs')
    ht = ht.annotate(csv_sample_ids=interval_ht.index(ht.locus).csv_sample_ids)

    ht = ht.annotate(
        carriers=ht.samples_num_alt['1'].union(ht.samples_num_alt['2']),
    )
    ht = ht.annotate(
        matching_samples=ht.carriers.intersection(ht.csv_sample_ids),
    )

    ht = ht.filter(ht.matching_samples.length() > 0)

    # Rare variant filter: keep if gnomAD popmax AF < 0.01 or missing (novel)
    gnomad_genomes_af = ht.gnomad_genomes.AF_POPMAX_OR_GLOBAL
    gnomad_exomes_af = ht.gnomad_exomes.AF_POPMAX_OR_GLOBAL
    ht = ht.filter(
        (hl.is_missing(gnomad_genomes_af) | (gnomad_genomes_af < MAX_POPMAX_AF))
        & (hl.is_missing(gnomad_exomes_af) | (gnomad_exomes_af < MAX_POPMAX_AF)),
    )

    # --- Step 5: Annotate BED fields ---
    # BED is 0-based half-open: start = pos - 1, end = pos - 1 + len(ref)
    ht = ht.annotate(
        bed_chrom=ht.locus.contig,
        bed_start=ht.locus.position - 1,
        bed_end=ht.locus.position - 1 + ht.alleles[0].length(),
    )

    ht = ht.drop('csv_sample_ids', 'carriers')

    n_variants = ht.count()
    logger.info(f'Found {n_variants} variants of interest')

    # --- Export full TSV with all annotations ---
    tsv_path = args.output.replace('.bed', '.tsv') if args.output.endswith('.bed') else args.output + '.tsv'
    ht_flat = ht.flatten()
    logger.info(f'Exporting full TSV to {tsv_path}')
    ht_flat.export(tsv_path, delimiter='\t')

    # --- Export minimal IGV-compatible BED ---
    bed_ht = ht.select(
        bed_chrom=ht.bed_chrom,
        bed_start=ht.bed_start,
        bed_end=ht.bed_end,
        name=hl.or_else(ht.mainTranscript.gene_symbol, 'intergenic')
        + '|'
        + hl.or_else(ht.mainTranscript.major_consequence, 'unknown'),
    )
    bed_ht = bed_ht.key_by()
    bed_ht = bed_ht.select('bed_chrom', 'bed_start', 'bed_end', 'name')
    logger.info(f'Exporting IGV BED to {args.output}')
    bed_ht.export(args.output, delimiter='\t', header=False)
    logger.info('Done.')


if __name__ == '__main__':
    main()
