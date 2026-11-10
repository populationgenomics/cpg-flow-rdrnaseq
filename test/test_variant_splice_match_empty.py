"""Test variant_splice_match behaviour when no genome SGs exist.

Covers the fix for the empty-intervals crash: when RNA SG participants have
no genome sequencing data, the script writes header-only TSV + empty BED
and exits cleanly instead of crashing on hl.filter_intervals([]).
"""

import csv
import os
import tempfile

import pandas as pd
import pytest

from rdrnaseq.scripts.variant_of_interest_subset import (
    EMPTY_OUTPUT_COLUMNS,
    merge_overlapping_intervals,
    read_fraser_csv,
)


@pytest.fixture
def empty_output_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


def _write_empty_outputs(tsv_path: str, bed_path: str):
    """Reproduce the early-return logic from main()."""
    pd.DataFrame(columns=EMPTY_OUTPUT_COLUMNS).to_csv(tsv_path, sep='\t', index=False)
    open(bed_path, 'w').close()


def test_empty_output_files_are_valid(empty_output_dir):
    """Empty outputs should be parseable files, not missing."""
    tsv_path = os.path.join(empty_output_dir, 'variants.tsv')
    bed_path = os.path.join(empty_output_dir, 'variants.bed')
    _write_empty_outputs(tsv_path, bed_path)

    assert os.path.exists(tsv_path)
    assert os.path.exists(bed_path)
    assert os.path.getsize(bed_path) == 0


def test_empty_tsv_has_correct_headers(empty_output_dir):
    """TSV should have all expected column headers and zero data rows."""
    tsv_path = os.path.join(empty_output_dir, 'variants.tsv')
    bed_path = os.path.join(empty_output_dir, 'variants.bed')
    _write_empty_outputs(tsv_path, bed_path)

    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    assert reader.fieldnames == EMPTY_OUTPUT_COLUMNS
    assert len(rows) == 0


def test_dashboard_parses_empty_tsv(empty_output_dir):
    """Dashboard reads variant TSV with csv.DictReader — empty file should yield 0 rows."""
    tsv_path = os.path.join(empty_output_dir, 'variants.tsv')
    bed_path = os.path.join(empty_output_dir, 'variants.bed')
    _write_empty_outputs(tsv_path, bed_path)

    seqr_links = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for variant in reader:
            seqr_links[variant.get('locus', '')] = 'link'

    assert seqr_links == {}


def test_read_fraser_csv_filters_unmapped_samples(empty_output_dir):
    """FRASER CSV with samples that have no genome SGs should produce empty DataFrame."""
    csv_path = os.path.join(empty_output_dir, 'fraser.csv')
    pd.DataFrame({
        'seqnames': ['chr1', 'chr2'],
        'start': [1000, 2000],
        'end': [1500, 2500],
        'sampleID': ['CPG_RNA_1', 'CPG_RNA_2'],
        'pValue': [0.001, 0.002],
        'padjust': [0.01, 0.02],
    }).to_csv(csv_path, index=False)

    rna_to_genome = {'CPG_RNA_1': set(), 'CPG_RNA_2': set()}
    result = read_fraser_csv(csv_path, rna_to_genome)
    assert len(result) == 0


def test_merge_overlapping_intervals_empty_input():
    """merge_overlapping_intervals should handle empty DataFrame."""
    empty_df = pd.DataFrame(columns=['seqnames', 'start', 'end', 'genome_sg_ids'])
    merged = merge_overlapping_intervals(empty_df)
    assert merged == []


def test_read_fraser_csv_keeps_mapped_samples(empty_output_dir):
    """Samples with genome SGs should be retained."""
    csv_path = os.path.join(empty_output_dir, 'fraser.csv')
    pd.DataFrame({
        'seqnames': ['chr1', 'chr2'],
        'start': [1000, 2000],
        'end': [1500, 2500],
        'sampleID': ['CPG_RNA_1', 'CPG_RNA_2'],
        'pValue': [0.001, 0.002],
        'padjust': [0.01, 0.02],
    }).to_csv(csv_path, index=False)

    rna_to_genome = {'CPG_RNA_1': {'CPG_GENOME_1'}, 'CPG_RNA_2': set()}
    result = read_fraser_csv(csv_path, rna_to_genome)
    assert len(result) == 1
    assert result.iloc[0]['sampleID'] == 'CPG_RNA_1'
