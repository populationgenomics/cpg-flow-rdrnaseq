"""Subset an AnnotateCohort MatrixTable to a set of sample column IDs."""

import argparse

import hail as hl

from cpg_utils.hail_batch import init_batch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='GCS path to source AnnotateCohort MT')
    parser.add_argument('--sgs', required=True, help='Path to newline-delimited genome SG IDs file')
    parser.add_argument('--output', required=True, help='GCS path to write subsetted MT')
    args = parser.parse_args()

    init_batch()

    with hl.current_backend().fs.open(args.sgs, 'r') as f:
        sg_ids = hl.literal({line.strip() for line in f if line.strip()})

    mt = hl.read_matrix_table(args.input)
    mt = mt.filter_cols(sg_ids.contains(mt.s))
    mt.write(args.output, overwrite=True)


if __name__ == '__main__':
    main()
