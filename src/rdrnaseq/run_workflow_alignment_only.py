#!/usr/bin/env python3

"""
This is the main entry point for the workflow.
This is a re-implementation of the canonical main.py file in production-pipelines.
The purpose of this script is to import all the Stages in the workflow (or at least the terminal workflow nodes)
and begin the CPG-Flow Stage discovery and graph construction process

This is re-implemented as a simpler form, only knowing how to build a single workflow, instead of choosing at runtime
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow

from rdrnaseq.stages import (
    Count,
    FastqScreen,
    GenerateSgQcReport,
    PicardRnaSeqMetrics,
    QcMultiQC,
    SamtoolsStats,
    Somalier,
    TrimAlignRNA,
)


def cli_main() -> None:
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    stages = [
        FastqScreen,
        TrimAlignRNA,
        Count,
        Somalier,
        SamtoolsStats,
        PicardRnaSeqMetrics,
        QcMultiQC,
        GenerateSgQcReport,
    ]
    run_workflow(name='rdrnaseq', stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
