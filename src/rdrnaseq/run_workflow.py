#!/usr/bin/env python3

"""
This is the main entry point for the workflow.
This is a re-implementation of the canonical main.py file in production-pipelines.
The purpose of this script is to import all the Stages in the workflow (or at least the terminal workflow nodes)
and begin the CPG-Flow Stage discovery and graph construction process

This is re-implemented as a simpler form, only knowing how to build a single workflow, instead of choosing at runtime
"""

from argparse import ArgumentParser

from cpg_flow.targets.types import Sex
from cpg_flow.workflow import run_workflow

from rdrnaseq.stages import (
    Count,
    Dashboard,
    FastqScreen,
    Fraser,
    GenerateSgQcReport,
    Outrider,
    PicardRnaSeqMetrics,
    QcMultiQC,
    SamtoolsStats,
    Somalier,
    TrimAlignRNA,
    VariantSpliceMatch,
)

_cpg_flow_sex_parse = Sex.parse


def _parse_sex_allowing_missing(input_sex: str | int | None) -> Sex:
    """
    Hotfix: cpg-flow reads pedigree sex as Sex.parse(str(ped_entry['sex'])), so a participant with no
    reported sex arrives as the string 'None' and raises. Treat a missing sex as UNKNOWN instead.
    """
    if input_sex is None or (isinstance(input_sex, str) and input_sex.strip().lower() in {'', 'none'}):
        return Sex.UNKNOWN
    return _cpg_flow_sex_parse(input_sex)


Sex.parse = staticmethod(_parse_sex_allowing_missing)


def cli_main() -> None:
    """
    CLI entrypoint - starts up the workflow
    """
    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    # Note - in production-pipelines the main.py script sets up layers of default configuration,
    # overlaid with workflow-specific configuration, and then runs the workflow.
    # If you want to re-use that model, this should be carried out before entering the workflow

    # Otherwise all configuration should be done by providing all relevant configs to analysis-runner
    # https://github.com/populationgenomics/team-docs/blob/main/cpg_utils_config.md#config-in-analysis-runner-jobs
    stages = [
        FastqScreen,
        TrimAlignRNA,
        Count,
        Somalier,
        SamtoolsStats,
        PicardRnaSeqMetrics,
        QcMultiQC,
        GenerateSgQcReport,
        Fraser,
        Outrider,
        VariantSpliceMatch,
        Dashboard,
    ]
    run_workflow(name='rdrnaseq', stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
