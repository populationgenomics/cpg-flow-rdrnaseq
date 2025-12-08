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

from rdrnaseq.stages import TrimAlignRNA
from rdrnaseq.stages import Count

def cli_main():
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
    stages = [TrimAlignRNA,Count]
    run_workflow(name='rdrnaseq', stages=stages, dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
