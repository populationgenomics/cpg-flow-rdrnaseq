"""
Run FastQ Screen to detect contamination from other organisms.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cpg_flow.resources import STANDARD
from cpg_utils import Path, config
from cpg_utils.hail_batch import command, get_batch

if TYPE_CHECKING:
    import hailtop.batch as hb
    from hailtop.batch.job import Job

    from cpg_flow.filetypes import FastqPairs


BT2_SUFFIXES = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2')


def _localize_genomes(b: hb.Batch) -> dict[str, hb.ResourceGroup]:
    """Localize bowtie2 index files for each configured genome via read_input_group."""
    genomes: dict[str, str] = config.config_retrieve(['references', 'fastq_screen_genomes'])
    localized = {}
    for name, gcs_prefix in genomes.items():
        files = {suffix.lstrip('.'): f'{gcs_prefix}{suffix}' for suffix in BT2_SUFFIXES}
        localized[name] = b.read_input_group(**files)
    return localized


def fastq_screen(
    input_fq_pairs: FastqPairs,
    output_txt: Path,
    output_html: Path,
    job_attrs: dict[str, str],
) -> list[Job]:
    """Screen R1 of the first FASTQ pair against a reference panel for contamination."""
    b = get_batch()
    nthreads = config.config_retrieve(['workflow', 'fastq_screen', 'nthreads'], 8)

    genomes = _localize_genomes(b)
    fq_resources = input_fq_pairs[0].as_resources(b)

    j = b.new_bash_job('FastqScreen', job_attrs | {'tool': 'fastq_screen'})
    j.image(config.config_retrieve(['images', 'fastq_screen']))
    STANDARD.set_resources(j=j, ncpu=nthreads, storage_gb=75)

    conf_lines = [f'THREADS\t{nthreads}']
    for name, rg in genomes.items():
        conf_lines.append(f'DATABASE\t{name}\t{rg}')
    write_conf = ' && '.join([f"echo '{line}' >> /tmp/fastq_screen.conf" for line in conf_lines])

    j.command(
        command(
            f"""\
            {write_conf}
            mkdir -p output
            fastq_screen --aligner bowtie2 --conf /tmp/fastq_screen.conf \
              --threads {nthreads} {fq_resources.r1} --outdir output/
            cp output/*_screen.txt {j.screen_txt}
            cp output/*_screen.html {j.screen_html}
            """,
            monitor_space=True,
        ),
    )

    b.write_output(j.screen_txt, str(output_txt))
    b.write_output(j.screen_html, str(output_html))

    return [j]
