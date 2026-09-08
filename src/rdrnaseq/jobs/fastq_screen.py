"""
Run FastQ Screen to detect contamination from other organisms.
"""

from hailtop.batch.job import Job

from cpg_flow.filetypes import FastqPairs
from cpg_flow.resources import STANDARD
from cpg_utils import Path, config
from cpg_utils.hail_batch import command, get_batch


def fastq_screen(
    input_fq_pairs: FastqPairs,
    output_txt: Path,
    output_html: Path,
    job_attrs: dict[str, str],
) -> list[Job]:
    """Screen R1 of the first FASTQ pair against a reference panel for contamination."""
    b = get_batch()
    nthreads = config.config_retrieve(['workflow', 'fastq_screen', 'nthreads'], 8)

    conf_file = b.read_input(config.config_retrieve(['references', 'fastq_screen_conf']))
    fq_resources = input_fq_pairs[0].as_resources(b)

    j = b.new_bash_job('FastqScreen', job_attrs | {'tool': 'fastq_screen'})
    j.image(config.config_retrieve(['images', 'fastq_screen']))
    STANDARD.set_resources(j=j, ncpu=nthreads, storage_gb=30)

    j.command(
        command(
            f"""\
            mkdir -p output
            fastq_screen --aligner bowtie2 --conf {conf_file} --threads {nthreads} {fq_resources.r1} --outdir output/
            cp output/*_screen.txt {j.screen_txt}
            cp output/*_screen.html {j.screen_html}
            """,
            monitor_space=True,
        ),
    )

    b.write_output(j.screen_txt, str(output_txt))
    b.write_output(j.screen_html, str(output_html))

    return [j]
