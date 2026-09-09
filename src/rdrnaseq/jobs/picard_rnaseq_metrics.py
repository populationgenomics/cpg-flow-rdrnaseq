"""
Run Picard CollectRnaSeqMetrics on a BAM file for RNA-specific QC.
"""

from hailtop.batch.job import Job

from cpg_flow.resources import STANDARD
from cpg_utils import Path, config
from cpg_utils.hail_batch import command, get_batch


def collect_rnaseq_metrics(
    input_bam: str | Path,
    output_metrics: Path,
    job_attrs: dict[str, str],
) -> Job:
    """Run Picard CollectRnaSeqMetrics and write the output."""
    b = get_batch()

    j = b.new_bash_job('PicardRnaSeqMetrics', job_attrs | {'tool': 'picard'})
    j.image(config.config_retrieve(['images', 'picard']))
    STANDARD.set_resources(j=j, ncpu=2, storage_gb=30)

    star_fasta = config.config_retrieve(['references', 'star', 'fasta'])
    star_dict = config.config_retrieve(['references', 'star', 'dict'], star_fasta.replace('.fa', '.dict'))
    reference = b.read_input_group(
        base=star_fasta,
        dict=star_dict,
        fai=f'{star_fasta}.fai',
    )
    bam_input = b.read_input_group(**{'bam': str(input_bam), 'bam.bai': f'{input_bam}.bai'})
    ref_flat = b.read_input(config.config_retrieve(['references', 'ref_flat']))

    rib_intervals_cmd = ''
    rib_intervals_path = config.config_retrieve(['references', 'ribosomal_intervals'], None)
    if rib_intervals_path:
        rib_intervals = b.read_input(rib_intervals_path)
        rib_intervals_cmd = f'-RIBOSOMAL_INTERVALS {rib_intervals}'

    j.command(
        command(
            f"""\
            picard -Xms1g -Xmx3g CollectRnaSeqMetrics \
              -I {bam_input.bam} \
              -O {j.metrics} \
              -R {reference.base} \
              -REF_FLAT {ref_flat} \
              {rib_intervals_cmd} \
              -STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND \
              -VALIDATION_STRINGENCY SILENT
            """,
            monitor_space=True,
        ),
    )

    b.write_output(j.metrics, str(output_metrics))
    return j
