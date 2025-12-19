"""
Align RNA-seq reads to the genome using STAR.

This module provides functions and classes to run the STAR aligner
within a Hail Batch workflow.
"""

import hailtop.batch as hb
from cpg_flow.filetypes import (
    BamPath,
    CramPath,
    FastqPair,
    FastqPairs,
)
from cpg_flow.resources import HIGHMEM, STANDARD
from cpg_utils import Path, config, to_path
from cpg_utils.config import image_path, reference_path
from cpg_utils.hail_batch import command, get_batch
from hailtop.batch.job import Job

from .bam_to_cram import bam_to_cram
from .markdups import markdup


class GCPStarReference:
    def __init__(self):
        genome_prefix = config.config_retrieve(['references', 'star', 'ref_dir'])

        # Ensure it's a string, strip trailing slash
        gp = str(genome_prefix).rstrip('/')

        # Dictionary comprehension is fine, but cleaner to list implicitly
        # The key names map to the expected attributes in the ResourceGroup
        files = {
            'chr_len': 'chrLength.txt',
            'chr_name_len': 'chrNameLength.txt',
            'chr_name': 'chrName.txt',
            'chr_start': 'chrStart.txt',
            'exon_ge_tr_info': 'exonGeTrInfo.tab',
            'exon_info': 'exonInfo.tab',
            'gene_info': 'geneInfo.tab',
            'genome': 'Genome',
            'genome_params': 'genomeParameters.txt',
            'sa': 'SA',
            'sa_idx': 'SAindex',
            'sjdb_info': 'sjdbInfo.txt',
            'sjdb_list_gtf': 'sjdbList.fromGTF.out.tab',
            'sjdb_list': 'sjdbList.out.tab',
            'transcript_info': 'transcriptInfo.tab',
        }

        # Prepend prefix
        self.genome_files = {k: f'{gp}/{v}' for k, v in files.items()}
        self.genome_res_group = get_batch().read_input_group(**self.genome_files)


def align(
    fastq_pairs: FastqPairs,
    sample_name: str,
    job_attrs: dict,
    output_bam: BamPath,
    output_cram: CramPath,
) -> list[Job] | None:
    """
    Align (potentially multiple) FASTQ pairs using STAR,
    merge the resulting BAMs (if necessary),
    then sort and index the resulting BAM.
    """

    b = get_batch()

    if not isinstance(fastq_pairs, FastqPairs):
        raise TypeError(f'fastq_pairs must be a FastqPairs object, not {type(fastq_pairs)}')
    if len(fastq_pairs) == 0:
        raise ValueError('fastq_pairs must contain at least one FastqPair')

    # Optimization: Instantiate Reference object once here, pass it down
    star_ref = GCPStarReference()

    merge = len(fastq_pairs) > 1
    jobs = []
    aligned_bams = []

    for _job_idx, fq_pair in enumerate(fastq_pairs, 1):
        if not isinstance(fq_pair, FastqPair):
            raise TypeError(f'fastq_pairs must contain FastqPair objects, not {type(fq_pair)}')
        j, bam = align_fq_pair(
            fastq_pair=fq_pair,
            sample_name=sample_name,
            star_ref=star_ref,  # Pass the instantiated ref
            job_attrs=job_attrs,
        )
        jobs.append(j)
        aligned_bams.append(bam)

    # Logic Flow:
    # 1. STAR produces Sorted BAMs.
    # 2. Merge (samtools merge) preserves sort order of inputs.
    # 3. Therefore, 'aligned_bam' below is already sorted.

    if merge:
        j, merged_bam = merge_bams(
            input_bams=aligned_bams,
            job_attrs=job_attrs,
        )
        jobs.append(j)
        aligned_bam = merged_bam
    else:
        aligned_bam = aligned_bams[0]

    # Optimization: Skip heavy re-sorting if we know STAR/Merge kept it sorted.
    # We pass 'assume_sorted=True' because STAR output is sorted and Merge preserves it.
    j, sorted_bam_group = sort_index_bam(
        input_bam=aligned_bam,
        job_attrs=job_attrs,
        assume_sorted=True,
    )
    jobs.append(j)

    # The output of sort_index_bam is a ResourceGroup containing .bam and .bam.bai
    # We extract the bam file for the next step.

    j, mkdup_bam = markdup(
        input_bam=sorted_bam_group,
        job_attrs=job_attrs,
        requested_nthreads=4,
    )
    jobs.append(j)
    out_bam = mkdup_bam

    # Output writing
    out_bam_path = to_path(output_bam.path)
    b.write_output(out_bam, str(out_bam_path.with_suffix('')))

    j, out_cram = bam_to_cram(
        input_bam=out_bam,
        job_attrs=job_attrs,
        requested_nthreads=4,
        reference_fasta_path=reference_path('broad/ref_fasta'),
    )
    jobs.append(j)
    out_cram_path = to_path(output_cram.path)
    b.write_output(out_cram, str(out_cram_path.with_suffix('')))

    return jobs


def align_fq_pair(
    fastq_pair: FastqPair,
    sample_name: str,
    star_ref: GCPStarReference,  # Received as argument
    job_attrs: dict,
) -> tuple[Job, hb.ResourceFile]:
    """
    Takes an input FastqPair object, and creates a job to align it using STAR.
    """
    b = get_batch()

    j = b.new_job(name='align_rna', attributes=job_attrs | {'tool': 'STAR'})
    j.image(image_path('star'))

    nthreads = 8

    # Optimize storage: 200GB is generous.
    # If possible, could check fastq size, but keeping safe default.
    res = HIGHMEM.set_resources(j=j, ncpu=nthreads, storage_gb=200)

    j.command(f"""
        STAR \\
        --runThreadN {(res.get_nthreads() - 1)} \\
        --genomeDir $(dirname {star_ref.genome_res_group.genome!s}) \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattrRGline ID:{sample_name} SN:{sample_name} \\
        --readFilesCommand zcat \\
        --readFilesIn {fastq_pair.r1} {fastq_pair.r2}
        mv Aligned.sortedByCoord.out.bam {j.output_bam}
    """)

    return j, j.output_bam


def merge_bams(
    input_bams: list[str | BamPath | Path],
    job_attrs: dict,
) -> tuple[Job, hb.ResourceFile]:
    """
    Merge a list of BAM files into a single BAM file.
    """
    b = get_batch()

    j = b.new_job(name='merge_bams', attributes=job_attrs | {'tool': 'samtools'})
    j.image(image_path('samtools'))

    nthreads = 4
    res = STANDARD.set_resources(j=j, ncpu=nthreads, storage_gb=50)

    # Samtools merge preserves input sort order.
    # Added -c -p usually helps with large merges to combine headers and comments properly.
    cmd = f'samtools merge -@ {res.get_nthreads() - 1} -o {j.merged_bam} {" ".join([str(b) for b in input_bams])}'
    j.command(command(cmd, monitor_space=True))
    return j, j.merged_bam


def sort_index_bam(
    input_bam: str | BamPath | Path,
    job_attrs: dict,
    assume_sorted: bool = False,
) -> tuple[Job, hb.ResourceGroup]:
    """
    Sort and index a BAM file.
    If assume_sorted is True, it skips sorting and only indexes.
    """
    b = get_batch()

    j = b.new_job(name='sort_index_bam', attributes=job_attrs | {'tool': 'samtools'})
    j.image(image_path('samtools'))

    nthreads = 4
    res = STANDARD.set_resources(j=j, ncpu=nthreads, storage_gb=50)

    j.declare_resource_group(
        sorted_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai',
        },
    )

    # FIXED: Replaced pipe-to-index with sequential commands.
    # 'samtools index' cannot reliably read from stdin in a pipe.

    if assume_sorted:
        # Efficiency: Just copy (or move) and index.
        # cp is fast. We cp to the resource group filename.
        cmd = f"""
        cp {input_bam} {j.sorted_bam['bam']} && \\
        samtools index -@ {res.get_nthreads() - 1} {j.sorted_bam['bam']} {j.sorted_bam['bam.bai']}
        """
    else:
        # Full sort and index
        cmd = f"""
        samtools sort -@ {res.get_nthreads() - 1} -o {j.sorted_bam['bam']} {input_bam} && \\
        samtools index -@ {res.get_nthreads() - 1} {j.sorted_bam['bam']} {j.sorted_bam['bam.bai']}
        """

    j.command(command(cmd, monitor_space=True))
    return j, j.sorted_bam
