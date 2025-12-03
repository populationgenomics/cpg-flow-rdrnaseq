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
from cpg_flow.utils import can_reuse
from cpg_utils import Path, to_path
from cpg_utils.config import image_path, reference_path
from cpg_utils.hail_batch import Batch, command
from hailtop.batch.job import Job

from rdrnaseq.jobs.bam_to_cram import bam_to_cram
from rdrnaseq.jobs.markdups import markdup


class STAR:
    """
    Construct a STAR command for aligning FASTQs.
    """

    def __init__(
        self,
        input_fastq_pair: FastqPair,
        sample_name: str,
        genome: hb.ResourceGroup,
        nthreads: int,
        read_group: dict[str, str] | None = None,
        bamout: bool = True,
        sort: bool = True,
    ):
        self.bamout = bamout
        self.sort = sort
        self.nthreads = nthreads

        # Define output filenames and types
        if self.bamout:
            if self.sort:
                self.outSAMtype = 'BAM SortedByCoordinate'
                self.outStd = 'BAM_SortedByCoordinate'
                self.output_filename = 'Aligned.sortedByCoord.out.bam'
            else:
                self.outSAMtype = 'BAM Unsorted'
                self.outStd = 'BAM_Unsorted'
                self.output_filename = 'Aligned.out.bam'
            self.output_extension = 'bam'
        else:
            self.outSAMtype = 'SAM'
            self.outStd = 'Log'
            self.output_filename = 'Aligned.out.sam'
            self.output_extension = 'sam'

        # Handle Read Group
        self.read_group = read_group or {}
        self.read_group.setdefault('ID', sample_name)
        self.read_group.setdefault('SM', sample_name)
        self.read_group_line = self._create_read_group_line()

        # Build Command
        # Note: input_fastq_pair.r1/r2 are typically typically Hail ResourceFiles,
        # casting to str gives their container path.
        self.command_parts = [
            'STAR',
            '--runThreadN',
            str(self.nthreads),
            '--genomeDir',
            str(genome.genome.dirname),
            '--outSAMtype',
            self.outSAMtype,
            '--outStd',
            self.outStd,
            '--outFileNamePrefix',
            'Aligned.',
            '--outSAMattrRGline',
            self.read_group_line,
            '--readFilesCommand',
            'zcat',
            '--readFilesIn',
            str(input_fastq_pair.r1),
            str(input_fastq_pair.r2),
        ]

    def __str__(self):
        return ' '.join(self.command_parts)

    def __repr__(self):
        return self.__str__()

    def _create_read_group_line(self) -> str:
        """
        Create a read group line for the STAR command.
        """
        parts = []
        for k, v in self.read_group.items():
            # Assign the final value using a conditional expression
            quoted_v = f'"{v}"' if ' ' in v and not (v.startswith('"') and v.endswith('"')) else v

            parts.append(f'{k}:{quoted_v}')
        return ' '.join(parts)


class GCPStarReference:
    def __init__(self, b: Batch, genome_prefix: str | Path):
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
        self.genome_res_group = b.read_input_group(**self.genome_files)


def align(
    b: Batch,
    fastq_pairs: FastqPairs,
    sample_name: str,
    genome_prefix: str | Path,
    mark_duplicates: bool = True,
    output_bam: BamPath | None = None,
    output_cram: CramPath | None = None,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> list[Job] | None:
    """
    Align (potentially multiple) FASTQ pairs using STAR,
    merge the resulting BAMs (if necessary),
    then sort and index the resulting BAM.
    """
    if output_cram and can_reuse(output_cram, overwrite):
        return None
    if not output_cram and output_bam and can_reuse(output_bam, overwrite):
        return None

    if not isinstance(fastq_pairs, FastqPairs):
        raise TypeError(f'fastq_pairs must be a FastqPairs object, not {type(fastq_pairs)}')
    if len(fastq_pairs) == 0:
        raise ValueError('fastq_pairs must contain at least one FastqPair')

    # Optimization: Instantiate Reference object once here, pass it down
    star_ref = GCPStarReference(b=b, genome_prefix=genome_prefix)

    merge = len(fastq_pairs) > 1
    jobs = []
    aligned_bams = []

    for job_idx, fq_pair in enumerate(fastq_pairs, 1):
        if not isinstance(fq_pair, FastqPair):
            raise TypeError(f'fastq_pairs must contain FastqPair objects, not {type(fq_pair)}')
        label = f'{extra_label} {job_idx}' if extra_label else f'{job_idx}'
        j, bam = align_fq_pair(
            b=b,
            fastq_pair=fq_pair,
            sample_name=sample_name,
            star_ref=star_ref,  # Pass the instantiated ref
            extra_label=label,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        jobs.append(j)
        aligned_bams.append(bam)

    # Logic Flow:
    # 1. STAR produces Sorted BAMs.
    # 2. Merge (samtools merge) preserves sort order of inputs.
    # 3. Therefore, 'aligned_bam' below is already sorted.

    if merge:
        j, merged_bam = merge_bams(
            b=b,
            input_bams=aligned_bams,
            extra_label=extra_label,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        jobs.append(j)
        aligned_bam = merged_bam
    else:
        aligned_bam = aligned_bams[0]

    # Optimization: Skip heavy re-sorting if we know STAR/Merge kept it sorted.
    # We pass 'assume_sorted=True' because STAR output is sorted and Merge preserves it.
    j, sorted_bam_group = sort_index_bam(
        b=b,
        input_bam=aligned_bam,
        extra_label=extra_label,
        job_attrs=job_attrs,
        requested_nthreads=requested_nthreads,
        assume_sorted=True,
    )
    jobs.append(j)

    # The output of sort_index_bam is a ResourceGroup containing .bam and .bam.bai
    # We extract the bam file for the next step.
    sorted_bam_file = sorted_bam_group['bam']

    if mark_duplicates:
        j, mkdup_bam = markdup(
            b=b,
            input_bam=sorted_bam_file,
            extra_label=extra_label,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        jobs.append(j)
        out_bam = mkdup_bam
    else:
        out_bam = sorted_bam_file

    # Output writing
    if output_bam:
        out_bam_path = to_path(output_bam.path)
        b.write_output(out_bam, str(out_bam_path.with_suffix('')))

    if output_cram:
        j, out_cram = bam_to_cram(
            b=b,
            input_bam=out_bam,
            extra_label=extra_label,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
            reference_fasta_path=reference_path('star/fasta'),
        )
        jobs.append(j)
        out_cram_path = to_path(output_cram.path)
        b.write_output(out_cram, str(out_cram_path.with_suffix('')))

    return jobs


def align_fq_pair(
    b: Batch,
    fastq_pair: FastqPair,
    sample_name: str,
    star_ref: GCPStarReference,  # Received as argument
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
    genome_prefix: str | Path | None = None,  # kept for backward compat, unused #noqa:ARG001
) -> tuple[Job, hb.ResourceFile]:
    """
    Takes an input FastqPair object, and creates a job to align it using STAR.
    """
    job_name = 'align_rna'
    if extra_label:
        job_name += f' {extra_label}'

    j_attrs = (job_attrs or {}) | {'label': job_name, 'tool': 'STAR'}
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('star'))

    nthreads = requested_nthreads or 8
    # Optimize storage: 200GB is generous.
    # If possible, could check fastq size, but keeping safe default.
    res = HIGHMEM.set_resources(j, ncpu=nthreads, storage_gb=200)

    # Use resources N-1 for STAR, leaving 1 for system overhead/zcat
    star_nthreads = max(1, res.get_nthreads() - 1)

    star = STAR(
        input_fastq_pair=fastq_pair,
        sample_name=sample_name,
        genome=star_ref.genome_res_group,
        nthreads=star_nthreads,
        bamout=True,
        sort=True,  # STAR will coordinate-sort
    )

    # Atomic move is safer
    cmd = command(f'{star!s} && mv {star.output_filename} {j.output_bam}', monitor_space=True)
    j.command(cmd)

    return j, j.output_bam


def merge_bams(
    b: Batch,
    input_bams: list[str | BamPath | Path],
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, hb.ResourceFile]:
    """
    Merge a list of BAM files into a single BAM file.
    """
    job_name = 'merge_bams'
    if extra_label:
        job_name += f' {extra_label}'

    j_attrs = (job_attrs or {}) | {'label': job_name, 'tool': 'samtools'}
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(j, ncpu=nthreads, storage_gb=50)

    # Samtools merge preserves input sort order.
    # Added -c -p usually helps with large merges to combine headers and comments properly.
    cmd = f'samtools merge -@ {res.get_nthreads() - 1} -o {j.merged_bam} {" ".join([str(b) for b in input_bams])}'
    j.command(command(cmd, monitor_space=True))
    return j, j.merged_bam


def sort_index_bam(
    b: Batch,
    input_bam: str | BamPath | Path,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
    assume_sorted: bool = False,
) -> tuple[Job, hb.ResourceGroup]:
    """
    Sort and index a BAM file.
    If assume_sorted is True, it skips sorting and only indexes.
    """
    job_name = 'sort_index_bam'
    if extra_label:
        job_name += f' {extra_label}'

    j_attrs = (job_attrs or {}) | {'label': job_name, 'tool': 'samtools'}
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(j, ncpu=nthreads, storage_gb=50)

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
