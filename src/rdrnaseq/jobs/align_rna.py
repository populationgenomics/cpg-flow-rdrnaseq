"""
Align RNA-seq reads to the genome using STAR.

This module provides functions and classes to run the STAR aligner
within a Hail Batch workflow.
"""

import re

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import image_path, reference_path
from cpg_utils.hail_batch import Batch, command
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
    FastqPair,
    FastqPairs,
)
from cpg_workflows.jobs.bam_to_cram import bam_to_cram
from cpg_workflows.jobs.markdups import markdup
from cpg_workflows.resources import HIGHMEM, STANDARD
from cpg_workflows.utils import can_reuse


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
        self.command = ['STAR']
        self.bamout = bamout
        self.sort = sort

        # Create read group line
        self.read_group = read_group or {}
        if 'ID' not in self.read_group:
            self.read_group['ID'] = sample_name
        if 'SM' not in self.read_group:
            self.read_group['SM'] = sample_name
        self.read_group_line = self._create_read_group_line()

        # Create outSAMtype and set the predictable output filename
        if self.bamout and self.sort:
            self.outSAMtype = 'BAM SortedByCoordinate'
            self.outStd = 'BAM_SortedByCoordinate'
            # This is the filename STAR will write to in the job
            self.output_filename = 'Aligned.sortedByCoord.out.bam'
        elif self.bamout:
            self.outSAMtype = 'BAM Unsorted'
            self.outStd = 'BAM_Unsorted'
            self.output_filename = 'Aligned.out.bam'
        else:
            self.outSAMtype = 'SAM'
            self.outStd = 'Log' # Ensures output goes to file (self.output_filename)
            self.output_filename = 'Aligned.out.sam'

        # Get output extension
        self.output_extension = 'bam' if self.bamout else 'sam'

        # Create command
        self.command.extend(
            [
                '--runThreadN',
                str(nthreads),
                '--genomeDir',
                # Cleaner way to get the resource directory path
                f'{str(genome.genome.dirname)}',
                '--outSAMtype',
                self.outSAMtype,
                '--outStd',
                # Only use 'Log' if we expect a file output, which is true here.
                # If STAR writes to STDOUT (e.g., self.outStd = 'BAM_Unsorted')
                # we handle the redirection in align_fq_pair below.
                self.outStd,
                '--outFileNamePrefix',
                # Tell STAR to use the prefix, so the output file is named correctly.
                # e.g., 'Aligned.sortedByCoord.out.bam'
                'Aligned.',
                '--outSAMattrRGline',
                self.read_group_line,
                '--readFilesCommand',
                'zcat',
                '--readFilesIn',
                str(input_fastq_pair.r1),
                str(input_fastq_pair.r2),
            ],
        )

    def __str__(self):
        # The STAR command itself. The file movement/redirection is handled
        # by the calling job function.
        return ' '.join(self.command)

    def __repr__(self):
        return self.__str__()

    def _create_read_group_line(self) -> str:
        """
        Create a read group line for the STAR command.
        """
        read_group_line = ''
        for k, v in self.read_group.items():
            v_quoted = re.sub(r'(^.*\s.*$)', r'"\1"', v)
            read_group_line += f'{k}:{v_quoted} '
        return read_group_line.strip()


class GCPStarReference:
    def __init__(self, b: Batch, genome_prefix: str | Path):
        gp = re.sub(r'/+$', '', genome_prefix)
        self.genome_files = {
            key: f'{gp}/{file}'
            for key, file in {
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
            }.items()
        }
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
    # Don't run if the output exists and can be reused
    if output_cram and can_reuse(output_cram, overwrite):
        return None
    if not output_cram and output_bam and can_reuse(output_bam, overwrite):
        return None

    if not isinstance(fastq_pairs, FastqPairs):
        raise TypeError(f'fastq_pairs must be a FastqPairs object, not {type(fastq_pairs)}')
    if len(fastq_pairs) == 0:
        raise ValueError('fastq_pairs must contain at least one FastqPair')

    merge = len(fastq_pairs) > 1

    jobs = []
    aligned_bams = []

    # baseline align jobs, no prior dependencies
    for job_idx, fq_pair in enumerate(fastq_pairs, 1):
        if not isinstance(fq_pair, FastqPair):
            raise TypeError(f'fastq_pairs must contain FastqPair objects, not {type(fq_pair)}')
        label = f'{extra_label} {job_idx}' if extra_label else f'{job_idx}'
        j, bam = align_fq_pair(
            b=b,
            fastq_pair=fq_pair,
            sample_name=sample_name,
            genome_prefix=genome_prefix,
            extra_label=label,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        jobs.append(j)
        aligned_bams.append(bam)

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

    j, sorted_bam = sort_index_bam(
        b=b,
        input_bam=aligned_bam,
        extra_label=extra_label,
        job_attrs=job_attrs,
        requested_nthreads=requested_nthreads,
    )

    # this job is always added
    jobs.append(j)

    # Mark duplicates
    if mark_duplicates:
        j, mkdup_bam = markdup(
            b=b,
            input_bam=sorted_bam,
            extra_label=extra_label,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        jobs.append(j)
        out_bam = mkdup_bam
    else:
        out_bam = sorted_bam

    if output_bam:
        out_bam_path = to_path(output_bam.path)
        b.write_output(out_bam, str(out_bam_path.with_suffix('')))

    # Convert to CRAM
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
    genome_prefix: str | Path,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, hb.ResourceFile]:
    """
    Takes an input FastqPair object, and creates a job to align it using STAR.
    """
    job_name = 'align_rna'
    if extra_label:
        job_name += f' {extra_label}'

    align_tool = 'STAR'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=align_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('star'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = HIGHMEM.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=200,  # TODO: make configurable
    )

    star_ref = GCPStarReference(b=b, genome_prefix=genome_prefix)

    # Use the number of threads STAR will use (N-1)
    star_nthreads = res.get_nthreads() - 1

    star = STAR(
        input_fastq_pair=fastq_pair,
        sample_name=sample_name,
        genome=star_ref.genome_res_group,
        nthreads=star_nthreads,
        # Keep bamout and sort flags if you want to allow them to be customized later
        bamout=True,
        sort=True,
        # Removed: output_path, stdout
    )

    # 1. Get the STAR command
    star_cmd = str(star)

    # 2. Add the file movement: move the predictable output file to the Hail Batch resource
    # NOTE: The check for 'BAM_SortedByCoordinate' or 'BAM_Unsorted' as outStd
    # means STAR outputs to a file, which is the desired behavior for Hail Batch.
    # We use 'mv' to rename the file.
    move_cmd = f'mv {star.output_filename} {j.output_bam}'

    # 3. Combine the full command
    cmd = command(f'{star_cmd} && {move_cmd}', monitor_space=True)

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

    merge_tool = 'samtools'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=merge_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,  # TODO: make configurable
    )

    cmd = f'samtools merge -@ {res.get_nthreads() - 1} -o {j.merged_bam} {" ".join([str(b) for b in input_bams])}'
    j.command(command(cmd, monitor_space=True))
    return j, j.merged_bam


def sort_index_bam(
    b: Batch,
    input_bam: str | BamPath | Path,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, hb.ResourceGroup]:
    """
    Sort and index a BAM file.
    """
    job_name = 'sort_index_bam'
    if extra_label:
        job_name += f' {extra_label}'

    sort_tool = 'samtools'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=sort_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,  # TODO: make configurable
    )

    j.declare_resource_group(
        sorted_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai',
        },
    )

    cmd = f'samtools sort -@ {res.get_nthreads() - 1} {input_bam} | tee {j.sorted_bam["bam"]} | samtools index -@ {res.get_nthreads() - 1} - {j.sorted_bam["bam.bai"]}'
    j.command(command(cmd, monitor_space=True))
    return j, j.sorted_bam