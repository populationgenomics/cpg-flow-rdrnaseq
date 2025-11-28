"""
Trim raw FASTQ reads using fastp or cutadapt
"""

from dataclasses import dataclass
from enum import Enum

from cpg_utils import config
from cpg_utils.hail_batch import Batch, command
from cpg_workflows.filetypes import FastqPair
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse
from cpg_workflows.workflow import SequencingGroup
from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job

# Default configuration constants
DEFAULT_MIN_LENGTH = 50
DEFAULT_STORAGE_GB = 50
DEFAULT_NTHREADS = 8
DEFAULT_QUALITY_TRIM = 20


class MissingFastqInputException(Exception):
    """Raise if alignment input is missing"""



class InvalidSequencingTypeException(Exception):
    """Raise if alignment type is not 'rna'"""


@dataclass(frozen=True)
class AdapterSequence:
    """
    A class to represent an adapter sequence.
    """

    sequence: str
    name: str | None = None


@dataclass(frozen=True)
class AdapterPair:
    """
    A class to represent a pair of adapter sequences.
    """

    r1: AdapterSequence
    r2: AdapterSequence


class AdapterPairs(Enum):
    """
    A class to represent a set of adapter pairs.
    """

    ILLUMINA_TRUSEQ = AdapterPair(
        r1=AdapterSequence(
            sequence='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
            name='TruSeq Adapter Index 1',
        ),
        r2=AdapterSequence(
            sequence='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
            name='TruSeq Adapter Index 2',
        ),
    )

    ILLUMINA_NEXTERA = AdapterPair(
        r1=AdapterSequence(
            sequence='CTGTCTCTTATACACATCT',
            name='Nextera Adapter Index 1',
        ),
        r2=AdapterSequence(
            sequence='CTGTCTCTTATACACATCT',
            name='Nextera Adapter Index 2',
        ),
    )


class Cutadapt:
    """
    Construct a cutadapt command for trimming FASTQs.
    """

    def __init__(
        self,
        input_fastq_pair: FastqPair,
        output_fastq_pair: FastqPair,
        adapter_type: str,
        paired: bool = True,
        min_length: int = 50,
        two_colour: bool = True,
        polyA: bool = False,
        quality_trim: int | None = None,
    ):
        try:
            adapters: AdapterPair = AdapterPairs[adapter_type].value
        except KeyError:
            raise ValueError(f'Invalid adapter type: {adapter_type}')

        self.command = [
            'cutadapt',
            *('-o', str(output_fastq_pair.r1)),
            *('-a', adapters.r1.sequence),
        ]
        if paired:
            self.command.extend(
                [
                    *('-p', str(output_fastq_pair.r2)),
                    *('-A', adapters.r2.sequence),
                ],
            )

        # FIX: Check for None to allow 0 as a valid value
        if quality_trim is not None:
            if two_colour:
                self.command.append(f'--nextseq-trim={quality_trim}')
            else:
                self.command.append(f'-q {quality_trim}')

        if min_length:
            self.command.append(f'--minimum-length={min_length}')
        if polyA:
            self.command.append('--poly-a')

        self.command.append(str(input_fastq_pair.r1))
        if paired:
            self.command.append(str(input_fastq_pair.r2))

    def __str__(self) -> str:
        return ' '.join(self.command)

    def __repr__(self) -> str:
        return str(self)


class Fastp:
    """
    Construct a fastp command for trimming FASTQs.
    """

    def __init__(
        self,
        input_fastq_pair: FastqPair,
        output_fastq_pair: FastqPair,
        adapter_type: str,
        paired: bool = True,
        min_length: int = 50,
        nthreads: int = 3,
        polyG: bool = True,
        polyX: bool = False,
    ):
        try:
            adapters: AdapterPair = AdapterPairs[adapter_type].value
        except KeyError:
            raise ValueError(f'Invalid adapter type: {adapter_type}')

        self.command = [
            'fastp',
            *('--in1', str(input_fastq_pair.r1)),
            *('--out1', str(output_fastq_pair.r1)),
            *('--length_required', str(min_length)),
            *('--adapter_sequence', adapters.r1.sequence),
            *('--thread', str(nthreads)),
        ]
        if paired:
            self.command.extend(
                [
                    *('--in2', str(input_fastq_pair.r2)),
                    *('--out2', str(output_fastq_pair.r2)),
                    *('--adapter_sequence_r2', adapters.r2.sequence),
                ],
            )
        if not polyG:
            self.command.append('--disable_trim_poly_g')
        if polyX:
            self.command.append('--trim_poly_x')

    def __str__(self) -> str:
        return ' '.join(self.command)

    def __repr__(self) -> str:
        return str(self)


def trim(
    b: Batch,
    sequencing_group: SequencingGroup,
    input_fq_pair: FastqPair,
    output_fq_pair: FastqPair | None = None,
    job_attrs: dict | None = None,
    extra_label: str | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> tuple[Job | None, FastqPair]:
    """
    Takes an input FastqPair object, and creates a job to trim the FASTQs using fastp.
    """
    # Validate inputs
    if not input_fq_pair.r1 or not input_fq_pair.r2:
        raise MissingFastqInputException("Both R1 and R2 FASTQ files must be provided")

    # Don't run if all output files exist and can be reused
    if output_fq_pair and can_reuse(output_fq_pair.r1, overwrite) and can_reuse(output_fq_pair.r2, overwrite):
        return None, output_fq_pair.as_resources(b)

    base_job_name = 'TrimFastqs'
    if extra_label:
        base_job_name += f' {extra_label}'

    if not config.config_retrieve(['workflow']['sequencing_type']) == 'transcriptome':
        raise InvalidSequencingTypeException(
            f"Invalid sequencing type '{config.config_retrieve(['workflow']['sequencing_type'])}'"
             f" for job type '{base_job_name}'; sequencing type must be 'transcriptome'",
        )

    # Get configuration with defaults
    trim_config = config.config_retrieve('trim', {})
    adapter_type = trim_config.get('adapter_type')
    if not adapter_type:
        raise ValueError('No adapter type specified in config file')

    min_length = trim_config.get('min_length', DEFAULT_MIN_LENGTH)
    storage_gb = trim_config.get('storage_gb', DEFAULT_STORAGE_GB)
    quality_trim = trim_config.get('quality_trim', DEFAULT_QUALITY_TRIM)

    trim_tool = trim_config.get('tool', 'fastp')

    trim_j_name = base_job_name
    trim_j_attrs = (job_attrs or {}) | dict(label=base_job_name, tool=trim_tool)
    trim_j = b.new_job(trim_j_name, trim_j_attrs)
    trim_j.image(image_path(trim_tool))

    # Set resource requirements with configurable values
    nthreads = requested_nthreads or trim_config.get('nthreads', DEFAULT_NTHREADS)
    res = STANDARD.set_resources(
        trim_j,
        ncpu=nthreads,
        storage_gb=storage_gb,
    )

    fastq_pair = input_fq_pair.as_resources(b)

    trim_j.declare_resource_group(output_r1={'fastq.gz': '{root}.fastq.gz'})
    trim_j.declare_resource_group(output_r2={'fastq.gz': '{root}.fastq.gz'})
    assert isinstance(trim_j.output_r1, ResourceGroup)
    assert isinstance(trim_j.output_r2, ResourceGroup)
    out_fqs = FastqPair(
        r1=trim_j.output_r1['fastq.gz'],
        r2=trim_j.output_r2['fastq.gz'],
    )

    # Create appropriate trimming command based on tool selection
    if trim_tool == 'cutadapt':
        trim_cmd = Cutadapt(
            input_fastq_pair=fastq_pair,
            output_fastq_pair=out_fqs,
            adapter_type=adapter_type,
            paired=True,
            min_length=min_length,
            quality_trim=quality_trim,
            two_colour=trim_config.get('two_colour', True),
            polyA=trim_config.get('polyA', False),
        )
    else:  # Default to fastp
        trim_cmd = Fastp(
            input_fastq_pair=fastq_pair,
            output_fastq_pair=out_fqs,
            adapter_type=adapter_type,
            paired=True,
            min_length=min_length,
            nthreads=res.get_nthreads(), # uses actual threads allocated
            polyG=trim_config.get('polyG', True),
            polyX=trim_config.get('polyX', False),
        )
    trim_j.command(command(str(trim_cmd), monitor_space=True))

    # Write output to file
    if output_fq_pair:
        b.write_output(out_fqs.r1, str(output_fq_pair.r1))
        b.write_output(out_fqs.r2, str(output_fq_pair.r2))

    return trim_j, out_fqs
