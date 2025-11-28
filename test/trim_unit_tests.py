import pytest
from unittest.mock import MagicMock, patch
from dataclasses import dataclass

# Import the classes and functions from your script
# Replace 'trim_script' with the actual name of your file
from rdrnaseq.jobs.trim import (
    Cutadapt,
    Fastp, 
    trim, 
    AdapterPairs, 
    InvalidSequencingTypeException,
    MissingFastqInputException
)

# --- Mocks for Data Classes ---
# We mock these to avoid needing the full cpg_workflows dependency tree for unit tests

@dataclass
class MockResource:
    path: str
    def __str__(self): return self.path

@dataclass
class MockFastqPair:
    r1: MockResource
    r2: MockResource

@pytest.fixture
def input_pair():
    return MockFastqPair(
        r1=MockResource("input_R1.fq.gz"), 
        r2=MockResource("input_R2.fq.gz")
    )

@pytest.fixture
def output_pair():
    return MockFastqPair(
        r1=MockResource("output_R1.fq.gz"), 
        r2=MockResource("output_R2.fq.gz")
    )

# --- Tests for Helper Classes ---

def test_adapter_enum():
    """Ensure Enum values are accessible and correct."""
    truseq = AdapterPairs['ILLUMINA_TRUSEQ'].value
    assert "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" == truseq.r1.sequence
    assert "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" == truseq.r2.sequence

# --- Tests for Command Generation (Cutadapt) ---

def test_cutadapt_command_generation_defaults(input_pair, output_pair):
    """Test Cutadapt command with standard defaults."""
    cmd_obj = Cutadapt(
        input_fastq_pair=input_pair,
        output_fastq_pair=output_pair,
        adapter_type='ILLUMINA_TRUSEQ',
        paired=True,
        min_length=50,
        quality_trim=20
    )
    cmd = str(cmd_obj)
    
    # Check basic structure
    assert "cutadapt" in cmd
    assert "-o output_R1.fq.gz" in cmd
    assert "-p output_R2.fq.gz" in cmd
    # Check adapters
    assert "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" in cmd
    assert "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" in cmd
    # Check flags
    assert "--nextseq-trim=20" in cmd  # Default two_colour=True
    assert "--minimum-length=50" in cmd
    assert "input_R1.fq.gz input_R2.fq.gz" in cmd

def test_cutadapt_custom_flags(input_pair, output_pair):
    """Test Cutadapt with specific flags toggled."""
    cmd_obj = Cutadapt(
        input_fastq_pair=input_pair,
        output_fastq_pair=output_pair,
        adapter_type='ILLUMINA_NEXTERA',
        two_colour=False, # Should use -q instead of --nextseq-trim
        polyA=True,
        quality_trim=30
    )
    cmd = str(cmd_obj)
    
    assert "-q 30" in cmd
    assert "--nextseq-trim" not in cmd
    assert "--poly-a" in cmd
    assert "CTGTCTCTTATACACATCT" in cmd  # Nextera sequence

def test_cutadapt_invalid_adapter(input_pair, output_pair):
    """Test that invalid adapter keys raise ValueError."""
    with pytest.raises(ValueError) as exc:
        Cutadapt(input_pair, output_pair, adapter_type='INVALID_ADAPTER')
    assert "Invalid adapter type" in str(exc.value)

# --- Tests for Command Generation (Fastp) ---

def test_fastp_command_generation(input_pair, output_pair):
    """Test Fastp command generation."""
    cmd_obj = Fastp(
        input_fastq_pair=input_pair,
        output_fastq_pair=output_pair,
        adapter_type='ILLUMINA_TRUSEQ',
        min_length=75,
        nthreads=4
    )
    cmd = str(cmd_obj)
    
    assert "fastp" in cmd
    assert "--in1 input_R1.fq.gz" in cmd
    assert "--out1 output_R1.fq.gz" in cmd
    assert "--length_required 75" in cmd
    assert "--thread 4" in cmd
    assert "--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" in cmd

def test_fastp_flags(input_pair, output_pair):
    """Test Fastp boolean flags."""
    cmd_obj = Fastp(
        input_fastq_pair=input_pair,
        output_fastq_pair=output_pair,
        adapter_type='ILLUMINA_TRUSEQ',
        polyG=False, # Should add disable flag
        polyX=True   # Should add trim flag
    )
    cmd = str(cmd_obj)
    
    assert "--disable_trim_poly_g" in cmd
    assert "--trim_poly_x" in cmd

# --- Tests for the Main Trim Function (Logic & Config) ---

@patch('rdrnaseq.jobs.trim.get_config')
@patch('rdrnaseq.jobs.trim.image_path')  # Mock image path lookup
def test_trim_workflow_success(mock_image, mock_get_config, input_pair):
    """
    Test the high-level trim function.
    Mocks the CPG Config and Batch object to ensure logic flow is correct.
    """
    # 1. Setup Config Mock
    mock_get_config.return_value = {
        'workflow': {'sequencing_type': 'transcriptome'},
        'trim': {
            'adapter_type': 'ILLUMINA_TRUSEQ',
            'tool': 'fastp',
            'min_length': 50
        }
    }
    mock_image.return_value = "fastp:latest"

    # 2. Setup Batch and Job Mocks
    mock_batch = MagicMock()
    mock_job = MagicMock()
    mock_batch.new_job.return_value = mock_job
    
    # Mock the ResourceGroup output behavior of the job
    mock_job.output_r1 = {'fastq.gz': 'out.r1.fq.gz'}
    mock_job.output_r2 = {'fastq.gz': 'out.r2.fq.gz'}

    # 3. Setup Sequencing Group Mock
    mock_sg = MagicMock()

    # 4. Run the function
    # We must mock input_pair.as_resources because it interacts with the Batch backend
    with patch.object(input_pair, 'as_resources', return_value=input_pair):
        job, out_pair = trim(mock_batch, mock_sg, input_pair)

    # 5. Assertions
    assert job is not None
    mock_batch.new_job.assert_called_with('TrimFastqs', {'label': 'TrimFastqs', 'tool': 'fastp'})
    
    # Verify command was set
    # command() is a CPG wrapper, we check if the job.command was called with a string containing fastp
    args, _ = mock_job.command.call_args
    assert "fastp" in str(args[0])

@patch('rdrnaseq.jobs.trim.get_config')
def test_trim_wrong_sequencing_type(mock_get_config, input_pair):
    """Test that the function raises an error if sequencing type is not transcriptome."""
    mock_get_config.return_value = {
        'workflow': {'sequencing_type': 'genome'}, # Wrong type
        'trim': {'adapter_type': 'ILLUMINA_TRUSEQ'}
    }
    
    mock_batch = MagicMock()
    mock_sg = MagicMock()

    with pytest.raises(InvalidSequencingTypeException):
        trim(mock_batch, mock_sg, input_pair)

def test_trim_missing_inputs():
    """Test that missing R1/R2 inputs raise an exception immediately."""
    bad_input = MockFastqPair(r1=None, r2=None)
    mock_batch = MagicMock()
    mock_sg = MagicMock()

    with pytest.raises(MissingFastqInputException):
        trim(mock_batch, mock_sg, bad_input)