"""Tests for the AnnotateCohort MT subset flow.

Tests the three-layer chain:
  Stage (VariantSpliceMatch.queue_jobs) → Job (subset_cohort_mt) → Script (subset_mt.py)

Script tests use a real local Hail backend. Job and stage tests mock Hail Batch
and Metamist but verify correct wiring and arguments.
"""

import sys
from unittest.mock import MagicMock, patch

import pytest

import hail as hl

from rdrnaseq.scripts.dashboard_utilities import build_rna_to_genome_map

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope='session', autouse=True)
def hail_init():
    hl.init(backend='local', quiet=True)


@pytest.fixture()
def test_mt(tmp_path):
    """Create a minimal MT with 4 genome sample columns and 5 rows."""
    mt_path = str(tmp_path / 'input.mt')
    sample_ids = ['CPGAAA', 'CPGBBB', 'CPGCCC', 'CPGDDD']
    mt = hl.utils.range_matrix_table(n_rows=5, n_cols=len(sample_ids))
    mt = mt.annotate_cols(s=hl.literal(sample_ids)[mt.col_idx])
    mt = mt.key_cols_by('s')
    mt.write(mt_path)
    return mt_path, sample_ids


@pytest.fixture()
def sg_ids_file(tmp_path):
    """Write a newline-delimited SG IDs file and return its path."""

    def _write(ids):
        path = tmp_path / 'sg_ids.txt'
        path.write_text('\n'.join(ids) + '\n')
        return str(path)

    return _write


FAKE_PEDIGREE_RESULT = {
    'project': {
        'sequencingGroups': [
            {
                'id': 'CPGRNA1',
                'sample': {
                    'participant': {
                        'externalId': 'P001',
                        'families': [{'externalId': 'FAM001'}],
                        'familyParticipants': [{'affected': 2}],
                        'samples': [
                            {'sequencingGroups': [{'id': 'CPGAAA', 'technology': 'short-read', 'type': 'genome'}]},
                        ],
                    },
                },
            },
            {
                'id': 'CPGRNA2',
                'sample': {
                    'participant': {
                        'externalId': 'P002',
                        'families': [{'externalId': 'FAM002'}],
                        'familyParticipants': [{'affected': 1}],
                        'samples': [
                            {'sequencingGroups': [{'id': 'CPGCCC', 'technology': 'short-read', 'type': 'genome'}]},
                        ],
                    },
                },
            },
        ],
    },
}


# ---------------------------------------------------------------------------
# Script tests: subset_mt.py with real local Hail
# ---------------------------------------------------------------------------


class TestSubsetMTScript:
    def _run_subset(self, mt_path, sg_file, output_path):
        with (
            patch('rdrnaseq.scripts.subset_mt.init_batch'),
            patch.object(sys, 'argv', ['subset_mt', '--input', mt_path, '--sgs', sg_file, '--output', output_path]),
        ):
            from rdrnaseq.scripts.subset_mt import main

            main()

    def test_keeps_requested_samples(self, test_mt, sg_ids_file, tmp_path):
        mt_path, _ = test_mt
        keep = ['CPGAAA', 'CPGCCC']
        output = str(tmp_path / 'subset.mt')

        self._run_subset(mt_path, sg_ids_file(keep), output)

        result = hl.read_matrix_table(output)
        assert result.count_cols() == 2
        assert set(result.s.collect()) == set(keep)

    def test_preserves_row_count(self, test_mt, sg_ids_file, tmp_path):
        mt_path, _ = test_mt
        output = str(tmp_path / 'one_sample.mt')

        self._run_subset(mt_path, sg_ids_file(['CPGBBB']), output)

        result = hl.read_matrix_table(output)
        assert result.count_rows() == 5
        assert result.count_cols() == 1

    def test_no_matching_ids_gives_zero_cols(self, test_mt, sg_ids_file, tmp_path):
        mt_path, _ = test_mt
        output = str(tmp_path / 'empty.mt')

        self._run_subset(mt_path, sg_ids_file(['CPGZZZ']), output)

        result = hl.read_matrix_table(output)
        assert result.count_cols() == 0

    def test_all_samples_kept(self, test_mt, sg_ids_file, tmp_path):
        mt_path, all_samples = test_mt
        output = str(tmp_path / 'all.mt')

        self._run_subset(mt_path, sg_ids_file(all_samples), output)

        result = hl.read_matrix_table(output)
        assert result.count_cols() == len(all_samples)
        assert set(result.s.collect()) == set(all_samples)

    def test_blank_lines_in_sg_file_ignored(self, test_mt, tmp_path):
        mt_path, _ = test_mt
        sg_file = tmp_path / 'sg_ids.txt'
        sg_file.write_text('\nCPGAAA\n\n\nCPGDDD\n\n')
        output = str(tmp_path / 'blanks.mt')

        self._run_subset(mt_path, str(sg_file), output)

        result = hl.read_matrix_table(output)
        assert set(result.s.collect()) == {'CPGAAA', 'CPGDDD'}


# ---------------------------------------------------------------------------
# Job tests: subset_cohort_mt and match_variants_and_splicing
# ---------------------------------------------------------------------------


class TestSubsetCohortMTJob:
    @patch('rdrnaseq.jobs.variant_splice_match.command', side_effect=lambda x: x)
    @patch('rdrnaseq.jobs.variant_splice_match.config')
    @patch('rdrnaseq.jobs.variant_splice_match.get_batch')
    def test_creates_job_with_correct_command(self, mock_get_batch, mock_config, _mock_cmd):
        mock_batch = MagicMock()
        mock_job = MagicMock()
        mock_batch.new_job.return_value = mock_job
        mock_batch.read_input.return_value = '/local/sg_ids.txt'
        mock_get_batch.return_value = mock_batch
        mock_config.config_retrieve.return_value = {'driver_image': 'test-image:latest'}

        from rdrnaseq.jobs.variant_splice_match import subset_cohort_mt

        j = subset_cohort_mt(
            source_mt_path='gs://bucket/input.mt',
            sg_id_file='gs://bucket/sg_ids.txt',
            output_mt_path='gs://bucket/output.mt',
            job_attrs={'cohort': 'COH123'},
        )

        assert j is mock_job
        mock_job.image.assert_called_once_with('test-image:latest')
        mock_job.storage.assert_called_once_with('10Gi')
        cmd = mock_job.command.call_args[0][0]
        assert 'rdrnaseq.scripts.subset_mt' in cmd
        assert 'gs://bucket/input.mt' in cmd
        assert 'gs://bucket/output.mt' in cmd

    @patch('rdrnaseq.jobs.variant_splice_match.command', side_effect=lambda x: x)
    @patch('rdrnaseq.jobs.variant_splice_match.config')
    @patch('rdrnaseq.jobs.variant_splice_match.get_batch')
    def test_localises_sg_id_file(self, mock_get_batch, mock_config, _mock_cmd):
        mock_batch = MagicMock()
        mock_batch.new_job.return_value = MagicMock()
        mock_batch.read_input.return_value = '/local/sg_ids.txt'
        mock_get_batch.return_value = mock_batch
        mock_config.config_retrieve.return_value = {'driver_image': 'img'}

        from rdrnaseq.jobs.variant_splice_match import subset_cohort_mt

        subset_cohort_mt('gs://in.mt', 'gs://bucket/ids.txt', 'gs://out.mt', {})

        mock_batch.read_input.assert_called_once_with('gs://bucket/ids.txt')
        cmd = mock_batch.new_job.return_value.command.call_args[0][0]
        assert '/local/sg_ids.txt' in cmd


class TestMatchVariantsAndSplicingJob:
    @patch('rdrnaseq.jobs.variant_splice_match.command', side_effect=lambda x: x)
    @patch('rdrnaseq.jobs.variant_splice_match.config')
    @patch('rdrnaseq.jobs.variant_splice_match.get_batch')
    def test_creates_job_pair_per_dataset(self, mock_get_batch, mock_config, _mock_cmd):
        mock_batch = MagicMock()
        mock_job = MagicMock()
        mock_reg_job = MagicMock()
        mock_batch.new_job.return_value = mock_job
        mock_batch.new_python_job.return_value = mock_reg_job
        mock_batch.read_input.return_value = '/local/fraser.csv'
        mock_get_batch.return_value = mock_batch
        mock_config.config_retrieve.return_value = {'driver_image': 'img'}

        from rdrnaseq.jobs.variant_splice_match import match_variants_and_splicing

        subset_job = MagicMock()
        jobs = match_variants_and_splicing(
            mt_path='gs://subset.mt',
            subset_job=subset_job,
            fraser_csv='gs://fraser.csv',
            sg_ids_by_dataset={'proj-a': ['CPGRNA1'], 'proj-b': ['CPGRNA2']},
            output_by_dataset={
                'proj-a': 'gs://out/a/voi',
                'coarse_proj-a': 'gs://tmp/a/coarse',
                'proj-b': 'gs://out/b/voi',
                'coarse_proj-b': 'gs://tmp/b/coarse',
            },
            cohort_id='COH123',
            job_attrs={},
            sequencing_type='transcriptome',
            cell_library_type='blood_polyRNA',
        )

        # 2 datasets × (1 match + 1 registration) = 4 jobs
        assert len(jobs) == 4
        # Match jobs depend on subset_job
        mock_job.depends_on.assert_called_with(subset_job)

    @patch('rdrnaseq.jobs.variant_splice_match.command', side_effect=lambda x: x)
    @patch('rdrnaseq.jobs.variant_splice_match.config')
    @patch('rdrnaseq.jobs.variant_splice_match.get_batch')
    def test_command_passes_rna_ids_and_mt(self, mock_get_batch, mock_config, _mock_cmd):
        mock_batch = MagicMock()
        mock_job = MagicMock()
        mock_batch.new_job.return_value = mock_job
        mock_batch.new_python_job.return_value = MagicMock()
        mock_batch.read_input.return_value = '/local/fraser.csv'
        mock_get_batch.return_value = mock_batch
        mock_config.config_retrieve.return_value = {'driver_image': 'img'}

        from rdrnaseq.jobs.variant_splice_match import match_variants_and_splicing

        match_variants_and_splicing(
            mt_path='gs://subset.mt',
            subset_job=MagicMock(),
            fraser_csv='gs://fraser.csv',
            sg_ids_by_dataset={'proj-a': ['CPGRNA1', 'CPGRNA2']},
            output_by_dataset={'proj-a': 'gs://out/voi', 'coarse_proj-a': 'gs://tmp/coarse'},
            cohort_id='COH123',
            job_attrs={},
            sequencing_type='transcriptome',
            cell_library_type='blood_polyRNA',
        )

        cmd = mock_job.command.call_args[0][0]
        assert 'gs://subset.mt' in cmd
        assert 'CPGRNA1 CPGRNA2' in cmd
        assert 'variant_of_interest_subset' in cmd

    @patch('rdrnaseq.jobs.variant_splice_match.command', side_effect=lambda x: x)
    @patch('rdrnaseq.jobs.variant_splice_match.config')
    @patch('rdrnaseq.jobs.variant_splice_match.get_batch')
    def test_writes_output_tsv_and_bed(self, mock_get_batch, mock_config, _mock_cmd):
        mock_batch = MagicMock()
        mock_job = MagicMock()
        mock_batch.new_job.return_value = mock_job
        mock_batch.new_python_job.return_value = MagicMock()
        mock_batch.read_input.return_value = '/local/fraser.csv'
        mock_get_batch.return_value = mock_batch
        mock_config.config_retrieve.return_value = {'driver_image': 'img'}

        from rdrnaseq.jobs.variant_splice_match import match_variants_and_splicing

        match_variants_and_splicing(
            mt_path='gs://subset.mt',
            subset_job=MagicMock(),
            fraser_csv='gs://fraser.csv',
            sg_ids_by_dataset={'proj-a': ['CPGRNA1']},
            output_by_dataset={'proj-a': 'gs://out/voi', 'coarse_proj-a': 'gs://tmp/coarse'},
            cohort_id='COH123',
            job_attrs={},
            sequencing_type='transcriptome',
            cell_library_type='blood_polyRNA',
        )

        write_calls = {c[0][1] for c in mock_batch.write_output.call_args_list}
        assert 'gs://out/voi.tsv' in write_calls
        assert 'gs://out/voi.bed' in write_calls


# ---------------------------------------------------------------------------
# Integration: pedigree mapping → SG ID file → subset MT (real Hail)
# ---------------------------------------------------------------------------


class TestStageToScriptIntegration:
    """End-to-end: simulate what the stage does, then actually run the subset script."""

    def test_pedigree_mapping_drives_subset(self, test_mt, tmp_path):
        mt_path, all_samples = test_mt

        # Stage step 1: parse pedigree → genome SG IDs
        rna_to_genome = build_rna_to_genome_map(FAKE_PEDIGREE_RESULT)
        assert rna_to_genome == {'CPGRNA1': {'CPGAAA'}, 'CPGRNA2': {'CPGCCC'}}

        all_genome_ids = sorted(set().union(*rna_to_genome.values()))
        assert all_genome_ids == ['CPGAAA', 'CPGCCC']

        # Stage step 2: write SG ID file
        sg_file = tmp_path / 'genome_sg_ids.txt'
        sg_file.write_text('\n'.join(all_genome_ids) + '\n')

        # Script step: subset MT
        output_path = str(tmp_path / 'subsetted.mt')
        with (
            patch('rdrnaseq.scripts.subset_mt.init_batch'),
            patch.object(
                sys, 'argv', ['subset_mt', '--input', mt_path, '--sgs', str(sg_file), '--output', output_path]
            ),
        ):
            from rdrnaseq.scripts.subset_mt import main

            main()

        result = hl.read_matrix_table(output_path)
        assert result.count_cols() == 2
        assert set(result.s.collect()) == {'CPGAAA', 'CPGCCC'}
        assert result.count_rows() == 5

    def test_empty_pedigree_gives_no_genome_ids(self):
        empty_result = {'project': {'sequencingGroups': []}}
        rna_to_genome = build_rna_to_genome_map(empty_result)
        assert rna_to_genome == {}


# ---------------------------------------------------------------------------
# resolve_annotate_cohort_mt tests
# ---------------------------------------------------------------------------


class TestResolveAnnotateCohortMT:
    @patch('rdrnaseq.jobs.variant_splice_match.config')
    def test_returns_config_override_when_set(self, mock_config):
        mock_config.config_retrieve.side_effect = lambda keys, **kw: (
            'gs://override/my.mt' if keys == ['variant_splice_match', 'mt_path'] else kw.get('default')
        )

        from rdrnaseq.jobs.variant_splice_match import resolve_annotate_cohort_mt

        assert resolve_annotate_cohort_mt('test-project') == 'gs://override/my.mt'

    @patch('rdrnaseq.jobs.variant_splice_match.query_for_latest_analysis')
    @patch('rdrnaseq.jobs.variant_splice_match.config')
    def test_queries_metamist_when_no_config(self, mock_config, mock_query):
        def _config_retrieve(keys, *args, **kw):
            if keys == ['variant_splice_match', 'mt_path']:
                return kw.get('default', args[0] if args else None)
            return False

        mock_config.config_retrieve.side_effect = _config_retrieve
        mock_query.return_value = 'gs://queried/annotate_cohort.mt'

        from rdrnaseq.jobs.variant_splice_match import resolve_annotate_cohort_mt

        result = resolve_annotate_cohort_mt('test-project')
        assert result == 'gs://queried/annotate_cohort.mt'
        mock_query.assert_called_once_with(
            dataset='seqr',
            analysis_type='matrixtable',
            sequencing_type='genome',
            long_read=False,
            stage_name='AnnotateCohort',
        )

    @patch('rdrnaseq.jobs.variant_splice_match.query_for_latest_analysis')
    @patch('rdrnaseq.jobs.variant_splice_match.config')
    def test_raises_when_no_mt_found(self, mock_config, mock_query):
        mock_config.config_retrieve.side_effect = lambda keys, *args, **kw: kw.get('default', args[0] if args else None)
        mock_query.return_value = None

        from rdrnaseq.jobs.variant_splice_match import resolve_annotate_cohort_mt

        with pytest.raises(RuntimeError, match='No AnnotateCohort MatrixTable found'):
            resolve_annotate_cohort_mt('test-project')
