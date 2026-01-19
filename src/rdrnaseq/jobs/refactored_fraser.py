import hailtop.batch as hb
import pandas as pd
from cpg_flow.filetypes import BamPath, CramPath
from cpg_flow.resources import HIGHMEM, STANDARD
from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, get_config, reference_path
from cpg_utils.hail_batch import command, get_batch
from hailtop.batch.job import Job

from rdrnaseq.jobs.bam_to_cram import cram_to_bam

# R Script Paths
R_INIT = 'RDrnaseq/fraser_init.R'
R_COUNT_SPLIT = 'RDrnaseq/fraser_count_split.R'
R_MERGE_SPLIT = 'RDrnaseq/fraser_merge_split.R'
R_COUNT_NON_SPLIT = 'RDrnaseq/fraser_count_non_split.R'
R_MERGE_NON_SPLIT = 'RDrnaseq/fraser_merge_non_split.R'
R_ANALYSIS = 'RDrnaseq/fraser_analysis.R'


def get_fraser_job(batch: hb.Batch, name: str, job_attrs: dict) -> Job:
    """Helper to create a standard FRASER job."""
    j = batch.new_job(name, attributes=job_attrs | {'tool': 'fraser'})
    j.image(config_retrieve(['images', 'fraser']))
    return j


def fraser_storage_required_gb(num_bams: int, base_storage_gb: int, per_bam_storage_gb: int) -> int:
    """Calculates disk space needed based on cohort size."""
    return base_storage_gb + num_bams * per_bam_storage_gb


# --- Orchestration Pipeline ---


def fraser_pipeline(
    input_bams_or_crams: list[tuple[str, BamPath, None] | tuple[str, CramPath, Path]],
    cohort_id: str,
    job_attrs: dict,
    output_fds_path: dict[str, Path],
    output_prefix: Path,
):
    b = get_batch()
    # Use output_prefix for intermediate steps and extract final paths from dict
    root = to_path(output_prefix / cohort_id)
    num_samples = len(input_bams_or_crams)

    # 0. Localize/Convert Inputs
    input_bams_localised: dict[str, hb.ResourceFile] = {}
    for sample_id, input_bam_or_cram, potential_bam_path in input_bams_or_crams:
        if isinstance(input_bam_or_cram, CramPath) and isinstance(potential_bam_path, Path):
            _, output_bam = cram_to_bam(
                input_cram=input_bam_or_cram.resource_group(b),
                output_bam=potential_bam_path,
                job_attrs=job_attrs,
                reference_fasta_path=reference_path('broad/ref_fasta'),
            )
            input_bams_localised[sample_id] = output_bam.bam
        elif isinstance(input_bam_or_cram, BamPath):
            # Localise BAM
            input_bams_localised[sample_id] = input_bam_or_cram.resource_group(b).bam

    # 1. Init
    init_fds_path = root / 'init' / 'fds-object.RDS'
    fds_res = fraser_init(b, input_bams_localised, cohort_id, job_attrs, init_fds_path)

    # 2. Count Split
    split_counts_res = {}
    for sid, localised_bam in input_bams_localised.items():
        out_p = root / 'split_counts' / f'{sid}.RDS'
        # Pass the specific BAM for this sample count
        split_counts_res[sid] = fraser_count_split_reads(b, fds_res, localised_bam, sid, cohort_id, job_attrs, out_p)

    # 3. Merge Split
    merge_split_paths = {
        'g_ranges_split': root / 'merge_split' / 'g_ranges_split_counts.RDS',
        'g_ranges_non_split': root / 'merge_split' / 'g_ranges_non_split_counts.RDS',
        'splice_site_coords': root / 'merge_split' / 'splice_site_coords.RDS',
    }
    merge_out_rg = fraser_merge_split_reads(b, fds_res, split_counts_res, cohort_id, job_attrs, merge_split_paths)

    # 4. Count Non-Split
    non_split_counts_res = {}
    for sid, localised_bam in input_bams_localised.items():
        out_p = root / 'non_split_counts' / f'{sid}.h5'
        non_split_counts_res[sid] = fraser_count_non_split_reads(
            b, fds_res, localised_bam, merge_out_rg.splice_site_coords, sid, cohort_id, job_attrs, out_p
        )

    # --- Step 5: Merge Non-Split ---
    # This fulfills the 'Rds_data' requirement from your Stage
    fds_tar_path = to_path(output_fds_path['Rds_data'])
    fds_tar_res = fraser_merge_non_split_reads(
        b,
        fds_res,
        non_split_counts_res,
        merge_out_rg.g_ranges_non_split,
        cohort_id,
        job_attrs,
        fds_tar_path,
        num_samples,
    )

    # --- Step 6: Analysis ---
    # This fulfills the 'seqr_data' requirement from your Stage
    analysis_paths = {
        'sig_results': root / 'results' / 'results.significant.csv',
        'all_results': to_path(output_fds_path['seqr_data']),
        'plots': root / 'results' / 'plots.tar.gz',
        'stats': root / 'results' / 'statistics_summary.txt',
    }
    fraser_analysis(b, fds_tar_res, cohort_id, job_attrs, analysis_paths, num_samples)


# --- Step Methods ---


def fraser_init(b: hb.Batch, input_bams: dict, cohort_id: str, job_attrs: dict, output_path: Path) -> hb.ResourceFile:
    if output_path.exists():
        return b.read_input(output_path)

    j = get_fraser_job(b, 'fraser_init', job_attrs)
    # Dynamic storage for Init
    storage = fraser_storage_required_gb(len(input_bams), base_storage_gb=50, per_bam_storage_gb=10)
    res = HIGHMEM.set_resources(j=j, ncpu=10, storage_gb=storage)

    df = pd.DataFrame([{'sample_id': s, 'bam_path': str(bam)} for s, bam in input_bams.items()])
    tmp_csv = to_path(output_path).parent / 'sample_map.csv'
    with tmp_csv.open('w') as f:
        df.to_csv(f, index=False)
    sample_map_input = b.read_input(str(tmp_csv))

    j.command(
        command(f"""
        mkdir -p /io/work
        Rscript {R_INIT} --cohort_id "{cohort_id}" --sample_map {sample_map_input} \\
            --work_dir "/io/work" --nthreads {res.get_nthreads()}
        mv /io/work/savedObjects/FRASER_{cohort_id}/fds-object.RDS {j.fds}
        rm -rf /io/work
    """)
    )
    b.write_output(j.fds, str(output_path))
    return j.fds


def fraser_count_split_reads(b, fds, bam, sample_id, cohort_id, job_attrs, output_path) -> hb.ResourceFile:
    if output_path.exists():
        return b.read_input(output_path)

    j = get_fraser_job(b, f'count_split_{sample_id}', job_attrs)
    res = STANDARD.set_resources(j=j, ncpu=4, storage_gb=20)
    j.command(
        command(f"""
        mkdir -p /io/work
        Rscript {R_COUNT_SPLIT} --fds_path {fds} --bam_path {bam} --cohort_id "{cohort_id}" \\
            --sample_id "{sample_id}" --work_dir "/io/work" --nthreads "{res.get_nthreads()}"
        mv /io/work/cache/splitCounts/splitCounts-{sample_id}.RDS {j.out}
        rm -rf /io/work
    """)
    )
    b.write_output(j.out, str(output_path))
    return j.out


def fraser_merge_split_reads(b, fds, split_counts, cohort_id, job_attrs, output_paths) -> hb.ResourceGroup:
    # Check if all outputs exist
    if all(p.exists() for p in output_paths.values()):
        return b.read_input_group(**{k: str(v) for k, v in output_paths.items()})

    j = get_fraser_job(b, 'fraser_merge_split', job_attrs)
    j.declare_resource_group(out={k: v.name for k, v in output_paths.items()})

    storage = fraser_storage_required_gb(len(split_counts), base_storage_gb=50, per_bam_storage_gb=10)
    res = HIGHMEM.set_resources(j=j, ncpu=10, storage_gb=storage)

    cache_path = '/io/work/cache/splitCounts'
    setup_cmds = [f'mkdir -p {cache_path}']
    for sid, r_file in split_counts.items():
        setup_cmds.append(f'ln -s {r_file} {cache_path}/splitCounts-{sid}.RDS')

    j.command(
        command(
            '\n'.join(setup_cmds)
            + f"""
        Rscript {R_MERGE_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --work_dir "/io/work" --nthreads "{res.get_nthreads()}"
        mv /io/work/g_ranges_split_counts.RDS {j.out.g_ranges_split}
        mv /io/work/g_ranges_non_split_counts.RDS {j.out.g_ranges_non_split}
        mv /io/work/splice_site_coords.RDS {j.out.splice_site_coords}
        rm -rf /io/work
    """
        )
    )
    for k, p in output_paths.items():
        b.write_output(j.out[k], str(p))
    return j.out


def fraser_count_non_split_reads(b, fds, bam, coords, sample_id, _cohort_id, job_attrs, output_path) -> hb.ResourceFile:
    if output_path.exists():
        return b.read_input(output_path)

    j = get_fraser_job(b, f'count_non_split_{sample_id}', job_attrs)
    res = STANDARD.set_resources(j=j, ncpu=4, storage_gb=20)
    j.command(
        command(f"""
        mkdir -p /io/work
        Rscript {R_COUNT_NON_SPLIT} --fds_path {fds} --bam_path {bam} --coords_path {coords} \\
            --sample_id "{sample_id}" --work_dir "/io/work" --nthreads "{res.get_nthreads()}"
        TARGET_FILE=$(find /io/work/cache/nonSplicedCounts/ -name "nonSplicedCounts-{sample_id}.*")
        mv "$TARGET_FILE" {j.out}
        rm -rf /io/work
    """)
    )
    b.write_output(j.out, str(output_path))
    return j.out


def fraser_merge_non_split_reads(
    b, fds, non_split_counts, filtered_ranges, cohort_id, job_attrs, output_path, num_samples
) -> hb.ResourceFile:
    if output_path.exists():
        return b.read_input(output_path)

    j = get_fraser_job(b, 'fraser_merge_non_split', job_attrs)
    storage = fraser_storage_required_gb(num_samples, base_storage_gb=50, per_bam_storage_gb=10)
    res = HIGHMEM.set_resources(j=j, ncpu=10, storage_gb=storage)

    setup_lines = [f'mkdir -p /io/work/cache/nonSplicedCounts/{cohort_id}']
    for sid, r in non_split_counts.items():
        setup_lines.append(f'EXT=$(echo {r} | sed "s/.*\\.//")')
        setup_lines.append(f'ln -s {r} /io/work/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sid}.$EXT')

    j.command(
        command(
            '\n'.join(setup_lines)
            + f"""
        Rscript {R_MERGE_NON_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --filtered_ranges_path {filtered_ranges} --work_dir "/io/work" --nthreads "{res.get_nthreads()}"
        tar -czvf {j.fds_tar} -C /io/work/savedObjects/ FRASER_{cohort_id}/
        rm -rf /io/work
    """
        )
    )
    b.write_output(j.fds_tar, str(output_path))
    return j.fds_tar


def fraser_analysis(b, fds_tar, cohort_id, job_attrs, output_paths, num_samples):
    if all(p.exists() for p in output_paths.values()):
        return None

    j = get_fraser_job(b, f'fraser_analysis_{cohort_id}', job_attrs)
    j.declare_resource_group(output={k: v.name for k, v in output_paths.items()})

    storage = fraser_storage_required_gb(num_samples, base_storage_gb=50, per_bam_storage_gb=10)
    res = HIGHMEM.set_resources(j=j, ncpu=10, storage_gb=storage)

    cfg = get_config().get('fraser', {})
    j.command(
        command(f"""
        mkdir -p /io/work
        tar -xvf {fds_tar} -C /io/work/
        Rscript {R_ANALYSIS} --fds_dir "/io/work" --cohort_id "{cohort_id}" \\
            --pval_cutoff {cfg.get('pval_cutoff', 0.05)} --delta_psi_cutoff {cfg.get('delta_psi_cutoff', 0.3)} \\
            --nthreads {res.get_nthreads()}
        tar -czvf {j.output.plots} -C /io/work/ plots/
        cp /io/work/results.significant.csv {j.output.sig_results}
        cp /io/work/results.all.csv {j.output.all_results}
        cp /io/work/statistics_summary.txt {j.output.stats}
        rm -rf /io/work
    """)
    )
    for k, p in output_paths.items():
        b.write_output(j.output[k], str(p))
    return j
