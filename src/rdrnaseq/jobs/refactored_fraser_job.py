import hailtop.batch as hb
from cpg_flow.resources import HIGHMEM, STANDARD
from cpg_utils.config import config_retrieve, get_config
from cpg_utils.hail_batch import command, get_batch
from hailtop.batch.job import Job, BashJob

from rdrnaseq.jobs.bam_to_cram import cram_to_bam


# Assuming scripts are in the same directory or available in the container
R_INIT = 'fraser_init.R'
R_COUNT_SPLIT = 'fraser_count_split.R'
R_MERGE_SPLIT = 'fraser_merge_split.R'
R_COUNT_NON_SPLIT = 'fraser_count_non_split.R'
R_MERGE_NON_SPLIT = 'fraser_merge_non_split.R'
R_ANALYSIS = 'fraser_analysis.R'


def fraser_storage_required_gb(num_bams: int, base_storage_gb: int, per_bam_storage_gb: int) -> int:
    return base_storage_gb + num_bams * per_bam_storage_gb


def get_fraser_job(batch: hb.Batch, name: str, job_attrs: dict) -> Job:
    """Helper to create a standard FRASER job with the right image and resources."""
    j = batch.new_job(name, attributes=job_attrs | {'tool': 'fraser'})
    j.image(config_retrieve(['images', 'fraser']))
    return j


def fraser_init(input_bams: dict[str, hb.ResourceFile], cohort_id: str, job_attrs: dict) -> tuple[Job, hb.ResourceFile]:
    b = get_batch()

    j = get_fraser_job(b, 'fraser_init', job_attrs)

    # Set resource requirements
    storage_required_gb = fraser_storage_required_gb(
        len(input_bams),
        config_retrieve(['workflow', 'fraser_init_storage'], 50),
        config_retrieve(['workflow', 'fraser_bam_single_storage_req'], 10),
    )

    res = HIGHMEM.set_resources(
        j=j,
        ncpu=config_retrieve(['workflow', 'fraser_init_cpu'], 10),
        storage_gb=storage_required_gb,
    )

    bam_files_r_str = ''
    sample_ids_r_str = ''
    sample_ids = []
    for sample_id, bam_file in input_bams.items():
        bam_files_r_str += f'"{bam_file}", '
        sample_ids_r_str += f'"{sample_id}", '
        sample_ids.append(sample_id)
    bam_files_r_str = bam_files_r_str[:-2]  # Remove trailing comma and space
    sample_ids_r_str = sample_ids_r_str[:-2]  # Remove trailing comma and space

    cmd = f"""
    Rscript {R_INIT} \\
        --cohort_id "{cohort_id}" \\
        --sample_ids "{sample_ids_r_str}" \\
        --bam_files "{bam_files_r_str}" \\
        --nthreads {res.get_nthreads()}

    mv output/savedObjects/FRASER_{cohort_id}/fds-object.RDS {j.fds}
    """
    j.command(command(cmd, monitor_space=True))
    return j, j.fds


def fraser_count_split_reads(
    fds: hb.ResourceFile,
    sample_id: str,
    bam: hb.ResourceFile,
    cohort_id: str,
    job_attrs: dict,
) -> tuple[Job, hb.ResourceFile]:
    b = get_batch()
    j = get_fraser_job(b, f'count_split_{sample_id}', job_attrs)

    res = STANDARD.set_resources(
        j=j,
        ncpu=config_retrieve(['workflow', 'fraser_count_job_cpu'], 4),
        storage_gb=20,
    )

    cmd = f"""
    Rscript {R_COUNT_SPLIT} \\
        --fds_path {fds} \\
        --cohort_id "{cohort_id}" \\
        --sample_id "{sample_id}" \\
        --nthreads "{res.get_nthreads()}"

    mv output/cache/splitCounts/splitCounts-{sample_id}.RDS {j.out}
    """
    j.command(command(cmd))
    return j, j.out


def fraser_merge_split_reads(
    fds: hb.ResourceFile, split_counts: dict[str, hb.ResourceFile], cohort_id: str, job_attrs: dict
) -> tuple[Job, hb.ResourceGroup]:
    """3. Merge split reads and extract splice site coordinates for non-split counting."""
    b = get_batch()
    j = get_fraser_job(b, 'fraser_merge_split', job_attrs)

    # Prepare resource group for multiple outputs
    j.BashJob.declare_resource_group(
        out={
            'g_ranges_split': 'g_ranges_split_counts.RDS',
            'g_ranges_non_split': 'g_ranges_non_split_counts.RDS',
            'splice_site_coords': 'splice_site_coords.RDS',
        }
    )

    # Link split counts into the expected cache directory
    # Adjusted to match FRASER's internal directory requirements
    cache_path = f'output/savedObjects/{cohort_id}/cache/splitCounts'

    storage_required_gb = fraser_storage_required_gb(
        len(split_counts),
        config_retrieve(['workflow', 'fraser_init_storage'], 50),
        config_retrieve(['workflow', 'fraser_bam_single_storage_req'], 10),
    )

    res = HIGHMEM.set_resources(
        j=j,
        ncpu=config_retrieve(['workflow', 'fraser_init_cpu'], 10),
        storage_gb=storage_required_gb,
    )

    setup_cache = '\n'.join(
        [
            f'mkdir -p {cache_path} && ln -s {res} {cache_path}/splitCounts-{sid}.RDS'
            for sid, res in split_counts.items()
        ]
    )
    cmd = f"""
    {setup_cache}

    Rscript {R_MERGE_SPLIT} \\
        --fds_path {fds} \\
        --cohort_id "{cohort_id}" \\
        --nthreads "{res.get_nthreads()}"

    mv g_ranges_split_counts.RDS {j.out.g_ranges_split}
    mv g_ranges_non_split_counts.RDS {j.out.g_ranges_non_split}
    mv splice_site_coords.RDS {j.out.splice_site_coords}
    """
    j.command(command(cmd))
    return j, j.out


def fraser_count_non_split_reads(
    fds: hb.ResourceFile, coords: hb.ResourceFile, sample_id: str, cohort_id: str, job_attrs: dict
) -> tuple[Job, hb.ResourceFile]:
    """4. Count non-split reads for a single sample using the merged splice site coordinates."""
    b = get_batch()
    j = get_fraser_job(b, f'count_non_split_{sample_id}', job_attrs)
    res = STANDARD.set_resources(
        j=j,
        ncpu=config_retrieve(['workflow', 'fraser_count_job_cpu'], 4),
        storage_gb=20,
    )

    cmd = f"""
    Rscript {R_COUNT_NON_SPLIT} \\
        --fds_path {fds} \\
        --cohort_id "{cohort_id}" \\
        --sample_id "{sample_id}" \\
        --coords_path {coords} \\
        --nthreads "{res.get_nthreads()}"

    # Note: FRASER 2.0 often produces .h5, fallback to .RDS if needed
    if [ -f "output/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sample_id}.h5" ]; then
        mv output/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sample_id}.h5 {j.out}
    else
        mv output/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sample_id}.RDS {j.out}
    fi
    """
    j.command(command(cmd))
    return j, j.out


def fraser_merge_non_split_reads(
    fds: hb.ResourceFile,
    non_split_counts: dict[str, hb.ResourceFile],
    filtered_ranges: hb.ResourceFile,
    cohort_id: str,
    job_attrs: dict,
) -> tuple[Job, hb.ResourceFile]:
    """5. Final merge of non-split counts and calculation of PSI/Jaccard values."""
    b = get_batch()
    j = get_fraser_job(b, 'fraser_merge_non_split', job_attrs)

    storage_required_gb = fraser_storage_required_gb(
        len(non_split_counts),
        config_retrieve(['workflow', 'fraser_init_storage'], 50),
        config_retrieve(['workflow', 'fraser_bam_single_storage_req'], 10),
    )

    res = HIGHMEM.set_resources(
        j=j,
        ncpu=config_retrieve(['workflow', 'fraser_merge_cpu'], 10),
        storage_gb=storage_required_gb,
    )

    setup_cache = '\n'.join(
        [
            f'mkdir -p output/cache/nonSplicedCounts/{cohort_id} && ln -s '
            f'{res} output/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sid}.h5'
            for sid, res in non_split_counts.items()
        ]
    )

    cmd = f"""
    {setup_cache}

    Rscript {R_MERGE_NON_SPLIT} \\
        --fds_path {fds} \\
        --cohort_id "{cohort_id}" \\
        --filtered_ranges_path {filtered_ranges} \\
        --nthreads "{res.get_nthreads()}"

    tar -czvf {j.fds_tar} output/savedObjects/FRASER_{cohort_id}/
    """
    j.command(command(cmd))
    return j, j.fds_tar


def fraser_analysis(fds_tar: hb.ResourceFile, cohort_id: str, job_attrs: dict) -> Job:
    """6. Statistical analysis, hyperparameter optimization, and result generation."""
    b = get_batch()
    j = get_fraser_job(b, f'fraser_analysis_{cohort_id}', job_attrs)

    storage_required_gb = fraser_storage_required_gb(
        len(fds_tar),
        config_retrieve(['workflow', 'fraser_init_storage'], 50),
        config_retrieve(['workflow', 'fraser_bam_single_storage_req'], 10),
    )

    res = HIGHMEM.set_resources(
        j=j,
        ncpu=config_retrieve(['workflow', 'fraser_main_cpu'], 10),
        storage_gb=storage_required_gb,
    )

    j.declare_resource_group(
        output={'sig_results': 'results.significant.csv', 'all_results': 'results.all.csv', 'plots': 'plots.tar.gz'}
    )

    cfg = get_config().get('fraser', {})
    pval = cfg.get('pval_cutoff', 0.05)
    dpsi = cfg.get('delta_psi_cutoff', 0.3)

    cmd = f"""
    tar -xvf {fds_tar}

    Rscript {R_ANALYSIS} \\
        --fds_dir "output" \\
        --cohort_id "{cohort_id}" \\
        --pval_cutoff {pval} \\
        --delta_psi_cutoff {dpsi} \\
        --nthreads {res.get_nthreads()}

    tar -czvf {j.output.plots} *.png
    cp results.significant.csv {j.output.sig_results}
    cp results.all.csv {j.output.all_results}
    """
    j.command(command(cmd))
    return j
