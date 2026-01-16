# ruff : noqa: PLR0912
# ruff : noqa: PLR0915
import hailtop.batch as hb
from cpg_flow.filetypes import (
    BamPath,
    CramPath,
)
from cpg_flow.resources import HIGHMEM, STANDARD
from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, get_config, reference_path
from cpg_utils.hail_batch import command, get_batch
from hailtop.batch.job import Job

from rdrnaseq.jobs.bam_to_cram import cram_to_bam

# Assuming scripts are in the same directory or available in the container
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


def fraser_pipeline(
    input_bams_or_crams: list[tuple[str, BamPath, None] | tuple[str, CramPath, Path]],
    cohort_id: str,
    job_attrs: dict,
    output_prefix: Path,
    output_fds_path: dict[str, hb.ResourceFile],
):
    """
    Main orchestration function with checkpointing.
    Checks if output exists in GCS before creating a job.
    """
    jobs: list[Job] = []
    b = get_batch()
    root = to_path(output_prefix / cohort_id)
    # --- Step 0: Check Bams ---
    input_bams_localised: dict[str, hb.ResourceFile] = {}
    for input_bam_or_cram_tuple in input_bams_or_crams:
        sample_id = input_bam_or_cram_tuple[0]
        input_bam_or_cram = input_bam_or_cram_tuple[1]
        potential_bam_path = input_bam_or_cram_tuple[2]
        if isinstance(input_bam_or_cram, CramPath) and isinstance(potential_bam_path, Path):
            j, output_bam = cram_to_bam(
                input_cram=input_bam_or_cram.resource_group(b),
                output_bam=potential_bam_path,
                job_attrs=job_attrs,
                reference_fasta_path=reference_path('broad/ref_fasta'),
            )
            if j and isinstance(j, Job):
                jobs.append(j)
            input_bams_localised[sample_id] = output_bam.bam
        elif isinstance(input_bam_or_cram, BamPath):
            # Localise BAM
            input_bams_localised[sample_id] = input_bam_or_cram.resource_group(b).bam
    # Use a generator to find the first element that is NOT a ResourceFile
    if any(not isinstance(f, hb.ResourceFile) for f in input_bams_localised.values()):
        raise TypeError('All elements in input_bams_localised must be instances of hb.ResourceFile.')
    input_bams = input_bams_localised
    # --- Step 1: Initialization ---
    init_fds_path = root / 'init' / 'fds-object.RDS'

    if init_fds_path.exists():
        print(f'Skipping Init: Found {init_fds_path}')
        fds_res = b.read_input(str(init_fds_path))
    else:
        j, fds_res = fraser_init(b, input_bams, cohort_id, job_attrs, output_fds_path)
        b.write_output(fds_res, str(init_fds_path))

    # --- Step 2: Count Split Reads (Parallel) ---
    split_counts_res = {}

    for seq_group_id in input_bams:
        sample_out_path = root / 'split_counts' / f'{seq_group_id}.RDS'

        if sample_out_path.exists():
            split_counts_res[seq_group_id] = b.read_input(str(sample_out_path))
        else:
            j, out = fraser_count_split_reads(b, fds_res, seq_group_id, cohort_id, job_attrs)
            b.write_output(out, str(sample_out_path))
            split_counts_res[seq_group_id] = out

    # --- Step 3: Merge Split Reads ---
    # We check one key file to determine if merge is done
    merge_split_paths = {
        'g_ranges_split': root / 'merge_split' / 'g_ranges_split_counts.RDS',
        'g_ranges_non_split': root / 'merge_split' / 'g_ranges_non_split_counts.RDS',
        'splice_site_coords': root / 'merge_split' / 'splice_site_coords.RDS',
    }

    if all(p.exists() for p in merge_split_paths.values()):
        print(f'Skipping Merge Split: Found all outputs in {root / "merge_split"}')
        merge_split_res = {k: b.read_input(str(p)) for k, p in merge_split_paths.items()}
        # ResourceGroup construction for convenience if needed, though we use dict here
        splice_site_coords_res = merge_split_res['splice_site_coords']
        filtered_ranges_res = merge_split_res['g_ranges_non_split']
    else:
        j, out_rg = fraser_merge_split_reads(b, fds_res, split_counts_res, cohort_id, job_attrs)

        # Write outputs
        b.write_output(out_rg.g_ranges_split, str(merge_split_paths['g_ranges_split']))
        b.write_output(out_rg.g_ranges_non_split, str(merge_split_paths['g_ranges_non_split']))
        b.write_output(out_rg.splice_site_coords, str(merge_split_paths['splice_site_coords']))

        splice_site_coords_res = out_rg.splice_site_coords
        filtered_ranges_res = out_rg.g_ranges_non_split

    # --- Step 4: Count Non-Split Reads (Parallel) ---
    non_split_counts_res = {}

    for sample_id in input_bams:
        # Check for both .h5 and .RDS as FRASER might produce either
        h5_path = root / 'non_split_counts' / f'{sample_id}.h5'
        rds_path = root / 'non_split_counts' / f'{sample_id}.RDS'

        if h5_path.exists():
            non_split_counts_res[sample_id] = b.read_input(str(h5_path))
        elif rds_path.exists():
            non_split_counts_res[sample_id] = b.read_input(str(rds_path))
        else:
            j, out = fraser_count_non_split_reads(b, fds_res, splice_site_coords_res, sample_id, cohort_id, job_attrs)
            # Default write to .h5 for consistency
            b.write_output(out, str(h5_path))
            non_split_counts_res[sample_id] = out

    # --- Step 5: Merge Non-Split Reads ---
    fds_tar_path = root / 'merge_non_split' / f'FRASER_{cohort_id}.tar.gz'

    if fds_tar_path.exists():
        print(f'Skipping Merge Non-Split: Found {fds_tar_path}')
        fds_tar_res = b.read_input(str(fds_tar_path))
    else:
        j, fds_tar_res = fraser_merge_non_split_reads(
            b, fds_res, non_split_counts_res, filtered_ranges_res, cohort_id, job_attrs
        )
        b.write_output(fds_tar_res, str(fds_tar_path))

    # --- Step 6: Analysis ---
    results_path = root / 'results'
    sig_csv_path = results_path / 'results.significant.csv'

    if sig_csv_path.exists():
        print(f'Skipping Analysis: Found {sig_csv_path}')
    else:
        j = fraser_analysis(b, fds_tar_res, cohort_id, job_attrs)

        # Write resource group outputs
        b.write_output(j.output.sig_results, str(sig_csv_path))
        b.write_output(j.output.all_results, str(results_path / 'results.all.csv'))
        b.write_output(j.output.plots, str(results_path / 'plots.tar.gz'))
        b.write_output(j.output.stats, str(results_path / 'statistics_summary.txt'))


def fraser_storage_required_gb(num_bams: int, base_storage_gb: int, per_bam_storage_gb: int) -> int:
    return base_storage_gb + num_bams * per_bam_storage_gb


def fraser_init(
    b: hb.Batch,
    input_bams: dict[str, hb.ResourceFile],
    cohort_id: str,
    job_attrs: dict,
    output_fds_path: dict[str, hb.ResourceFile],
) -> tuple[Job, hb.ResourceFile]:
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

    # collect all the IDs in the loop
    bam_ids = []
    sample_ids = []

    for sample_id, bam_file in input_bams.items():
        bam_ids.append(bam_file)
        sample_ids.append(sample_id)

    # generate commaspace-delimited lists of quoted IDs
    bam_files_r_str = ', '.join(f'"{bam}"' for bam in bam_ids)
    sample_ids_r_str = ', '.join(f'"{sam}"' for sam in sample_ids)

    cmd = f"""
    Rscript {R_INIT} \\
        --cohort_id "{cohort_id}" \\
        --sample_ids "{sample_ids_r_str}" \\
        --bam_files "{bam_files_r_str}" \\
        --nthreads {res.get_nthreads()}

    mv output/savedObjects/FRASER_{cohort_id}/fds-object.RDS {j.fds}
    """
    j.command(command(cmd, monitor_space=True))
    b.write_output(j.fds, output_fds_path['temp_data'])

    return j, j.fds


def fraser_count_split_reads(
    b: hb.Batch, fds: hb.ResourceFile, sample_id: str, cohort_id: str, job_attrs: dict
) -> tuple[Job, hb.ResourceFile]:
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
    b: hb.Batch, fds: hb.ResourceFile, split_counts: dict[str, hb.ResourceFile], cohort_id: str, job_attrs: dict
) -> tuple[Job, hb.ResourceGroup]:
    j = get_fraser_job(b, 'fraser_merge_split', job_attrs)

    # Prepare resource group for multiple outputs
    j.declare_resource_group(
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

    for sid, res in split_counts.items():
        j.command(f'mkdir -p {cache_path} && ln -s {res} {cache_path}/splitCounts-{sid}.RDS')

    cmd = f"""
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
    b: hb.Batch, fds: hb.ResourceFile, coords: hb.ResourceFile, sample_id: str, cohort_id: str, job_attrs: dict
) -> tuple[Job, hb.ResourceFile]:
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
    b: hb.Batch,
    fds: hb.ResourceFile,
    non_split_counts: dict[str, hb.ResourceFile],
    filtered_ranges: hb.ResourceFile,
    cohort_id: str,
    job_attrs: dict,
) -> tuple[Job, hb.ResourceFile]:
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


def fraser_analysis(b: hb.Batch, fds_tar: hb.ResourceFile, cohort_id: str, job_attrs: dict) -> Job:
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
        output={
            'sig_results': 'results.significant.csv',
            'all_results': 'results.all.csv',
            'plots': 'plots.tar.gz',
        }
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

    tar -czvf {j.output.plots} plots/
    cp results.significant.csv {j.output.sig_results}
    cp results.all.csv {j.output.all_results}
    cp statistics_summary.txt {j.output.stats}
    """
    j.command(command(cmd))
    return j
