import hailtop.batch as hb
import pandas as pd
from cpg_flow.utils import exists
from cpg_flow.resources import HIGHMEM, STANDARD
from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, get_config
from cpg_utils.hail_batch import command, get_batch
from hailtop.batch.job import Job

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
    input_bams: list[tuple[str, str]],
    cohort_id: str,
    job_attrs: dict,
    output_fds_path: dict[str, Path],
    output_prefix: str,
) -> list[Job]:
    b = get_batch()
    # Use output_prefix for intermediate steps and extract final paths from dict
    root = to_path(output_prefix)
    all_jobs = []

    # 0. Localize/Convert Inputs
    input_bams_localised: dict[str, hb.ResourceFile] = {
        sample_id: b.read_input_group(
            **{
                'bam': input_bam,
                'bam.bai': f'{input_bam}.bai',
            }
        ).bam
        for sample_id, input_bam in input_bams
    }

    # 1. Init
    init_fds_path = root / 'init' / 'fds-object.RDS'
    fds_res, j_init = fraser_init(b, input_bams_localised, cohort_id, job_attrs, init_fds_path)
    all_jobs.append(j_init)

    # 2. Count Split
    split_counts_res = {}
    for sid, bam_rg in input_bams_localised.items():
        out_p = root / 'split_counts' / f'{sid}.RDS'
        res, j = fraser_count_split_reads(b, fds_res, bam_rg, sid, cohort_id, job_attrs, out_p)
        split_counts_res[sid] = res
        all_jobs.append(j)

    # 3. Merge Split
    merge_split_paths = {
        'g_ranges_split': root / 'merge_split' / 'g_ranges_split_counts.RDS',
        'g_ranges_non_split': root / 'merge_split' / 'g_ranges_non_split_counts.RDS',
        'splice_site_coords': root / 'merge_split' / 'splice_site_coords.RDS',
    }
    merge_out_rg, j_merge_split = fraser_merge_split_reads(
        b, fds_res, split_counts_res, cohort_id, job_attrs, merge_split_paths
    )
    all_jobs.append(j_merge_split)

    # 4. Count Non-Split
    non_split_counts_res = {}
    for sid, bam_rg in input_bams_localised.items():
        out_p = root / 'non_split_counts' / f'{sid}.h5'
        res, j = fraser_count_non_split_reads(
            b, fds_res, bam_rg, merge_out_rg.splice_site_coords, sid, job_attrs, out_p
        )
        non_split_counts_res[sid] = res
        all_jobs.append(j)

    # 5. Merge Non-Split
    fds_tar_path = to_path(output_fds_path['Rds_data'])
    fds_tar_res, j_merge_non = fraser_merge_non_split_reads(
        b,
        fds_res,
        non_split_counts_res,
        merge_out_rg.g_ranges_non_split,
        cohort_id,
        job_attrs,
        fds_tar_path,
        len(input_bams_localised),
    )
    all_jobs.append(j_merge_non)

    # 6. Analysis
    analysis_paths = {
        'sig_results': root / 'results' / 'results.significant.csv',
        'all_results': to_path(output_fds_path['seqr_data']),
        'plots': root / 'results' / 'plots.tar.gz',
        'stats': root / 'results' / 'statistics_summary.txt',
    }
    j_analysis = fraser_analysis(b, fds_tar_res, cohort_id, job_attrs, analysis_paths, len(input_bams_localised))
    all_jobs.append(j_analysis)

    return all_jobs


def fraser_init(b, input_bams, cohort_id, job_attrs, output_path) -> tuple[hb.ResourceFile, Job|None]:
    if exists(output_path):
        return b.read_output(str(output_path)), None

    j = get_fraser_job(b, 'fraser_init', job_attrs)

    storage = fraser_storage_required_gb(len(input_bams), 50, 10)
    res = HIGHMEM.set_resources(j=j, ncpu=10, storage_gb=storage)

    sample_df = pd.DataFrame([{'sample_id': s, 'bam_path': str(bam.bam)} for s, bam in input_bams.items()])
    tmp_csv = to_path(output_path).parent / 'sample_map.csv'
    with tmp_csv.open('w') as f:
        sample_df.to_csv(f, index=False)

    j.command(
        command(f"""
        Rscript {R_INIT} --cohort_id "{cohort_id}" --sample_map {b.read_input(str(tmp_csv))} \\
            --work_dir "/io/work" --nthreads {res.get_nthreads()}
        mv /io/work/savedObjects/FRASER_{cohort_id}/fds-object.RDS {j.fds}
    """)
    )
    b.write_output(j.fds, str(output_path))
    return j.fds, j


def fraser_count_split_reads(b, fds, bam_rg, sample_id, cohort_id, job_attrs, output_path):
    if exists(output_path):
        return b.read_output(str(output_path)), None

    j = get_fraser_job(b, f'count_split_{sample_id}', job_attrs)
    res = STANDARD.set_resources(j=j, ncpu=4, storage_gb=20)
    j.command(
        command(f"""
        ln -s {bam_rg.bam} /io/input.bam
        ln -s {bam_rg.bai} /io/input.bam.bai
        Rscript {R_COUNT_SPLIT} --fds_path {fds} --bam_path /io/input.bam --cohort_id "{cohort_id}" \\
            --sample_id "{sample_id}" --work_dir "/io/work" --nthreads "{res.get_nthreads()}"
        mv /io/work/cache/splitCounts/splitCounts-{sample_id}.RDS {j.out}
    """)
    )
    b.write_output(j.out, str(output_path))
    return j.out, j


def fraser_merge_split_reads(b, fds, split_counts, cohort_id, job_attrs, output_paths):
    if all(exists(outputs) for outputs in output_paths.values()):
        return b.read_output(str(output_paths)), None

    j = get_fraser_job(b, 'fraser_merge_split', job_attrs)
    j.declare_resource_group(out={k: v.name for k, v in output_paths.items()})
    storage = fraser_storage_required_gb(len(split_counts), 50, 10)
    res = HIGHMEM.set_resources(j=j, ncpu=10, storage_gb=storage)

    setup = ['mkdir -p /io/work/cache/splitCounts']
    setup += [f'ln -s {r} /io/work/cache/splitCounts/splitCounts-{sid}.RDS' for sid, r in split_counts.items()]

    j.command(
        command(
            '\n'.join(setup)
            + f"""
        Rscript {R_MERGE_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --work_dir "/io/work" --nthreads "{res.get_nthreads()}"
        mv /io/work/g_ranges_split_counts.RDS {j.out.g_ranges_split}
        mv /io/work/g_ranges_non_split_counts.RDS {j.out.g_ranges_non_split}
        mv /io/work/splice_site_coords.RDS {j.out.splice_site_coords}
    """
        )
    )
    for k, p in output_paths.items():
        b.write_output(j.out[k], str(p))
    return j.out, j


def fraser_count_non_split_reads(b, fds, bam_rg, coords, sample_id, job_attrs, output_path):
    if exists(output_path):
        return b.read_output(str(output_path)), None

    j = get_fraser_job(b, f'count_non_split_{sample_id}', job_attrs)
    res = STANDARD.set_resources(j=j, ncpu=4, storage_gb=20)
    j.command(
        command(f"""
        ln -s {bam_rg.bam} /io/input.bam
        ln -s {bam_rg.bai} /io/input.bam.bai
        Rscript {R_COUNT_NON_SPLIT} --fds_path {fds} --bam_path /io/input.bam --coords_path {coords} \\
            --sample_id "{sample_id}" --work_dir "/io/work" --nthreads "{res.get_nthreads()}"
        find /io/work/cache/nonSplicedCounts/ -name "nonSplicedCounts-{sample_id}.*" -exec mv {{}} {j.out} \\;
    """)
    )
    b.write_output(j.out, str(output_path))
    return j.out, j


def fraser_merge_non_split_reads(
    b, fds, non_split_counts, filtered_ranges, cohort_id, job_attrs, output_path, num_samples
):
    if exists(output_path):
        return b.read_output(str(output_path)), None
    j = get_fraser_job(b, 'fraser_merge_non_split', job_attrs)
    storage = fraser_storage_required_gb(num_samples, base_storage_gb=50, per_bam_storage_gb=10)
    res = HIGHMEM.set_resources(j=j, ncpu=10, storage_gb=storage)

    setup = [f'mkdir -p /io/work/cache/nonSplicedCounts/{cohort_id}']
    for sid, r in non_split_counts.items():
        setup.append(f'EXT=$(echo {r} | sed "s/.*\\.//")')
        setup.append(f'ln -s {r} /io/work/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sid}.$EXT')

    j.command(
        command(
            '\n'.join(setup)
            + f"""
        Rscript {R_MERGE_NON_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --filtered_ranges_path {filtered_ranges} --work_dir "/io/work" --nthreads "{res.get_nthreads()}"
        tar -czf {j.fds_tar} -C /io/work/savedObjects/ FRASER_{cohort_id}/
    """
        )
    )
    b.write_output(j.fds_tar, str(output_path))
    return j.fds_tar, j


def fraser_analysis(b, fds_tar, cohort_id, job_attrs, output_paths, num_samples):
    if all(exists(x) for x in output_paths.values()):
        return None

    j = get_fraser_job(b, 'fraser_analysis', job_attrs)
    j.declare_resource_group(out={k: v.name for k, v in output_paths.items()})
    storage = fraser_storage_required_gb(num_samples, 50, 10)
    res = HIGHMEM.set_resources(j=j, ncpu=10, storage_gb=storage)

    cfg = get_config().get('fraser', {})
    j.command(
        command(f"""
        mkdir -p /io/work
        tar -xf {fds_tar} -C /io/work/
        Rscript {R_ANALYSIS} --fds_dir "/io/work" --cohort_id "{cohort_id}" \\
            --pval_cutoff {cfg.get('pval_cutoff', 0.05)} --delta_psi_cutoff {cfg.get('delta_psi_cutoff', 0.3)} \\
            --nthreads {res.get_nthreads()}
        tar -czf {j.out.plots} -C /io/work/ plots/
        cp /io/work/results.significant.csv {j.out.sig_results}
        cp /io/work/results.all.csv {j.out.all_results}
        cp /io/work/statistics_summary.txt {j.out.stats}
    """)
    )
    for k, p in output_paths.items():
        b.write_output(j.out[k], str(p))
    return j
