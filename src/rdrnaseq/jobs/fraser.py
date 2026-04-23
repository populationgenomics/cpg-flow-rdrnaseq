# ruff: noqa: PLR0912
# ruff: noqa: PLR0915
import hailtop.batch as hb
from cpg_flow.resources import HIGHMEM, STANDARD, MachineType
from cpg_flow.utils import exists
from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, get_config
from cpg_utils.hail_batch import command, get_batch
from hailtop.batch.job import Job
from loguru import logger

R_INIT = 'RDrnaseq/fraser_init.R'
R_COUNT_SPLIT = 'RDrnaseq/fraser_count_split.R'
R_MERGE_SPLIT = 'RDrnaseq/fraser_merge_split.R'
R_COUNT_NON_SPLIT = 'RDrnaseq/fraser_count_non_split.R'
R_MERGE_NON_SPLIT = 'RDrnaseq/fraser_merge_non_split.R'
R_JOIN_COUNTS = 'RDrnaseq/fraser_join_counts.R'
R_ANALYSIS = 'RDrnaseq/fraser_analysis.R'

BASE_STORAGE_GB_COHORT = config_retrieve(['cohort_job_resources', 'base_storage_gb'], 100)
PER_BAM_STORAGE_COHORT = config_retrieve(['cohort_job_resources', 'per_bam_storage'], 10)
NCPU_COHORT = config_retrieve(['cohort_job_resources', 'ncpu'], 10)

BASE_STORAGE_GB_SAMPLE = config_retrieve(['sample_job_resources', 'base_storage_gb'], 100)
PER_BAM_STORAGE_SAMPLE = config_retrieve(['sample_job_resources', 'per_bam_storage'], 10)
NCPU_SAMPLE = config_retrieve(['sample_job_resources', 'ncpu'], 16)


def fraser_storage_required_gb(num_bams: int, base_storage_gb: int, per_bam_storage_gb: int) -> int:
    """Calculate disk space needed based on cohort size."""
    return base_storage_gb + num_bams * per_bam_storage_gb


def get_fraser_job(
    batch: hb.Batch,
    name: str,
    job_attrs: dict,
    n_samples: int,
    base_storage_gb: int,
    per_bam_storage: int,
    ncpu: int,
    machine_required: MachineType,  # would have to import this abstract class
) -> tuple[Job, int]:
    """Create a standard FRASER job with common configuration."""
    storage = fraser_storage_required_gb(n_samples, base_storage_gb, per_bam_storage)
    j = batch.new_job(name, attributes=job_attrs | {'tool': 'fraser'})
    j.image(config_retrieve(['images', 'fraser']))
    j.command('export HDF5_USE_FILE_LOCKING=FALSE')

    res = machine_required.set_resources(j=j, ncpu=ncpu, storage_gb=storage)

    return j, res.get_nthreads()


def fraser_pipeline(
    input_bams: list[tuple[str, str]],
    cohort_id: str,
    job_attrs: dict,
    output_fds_path: dict[str, Path],
    output_prefix: str,
) -> list[Job]:
    b = get_batch()
    root = to_path(output_prefix)

    all_jobs = []

    # Localize inputs as ResourceGroups
    input_bams_localised: dict[str, hb.ResourceGroup] = {
        sample_id: b.read_input_group(
            **{
                'bam': input_bam,
                'bam.bai': f'{input_bam}.bai',
            }
        )
        for sample_id, input_bam in input_bams
    }

    # Init
    init_fds_path = root / 'init' / 'fds-object.RDS'
    logger.info('Planning Fraser init job')
    fds_init_res, j_init = fraser_init(b, input_bams_localised, cohort_id, job_attrs, init_fds_path)
    if j_init:
        all_jobs.append(j_init)

    # Count Split
    split_counts_res = {}
    count_split_jobs = []
    logger.info('Planning Fraser count-split jobs')
    for sid, bam_rg in input_bams_localised.items():
        out_p = root / 'split_counts' / f'{sid}.RDS'
        res, j = fraser_count_split_reads(b, fds_init_res, bam_rg, sid, cohort_id, job_attrs, out_p)
        split_counts_res[sid] = res
        if j:
            count_split_jobs.append(j)

    logger.info(f'{len(count_split_jobs)} count-split jobs found')

    if all_jobs and count_split_jobs:
        for job in count_split_jobs:
            job.depends_on(all_jobs[-1])
    if count_split_jobs:
        all_jobs.extend(count_split_jobs)

    # Merge Split
    merge_split_paths = {
        'fds_object': root / 'merge_split' / 'fds-object.RDS',
        'g_ranges_split': root / 'merge_split' / 'g_ranges_split_counts.RDS',
        'g_ranges_non_split': root / 'merge_split' / 'g_ranges_non_split_counts.RDS',
        'splice_site_coords': root / 'merge_split' / 'splice_site_coords.RDS',
    }

    ref_sid = list(input_bams_localised.keys())[0]  # noqa: RUF015

    merge_split_res, j_merge_split = fraser_merge_split_reads(
        b,
        fds_init_res,
        split_counts_res,
        cohort_id,
        job_attrs,
        merge_split_paths,
        reference_bam_rg=input_bams_localised[ref_sid],  # Pass the ResourceGroup here
    )
    if j_merge_split and all_jobs:
        j_merge_split.depends_on(*all_jobs)
    if j_merge_split:
        all_jobs.append(j_merge_split)

    # Count Non-Split
    non_split_counts_res = {}
    count_non_split_jobs = []
    logger.info('Planning Count non-split job')
    for sid, bam_rg in input_bams_localised.items():
        out_p = root / 'non_split_counts' / f'{sid}.h5'
        res, j = fraser_count_non_split_reads(
            b, merge_split_res.fds_object, bam_rg, merge_split_res.splice_site_coords, sid, cohort_id, job_attrs, out_p
        )
        non_split_counts_res[sid] = res
        if j:
            count_non_split_jobs.append(j)

    logger.info(f'{len(count_non_split_jobs)} count-non-split jobs found')

    if all_jobs and count_non_split_jobs:
        for job in count_non_split_jobs:
            job.depends_on(all_jobs[-1])

    if count_non_split_jobs:
        all_jobs.extend(count_non_split_jobs)

    # Merge Non-Split
    fds_tar_path = to_path(output_fds_path['Rds_data'])
    logger.info('Planning Merge non-split job')
    fds_tar_res, j_merge_non = fraser_merge_non_split_reads(
        b,
        merge_split_res.fds_object,
        non_split_counts_res,
        merge_split_res.splice_site_coords,
        cohort_id,
        job_attrs,
        fds_tar_path,
        len(input_bams_localised),
        reference_bam_rg=input_bams_localised[ref_sid],
    )

    if j_merge_non and all_jobs:
        j_merge_non.depends_on(*all_jobs)
    if j_merge_non:
        all_jobs.append(j_merge_non)

    # Join counts
    joined_tar_path = root / 'joined' / 'fds_joined.tar.gz'
    fds_joined_res, j_join = fraser_join_counts(
        b,
        merge_split_res=merge_split_res,
        merge_non_split_tar=fds_tar_res,
        cohort_id=cohort_id,
        job_attrs=job_attrs,
        output_path=joined_tar_path,
        num_samples=len(input_bams_localised),
    )

    if j_join and all_jobs:
        j_join.depends_on(*all_jobs)
    if j_join:
        all_jobs.append(j_join)

    # Analysis
    analysis_paths = {
        'all_results': to_path(output_fds_path['seqr_data']),
        'significant_results': to_path(output_fds_path['sig_results']),
        'plots': root / 'results' / 'plots.tar.gz',
        'final_fds': root / 'results' / 'final_object.tar.gz',
    }
    logger.info('Planning Fraser analysis job')
    j_analysis = fraser_analysis(b, fds_joined_res, cohort_id, job_attrs, analysis_paths, len(input_bams_localised))

    if j_analysis and all_jobs:
        j_analysis.depends_on(*all_jobs)

    if j_analysis:
        all_jobs.append(j_analysis)

    return all_jobs


def fraser_init(b, input_bams, cohort_id, job_attrs, output_path) -> tuple[hb.Resource, Job | None]:
    if exists(output_path):
        return b.read_input(output_path), None

    j, threads = get_fraser_job(
        b,
        'fraser_init',
        job_attrs,
        n_samples=len(input_bams),
        base_storage_gb=BASE_STORAGE_GB_COHORT,
        per_bam_storage=PER_BAM_STORAGE_COHORT,
        ncpu=NCPU_COHORT,
        machine_required=HIGHMEM,
    )
    j.command('mkdir -p /io/batch/input_bams /io/work')
    j.command('echo "sample_id,bam_path" > /io/work/sample_map.csv')
    for sid, bam_rg in input_bams.items():
        j.command(f'ln -s {bam_rg.bam} /io/batch/input_bams/{sid}.bam')
        j.command(f'ln -s {bam_rg["bam.bai"]} /io/batch/input_bams/{sid}.bam.bai')
        j.command(f'echo "{sid},/io/batch/input_bams/{sid}.bam" >> /io/work/sample_map.csv')

    j.command(
        command(f"""
        Rscript {R_INIT} --cohort_id "{cohort_id}" --sample_map "/io/work/sample_map.csv" \\
            --work_dir "/io/work" --nthreads {threads}
        mv /io/work/savedObjects/FRASER_{cohort_id}/fds-object.RDS {j.fds}
    """)
    )
    b.write_output(j.fds, str(output_path))
    return j.fds, j


def fraser_count_split_reads(b, fds, bam_rg, sample_id, cohort_id, job_attrs, output_path):
    if exists(output_path):
        return b.read_input(output_path), None

    j, threads = get_fraser_job(
        b,
        f'count_split_{sample_id}',
        job_attrs,
        n_samples=1,
        base_storage_gb=BASE_STORAGE_GB_SAMPLE,
        per_bam_storage=PER_BAM_STORAGE_SAMPLE,
        ncpu=NCPU_SAMPLE,
        machine_required=HIGHMEM,
    )

    j.command('mkdir -p /io/batch/input_bams')
    j.command(f'ln -s {bam_rg.bam} /io/batch/input_bams/{sample_id}.bam')
    j.command(f'ln -s {bam_rg["bam.bai"]} /io/batch/input_bams/{sample_id}.bam.bai')

    j.command(
        command(f"""
        Rscript {R_COUNT_SPLIT} --fds_path {fds} --bam_path "/io/batch/input_bams/{sample_id}.bam" \\
            --cohort_id "{cohort_id}" --sample_id "{sample_id}" --work_dir "/io/work" \\
            --nthreads "{threads}"
        mv /io/work/cache/splitCounts/splitCounts-{sample_id}.RDS {j.out}
    """)
    )
    b.write_output(j.out, str(output_path))
    return j.out, j


def fraser_merge_split_reads(
    b,
    fds,
    split_counts,
    cohort_id,
    job_attrs,
    output_paths,
    reference_bam_rg=None,
):
    """Merge split read counts and create GRanges objects for split and non-split counts.
    Also localizing one reference bam to pass fraser metadata check even though the bam is not used in the merge step"""
    all_output_paths = output_paths | {'split_counts_tar': output_paths['fds_object'].parent / 'split_counts.tar.gz'}

    if all(exists(p) for p in all_output_paths.values()):
        return b.read_input_group(**all_output_paths), None

    j, threads = get_fraser_job(
        b,
        'fraser_merge_split',
        job_attrs,
        n_samples=len(split_counts),
        base_storage_gb=BASE_STORAGE_GB_COHORT,
        per_bam_storage=PER_BAM_STORAGE_COHORT,
        ncpu=NCPU_COHORT,
        machine_required=HIGHMEM,
    )
    j.declare_resource_group(out={k: v.name for k, v in all_output_paths.items()})

    setup = ['mkdir -p /io/work/cache/splitCounts /io/work/savedObjects /io/batch/input_bams']

    for sid, r in split_counts.items():
        setup.append(f'ln -s {r} /io/work/cache/splitCounts/splitCounts-{sid}.RDS')

    if reference_bam_rg:
        setup.append(f'ln -s {reference_bam_rg.bam} /io/batch/input_bams/reference.bam')
        setup.append(f'ln -s {reference_bam_rg["bam.bai"]} /io/batch/input_bams/reference.bam.bai')

    j.command(
        command(
            'ulimit -n 4096\n'
            + '\n'.join(setup)
            + f"""
        Rscript {R_MERGE_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --work_dir "/io/work" --nthreads "{threads}"

    mv /io/work/savedObjects/FRASER_{cohort_id}/fds-object.RDS {j.out.fds_object}
    mv /io/work/g_ranges_split_counts.RDS {j.out.g_ranges_split}
    mv /io/work/g_ranges_non_split_counts.RDS {j.out.g_ranges_non_split}
    mv /io/work/splice_site_coords.RDS {j.out.splice_site_coords}

        tar -czf {j.out.split_counts_tar} -C /io/work/savedObjects/FRASER_{cohort_id}/ splitCounts/
"""
        )
    )

    for k, p in all_output_paths.items():
        b.write_output(j.out[k], str(p))
    return j.out, j


def fraser_count_non_split_reads(b, fds, bam_rg, coords, sample_id, cohort_id, job_attrs, output_path):
    if exists(output_path):
        return b.read_input(output_path), None

    j, threads = get_fraser_job(
        b,
        f'count_non_split_{sample_id}',
        job_attrs,
        n_samples=1,
        base_storage_gb=BASE_STORAGE_GB_SAMPLE,
        per_bam_storage=PER_BAM_STORAGE_SAMPLE,
        ncpu=NCPU_SAMPLE,
        machine_required=HIGHMEM,
    )

    j.command('mkdir -p /io/batch/input_bams')
    j.command(f'ln -s {bam_rg.bam} /io/batch/input_bams/{sample_id}.bam')
    j.command(f'ln -s {bam_rg["bam.bai"]} /io/batch/input_bams/{sample_id}.bam.bai')

    j.command(
        command(f"""
        Rscript {R_COUNT_NON_SPLIT} --fds_path {fds} --bam_path "/io/batch/input_bams/{sample_id}.bam" \\
            --coords_path {coords} --sample_id "{sample_id}" --cohort_id "{cohort_id}" \\
            --work_dir "/io/work" --nthreads "{threads}"

        if [ -f /io/work/cache/nonSplicedCounts/FRASER_{cohort_id}/nonSplicedCounts-{sample_id}.h5 ]; then
            echo "Found H5 file, moving to output"
            mv /io/work/cache/nonSplicedCounts/FRASER_{cohort_id}/nonSplicedCounts-{sample_id}.h5 {j.out}
        else
            echo "ERROR: No output file found for {sample_id}"
            find /io/work/cache/nonSplicedCounts/ -type f
            exit 1
        fi
    """)
    )
    b.write_output(j.out, str(output_path))
    return j.out, j


def fraser_merge_non_split_reads(
    b, fds, non_split_counts, filtered_ranges, cohort_id, job_attrs, output_path, num_samples, reference_bam_rg
):
    if exists(output_path):
        return b.read_input(output_path), None

    j, threads = get_fraser_job(
        b,
        'fraser_merge_non_split',
        job_attrs,
        n_samples=num_samples,
        base_storage_gb=BASE_STORAGE_GB_COHORT,
        per_bam_storage=PER_BAM_STORAGE_COHORT,
        ncpu=NCPU_COHORT,
        machine_required=HIGHMEM,
    )
    fds_name = f'FRASER_{cohort_id}'
    cache_h5_path = '/io/work/cache/nonSplicedCounts'

    setup = [f'mkdir -p {cache_h5_path} /io/batch/input_bams', f'mkdir -p /io/work/savedObjects/{fds_name}']

    if not non_split_counts:
        raise ValueError('non_split_counts dictionary is empty')

    for sid, r in non_split_counts.items():
        setup.append(f'ln -s {r} {cache_h5_path}/nonSplicedCounts-{sid}.h5')

    if reference_bam_rg:
        setup.append(f'ln -s {reference_bam_rg.bam} /io/batch/input_bams/reference.bam')
        setup.append(f'ln -s {reference_bam_rg["bam.bai"]} /io/batch/input_bams/reference.bam.bai')

    j.command(
        command(
            'ulimit -n 4096\n'
            + '\n'.join(setup)
            + f"""

        Rscript {R_MERGE_NON_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --filtered_ranges_path {filtered_ranges} --work_dir "/io/work" --nthreads "{threads}"

        tar -cvzf {j.fds_tar} -C /io/work/savedObjects/ {fds_name}/
    """
        )
    )
    b.write_output(j.fds_tar, str(output_path))
    return j.fds_tar, j


def fraser_join_counts(
    b: hb.Batch,
    merge_split_res: hb.ResourceGroup,
    merge_non_split_tar: hb.ResourceFile,
    cohort_id: str,
    job_attrs: dict,
    output_path: Path,
    num_samples: int,
) -> tuple[hb.Resource, Job | None]:
    if exists(output_path):
        return b.read_input(output_path), None

    j, threads = get_fraser_job(
        b,
        'fraser_join_counts',
        job_attrs,
        n_samples=num_samples,
        base_storage_gb=BASE_STORAGE_GB_COHORT,
        per_bam_storage=PER_BAM_STORAGE_COHORT,
        ncpu=NCPU_COHORT,
        machine_required=HIGHMEM,
    )

    fds_name = f'FRASER_{cohort_id}'
    work_dir = '/io/work'

    setup = [
        f'mkdir -p {work_dir}/savedObjects/{fds_name}',
        f'mkdir -p {work_dir}/savedObjects',
        f'tar -xzf {merge_non_split_tar} -C {work_dir}/savedObjects/',
        f'tar -xzf {merge_split_res.split_counts_tar} -C {work_dir}/savedObjects/{fds_name}/',
        f'cp {merge_split_res.splice_site_coords} {work_dir}/splice_site_coords.RDS',
    ]

    cmd = f"""
ulimit -n 4096
{chr(10).join(setup)}

Rscript {R_JOIN_COUNTS} --fds_path "{work_dir}/savedObjects/{fds_name}/fds-object.RDS" \\
    --cohort_id "{cohort_id}" \\
    --work_dir "{work_dir}" \\
    --nthreads "{threads}"

echo "=== R script completed, creating tar archive ==="
tar -cvzf {j.fds_tar} -C {work_dir}/savedObjects/ {fds_name}/
"""

    j.command(command(cmd))
    b.write_output(j.fds_tar, str(output_path))
    return j.fds_tar, j


def fraser_analysis(b, fds_tar, cohort_id, job_attrs, output_paths, num_samples) -> Job | None:
    if all(exists(x) for x in output_paths.values()):
        return None

    j, threads = get_fraser_job(
        b,
        'fraser_analysis',
        job_attrs,
        n_samples=num_samples,
        base_storage_gb=BASE_STORAGE_GB_COHORT,
        per_bam_storage=PER_BAM_STORAGE_COHORT,
        ncpu=NCPU_COHORT,
        machine_required=HIGHMEM,
    )
    j.declare_resource_group(out={k: v.name for k, v in output_paths.items()})
    cfg = get_config().get('fraser', {})
    z_cutoff_arg = f'--z_cutoff {cfg["z_cutoff"]}' if 'z_cutoff' in cfg else ''

    j.command(
        command(f"""
        mkdir -p /io/work/savedObjects
        tar -xzf {fds_tar} -C /io/work/savedObjects/

        Rscript {R_ANALYSIS} --fds_dir "/io/work/savedObjects" --cohort_id "{cohort_id}" \\
            --pval_cutoff {cfg.get('pval_cutoff', 0.05)} \\
            --delta_psi_cutoff {cfg.get('delta_psi_cutoff', 0.3)} \\
            --min_count {cfg.get('min_count', 5)} \\
            --nthreads {threads} {z_cutoff_arg}

        tar -czvf {j.out.plots} qc_plots
        cp {cohort_id}.significant.csv {j.out.significant_results}
        cp {cohort_id}.all_results.csv.gz {j.out.all_results}
        tar -czvf {j.out.final_fds} -C savedObjects {cohort_id}_final/
    """)
    )

    for k, p in output_paths.items():
        b.write_output(j.out[k], str(p))
    return j
