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

def fraser_storage_required_gb(num_bams: int, base_storage_gb: int, per_bam_storage_gb: int) -> int:
    return base_storage_gb + num_bams * per_bam_storage_gb

# --- Orchestration Pipeline ---

def fraser_pipeline(
    input_bams_or_crams: list[tuple[str, BamPath, None] | tuple[str, CramPath, Path]],
    cohort_id: str,
    job_attrs: dict,
    output_prefix: Path,
):
    b = get_batch()
    root = to_path(output_prefix / cohort_id)

    # 0. Localize/Convert Inputs
    input_bams_localised: dict[str, hb.ResourceFile] = {}
    for sample_id, input_bam_or_cram, potential_bam_path in input_bams_or_crams:
        if isinstance(input_bam_or_cram, CramPath) and isinstance(potential_bam_path, Path):
            j, output_bam = cram_to_bam(
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
    for sid in input_bams_localised:
        out_p = root / 'split_counts' / f'{sid}.RDS'
        split_counts_res[sid] = fraser_count_split_reads(b, fds_res, sid, cohort_id, job_attrs, out_p)

    # 3. Merge Split
    merge_split_paths = {
        'g_ranges_split': root / 'merge_split' / 'g_ranges_split_counts.RDS',
        'g_ranges_non_split': root / 'merge_split' / 'g_ranges_non_split_counts.RDS',
        'splice_site_coords': root / 'merge_split' / 'splice_site_coords.RDS',
    }
    merge_out_rg = fraser_merge_split_reads(b, fds_res, split_counts_res, cohort_id, job_attrs, merge_split_paths)

    # 4. Count Non-Split
    non_split_counts_res = {}
    for sid in input_bams_localised:
        # Note: FRASER 2.0 uses .h5 primarily
        out_p = root / 'non_split_counts' / f'{sid}.h5'
        non_split_counts_res[sid] = fraser_count_non_split_reads(
            b, fds_res, merge_out_rg.splice_site_coords, sid, cohort_id, job_attrs, out_p
        )

    # 5. Merge Non-Split
    fds_tar_path = root / 'merge_non_split' / f'FRASER_{cohort_id}.tar.gz'
    fds_tar_res = fraser_merge_non_split_reads(
        b, fds_res, non_split_counts_res, merge_out_rg.g_ranges_non_split, cohort_id, job_attrs, fds_tar_path
    )

    # 6. Analysis
    analysis_paths = {
        'sig_results': root / 'results' / 'results.significant.csv',
        'all_results': root / 'results' / 'results.all.csv',
        'plots': root / 'results' / 'plots.tar.gz',
        'stats': root / 'results' / 'statistics_summary.txt'
    }
    fraser_analysis(b, fds_tar_res, cohort_id, job_attrs, analysis_paths)


# --- Step Methods ---

def fraser_init(b, input_bams, cohort_id, job_attrs, output_path) -> hb.ResourceFile:
    if output_path.exists():
        return b.read_input(output_path)

    j = get_fraser_job(b, 'fraser_init', job_attrs)
    res = HIGHMEM.set_resources(j, ncpu=10, storage_gb=100)

    bam_files_r = ', '.join(f'"{bam}"' for bam in input_bams.values())
    sample_ids_r = ', '.join(f'"{sam}"' for sam in input_bams.keys())

    j.command(command(f"""
        Rscript {R_INIT} --cohort_id "{cohort_id}" --sample_ids "{sample_ids_r}" \\
            --bam_files "{bam_files_r}" --nthreads {res.get_nthreads()}
        mv output/savedObjects/FRASER_{cohort_id}/fds-object.RDS {j.fds}
    """))
    b.write_output(j.fds, output_path)
    return j.fds

def fraser_count_split_reads(b, fds, sample_id, cohort_id, job_attrs, output_path) -> hb.ResourceFile:
    if output_path.exists():
        return b.read_input(output_path)

    j = get_fraser_job(b, f'count_split_{sample_id}', job_attrs)
    res = STANDARD.set_resources(j, ncpu=4, storage_gb=20)
    j.command(command(f"""
        Rscript {R_COUNT_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --sample_id "{sample_id}" --nthreads "{res.get_nthreads()}"
        mv output/cache/splitCounts/splitCounts-{sample_id}.RDS {j.out}
    """))
    b.write_output(j.out, str(output_path))
    return j.out

def fraser_merge_split_reads(b, fds, split_counts, cohort_id, job_attrs, output_paths) -> hb.ResourceGroup:
    # Check if all outputs exist
    if all(p.exists() for p in output_paths.values()):
        return b.read_input_group(**{k: str(v) for k, v in output_paths.items()})

    j = get_fraser_job(b, 'fraser_merge_split', job_attrs)
    j.declare_resource_group(out={k: v.name for k, v in output_paths.items()})
    res = HIGHMEM.set_resources(j, ncpu=10, storage_gb=50)

    cache_path = f'output/savedObjects/{cohort_id}/cache/splitCounts'
    setup_cmds = [f'mkdir -p {cache_path}']
    for sid, r_file in split_counts.items():
        setup_cmds.append(f'ln -s {r_file} {cache_path}/splitCounts-{sid}.RDS')

    j.command(command("\n".join(setup_cmds) + f"""
        Rscript {R_MERGE_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" --nthreads "{res.get_nthreads()}"
        mv g_ranges_split_counts.RDS {j.out.g_ranges_split}
        mv g_ranges_non_split_counts.RDS {j.out.g_ranges_non_split}
        mv splice_site_coords.RDS {j.out.splice_site_coords}
    """))

    for k, p in output_paths.items():
        b.write_output(j.out[k], str(p))
    return j.out

def fraser_count_non_split_reads(b, fds, coords, sample_id, cohort_id, job_attrs, output_path) -> hb.ResourceFile:
    if output_path.exists():
        return b.read_input(output_path)

    j = get_fraser_job(b, f'count_non_split_{sample_id}', job_attrs)
    res = STANDARD.set_resources(j, ncpu=4, storage_gb=20)
    j.command(command(f"""
        Rscript {R_COUNT_NON_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --sample_id "{sample_id}" --coords_path {coords} --nthreads "{res.get_nthreads()}"
        # Fallback for .h5 vs .RDS
        if [ -f "output/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sample_id}.h5" ]; then
            mv output/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sample_id}.h5 {j.out}
        else
            mv output/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sample_id}.RDS {j.out}
        fi
    """))
    b.write_output(j.out, str(output_path))
    return j.out

def fraser_merge_non_split_reads(b, fds, non_split_counts, filtered_ranges, cohort_id, job_attrs, output_path) -> hb.ResourceFile:
    if output_path.exists():
        return b.read_input(output_path)

    j = get_fraser_job(b, 'fraser_merge_non_split', job_attrs)
    res = HIGHMEM.set_resources(j, ncpu=10, storage_gb=50)

    setup_cache = '\n'.join([
        f'mkdir -p output/cache/nonSplicedCounts/{cohort_id} && ln -s {r} output/cache/nonSplicedCounts/{cohort_id}/nonSplicedCounts-{sid}.h5'
        for sid, r in non_split_counts.items()
    ])

    j.command(command(f"""
        {setup_cache}
        Rscript {R_MERGE_NON_SPLIT} --fds_path {fds} --cohort_id "{cohort_id}" \\
            --filtered_ranges_path {filtered_ranges} --nthreads "{res.get_nthreads()}"
        tar -czvf {j.fds_tar} output/savedObjects/FRASER_{cohort_id}/
    """))
    b.write_output(j.fds_tar, str(output_path))
    return j.fds_tar

def fraser_analysis(b, fds_tar, cohort_id, job_attrs, output_paths):
    if all(p.exists() for p in output_paths.values()):
        return None

    j = get_fraser_job(b, f'fraser_analysis_{cohort_id}', job_attrs)
    j.declare_resource_group(output={k: v.name for k, v in output_paths.items()})
    res = HIGHMEM.set_resources(j, ncpu=10, storage_gb=50)

    cfg = get_config().get('fraser', {})
    j.command(command(f"""
        tar -xvf {fds_tar}
        Rscript {R_ANALYSIS} --fds_dir "output" --cohort_id "{cohort_id}" \\
            --pval_cutoff {cfg.get('pval_cutoff', 0.05)} --delta_psi_cutoff {cfg.get('delta_psi_cutoff', 0.3)} \\
            --nthreads {res.get_nthreads()}
        tar -czvf {j.output.plots} plots/
        cp results.significant.csv {j.output.sig_results}
        cp results.all.csv {j.output.all_results}
        cp statistics_summary.txt {j.output.stats}
    """))

    for k, p in output_paths.items():
        b.write_output(j.output[k], str(p))
    return j
