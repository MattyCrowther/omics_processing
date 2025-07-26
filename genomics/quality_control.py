from pathlib import Path
import sys
import os

sys.path.append("../")

from qc.alignment import (
    flagstat_summary,
    alignment_summary,
    coverage_metrics,
    size_distribution,
    quality_depth,
    gc_bias,
    qualimap_bam,
    coverage_depth_distribution
)

def perform_qc(bam_file, reference_fasta, base_dir=Path(".")):
    """
    Run post-alignment QC metrics on a BAM file using reference genome.
    All outputs are saved into a specified QC directory.

    Parameters:
        bam_file (str or Path): Path to deduplicated, coordinate-sorted BAM file.
        reference_fasta (str or Path): Path to reference FASTA used in alignment.
        qc_dir (str or Path): Directory to write all QC output files. Defaults to 'qc/'.

    Returns:
        list of Path: All QC result file paths.
    """
    qc_dir = base_dir / "qc"
    bam_file = Path(bam_file)
    reference_fasta = Path(reference_fasta)
    qc_dir = Path(qc_dir)
    qc_dir.mkdir(exist_ok=True)

    print("\nRunning post-alignment QC...")

    alignment_dir = qc_dir / "alignment"
    coverage_dir = qc_dir / "coverage"
    insert_dir = qc_dir / "insert_size"
    stats_dir = qc_dir / "stats"
    bias_dir = qc_dir / "bias"
    qualimap_dir = qc_dir / "qualimap"

    # Ensure subdirectories exist
    for subdir in [alignment_dir, coverage_dir, insert_dir, stats_dir, bias_dir, qualimap_dir]:
        subdir.mkdir(parents=True, exist_ok=True)

    # Output paths
    flagstat_out = alignment_dir / f"flagstat_{bam_file.stem}.txt"
    align_metrics = alignment_dir / "alignment_metrics.txt"
    
    coverage_prefix = coverage_dir / bam_file.stem
    per_base_file = coverage_prefix.with_name(f"{coverage_prefix.name}.per-base.bed.gz")
    depth_plot = coverage_dir / f"{bam_file.stem}_depth_hist.png"

    insert_metrics = insert_dir / "insert_size_metrics.txt"

    stats_out = stats_dir / "samtools_stats.txt"
    gc_metrics = bias_dir / "gc_bias_metrics.txt"

    flagstat_summary(bam_file, flagstat_out)
    alignment_summary(reference_fasta, bam_file, align_metrics)
    coverage_metrics(bam_file, prefix=coverage_prefix)
    size_distribution(bam_file, output_filename=insert_metrics)
    quality_depth(bam_file, stats_out)
    gc_bias(bam_file, reference_fasta, gc_metrics)
    qualimap_bam(bam_file, qualimap_dir)

    # Optional: depth plot
    per_base_file = coverage_prefix.with_name(f"{coverage_prefix.name}.per-base.bed.gz")
    depth_plot = qc_dir / f"{bam_file.stem}_depth_hist.png"
    coverage_depth_distribution(per_base_file, depth_plot)

    return qc_dir
