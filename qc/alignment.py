import subprocess
from pathlib import Path
import pandas as pd 
import matplotlib.pyplot as plt  

def fast_qc(reads, qc_dir):
    """
    Run FastQC on a list of FASTQ files and save results to qc_dir.

    Parameters:
        reads (list of str or Path): List of input FASTQ files.
        qc_dir (str or Path): Directory to save FastQC outputs.

    Returns:
        Path: Path to the output directory.
    """
    reads = [Path(r) for r in reads]
    qc_dir = Path(qc_dir)
    qc_dir.mkdir(exist_ok=True)

    subprocess.run(["fastqc", "-o", str(qc_dir), 
                    *[str(f) for f in reads]], 
                    check=True)
    return qc_dir

def multi_qc(target_dir, output_dir=None):
    """
    Run MultiQC on a target directory containing QC reports.

    Parameters:
        target_dir (str or Path): Directory containing input QC reports.
        output_dir (str or Path, optional): Output directory for MultiQC report.
            If None, MultiQC writes to the current working directory.

    Returns:
        Path: Path to the directory containing the MultiQC report.
    """
    target_dir = Path(target_dir)
    cmd = ["multiqc", str(target_dir)]

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        cmd.extend(["-o", str(output_dir)])
    else:
        output_dir = Path.cwd()
    subprocess.run(cmd, check=True)
    return output_dir

def flagstat_summary(bam_filename, output_file=None):
    """
    Run samtools flagstat on a BAM file.

    Parameters:
        bam_filename (str or Path): Path to input BAM file (coordinate-sorted, duplicate-marked alignments).
        output_file (str or Path, optional): Output file path. If not provided, defaults to 'flagstat_{basename}.txt'.

    Returns:
        Path: Path to the generated flagstat summary file.
    """
    bam_path = Path(bam_filename)
    if output_file:
        output_file = Path(output_file)  
    else: 
        Path(f"flagstat_{bam_path.stem}.txt")

    with open(output_file, "w") as out:
        subprocess.run(["samtools", "flagstat", str(bam_path)], 
                       stdout=out, check=True)

    return output_file

def alignment_summary(reference_fasta, bam_filename, 
                      output_filename=None):
    """
    Run Picard CollectAlignmentSummaryMetrics on a BAM file.

    Parameters:
        reference_fasta (str or Path): Path to reference FASTA file used for alignment.
        bam_filename (str or Path): Path to input BAM file (coordinate-sorted, duplicate-marked alignments).
        output_filename (str or Path, optional): Output file path. Defaults to 'alignment_metrics.txt'.

    Returns:
        Path: Path to the generated alignment metrics file.
    """
    reference_fasta = Path(reference_fasta)
    bam_filename = Path(bam_filename)

    if output_filename:
        output_filename = Path(output_filename)
    else:
        output_filename = Path("alignment_metrics.txt")

    subprocess.run([
        "picard",
        "CollectAlignmentSummaryMetrics",
        f"R={reference_fasta}",
        f"I={bam_filename}",
        f"O={output_filename}"
    ], check=True)

    return output_filename

def coverage_metrics(bam_filename, prefix=None):
    """
    Run mosdepth to compute coverage metrics for a BAM file.

    Parameters:
        bam_filename (str or Path): Path to input BAM file (coordinate-sorted, duplicate-marked alignments).
        prefix (str or Path, optional): Output prefix for mosdepth files (default: uses BAM basename without extension).

    Returns:
        dict: Paths to key mosdepth output files.
    """
    bam_filename = Path(bam_filename)

    if prefix:
        prefix = Path(prefix)
    else:
        prefix = Path(bam_filename.stem)

    subprocess.run([
        "mosdepth", "--threads", "4", str(prefix), str(bam_filename)
    ], check=True)

    outputs = {
        "summary": prefix.with_suffix(".mosdepth.summary.txt"),
        "regions": prefix.with_name(f"{prefix.name}.regions.bed.gz"),
        "per_base": prefix.with_name(f"{prefix.name}.per-base.bed.gz")
    }

    return outputs

def size_distribution(bam_filename, output_filename=None):
    """
    Run Picard CollectInsertSizeMetrics to compute fragment size distribution.

    Parameters:
        bam_filename (str or Path): Path to input BAM file (coordinate-sorted, duplicate-marked alignments).
        output_filename (str or Path, optional): Text output file for metrics. Defaults to 'insert_size_metrics.txt'.

    Returns:
        dict: Paths to the metrics text file and PDF histogram.
    """
    bam_filename = Path(bam_filename)

    if output_filename:
        output_filename = Path(output_filename)
    else:
        output_filename = Path("insert_size_metrics.txt")

    pdf_output = output_filename.with_suffix(".pdf")

    subprocess.run([
        "picard", "CollectInsertSizeMetrics",
        f"I={bam_filename}",
        f"O={output_filename}",
        f"H={pdf_output}",
        "M=0.5"
    ], check=True)

    return output_filename, pdf_output

def quality_depth(bam_filename, output_filename=None):
    """
    Run samtools stats to gather alignment quality and depth metrics.

    Parameters:
        bam_filename (str or Path): Path to input BAM file (coordinate-sorted, duplicate-marked alignments).
        output_filename (str or Path, optional): Output text file for samtools stats.
            Defaults to 'samtools_stats.txt'.

    Returns:
        Path: Path to the generated stats file.
    """
    bam_filename = Path(bam_filename)

    if output_filename:
        output_filename = Path(output_filename)
    else:
        output_filename = Path("samtools_stats.txt")

    with open(output_filename, "w") as out:
        subprocess.run(["samtools", "stats", str(bam_filename)],
                       stdout=out, check=True)

    return output_filename

def gc_bias(bam_filename, reference_filename, output_filename=None):
    """
    Run Picard CollectGcBiasMetrics to assess GC bias in aligned reads.

    Parameters:
        bam_filename (str or Path): Path to input BAM file (coordinate-sorted, duplicate-marked alignments).
        reference_filename (str or Path): Path to reference genome FASTA used for alignment.
        output_filename (str or Path, optional): Output metrics text file.
            Defaults to 'gc_bias_metrics.txt'.

    Returns:
        list: [metrics_file, chart_file, summary_file]
    """
    bam_filename = Path(bam_filename)
    reference_filename = Path(reference_filename)

    if output_filename:
        output_filename = Path(output_filename)
    else:
        output_filename = Path("gc_bias_metrics.txt")

    base_stem = output_filename.with_suffix("").name
    chart_filename = Path(f"{base_stem}.pdf")
    summary_filename = Path(f"{base_stem}.txt")

    subprocess.run([
        "picard", "CollectGcBiasMetrics",
        f"I={bam_filename}",
        f"O={output_filename}",
        f"CHART={chart_filename}",
        f"S={summary_filename}",
        f"R={reference_filename}"
    ], check=True)

    return output_filename, chart_filename, summary_filename
    
def qualimap_bam(bam_filename, out_dir=None):
    """
    Run Qualimap bamqc to generate a QC report on aligned reads.

    Parameters:
        bam_filename (str or Path): Path to input BAM file (coordinate-sorted, duplicate-marked alignments).
        out_dir (str or Path, optional): Output directory for the Qualimap report.
            Defaults to 'qualimap_report'.

    Returns:
        Path: Path to the output report directory.
    """
    bam_filename = Path(bam_filename)

    if out_dir:
        out_dir = Path(out_dir)
    else:
        out_dir = Path("qualimap_report")

    subprocess.run([
        "qualimap", "bamqc",
        "-bam", str(bam_filename),
        "-outdir", str(out_dir),
        "-nt", "4"
    ], check=True)

    return out_dir

def coverage_depth_distribution(sample_file, output_file=None):
    """
    Plot and save a histogram of per-base coverage depth from a mosdepth BED file.

    Parameters:
        sample_file (str or Path): Path to mosdepth .per-base.bed.gz file.
        output_file (str or Path, optional): Output path for the saved plot.
            Defaults to '{sample_file stem}_depth_hist.png'.

    Returns:
        Path: Path to the saved plot image.
    """
    sample_file = Path(sample_file)

    if output_file:
        output_file = Path(output_file)
    else:
        output_file = sample_file.with_name(f"{sample_file.stem}_depth_hist.png")

    df = pd.read_csv(sample_file, sep='\t', header=None,
                     names=["chrom", "start", "end", "depth"])

    plt.figure()
    plt.hist(df["depth"], bins=100, range=(0, 100))
    plt.xlabel("Depth")
    plt.ylabel("Number of Bases")
    plt.title("Coverage Depth Distribution")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    return output_file