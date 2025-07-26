import subprocess

def cut_adapt(reads, trimmed_dir, adapter="AGATCGGAAGAGC"):
    """Trim reads using cutadapt and output to trimmed_dir."""
    trimmed_dir.mkdir(exist_ok=True)
    if len(reads) == 2:
        # Paired-end
        out1 = trimmed_dir / f"{reads[0].stem}_trimmed.fastq.gz"
        out2 = trimmed_dir / f"{reads[1].stem}_trimmed.fastq.gz"
        cmd = [
            "cutadapt", "-a", adapter, "-A", adapter,
            "-o", str(out1), "-p", str(out2),
            str(reads[0]), str(reads[1]),
            "--quality-cutoff", "20", "--minimum-length", "30"
        ]
        outputs = [out1, out2]
    else:
        # Single-end
        out1 = trimmed_dir / f"{reads[0].stem}_trimmed.fastq.gz"
        cmd = [
            "cutadapt", "-a", adapter,
            "-o", str(out1), str(reads[0]),
            "--quality-cutoff", "20", "--minimum-length", "30"
        ]
        outputs = [out1]

    print("\nTrimming reads with cutadapt...")
    subprocess.run(cmd, check=True)
    print("Trimming complete.")
    return outputs

def fastp_trim(read1, read2, out1, out2):
    print(f"\n⚡ Trimming reads with fastp...")
    subprocess.run([
        "fastp",
        "-i", read1, "-I", read2,
        "-o", out1, "-O", out2,
        "--detect_adapter_for_pe",
        "--thread", "4",
        "--html", "fastp_report.html",
        "--json", "fastp_report.json"
    ], check=True)
    print("fastp trimming + QC complete.")