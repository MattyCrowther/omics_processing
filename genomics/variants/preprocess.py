import subprocess
from pathlib import Path

def add_read_groups(
    input_bam: Path,
    output_bam: Path,
    rgid: str = "1",
    rglb: str = "lib1",
    rgpl: str = "illumina",
    rgpu: str = "unit1",
    rgsm: str = "sample",
):
    """
    Adds or replaces read group information in a BAM file using Picard.

    Parameters:
        input_bam (Path): Input BAM file.
        output_bam (Path): Output BAM file with read groups.
        rgid (str): Read Group ID.
        rglb (str): Read Group Library.
        rgpl (str): Platform (e.g., illumina).
        rgpu (str): Platform unit.
        rgsm (str): Sample name.
    """
    cmd = [
        "picard",
        "AddOrReplaceReadGroups",
        f"I={input_bam}",
        f"O={output_bam}",
        f"RGID={rgid}",
        f"RGLB={rglb}",
        f"RGPL={rgpl}",
        f"RGPU={rgpu}",
        f"RGSM={rgsm}"
    ]

    subprocess.run(cmd, check=True)


def generate_bqsr_table(
    input_bam: Path,
    reference_fasta: Path,
    known_sites_vcf: Path,
    output_table: Path
):
    """
    Runs GATK BaseRecalibrator to create a BQSR table.

    Parameters:
        input_bam (Path): BAM file with read groups.
        reference_fasta (Path): Reference genome FASTA.
        known_sites_vcf (Path): Known variant sites (VCF).
        output_table (Path): Output recalibration table.
    """
    cmd = [
        "gatk",
        "BaseRecalibrator",
        "-I", str(input_bam),
        "-R", str(reference_fasta),
        "--known-sites", str(known_sites_vcf),
        "-O", str(output_table)
    ]
    subprocess.run(cmd, check=True)


def apply_bqsr(
    input_bam: Path,
    reference_fasta: Path,
    bqsr_table: Path,
    output_bam: Path
):
    """
    Applies BQSR using a precomputed recalibration table.

    Parameters:
        input_bam (Path): BAM file with read groups.
        reference_fasta (Path): Reference genome FASTA.
        bqsr_table (Path): Recalibration table from BaseRecalibrator.
        output_bam (Path): Recalibrated output BAM file.
    """
    cmd = [
        "gatk",
        "ApplyBQSR",
        "-I", str(input_bam),
        "-R", str(reference_fasta),
        "--bqsr-recal-file", str(bqsr_table),
        "-O", str(output_bam)
    ]
    subprocess.run(cmd, check=True)
