import subprocess
from pathlib import Path
from typing import List


def freebayes(bam_path, reference_fasta, 
              output_vcf, extra_args=None):
    """
    Run FreeBayes to call variants from a BAM file.

    Parameters:
    - bam_path (str or Path): Path to the deduplicated, sorted BAM file
    - reference_fasta (str or Path): Path to the reference genome in FASTA format
    - output_vcf (str or Path): Output path for the raw VCF file
    - extra_args (list): Optional list of additional arguments to pass to FreeBayes
    """
    bam_path = str(bam_path)
    reference_fasta = str(reference_fasta)
    output_vcf = Path(output_vcf)

    cmd = [
        "freebayes",
        "-f", reference_fasta,
        bam_path
    ]
    if extra_args:
        cmd.extend(extra_args)
    
    output_vcf.parent.mkdir(exist_ok=True)
    with output_vcf.open("w") as out_vcf:
        subprocess.run(cmd, stdout=out_vcf, check=True)


def haplotype_caller(bam_path: str, reference_fasta: str, 
                     output_gvcf: str):
    """
    Generate a gVCF file for a single sample using GATK HaplotypeCaller.
    """
    subprocess.run([
        "gatk", "HaplotypeCaller",
        "-R", reference_fasta,
        "-I", bam_path,
        "-O", output_gvcf,
        "-ERC", "GVCF"
    ], check=True)


def combine_gvcfs(reference_fasta: str, gvcf_paths: List[str], 
                  output_path: str):
    """
    Combine multiple gVCFs into a single gVCF (for ≤50 samples).
    """
    cmd = [
        "gatk", "CombineGVCFs",
        "-R", reference_fasta,
    ]
    for path in gvcf_paths:
        cmd += ["--variant", path]
    cmd += ["-O", output_path]
    subprocess.run(cmd, check=True)


def import_gvcfs_to_db(samples_list_file: str, 
                       intervals_file: str, db_path: str, 
                       threads: int = 4):
    """
    Import gVCFs into a GenomicsDB for joint genotyping in large cohorts.
    """
    subprocess.run([
        "gatk", "GenomicsDBImport",
        "--genomicsdb-workspace-path", db_path,
        "--sample-name-map", samples_list_file,
        "--L", intervals_file,
        "--reader-threads", str(threads)
    ], check=True)


def genotype_gvcfs(reference_fasta: str, input_vcf_or_db: str, 
                   output_vcf: str, is_db: bool = False):
    """
    Perform joint genotyping on combined gVCF or a GenomicsDBImport database.
    """
    input_path = f"gendb://{input_vcf_or_db}" if is_db else input_vcf_or_db
    subprocess.run([
        "gatk", "GenotypeGVCFs",
        "-R", reference_fasta,
        "-V", input_path,
        "-O", output_vcf
    ], check=True)


def run_manta(bam_file, reference_fa, output_dir="manta_sv"):
    """Configure and run Manta structural variant caller."""
    output_dir = Path(output_dir)
    
    subprocess.run([
        "configManta.py",
        "--bam", str(bam_file),
        "--referenceFasta", str(reference_fa),
        "--runDir", str(output_dir)
    ], check=True)

    subprocess.run([
        str(output_dir / "runWorkflow.py"),
        "-m", "local",
        "-j", "4"
    ], check=True)

    return output_dir / "results" / "variants" / "diploidSV.vcf.gz"
