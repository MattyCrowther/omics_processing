import subprocess
from pathlib import Path

def freebayes(bam_path, reference_fasta, output_vcf, extra_args=None):
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

    with output_vcf.open("w") as out_vcf:
        subprocess.run(cmd, stdout=out_vcf, check=True)
