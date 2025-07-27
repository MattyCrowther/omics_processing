import subprocess
def run_quast(contigs_fasta, output_dir="quast_output", reference_fasta=None, reference_gff=None):
    """
    Evaluate genome assembly quality using QUAST.
    
    Parameters:
        contigs_fasta (str): Path to contigs FASTA file
        output_dir (str): Output directory for QUAST report
        reference_fasta (str): Optional reference genome FASTA
        reference_gff (str): Optional reference annotation GFF
    """
    cmd = [
        "quast",
        contigs_fasta,
        "-o", output_dir
    ]
    
    if reference_fasta:
        cmd += ["-r", reference_fasta]
    if reference_gff:
        cmd += ["-g", reference_gff]
    
    subprocess.run(cmd)