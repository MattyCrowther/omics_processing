import subprocess

def run_spades(read1, read2, output_dir="spades_output", threads=4, memory=16):
    """
    Assemble genome de novo using SPAdes.
    
    Parameters:
        read1 (str): Path to trimmed read 1 (FASTQ.gz)
        read2 (str): Path to trimmed read 2 (FASTQ.gz)
        output_dir (str): Output directory for SPAdes results
        threads (int): Number of threads to use
        memory (int): Max memory (in GB)
    """
    subprocess.run([
        "spades.py",
        "-1", read1,
        "-2", read2,
        "-o", output_dir,
        "--threads", str(threads),
        "--memory", str(memory)
    ])
    