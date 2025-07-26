import subprocess
from pathlib import Path

def alignment_index(fasta_path: Path):
    """Run BWA-MEM2 index on a FASTA file if index files are missing."""
    fasta_path = Path(fasta_path)
    expected_extensions = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".sa"]
    index_files = [fasta_path.with_suffix(fasta_path.suffix + ext) for ext in expected_extensions]

    if all(f.exists() for f in index_files):
        print(f"BWA-MEM2 index already exists for: {fasta_path}")
        return

    print(f"Indexing reference: {fasta_path}")
    subprocess.run(["bwa-mem2", "index", str(fasta_path)], check=True)


def align_ends(reference: str, sequence1: str, sequence2: str, out_filename=None) -> str:
    sam_file = Path(out_filename) if out_filename else Path("aln.sam")
    print(f"Aligning reads → {sam_file}")
    with sam_file.open("w") as out:
        subprocess.run([
            "bwa-mem2", "mem", "-t", "4", reference,
            sequence1, sequence2
        ], stdout=out, check=True)
    return str(sam_file)
