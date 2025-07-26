import subprocess
from pysradb.sraweb import SRAweb
from pathlib import Path
import urllib.request 

def sra_metadata(identifier):
    """Query and print study and experiment metadata for a given SRA run."""
    db = SRAweb()
    df = db.sra_metadata(identifier, detailed=True)
    print(df[['study_accession', 'experiment_accession']])
    return df

def download_sra(identifier, raw_dir):
    """Download .sra file and extract compressed FASTQ files to raw_dir."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nDownloading data for {identifier}...")
    subprocess.run(["prefetch", identifier], check=True)
    subprocess.run(["fasterq-dump", identifier, "--split-files", "-O", str(raw_dir)], check=True)
    print("FASTQ download and extraction complete.")

def download_reference_genome(url: str, filename: str, target_dir="ref") -> Path:
    """Download the reference genome FASTA to a target directory."""
    ref_dir = Path(target_dir)
    ref_dir.mkdir(exist_ok=True)
    dest_path = ref_dir / filename
    if not dest_path.exists():
        print(f"Downloading: {filename}")
        urllib.request.urlretrieve(url, dest_path)
    else:
        print(f"Found existing reference genome: {dest_path}")
    return dest_path