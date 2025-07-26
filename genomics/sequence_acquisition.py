from pathlib import Path
import sys

sys.path.append("../")

from downloader import download_sra
from downloader import sra_metadata
from qc.alignment import fast_qc
from qc.alignment import multi_qc
from read_trim import cut_adapt


RAW_DIR = Path("raw")
TRIMMED_DIR = Path("trimmed")
QC_DIR = Path("qc")

def detect_fastq_files(identifier, directory):
    """Return list of FASTQ files for a given identifier in the target directory."""
    fastq_files = list(directory.glob(f"{identifier}*.fastq*"))
    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found for {identifier} in {directory}")
    return sorted(fastq_files)

def get_sequence(identifier, base_dir=Path("."), do_trimming=True):
    """Sequence acquisition for a given SRA accession."""
    print(f"\nStarting processing for: {identifier}")
    
    raw_dir = base_dir / RAW_DIR
    trimmed_dir = base_dir / TRIMMED_DIR
    qc_dir = base_dir / QC_DIR

    sra_metadata(identifier)
    download_sra(identifier, raw_dir)
    raw_fastqs = detect_fastq_files(identifier, raw_dir)

    fast_qc(raw_fastqs, qc_dir)

    if do_trimming:
        trimmed_fastqs = cut_adapt(raw_fastqs, trimmed_dir)
        fast_qc(trimmed_fastqs, qc_dir)

    multi_qc(qc_dir)
    
    if do_trimming:
        return trimmed_fastqs
    return raw_fastqs
    

