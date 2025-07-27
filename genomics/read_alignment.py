from pathlib import Path
import sys
import os

sys.path.append("../")

from downloader import download_reference_genome
from alignment import alignment_index
from utils import decompress_gzip
from alignment import align_ends
from converter import SAM_to_BAM
from converter import sort_bam
from converter import index_bam
from converter import mark_duplicates


    
def reference_genome(ref_url, base_dir=Path(".")):
    """Download and prepare reference genome and annotation under base_dir/ref/."""
    ref_dir = base_dir / "ref"
    
    ref_filename = os.path.basename(ref_url)
    print(f"Downloading reference: {ref_filename}")
    gz_path = download_reference_genome(ref_url, ref_filename, target_dir=ref_dir)
    fa_path = decompress_gzip(gz_path)
    print(f"Indexing: {fa_path}")
    alignment_index(fa_path)
    return Path(fa_path)

def align_reads(fa_path, read1, read2,base_dir=Path(".")):
    sam_filename = base_dir/"aln.sam"
    sam_file = align_ends(fa_path, read1, read2, sam_filename)
    bam_file = SAM_to_BAM(sam_file)
    sorted_bam_file = sort_bam(bam_file)
    index_bam(sorted_bam_file)
    dedup_bam = mark_duplicates(sorted_bam_file)
    return dedup_bam




