from pathlib import Path
import subprocess

def SAM_to_BAM(sam_file):
    bam_file = Path(sam_file).with_suffix(".bam")
    print(f"Converting SAM to BAM: {sam_file} → {bam_file}")

    if not Path(sam_file).exists():
        raise FileNotFoundError(f"SAM file not found: {sam_file}")
    if Path(sam_file).stat().st_size == 0:
        raise ValueError(f"SAM file is empty: {sam_file}")

    subprocess.run(["samtools", "view", "-b", str(sam_file), "-o", str(bam_file)], check=True)

    if not bam_file.exists():
        raise FileNotFoundError(f"BAM file not created: {bam_file}")

    # cleanup
    Path(sam_file).unlink()
    print(f"Deleted intermediate SAM: {sam_file}")

    return str(bam_file)


def sort_bam(bam_file):
    sorted_bam_file = Path(bam_file).with_name(Path(bam_file).stem + "_sorted.bam")
    print(f"Sorting BAM: {bam_file} → {sorted_bam_file}")
    subprocess.run(["samtools", "sort", "-o", str(sorted_bam_file), str(bam_file)], check=True)

    # cleanup
    Path(bam_file).unlink()
    print(f"Deleted unsorted BAM: {bam_file}")

    return str(sorted_bam_file)


def index_bam(bam_file):
    print(f"Indexing BAM: {bam_file}")
    subprocess.run(["samtools", "index", str(bam_file)], check=True)

def mark_duplicates(input_bam):
    input_bam = Path(input_bam)
    name_sorted = input_bam.with_name(input_bam.stem + "_namesort.bam")
    fixmate_bam = input_bam.with_name(input_bam.stem + "_fixmate.bam")
    coord_sorted = input_bam.with_name(input_bam.stem + "_coord.bam")
    dedup_bam = input_bam.with_name(input_bam.stem + "_dedup.bam")
    coord_bai = coord_sorted.with_suffix(".bam.bai")  # <- .bai of coord-sorted BAM

    print(f"Name-sorting BAM: {input_bam} → {name_sorted}")
    subprocess.run(["samtools", "sort", "-n", "-o", str(name_sorted), str(input_bam)], check=True)

    print(f"Fixing mates: {name_sorted} → {fixmate_bam}")
    subprocess.run(["samtools", "fixmate", "-m", str(name_sorted), str(fixmate_bam)], check=True)

    print(f"Coordinate-sorting fixed BAM: {fixmate_bam} → {coord_sorted}")
    subprocess.run(["samtools", "sort", "-o", str(coord_sorted), str(fixmate_bam)], check=True)

    print(f"Marking duplicates: {coord_sorted} → {dedup_bam}")
    subprocess.run(["samtools", "markdup", "-r", str(coord_sorted), str(dedup_bam)], check=True)

    print(f"Indexing final BAM: {dedup_bam}")
    subprocess.run(["samtools", "index", str(dedup_bam)], check=True)

    # Cleanup
    for f in [name_sorted, fixmate_bam, coord_sorted, input_bam,coord_bai]:
        if f.exists():
            f.unlink()
            print(f"Deleted intermediate: {f}")

    return str(dedup_bam)