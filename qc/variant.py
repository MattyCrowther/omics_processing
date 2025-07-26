import subprocess
from pathlib import Path
import matplotlib.pyplot as plt
import gzip

def stats(vcf_file, output_file=None):
    """
    Run bcftools stats on a VCF file.

    Parameters:
        vcf_file (str or Path): Path to input VCF.
        output_file (str or Path, optional): Output stats file. Defaults to vcf_file.stem + ".bcftools_stats.txt".

    Returns:
        Path: Path to stats file.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    vcf_file = Path(vcf_file)
    if output_file is None:
        output_file = vcf_file.with_suffix(".bcftools_stats.txt")
    else:
        output_file = Path(output_file)

    with open(output_file, "w") as out:
        subprocess.run(["bcftools", "stats", str(vcf_file)], stdout=out, check=True)

    return output_file

def validate(vcf_path):
    """
    Run vcf-validator to check VCF format validity.

    Parameters:
        vcf_path (str or Path): Path to the VCF file.

    Returns:
        bool: True if no validation errors, False otherwise.
    """
    vcf_path = Path(vcf_path)

    try:
        subprocess.run(["vcf-validator", str(vcf_path)],
                       check=True)
        print(f"[✓] VCF validation passed: {vcf_path.name}")
        return True
    except subprocess.CalledProcessError:
        print(f"[✗] VCF validation failed: {vcf_path.name}")
        return False
    

def qual_distribution(vcf_file, output_file=None, bins=100):
    """
    Plot histogram of QUAL scores from a VCF file.

    Parameters:
        vcf_file (str or Path): VCF file to parse.
        output_file (str or Path, optional): Path to save histogram PNG.
        bins (int): Number of histogram bins.

    Returns:
        Path: Path to saved plot.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    vcf_file = Path(vcf_file)
    qual_scores = []

    opener = gzip.open if vcf_file.suffix == ".gz" else open
    with opener(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            try:
                qual = float(parts[5])
                qual_scores.append(qual)
            except ValueError:
                continue

    plt.figure()
    plt.hist(qual_scores, bins=bins)
    plt.xlabel("QUAL Score")
    plt.ylabel("Number of Variants")
    plt.title("VCF Variant Quality Score Distribution")
    plt.tight_layout()

    if output_file is None:
        output_file = vcf_file.with_name(f"{vcf_file.stem}_qual_hist.png")
    else:
        output_file = Path(output_file)

    plt.savefig(output_file)
    plt.close()
    return output_file

def count_variant_types(vcf_file, output_file=None):
    """
    Count number of SNPs and indels in a VCF file and optionally write results.

    Parameters:
        vcf_file (str or Path): Path to input VCF (.vcf or .vcf.gz).
        output_file (str or Path, optional): Path to write summary results as text.

    Returns:
        dict: {'SNPs': int, 'Indels': int}
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    vcf_file = Path(vcf_file)
    snps = 0
    indels = 0

    opener = gzip.open if vcf_file.suffix == ".gz" else open

    with opener(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            ref = fields[3]
            alt = fields[4]

            for allele in alt.split(","):
                if len(ref) == 1 and len(allele) == 1:
                    snps += 1
                else:
                    indels += 1

    result = {"SNPs": snps, "Indels": indels}

    if output_file:
        output_file = Path(output_file)
        with open(output_file, "w") as out:
            for k, v in result.items():
                out.write(f"{k}: {v}\n")
        print(f"[✓] Variant type counts written to: {output_file.name}")

    return result



