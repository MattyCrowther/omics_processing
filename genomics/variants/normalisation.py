import subprocess
from pathlib import Path


def normalize(
    input_vcf,
    reference_fa,
    output_vcf=None,
    split_multiallelics=True,
    check_ref=True,
    validate_only=False,
):
    """
    Normalize or validate variants in a VCF file using bcftools norm:
    - Left-align indels
    - Optionally split multiallelics
    - Optionally check REF allele against reference
    - Optionally run in validation-only mode (no left-align or split)

    Parameters:
    - input_vcf (str or Path): Input VCF path (.vcf or .vcf.gz)
    - reference_fa (str or Path): Reference genome in FASTA format
    - output_vcf (str or Path): Output path (.vcf.gz or .vcf); if None, auto-generated
    - split_multiallelics (bool): If True, split into biallelics (unless validate_only)
    - check_ref (bool): If True, enforce REF allele check against FASTA
    - validate_only (bool): If True, only validate REF/ALT against reference without modifying VCF

    Returns:
    - Path to output VCF (validated or normalized)
    """
    input_vcf = Path(input_vcf)
    output_vcf = (
        Path(output_vcf)
        if output_vcf
        else input_vcf.with_name(
            input_vcf.stem + ("_validated.vcf" if validate_only else "_norm.vcf.gz")
        )
    )

    print(
        f"[normalize] {'Validating' if validate_only else 'Normalizing'}: "
        f"{input_vcf.name} → {output_vcf.name}"
    )

    cmd = [
        "bcftools", "norm",
        "-f", str(reference_fa),
        "-o", str(output_vcf),
        "-O", "v" if validate_only else "z",  # uncompressed VCF if just validating
    ]

    if not validate_only and split_multiallelics:
        cmd += ["-m", "-any"]

    if check_ref:
        cmd += ["-c", "s"]

    cmd.append(str(input_vcf))

    subprocess.run(cmd, check=True)
    return output_vcf


def index(vcf_gz):
    """
    Index a bgzipped VCF using tabix.

    Parameters:
    - vcf_gz (str or Path): Path to .vcf.gz file

    Returns:
    - Path to created index file (.tbi)
    """
    vcf_gz = Path(vcf_gz)
    print(f"[index] Indexing VCF: {vcf_gz.name}")
    subprocess.run(["tabix", "-p", "vcf", str(vcf_gz)], check=True)
    return Path(str(vcf_gz) + ".tbi")