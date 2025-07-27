import subprocess
from pathlib import Path

def assign_rsid(input_vcf, dbsnp_vcf, output_vcf=None):
    """
    Annotates VCF records with rsIDs from a known dbSNP-style reference.

    Parameters:
        input_vcf (str or Path): Path to the normalized input VCF (e.g., variants_norm.vcf.gz)
        dbsnp_vcf (str or Path): Path to a VCF file with known variant IDs (e.g., dbSNP or Ensembl VCF)
        output_vcf (str or Path, optional): Output path for the annotated VCF.
                                            If None, a default name is generated.
    """
    input_vcf = Path(input_vcf)
    dbsnp_vcf = Path(dbsnp_vcf)

    if output_vcf is None:
        output_vcf = input_vcf.with_name(input_vcf.stem + "_rsid.vcf")
    else:
        output_vcf = Path(output_vcf)

    print(f"[assign_rsid] Annotating with rsIDs using: {dbsnp_vcf}")
    
    subprocess.run([
        "bcftools", "annotate",
        "-a", str(dbsnp_vcf),
        "-c", "ID",
        "-o", str(output_vcf),
        "-O", "v",
        str(input_vcf)
    ], check=True)

    return output_vcf

def tidy_fields(input_vcf, fields_to_remove=None, output_vcf=None):
    """
    Removes specified INFO and FORMAT fields from a VCF file using bcftools annotate.

    Parameters:
        input_vcf (str or Path): Path to the VCF file to clean (e.g., variants_validated.vcf)
        fields_to_remove (list of str): List of field names to remove (e.g., ["INFO/OLD_TAG", "FORMAT/UNUSED_TAG"])
        output_vcf (str or Path, optional): Output path for the cleaned VCF.
                                            If None, appends '_tidy.vcf' to input filename.

    Returns:
        Path to the cleaned VCF file
    """
    input_vcf = Path(input_vcf)
    fields_to_remove = fields_to_remove or []
    output_vcf = Path(output_vcf) if output_vcf else input_vcf.with_name(input_vcf.stem + "_tidy.vcf")

    if not fields_to_remove:
        print("[tidy_vcf_fields] No fields specified for removal. Skipping.")
        return input_vcf

    remove_arg = ",".join(fields_to_remove)
    print(f"[tidy_vcf_fields] Removing fields: {remove_arg}")

    subprocess.run([
        "bcftools", "annotate",
        "--remove", remove_arg,
        "-o", str(output_vcf),
        "-O", "v",
        str(input_vcf)
    ], check=True)

    return output_vcf


def sort(input_vcf, output_vcf=None):
    """
    Sorts a VCF by chromosomal position, compresses with bgzip, and indexes it with tabix.

    Parameters:
        input_vcf (str or Path): Path to the input VCF (.vcf or .vcf.gz)
        output_vcf (str or Path, optional): Output path for the sorted/compressed VCF (.vcf.gz).
                                            If None, generates '<stem>_final.vcf.gz'

    Returns:
        Path to the compressed and indexed VCF file
    """
    input_vcf = Path(input_vcf)
    if output_vcf is None:
        base = input_vcf
        if base.suffix == ".gz" and base.with_suffix("").suffix == ".vcf":
            base = base.with_suffix("")  # Remove .gz
            base = base.with_suffix("")  # Remove .vcf
        else:
            base = base.with_suffix("")  # Remove single suffix
        output_vcf = base.with_name(base.name + "_sorted.vcf.gz")
    
    print(f"[sort_and_index_vcf] Sorting and compressing: {input_vcf.name} → {output_vcf.name}")

    subprocess.run([
        "bcftools", "sort",
        str(input_vcf),
        "-Oz",
        "-o", str(output_vcf)
    ], check=True)

    return output_vcf
