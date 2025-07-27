from pathlib import Path
import subprocess


def apply_filter(input_vcf, output_vcf, expression, label=None, include=True, fn=None):
    """
    Apply a bcftools filter to a VCF file.

    Parameters:
    - input_vcf (str or Path): Input VCF path
    - output_vcf (str or Path): Output VCF path
    - expression (str): bcftools filter expression
    - label (str): FILTER tag label to assign to failing variants (optional)
    - include (bool): Use '-i' to include matching variants, or '-e' to exclude them
    - fn (str): Optional name of the calling function, for logging/debugging

    Returns:
    - Path to the filtered output VCF
    """
    input_vcf = str(input_vcf)
    output_vcf = Path(output_vcf)

    if fn:
        print(f"[{fn}] Applying filter: {'-i' if include else '-e'} \"{expression}\"")

    cmd = [
        "bcftools",
        "filter",
        "-o",
        str(output_vcf),
        "-O",
        "v",
        "-i" if include else "-e",
        expression,
    ]

    if label and not include:
        cmd += ["-s", label]

    cmd.append(input_vcf)

    subprocess.run(cmd, check=True)
    return output_vcf


def quality_and_depth(input_vcf, output_vcf=None, qual_thresh=20, dp_thresh=10):
    """
    Filter variants based on quality and depth thresholds.
    If output_vcf is not provided, appends '_qualdepth.vcf' to the input filename.
    """
    input_vcf = Path(input_vcf)
    output_vcf = (
        Path(output_vcf)
        if output_vcf
        else input_vcf.with_name(input_vcf.stem + "_qualdepth.vcf")
    )

    return apply_filter(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        expression=f"QUAL>{qual_thresh} && INFO/DP>{dp_thresh}",
        include=True,
        fn="filter_by_qual_and_depth",
    )


def label_low_quality(input_vcf, output_vcf=None, qual_thresh=20, dp_thresh=10):
    """Label low-quality variants with a FILTER tag instead of removing them."""
    input_vcf = Path(input_vcf)
    output_vcf = (
        Path(output_vcf)
        if output_vcf
        else input_vcf.with_name(input_vcf.stem + "_labeled.vcf")
    )

    return apply_filter(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        expression=f"QUAL<={qual_thresh} || INFO/DP<={dp_thresh}",
        label="LowQual",
        include=False,
        fn="label_low_quality",
    )


def max_depth_filter(input_vcf, output_vcf=None, max_dp=500):
    """Exclude variants with depth above max_dp (e.g., PCR artifacts or repeats)."""
    input_vcf = Path(input_vcf)
    output_vcf = (
        Path(output_vcf)
        if output_vcf
        else input_vcf.with_name(input_vcf.stem + f"_dp{max_dp}.vcf")
    )

    return apply_filter(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        expression=f"INFO/DP<{max_dp}",
        include=True,
        fn="add_max_depth_filter",
    )


def low_af_and_mq(input_vcf, output_vcf=None, af_thresh=0.2, mq_thresh=40):
    """Label variants with low allele frequency or low mapping quality (if present)."""
    input_vcf = Path(input_vcf)
    output_vcf = (
        Path(output_vcf)
        if output_vcf
        else input_vcf.with_name(input_vcf.stem + "_afmq.vcf")
    )

    # Check if MQ is in header
    has_mq = (
        "MQ,"
        in subprocess.run(
            ["bcftools", "view", "-h", input_vcf], capture_output=True, text=True
        ).stdout
    )

    expr = f"AF<{af_thresh}" + (f" || MQ<{mq_thresh}" if has_mq else "")
    label = "LowAF" if not has_mq else "LowAF_MQ"

    return apply_filter(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        expression=expr,
        label=label,
        include=False,
        fn="filter_low_af_and_mq",
    )


def strand_bias(input_vcf, output_vcf=None):
    """Label strand-biased variants not supported on both DNA strands."""
    input_vcf = Path(input_vcf)
    output_vcf = (
        Path(output_vcf)
        if output_vcf
        else input_vcf.with_name(input_vcf.stem + "_strand.vcf")
    )

    return apply_filter(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        expression="SAF=0 || SAR=0",
        label="StrandBias",
        include=False,
        fn="filter_strand_bias",
    )


def strict_high_confidence(
    input_vcf,
    output_vcf=None,
    qual_thresh=30,
    dp_thresh=15,
    af_thresh=0.3,
    mq_thresh=50,
):
    """
    Apply a strict high-confidence filter combining quality, depth, allele frequency,
    and optionally mapping quality (MQ, if present).
    """
    input_vcf = Path(input_vcf)
    output_vcf = (
        Path(output_vcf)
        if output_vcf
        else input_vcf.with_name(input_vcf.stem + "_strict.vcf")
    )

    # Check if MQ is present in the VCF header
    has_mq = (
        "MQ,"
        in subprocess.run(
            ["bcftools", "view", "-h", input_vcf], capture_output=True, text=True
        ).stdout
    )

    expr = f"QUAL<{qual_thresh} || INFO/DP<{dp_thresh} || AF<{af_thresh}"
    if has_mq:
        expr += f" || MQ<{mq_thresh}"

    return apply_filter(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        expression=expr,
        label="Strict",
        include=False,
        fn="filter_strict_high_confidence",
    )


def sample_coverage(input_vcf, output_vcf=None, sample_idx=0, min_dp=10):
    """Filter variants with low depth in a specific sample (e.g. sample 0)."""
    input_vcf = Path(input_vcf)
    output_vcf = (
        Path(output_vcf)
        if output_vcf
        else input_vcf.with_name(input_vcf.stem + f"_sample{sample_idx}_dp.vcf")
    )

    return apply_filter(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        expression=f"FORMAT/DP[{sample_idx}]<{min_dp}",
        label="LowSampleDepth",
        include=False,
        fn="filter_sample_coverage",
    )


def separate_snps(input_vcf, output_vcf):
    """Extract only SNPs."""
    subprocess.run(
        ["bcftools", "view", "-v", "snps", input_vcf, "-o", output_vcf, "-O", "v"],
        check=True,
    )


def separate_indels(input_vcf, output_vcf):
    """Extract only indels."""
    subprocess.run(
        ["bcftools", "view", "-v", "indels", input_vcf, "-o", output_vcf, "-O", "v"],
        check=True,
    )
