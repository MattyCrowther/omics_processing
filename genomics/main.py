from pathlib import Path
from sequence_acquisition import get_sequence
from read_alignment import reference_genome
from read_alignment import align_reads
from quality_control import perform_qc
from genomics.variants.callers import freebayes
from genomics.variants.qc import validate
from genomics.variants.qc import stats
from genomics.variants.qc import count_variant_types
from genomics.variants.qc import qual_distribution
from genomics.variants.qc import qual_distribution

from genomics.variants.filters import quality_and_depth
from genomics.variants.filters import low_af_and_mq
from genomics.variants.filters import strand_bias
from genomics.variants.filters import max_depth_filter
from genomics.variants.filters import strict_high_confidence
from genomics.variants.filters import sample_coverage

from genomics.variants.normalisation import normalize
from genomics.variants.normalisation import index

from genomics.variants.annotate import tidy_fields
from genomics.variants.annotate import sort

# READS
base_dir = Path("yeast")


read_1,read_2 = get_sequence("SRR34533466",base_dir=base_dir)

# GENOME ALIGN
ref_url = ("ftp://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz")
fa_path = reference_genome(ref_url,base_dir=base_dir)
final_bam = align_reads(fa_path,read_1,read_2,base_dir=base_dir)

# QC
#perform_qc(final_bam,fa_path)
vcf_dir = base_dir / "vcf"
output_vcf = vcf_dir / "variants_raw.vcf"
v_qc_out = vcf_dir / "qc" / "variants"

freebayes(final_bam, fa_path, output_vcf=output_vcf)
# Filter by QUAL and INFO/DP 
vcf1 = quality_and_depth(output_vcf)
vcf2 = low_af_and_mq(vcf1,af_thresh=0.8)
vcf3 = strand_bias(vcf2)
norm_vcf = normalize(vcf3, fa_path)
tidied_vcf = tidy_fields(norm_vcf)
sorted_vcf = sort(tidied_vcf)
index(sorted_vcf)
final_vcf = sorted_vcf
# Perform QC
validate(final_vcf)
stats(final_vcf, v_qc_out / "stats.txt")
rv = count_variant_types(final_vcf, v_qc_out / "variant_type_count.txt")
qual_distribution(final_vcf, v_qc_out / "qual_distribution.png")


intermediate_files = [
    vcf1, vcf2, vcf3, norm_vcf, tidied_vcf, sorted_vcf,
    vcf1.with_suffix(".vcf.csi"), vcf1.with_suffix(".vcf.idx")
]
for f in intermediate_files:
    if Path(f).exists():
        Path(f).unlink()