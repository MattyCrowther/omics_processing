from pathlib import Path
from sequence_acquisition import get_sequence
from read_alignment import reference_genome
from read_alignment import align_reads
from quality_control import perform_qc
from genomics.variants.callers import freebayes
from qc.variant import validate
from qc.variant import stats
from qc.variant import count_variant_types
from qc.variant import qual_distribution
# READS
base_dir = Path("yeast")

'''
read_1,read_2 = get_sequence("SRR34533466",base_dir=base_dir)

# GENOME ALIGN
ref_url = ("ftp://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz")
anno_url = ("ftp://ftp.ensemblgenomes.org/pub/fungi/release-56/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.56.gtf.gz")
fa_path,anno_path = reference_genome(ref_url,anno_url,base_dir=base_dir)

final_bam = align_reads(fa_path,read_1,read_2,base_dir=base_dir)

# QC
perform_qc(final_bam,fa_path)
'''

'''
final_bam = "yeast/aln_sorted_dedup.bam"
fa_path = "yeast/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
print(final_bam,fa_path)
'''
output_vcf= base_dir/"variants_raw.vcf"
#freebayes(final_bam,fa_path,output_vcf=output_vcf)

v_qc_out = Path("yeast/qc/variants")
validate(output_vcf)
stats(output_vcf,v_qc_out/"stats.txt")
rv = count_variant_types(output_vcf,v_qc_out/"variant_type_count.txt")
qual_distribution(output_vcf,v_qc_out/"qual_distribution.png")