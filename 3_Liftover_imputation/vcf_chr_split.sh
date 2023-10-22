#!/bin/bash
#SBATCH -o chr_split.out
#SBATCH -J chr_split_CSY
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -p short
#SBATCH -t 12:00:00

module load gcc/6.2.0
module load htslib/1.10.2

tabix /n/data1/hms/dbmi/zaklab/csy/vqsr_neer/neer_snp_vqsr_PASS_only_trim_GQ20.recode.no_multi_alt.vcf.gz

VCFGZ="/n/data1/hms/dbmi/zaklab/csy/vqsr_neer/neer_snp_vqsr_PASS_only_trim_GQ20.recode.no_multi_alt.vcf.gz"

#tabix --list-chroms $VCFGZ > chromosomes.txt  # save all the chromosome names into a file

while IFS= read -r line; do
  tabix -h $VCFGZ $line > $line.vcf;
  bgzip -@ 8 $line.vcf;
  rm $line.vcf;
  tabix $line.vcf.gz
done < chromosomes.txt  # make an individual vcf for each chromosome

## The output file is in ./OUTPUT/neer_split_no_multi_alt/ ##

## MOVE ALL VCF FILES OF AUTOSOMAL CHROMOSOMES INTO ./OUTPUT/neer_split_no_multi_alt/imputation_input/ ##


