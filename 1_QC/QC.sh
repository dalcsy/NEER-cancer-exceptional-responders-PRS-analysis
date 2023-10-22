#!/bin/bash
#SBATCH -o bcftools_out
#SBATCH -J annorvar_CSY
#SBATCH --mem 80g
#SBATCH -c 20
#SBATCH -t 0-12:00:00
#SBATCH -p short

module load gcc/6.2.0
module load bcftools
module load htslib
module load vcftools

## To work directory
cd /n/data1/hms/dbmi/zaklab/csy/vqsr_neer/

## after running vqsr, it only add PASS to column 7. Filter with PASS.
vcftools --gzvcf snp.recalibrated.vcf.gz --remove-filtered-all --recode --stdout | gzip -c > neer_snp_vqsr_PASS_only.vcf.gz
## index
tabix neer_snp_vqsr_PASS_only.vcf.gz

## Only include GQ>20
vcftools --gzvcf neer_snp_vqsr_PASS_only.vcf.gz --out neer_snp_vqsr_PASS_only_trim_GQ20 --recode --minGQ 20

## If one row include multiple variants, split it.
bcftools norm -m - -o neer_snp_vqsr_PASS_only_trim_GQ20.recode.no_multi_alt.vcf.gz -O z --threads 80 neer_snp_vqsr_PASS_only_trim_GQ20.recode.vcf



