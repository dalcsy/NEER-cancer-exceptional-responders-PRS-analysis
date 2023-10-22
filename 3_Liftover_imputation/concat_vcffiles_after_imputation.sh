#!/bin/bash
#SBATCH -o imputationQC
#SBATCH -J imputationQC
#SBATCH --mem 32g
#SBATCH -c 8
#SBATCH -t 0-4:00:00
#SBATCH -p short

module load gcc/6.2.0
module load bcftools

cd ../OUTPUT/neer_imputed_vqsr_no_multi_alt_GQ20/impute_trimmed/

bcftools concat chr*.dose.vcf.gz_trim.vcf.gz -o merged.vcf.gz -O z --threads 8
bcftools sort merged.vcf.gz -O z -o neer_imputed_trim_sorted.vcf.gz


