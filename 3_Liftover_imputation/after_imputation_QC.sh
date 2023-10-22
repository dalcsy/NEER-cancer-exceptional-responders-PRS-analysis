#!/bin/bash
#SBATCH -o imputationQC
#SBATCH -J imputationQC
#SBATCH --mem 32g
#SBATCH -c 8
#SBATCH -t 0-4:00:00
#SBATCH -p short

module load gcc/6.2.0
module load bcftools

cd ../OUTPUT/neer_imputed_vqsr_no_multi_alt_GQ20/
for file in *vcf.gz;do bcftools view -i 'INFO/R2[0] >0.8'  -O z --threads 8 ${file} -o ./impute_trimmed/${file}_trim.vcf.gz; done
