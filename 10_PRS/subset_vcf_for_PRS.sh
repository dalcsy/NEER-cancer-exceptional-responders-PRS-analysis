#!/bin/bash
#SBATCH -o vcf_subset
#SBATCH -J vcf_subset
#SBATCH --mem 32g
#SBATCH -c 8
#SBATCH -t 0-6:00:00
#SBATCH -p short


module load gcc/6.2.0
module load bcftools



# This is just an example for the IBD model. It was done for all models except the IBD model. Bed files were generated in the R file.
# The input files were also in the DATA directory
bcftools view -R T1D_PRS_weight.bed --threads 24 -O z -o PCAWG_T1D_snps.vcf.gz /n/data1/hms/dbmi/zaklab/csy/PCAWG/combined_germline_files/pcawg_hs37d5/pcawg8.snps.indels.svs.phased.merged.v2.controlled.hs37d5.no_multi_alt.vcf.gz
bcftools view -R T1D_PRS_weight.bed --threads 24 -O z -o NEER_T1D_snps.vcf.gz /n/data1/hms/dbmi/zaklab/csy/neer_imputed/vqsr_no_multi_alt_GQ20/impute_trimmed/neer_imputed_trim_sorted_GQ20.vcf.gz

