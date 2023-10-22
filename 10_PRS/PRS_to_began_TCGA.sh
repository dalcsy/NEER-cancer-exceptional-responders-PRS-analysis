#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=64G
#SBATCH --job-name=begin_prs_csy
#SBATCH -o begin_prs_csy.out
#SBATCH --mail-user=siyuanchen@hsph.harvard.edu
#SBATCH --mail-type=FAIL
#SBATCH -a 1-22

ln -s /n/app/openblas/0.2.19/lib/liblapack.so liblapack.so.3
export LD_LIBRARY_PATH=$PWD
module load gcc/6.2.0
module load bcftools
module load htslib

# Remove InDels (we don't use INDELs)
bcftools view -v snps -O z -o chr$SLURM_ARRAY_TASK_ID.snps.vcf.gz /n/data1/hms/dbmi/zaklab/csy/TCGA_chr_split/imputation_input/${SLURM_ARRAY_TASK_ID}.vcf.gz

# Split multi-allelic sites into biallelic sites
bcftools norm -m - chr$SLURM_ARRAY_TASK_ID.snps.vcf.gz|gzip -c > chr$SLURM_ARRAY_TASK_ID.snps.biallelic.vcf.gz

# If your data are phased, remove all phasing information. Qctool PRS
# calculation has a compatibility issue with phased genome.
# 0|0 => 0/0; 0|1 => 0/1; 1|0 => 1/0; 1|1 => 1/1

zcat chr$SLURM_ARRAY_TASK_ID.snps.biallelic.vcf.gz | sed 's/0|0/0\/0/g' | sed 's/0|1/0\/1/g' | sed 's/1|0/1\/0/g' | sed 's/1|1/1\/1/g' | gzip -c > chr$SLURM_ARRAY_TASK_ID.snps.biallelic.unphased.vcf.gz

# vcf to bgen
/n/data1/hms/dbmi/zaklab/csy/PRS/NPS/qctools/qctool -g chr$SLURM_ARRAY_TASK_ID.snps.biallelic.unphased.vcf.gz -og chrom$SLURM_ARRAY_TASK_ID.neer.bgen -os chrom$SLURM_ARRAY_TASK_ID.tcga.sample -vcf-genotype-field GT


