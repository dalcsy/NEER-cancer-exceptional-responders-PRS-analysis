#!/bin/bash
#SBATCH -o icgc_neer_annovar_out
#SBATCH -J neer_annorvar_CSY
#SBATCH --mem 32g
#SBATCH -c 2
#SBATCH -t 0-12:00:00
#SBATCH -p short

module load gcc/6.2.0
module load annovar/20170601
module load bcftools


perl /n/app/annovar/20170601/table_annovar.pl --vcfinput --remove --protocol refGene,avsnp150,exac03,gnomad30_genome,dbnsfp41a,clinvar_20210123  --operation g,f,f,f,f,f --nastring . neer_snp_vqsr_PASS_only_trim_GQ20.recode.no_multi_alt.vcf.gz /n/data1/hms/dbmi/zaklab/csy/humandb/ -buildver hg38 -out neer_snp_vqsr_PASS_only_trim_GQ20.recode.no_multi_alt_anno


