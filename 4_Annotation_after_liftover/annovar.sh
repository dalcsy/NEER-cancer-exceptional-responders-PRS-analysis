#!/bin/bash
#SBATCH -o neer_annovar_out
#SBATCH -J neer_annorvar_CSY
#SBATCH --mem 32g
#SBATCH -c 2
#SBATCH -t 0-12:00:00
#SBATCH -p short

module load gcc/6.2.0
module load annovar/20170601
module load htslib


# The result file from this step is not the exonic variants we want. The GENCODE database include many non-coding transcripts and some to-be-validated genes.
# The main purpose is to cut the size of dataset.
bcftools view -R /n/data1/hms/dbmi/zaklab/csy/ref/genecode_annotation_files/hg19/gencode.v19.annotation_exon_splice_sort_noverlap_nochr.bed -O z --threads 24 -o neer_imputed_trim_sorted_GQ20_exonic.vcf.gz neer_imputed_trim_sorted_GQ20.vcf.gz

# Annotation
perl /n/app/annovar/20170601/table_annovar.pl --remove --protocol refGene,avsnp150,exac03,gnomad211_genome,dbnsfp41a,clinvar_20210123 --operation g,f,f,f,f,f --nastring . -vcfinput neer_imputed_trim_sorted_GQ20_exonic.vcf.gz /n/data1/hms/dbmi/zaklab/csy/humandb/hg19/ -buildver hg19 -out neer_imputed_trim_sorted_GQ20_exonic_anno_vcf



