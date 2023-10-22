


# LD annotation
# Example for T1D
# 4 significant models could be done based on the 4 input vcf datafiles in DATA/LDannotation
cd ../DATA/LDannotation
# Grep SNP list from vcf files to fulfill input requirement
cat NEER_T1D_snps_anno.vcf.gz.hg19_multianno.vcf |grep -v "^#"|cut -f 1 > ../../OUTPUT/LDannotation/NEER_T1D_SNPS1.txt
cat NEER_T1D_snps_anno.vcf.gz.hg19_multianno.vcf |grep -v "^#"|cut -f 2 > ../../OUTPUT/LDannotation/NEER_T1D_SNPS2.txt
paste ../../OUTPUT/LDannotation/NEER_T1D_SNPS1.txt ../../OUTPUT/LDannotation/NEER_T1D_SNPS2.txt -d _ > ../../OUTPUT/LDannotation/NEER_T1D_SNPS3.txt
paste ../../OUTPUT/LDannotation/NEER_T1D_SNPS1.txt ../../OUTPUT/LDannotation/NEER_T1D_SNPS2.txt ../../OUTPUT/LDannotation/NEER_T1D_SNPS3.txt > ../../OUTPUT/LDannotation/NEER_T1D_SNPS.txt


# home address: /home/sic539
# LD-annot could be found at https://github.com/ArnaudDroitLab/LD-annot

python3 ~/LDanno/LD-annot/LD-annot0.4.py NEER_T1D_snps_anno.vcf.gz.hg19_multianno.vcf gencode.v19.annotation_nochr.gff3 ../../OUTPUT/LDannotation/NEER_T1D_SNPS.txt gene 0.9 ../../OUTPUT/LDannotation/T1D_SNPs_LD_genes_anno.txt

