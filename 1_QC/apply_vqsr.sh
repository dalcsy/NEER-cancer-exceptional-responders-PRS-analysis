#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=64G
#SBATCH --job-name=vqsr_csy
#SBATCH -o vqsr_csy.out
#SBATCH --mail-user=siyuanchen@hsph.harvard.edu
#SBATCH --mail-type=FAIL

module load gcc/6.2.0
module load bcftools htslib python/2.7.12-ucs4 gatk

VCFPATH=/n/data1/hms/dbmi/zaklab/csy/vqsr_neer
GATKPATH=/n/data1/hms/dbmi/zaklab/gatk
VCFFILE=/n/data1/hms/dbmi/zaklab/NEER_VCF/unannotated.vcf.gz

# create the ExcessHet filter
#gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
#    -V $VCFFILE \
#    --filter-expression "ExcessHet > 54.69" \
#    --filter-name ExcessHet \
#    -O $VCFPATH/neer_excesshet.vcf.gz 

# remove per-sample columns
#gatk MakeSitesOnlyVcf \
#        -I $VCFPATH/neer_excesshet.vcf.gz \
#        -O $VCFPATH/neer_sitesonly.vcf.gz

# calculate VQSLOD tranches for indels
#gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
#    -V $VCFPATH/neer_sitesonly.vcf.gz \
#    --trust-all-polymorphic \
#    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
#    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL \
#    --max-gaussians 4 \
#    -resource:mills,known=false,training=true,truth=true,prior=12 $GATKPATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
#    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 $GATKPATH/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
#    -resource:dbsnp,known=true,training=false,truth=false,prior=2 $GATKPATH/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
#    -O $VCFPATH/neer_indels.recal \
#    --tranches-file $VCFPATH/neer_indels.tranches

# calculate VQSLOD tranches for SNPs
gatk --java-options "-Xmx32g -Xms32g" VariantRecalibrator \
    -V $VCFPATH/neer_sitesonly.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -mode SNP \
    --max-gaussians 6 \
    -resource:hapmap,known=false,training=true,truth=true,prior=15 $GATKPATH/hapmap_3.3.hg38.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12 $GATKPATH/1000G_omni2.5.hg38.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10 $GATKPATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7 $GATKPATH/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -O $VCFPATH/neer_snps.recal \
    --tranches-file $VCFPATH/neer_snps.tranches

gatk IndexFeatureFile \
      -I neer_snps.recal



# fill in the FILTER coluumn for indels
#gatk --java-options "-Xmx5g -Xms5g" \
#    ApplyVQSR \
#    -V $VCFPATH/neer_excesshet.vcf.gz \
#    --recal-file $VCFPATH/neer_indels.recal \
#    --tranches-file $VCFPATH/neer_indels.tranches \
#    --truth-sensitivity-filter-level 99.7 \
#    --create-output-variant-index true \
#    -mode INDEL \
#    -O $VCFPATH/indel.recalibrated.vcf.gz

# fill in the FILTER column for SNPs
gatk --java-options "-Xmx32g -Xms32g" \
    ApplyVQSR \
    -V $VCFPATH/indel.recalibrated.vcf.gz \
    --recal-file $VCFPATH/neer_snps.recal \
    --tranches-file $VCFPATH/neer_snps.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode SNP \
    -O $VCFPATH/snp.recalibrated.vcf.gz



