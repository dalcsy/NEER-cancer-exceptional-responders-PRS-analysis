setwd("/n/data1/hms/dbmi/zaklab/csy/NEER_germline_project/12_PRS_GSEA/SCRIPTS")


### Aggregate SNPs to genes ###
### Gene level score = sum(|beta| * N), where N = number of alleles in NEER ###
calc_variant_importance <- function(GT_df){
  
  # Calculate importance of variant when ref allele being the effect allele
  
  # If we are calculating hypothyroidism PRS, use GT_df$effect
  # If we are calculating Psoriasis PRS, use GT_df$effect_allele
  
  GT_df$REF_effect <- GT_df$REF==GT_df$effect
  GT_df$ALT_effect <- GT_df$ALT==GT_df$effect
  
  
  
  GT_df_ref <- GT_df[GT_df$REF_effect==T,]
  if(nrow(GT_df_ref)!=0){
    GT_df_ref[is.na(GT_df_ref)] <- 0
    GT_df_ref[GT_df_ref=="0|0"|GT_df_ref=="0/0"] <- 2
    GT_df_ref[GT_df_ref=="0|1"|GT_df_ref=="0/1"|GT_df_ref=="1|0"|GT_df_ref=="1/0"] <- 1
    GT_df_ref[GT_df_ref=="1|1"|GT_df_ref=="1/1"] <- 0
    
    GT_df_ref_num <- as.data.frame(lapply(GT_df_ref[,grep("er",colnames(GT_df_ref))],as.numeric))
    
    # If we are calculating hypothyroidism PRS, use GT_df_ref$b.cond
    # If we are calculating IBD PRS, use GT_df_alt$additive_beta
    # If we are calculating Psoriasis PRS, use GT_df_alt$effect_weight
    GT_df_ref$importance <- abs(GT_df_ref$b.cond) *  rowSums(as.matrix(GT_df_ref_num))
  }
  print("REF allele importance was calculated.")
  
  # Calculate importance of variant when alt allele being the effect allele
  
  
  GT_df_alt <- GT_df[GT_df$ALT_effect==T,]
  if(nrow(GT_df_alt)!=0){
    GT_df_alt[GT_df_alt=="0|0"|GT_df_alt=="0/0"] <- 0
    GT_df_alt[GT_df_alt=="0|1"|GT_df_alt=="0/1"|GT_df_alt=="1|0"|GT_df_alt=="1/0"] <- 1
    GT_df_alt[GT_df_alt=="1|1"|GT_df_alt=="1/1"] <- 2
    
    GT_df_alt_num <- as.data.frame(lapply(GT_df_alt[,grep("_sample",colnames(GT_df_alt))],as.numeric))
    # If we are calculating hypothyroidism PRS, use GT_df_alt$b.cond
    # If we are calculating IBD PRS, use GT_df_alt$additive_beta
    # If we are calculating Psoriasis PRS, use GT_df_alt$effect_weight
    
    GT_df_alt$importance <- abs(GT_df_alt$b.cond) *  rowSums(as.matrix(GT_df_alt_num))
  }
  
  print("ALT allele importance was calculated.")
  
  if(nrow(GT_df_ref)==0){
    return(GT_df_alt)
  }else if(nrow(GT_df_alt)==0){
    return(GT_df_ref)
  }else{
    
    return(rbind(GT_df_ref,GT_df_alt))
  }
}

calc_gene_importance <- function(variant_importance_df){
  library(dplyr)
  gene_importance <- variant_importance_df %>% group_by(annotation)%>% summarise(gene_importance = sum(importance))
  gene_importance$annotation <- str_split(gene_importance$annotation,pattern = ";")
  gene_importance$gene_name <- unlist(lapply(gene_importance$annotation,function(x){x[[6]]}))
  gene_importance$gene_id <- unlist(lapply(gene_importance$annotation,function(x){x[[2]]}))
  gene_importance$gene_name <- str_replace(gene_importance$gene_name,"gene_name=","")
  gene_importance$gene_id <- substr(gene_importance$gene_id,9,23)
  
  return(gene_importance[,2:4])
}


### GSEA data preparation ###
### GET GENE IMPORTANCE FILE ###

# Hypothyroidism gene importance
hypoT_weight_matrix <- read.delim("../DATA/GSEA/HypoT/NEER_PRS_df.txt",stringsAsFactors = F)
hypoT_gene_annotation <- read.delim("../OUTPUT/LDannotation/HypoT_SNPs_LD_genes_anno.txt",header=T,stringsAsFactors = F)

# get variant position column
library(stringr)
hypoT_weight_matrix$variant_pos <- paste(hypoT_weight_matrix$hg19.chrom,hypoT_weight_matrix$hg19.pos,sep="_")
hypoT_gene_annotation$variant_pos <- str_split(hypoT_gene_annotation$SNP,pattern = "_")
hypoT_gene_annotation$variant_pos <- lapply(hypoT_gene_annotation$variant_pos,function(x){paste(x[[1]],x[[2]],sep="_")})

# merge data
hypoT <- merge(hypoT_weight_matrix,hypoT_gene_annotation,by="variant_pos")

hypoT_variant_importance <- calc_variant_importance(hypoT)
hypoT_gene_importance <- calc_gene_importance(hypoT_variant_importance)


write.table(hypoT_gene_importance,"../OUTPUT/GSEA/hypoT_gene_importance.txt",sep="\t",row.names=F,col.names=T,quote=F)


# IBD gene importance
IBD_weight_matrix <- read.delim("../DATA/GSEA/IBD/IBD_PRS_df.txt",stringsAsFactors = F)
IBD_gene_annotation <- read.delim("../OUTPUT/LDannotation/IBD_SNPs_LD_genes_anno.txt",header=T,stringsAsFactors = F)

# get variant position column
library(stringr)
IBD_gene_annotation$variant_pos <- str_split(IBD_gene_annotation$SNP,pattern = "_")
IBD_gene_annotation$variant_pos <- lapply(IBD_gene_annotation$variant_pos,function(x){paste(x[[1]],x[[2]],sep="_")})


# merge data
IBD <- merge(IBD_weight_matrix,IBD_gene_annotation,by.x="chr_pos",by.y="variant_pos")
IBD$effect <- IBD$ALT
IBD_variant_importance <- calc_variant_importance(IBD)
IBD_gene_importance <- calc_gene_importance(IBD_variant_importance)


write.table(IBD_gene_importance,"../OUTPUT/GSEA/IBD_gene_importance.txt",sep="\t",row.names=F,col.names=T,quote=F)


# Psoriasis gene importance
Psoriasis_weight_matrix <- read.delim("../DATA/GSEA/Psoriasis/NEER_PRS_df.txt",stringsAsFactors = F)
Psoriasis_gene_annotation <- read.delim("../OUTPUT/LDannotation/Psoriasis_SNPs_LD_genes_anno.txt",header=T,stringsAsFactors = F)

# get variant position column
library(stringr)
Psoriasis_gene_annotation$variant_pos <- str_split(Psoriasis_gene_annotation$SNP,pattern = "_")
Psoriasis_gene_annotation$variant_pos <- lapply(Psoriasis_gene_annotation$variant_pos,function(x){paste(x[[1]],x[[2]],sep="_")})
Psoriasis_weight_matrix$variant_pos <- paste(Psoriasis_weight_matrix$CHROM,Psoriasis_weight_matrix$POS,sep="_")

# merge data
Psoriasis <- merge(Psoriasis_weight_matrix,Psoriasis_gene_annotation,by="variant_pos")
Psoriasis_variant_importance <- calc_variant_importance(Psoriasis)
Psoriasis_gene_importance <- calc_gene_importance(Psoriasis_variant_importance)


write.table(Psoriasis_gene_importance,"../OUTPUT/GSEA/Psoriasis_gene_importance.txt",sep="\t",row.names=F,col.names=T,quote=F)




### GSEA ###
### RUN GSEA ANALYSIS ###

library(fgsea)
library(dplyr)
# order by importance, generate data for GSEA
get_gsea_df <- function(importance_df){
  importance_df_new <- importance_df %>% group_by(gene_name)%>% summarise(gene_importance = sum(gene_importance))
  gsea_df <- importance_df_new[order(importance_df_new$gene_importance,decreasing=T),]
  gsea_df <- na.omit(gsea_df)
  ranks <- gsea_df$gene_importance
  names(ranks)<-gsea_df$gene_name
  return(ranks)
}


# gsea pathways loading
GSEA_HPO <- gmtPathways("../DATA/ref/c5.hpo.v7.5.symbols.gmt")
GSEA_GOBP <- gmtPathways("../DATA/ref/c5.go.bp.v7.5.symbols.gmt")
GSEA_KEGG <- gmtPathways("../DATA/ref/c2.cp.kegg.v7.5.symbols.gmt")
GSEA_ONCO <- gmtPathways("../DATA/ref/c6.all.v7.5.oncogenic_symbols.gmt")
GSEA_IMMUNE <- gmtPathways("../DATA/ref/c7.all.v7.5.immunologic_signature_symbols.gmt")


# Hypothyroidism GSEA
hypoT_importance <- read.delim("../OUTPUT/GSEA/hypoT_gene_importance.txt",header=T,stringsAsFactors = F)
hypoT_gsea_df <- get_gsea_df(hypoT_importance)

# HPO
HPOgseaRes <- fgsea(GSEA_HPO, hypoT_gsea_df, minSize=15, maxSize = 500, nperm=1000)
HPOgseaRes$leadingEdge <- as.character(HPOgseaRes$leadingEdge)
write.table(HPOgseaRes,"../OUTPUT/GSEA/GSEA_output/HypoT_HPO.txt",row.names=F,col.names=T,quote=F,sep="\t")


head(HPOgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(GSEA_HPO[["HP_ABNORMAL_CELL_MORPHOLOGY"]], hypoT_gsea_df)+
  ggtitle("GSEA of HP_ABNORMAL_CELL_MORPHOLOGY")

# GOBP
GOBPgseaRes <- fgsea(GSEA_GOBP, hypoT_gsea_df, minSize=15, maxSize = 500, nperm=1000)
GOBPgseaRes$leadingEdge <- as.character(GOBPgseaRes$leadingEdge)
write.table(GOBPgseaRes,"../OUTPUT/GSEA/GSEA_output/HypoT_GOBP.txt",row.names=F,col.names=T,quote=F,sep="\t")


head(GOBPgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(GSEA_GOBP[["GOBP_NEGATIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY"]], hypoT_gsea_df)+
  ggtitle("GOBP_NEGATIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY")

# KEGG
KEGGgseaRes <- fgsea(GSEA_KEGG, hypoT_gsea_df, minSize=15, maxSize = 500, nperm=1000)
KEGGgseaRes$leadingEdge <- as.character(KEGGgseaRes$leadingEdge)
write.table(KEGGgseaRes,"../OUTPUT/GSEA/GSEA_output/HypoT_KEGG.txt",row.names=F,col.names=T,quote=F,sep="\t")


head(KEGGgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(GSEA_KEGG[["KEGG_AUTOIMMUNE_THYROID_DISEASE"]], hypoT_gsea_df)+
  ggtitle("KEGG_AUTOIMMUNE_THYROID_DISEASE")

# ONCO
ONCOgseaRes <- fgsea(GSEA_ONCO, hypoT_gsea_df, minSize=15, maxSize = 500, nperm=1000)
ONCOgseaRes$leadingEdge <- as.character(ONCOgseaRes$leadingEdge)
write.table(ONCOgseaRes,"../OUTPUT/GSEA/GSEA_output/HypoT_ONCO.txt",row.names=F,col.names=T,quote=F,sep="\t")



head(ONCOgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(GSEA_ONCO[["PDGF_ERK_DN.V1_UP"]], hypoT_gsea_df)+
  ggtitle("PDGF_ERK_DN.V1_UP")

# IMMUNE
IMMUNEgseaRes <- fgsea(GSEA_IMMUNE, hypoT_gsea_df, minSize=15, maxSize = 500, nperm=1000)
IMMUNEgseaRes$leadingEdge <- as.character(IMMUNEgseaRes$leadingEdge)
write.table(IMMUNEgseaRes,"../OUTPUT/GSEA/GSEA_output/HypoT_IMMUNE.txt",row.names=F,col.names=T,quote=F,sep="\t")


head(IMMUNEgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(GSEA_IMMUNE[["SOBOLEV_T_CELL_PANDEMRIX_AGE_18_64YO_7DY_DN"]], hypoT_gsea_df)+
  ggtitle("SOBOLEV_T_CELL_PANDEMRIX_AGE_18_64YO_7DY_DN")


# IBD GSEA
IBD_importance <- read.delim("../OUTPUT/GSEA/IBD_gene_importance.txt",header=T,stringsAsFactors = F)
IBD_gsea_df <- get_gsea_df(IBD_importance)

# HPO
HPOgseaRes <- fgsea(GSEA_HPO, IBD_gsea_df, minSize=15, maxSize = 500, nperm=1000)
HPOgseaRes$leadingEdge <- as.character(HPOgseaRes$leadingEdge)
write.table(HPOgseaRes,"../OUTPUT/GSEA/GSEA_output/IBD_HPO.txt",row.names=F,col.names=T,quote=F,sep="\t")


head(HPOgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_HPO[["HP_ABNORMAL_CELL_MORPHOLOGY"]], IBD_gsea_df)+
#  ggtitle("GSEA of HP_ABNORMAL_CELL_MORPHOLOGY")

# GOBP
GOBPgseaRes <- fgsea(GSEA_GOBP, IBD_gsea_df, minSize=15, maxSize = 500, nperm=1000)
GOBPgseaRes$leadingEdge <- as.character(GOBPgseaRes$leadingEdge)
write.table(GOBPgseaRes,"../OUTPUT/GSEA/GSEA_output/IBD_GOBP.txt",row.names=F,col.names=T,quote=F,sep="\t")



head(GOBPgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_GOBP[["GOBP_NEGATIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY"]], IBD_gsea_df)+
#  ggtitle("GOBP_NEGATIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY")

# KEGG
KEGGgseaRes <- fgsea(GSEA_KEGG, IBD_gsea_df, minSize=15, maxSize = 500, nperm=1000)
KEGGgseaRes$leadingEdge <- as.character(KEGGgseaRes$leadingEdge)
write.table(KEGGgseaRes,"../OUTPUT/GSEA/GSEA_output/IBD_KEGG.txt",row.names=F,col.names=T,quote=F,sep="\t")



head(KEGGgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_KEGG[["KEGG_AUTOIMMUNE_THYROID_DISEASE"]], IBD_gsea_df)+
#  ggtitle("KEGG_AUTOIMMUNE_THYROID_DISEASE")

# ONCO
ONCOgseaRes <- fgsea(GSEA_ONCO, IBD_gsea_df, minSize=15, maxSize = 500, nperm=1000)
ONCOgseaRes$leadingEdge <- as.character(ONCOgseaRes$leadingEdge)
write.table(ONCOgseaRes,"../OUTPUT/GSEA/GSEA_output/IBD_ONCO.txt",row.names=F,col.names=T,quote=F,sep="\t")

head(ONCOgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_ONCO[["PDGF_ERK_DN.V1_UP"]], IBD_gsea_df)+
#  ggtitle("PDGF_ERK_DN.V1_UP")

# IMMUNE
IMMUNEgseaRes <- fgsea(GSEA_IMMUNE, IBD_gsea_df, minSize=15, maxSize = 500, nperm=1000)
IMMUNEgseaRes$leadingEdge <- as.character(IMMUNEgseaRes$leadingEdge)
write.table(IMMUNEgseaRes,"../OUTPUT/GSEA/GSEA_output/IBD_IMMUNE.txt",row.names=F,col.names=T,quote=F,sep="\t")

head(IMMUNEgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_IMMUNE[["SOBOLEV_T_CELL_PANDEMRIX_AGE_18_64YO_7DY_DN"]], IBD_gsea_df)+
#  ggtitle("SOBOLEV_T_CELL_PANDEMRIX_AGE_18_64YO_7DY_DN")


# Psoriasis GSEA
Psoriasis_importance <- read.delim("../OUTPUT/GSEA/Psoriasis_gene_importance.txt",header=T,stringsAsFactors = F)
Psoriasis_gsea_df <- get_gsea_df(Psoriasis_importance)

# HPO
HPOgseaRes <- fgsea(GSEA_HPO, Psoriasis_gsea_df, minSize=15, maxSize = 500, nperm=1000)
HPOgseaRes$leadingEdge <- as.character(HPOgseaRes$leadingEdge)
write.table(HPOgseaRes,"../OUTPUT/GSEA/GSEA_output/Psoriasis_HPO.txt",row.names=F,col.names=T,quote=F,sep="\t")

head(HPOgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_HPO[["HP_ABNORMAL_CELL_MORPHOLOGY"]], Psoriasis_gsea_df)+
#  ggtitle("GSEA of HP_ABNORMAL_CELL_MORPHOLOGY")

# GOBP
GOBPgseaRes <- fgsea(GSEA_GOBP, Psoriasis_gsea_df, minSize=15, maxSize = 500, nperm=1000)
GOBPgseaRes$leadingEdge <- as.character(GOBPgseaRes$leadingEdge)
write.table(GOBPgseaRes,"../OUTPUT/GSEA/GSEA_output/Psoriasis_GOBP.txt",row.names=F,col.names=T,quote=F,sep="\t")


head(GOBPgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_GOBP[["GOBP_NEGATIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY"]], Psoriasis_gsea_df)+
#  ggtitle("GOBP_NEGATIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY")

# KEGG
KEGGgseaRes <- fgsea(GSEA_KEGG, Psoriasis_gsea_df, minSize=15, maxSize = 500, nperm=1000)
KEGGgseaRes$leadingEdge <- as.character(KEGGgseaRes$leadingEdge)
write.table(KEGGgseaRes,"../OUTPUT/GSEA/GSEA_output/Psoriasis_KEGG.txt",row.names=F,col.names=T,quote=F,sep="\t")

head(KEGGgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_KEGG[["KEGG_AUTOIMMUNE_THYROID_DISEASE"]], Psoriasis_gsea_df)+
#  ggtitle("KEGG_AUTOIMMUNE_THYROID_DISEASE")

# ONCO
ONCOgseaRes <- fgsea(GSEA_ONCO, Psoriasis_gsea_df, minSize=15, maxSize = 500, nperm=1000)
ONCOgseaRes$leadingEdge <- as.character(ONCOgseaRes$leadingEdge)
write.table(ONCOgseaRes,"../OUTPUT/GSEA/GSEA_output/Psoriasis_ONCO.txt",row.names=F,col.names=T,quote=F,sep="\t")


head(ONCOgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_ONCO[["PDGF_ERK_DN.V1_UP"]], Psoriasis_gsea_df)+
#  ggtitle("PDGF_ERK_DN.V1_UP")

# IMMUNE
IMMUNEgseaRes <- fgsea(GSEA_IMMUNE, Psoriasis_gsea_df, minSize=15, maxSize = 500, nperm=1000)
IMMUNEgseaRes$leadingEdge <- as.character(IMMUNEgseaRes$leadingEdge)
write.table(IMMUNEgseaRes,"../OUTPUT/GSEA/GSEA_output/Psoriasis_IMMUNE.txt",row.names=F,col.names=T,quote=F,sep="\t")




head(IMMUNEgseaRes[order(padj, -abs(NES)), ], n=10)
#plotEnrichment(GSEA_IMMUNE[["SOBOLEV_T_CELL_PANDEMRIX_AGE_18_64YO_7DY_DN"]], Psoriasis_gsea_df)+
#  ggtitle("SOBOLEV_T_CELL_PANDEMRIX_AGE_18_64YO_7DY_DN")