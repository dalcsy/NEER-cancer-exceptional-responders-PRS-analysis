# visualization of PRS analysis
# set work direction
setwd("E:\\harvard\\Lab\\PRS_visualization")


neer_results <- read.delim("./results/NEER_all_results.txt",header=T)
pcawg_results <- read.delim("./results/PCAWG_all_results.txt",header=T)

# Only EUR patients
library(stringr)
#pcawg_supp <- read.csv("E:\\harvard\\Lab\\paper\\PCAWG\\Supplementary Table 1.csv",header=T)
#eur_pcawg_id <- pcawg_supp[pcawg_supp$ancestry_primary=="EUR",]$normal_specimen_aliquot_id
pcawg_results$ID <- str_replace(pcawg_results$ID,pattern = "X","")
pcawg_results$ID <- str_replace_all(pcawg_results$ID,pattern = "[.]","-")
#pcawg_prs_eur <- pcawg_results[pcawg_results$ID%in%eur_pcawg_id,]

# remove complete remission ....
#EUR_samples <- read.delim("EUR_lowsurvival_withNA_pcawg.txt",stringsAsFactors = F,header=F)
#pcawg_prs_eur <- pcawg_prs_eur[pcawg_prs_eur$ID%in%EUR_samples$V1,]



# immune checkpoint inhibitors "er_1661_sample","er_1671_sample","er_1957_sample","er_2311_sample"
#neer_results <- neer_results[neer_results$ID %in%c("er_2605_sample","er_3238_sample")==F,]

# matched cancer types
#EUR_samples <- read.delim("E:\\harvard\\Lab\\PCAWG\\comparison_separate_results\\binomial_test\\EUR_lowsurvival_withNA_matched_cancertype_pcawg.txt",stringsAsFactors = F,header=F)
#pcawg_prs_eur <- pcawg_prs_eur[pcawg_prs_eur$ID%in%EUR_samples$V1,]




### Use harmonized NEER and PCAWG dataset from Maria!!!
neer_clinical <- read.csv("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Maria_harmonized\\neer.csv",header=T)
pcawg_clinical <- read.csv("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Maria_harmonized\\pcawg_matched.csv",header=T)

neer_results <- neer_results[neer_results$ID%in%neer_clinical$ID_match_gen,]

library(stringr)
pcawg_results$ID <- str_replace(pcawg_results$ID,pattern = "X","")
pcawg_results$ID <- str_replace_all(pcawg_results$ID,pattern = "[.]","-")
eur_pcawg_id <- pcawg_clinical$normal_specimen_aliquot_id
pcawg_prs_eur <- pcawg_results[pcawg_results$ID%in%eur_pcawg_id,]


# get normalizaed data based on PCAWG distribution
for(phen in c(colnames(neer_results)[2:8])){
  mean_pcawg <- mean(pcawg_prs_eur[,phen])
  std_pcawg <- sd(pcawg_prs_eur[,phen])
  pcawg_prs_eur[,paste(phen,"norm",sep="_")] <- (pcawg_prs_eur[,phen]-mean_pcawg)/std_pcawg
  neer_results[,paste(phen,"norm",sep="_")] <- (neer_results[,phen]-mean_pcawg)/std_pcawg
}

# draw a joy plot
library(ggridges)

df_joy <- data.frame(Z_score = c(neer_results$Hypothyroidism_norm,neer_results$T1D_norm,
                                 neer_results$Multiple_sclerosis_norm,neer_results$Psoriasis_norm,neer_results$IBD_norm,
                                 neer_results$Rheumatoid_arthritis_norm,neer_results$Celiac_norm,pcawg_prs_eur$Hypothyroidism_norm,
                                 pcawg_prs_eur$T1D_norm,pcawg_prs_eur$Multiple_sclerosis_norm,pcawg_prs_eur$Psoriasis_norm,
                                 pcawg_prs_eur$IBD_norm,pcawg_prs_eur$Rheumatoid_arthritis_norm,pcawg_prs_eur$Celiac_norm
), 

Variables = c(rep("Hypothyroidism",nrow(neer_results))
              ,rep("T1D",nrow(neer_results)),rep("Multiple_sclerosis",nrow(neer_results))
              ,rep("Psoriasis",nrow(neer_results)),rep("IBD",nrow(neer_results)),rep("Rheumatoid_arthritis",nrow(neer_results))
              ,rep("Celiac",nrow(neer_results)),
              rep("Hypothyroidism",nrow(pcawg_prs_eur))
              ,rep("T1D",nrow(pcawg_prs_eur)),rep("Multiple_sclerosis",nrow(pcawg_prs_eur))
              ,rep("Psoriasis",nrow(pcawg_prs_eur)),rep("IBD",nrow(pcawg_prs_eur)),rep("Rheumatoid_arthritis",nrow(pcawg_prs_eur)),rep("Celiac",nrow(pcawg_prs_eur))
),Cohorts=c(rep("NEER",nrow(neer_results)*7),rep("PCAWG",nrow(pcawg_prs_eur)*7)))

library(ggplot2)
PRSjoyplot <- ggplot(df_joy, aes(x=Z_score, y=Variables,height=..density..,fill=Cohorts,col=Cohorts)) +
  geom_density_ridges(alpha=0.2,scale=1.3)+
  ylab('')+
  theme(text = element_text(size=20))+
  theme_classic()+theme(text = element_text(size = 9),
                        axis.text = element_text(size = 9),
                        axis.title = element_text(size = 9))


PRSjoyplot
ggsave("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Fig1_PRS_joyplot.pdf",PRSjoyplot,dpi=600,width = 13,height = 15,units = "cm",device = "pdf")


# Breast cancer
neer_results_breast <- neer_results[neer_results$ID%in%neer_clinical[neer_clinical$histology_tier2=="Breast",]$ID_match_gen,]
pcawg_prs_eur_breast <- pcawg_prs_eur[pcawg_prs_eur$ID%in%pcawg_clinical[pcawg_clinical$histology_tier2=="Breast",]$normal_specimen_aliquot_id,]


# draw a joy plot
library(ggridges)

df_joy <- data.frame(Z_score = c(neer_results_breast$Hypothyroidism_norm,neer_results_breast$T1D_norm,
                                 neer_results_breast$Multiple_sclerosis_norm,neer_results_breast$Psoriasis_norm,neer_results_breast$IBD_norm,
                                 neer_results_breast$Rheumatoid_arthritis_norm,neer_results_breast$Celiac_norm,pcawg_prs_eur_breast$Hypothyroidism_norm,
                                 pcawg_prs_eur_breast$T1D_norm,pcawg_prs_eur_breast$Multiple_sclerosis_norm,pcawg_prs_eur_breast$Psoriasis_norm,
                                 pcawg_prs_eur_breast$IBD_norm,pcawg_prs_eur_breast$Rheumatoid_arthritis_norm,pcawg_prs_eur_breast$Celiac_norm
), 

Variables = c(rep("Hypothyroidism",nrow(neer_results_breast))
              ,rep("T1D",nrow(neer_results_breast)),rep("Multiple_sclerosis",nrow(neer_results_breast))
              ,rep("Psoriasis",nrow(neer_results_breast)),rep("IBD",nrow(neer_results_breast)),rep("Rheumatoid_arthritis",nrow(neer_results_breast))
              ,rep("Celiac",nrow(neer_results_breast)),
              rep("Hypothyroidism",nrow(pcawg_prs_eur_breast))
              ,rep("T1D",nrow(pcawg_prs_eur_breast)),rep("Multiple_sclerosis",nrow(pcawg_prs_eur_breast))
              ,rep("Psoriasis",nrow(pcawg_prs_eur_breast)),rep("IBD",nrow(pcawg_prs_eur_breast)),rep("Rheumatoid_arthritis",nrow(pcawg_prs_eur_breast)),rep("Celiac",nrow(pcawg_prs_eur_breast))
),Cohorts=c(rep("NEER",nrow(neer_results_breast)*7),rep("PCAWG",nrow(pcawg_prs_eur_breast)*7)))

library(ggplot2)
PRSjoyplot_breast <- ggplot(df_joy, aes(x=Z_score, y=Variables,height=..density..,fill=Cohorts,col=Cohorts)) +
  geom_density_ridges(alpha=0.2,scale=1.3)+
  ylab('')+
  theme(text = element_text(size=20))+
  theme_classic()+theme(text = element_text(size = 9),
                        axis.text = element_text(size = 9),
                        axis.title = element_text(size = 9))+ggtitle("Breast cancer")


PRSjoyplot_breast




# pancreas cancer
neer_results_pancreas <- neer_results[neer_results$ID%in%neer_clinical[neer_clinical$histology_tier2=="Pancreas",]$ID_match_gen,]
pcawg_prs_eur_pancreas <- pcawg_prs_eur[pcawg_prs_eur$ID%in%pcawg_clinical[pcawg_clinical$histology_tier2=="Pancreas",]$normal_specimen_aliquot_id,]


# draw a joy plot
library(ggridges)

df_joy <- data.frame(Z_score = c(neer_results_pancreas$Hypothyroidism_norm,neer_results_pancreas$T1D_norm,
                                 neer_results_pancreas$Multiple_sclerosis_norm,neer_results_pancreas$Psoriasis_norm,neer_results_pancreas$IBD_norm,
                                 neer_results_pancreas$Rheumatoid_arthritis_norm,neer_results_pancreas$Celiac_norm,pcawg_prs_eur_pancreas$Hypothyroidism_norm,
                                 pcawg_prs_eur_pancreas$T1D_norm,pcawg_prs_eur_pancreas$Multiple_sclerosis_norm,pcawg_prs_eur_pancreas$Psoriasis_norm,
                                 pcawg_prs_eur_pancreas$IBD_norm,pcawg_prs_eur_pancreas$Rheumatoid_arthritis_norm,pcawg_prs_eur_pancreas$Celiac_norm
), 

Variables = c(rep("Hypothyroidism",nrow(neer_results_pancreas))
              ,rep("T1D",nrow(neer_results_pancreas)),rep("Multiple_sclerosis",nrow(neer_results_pancreas))
              ,rep("Psoriasis",nrow(neer_results_pancreas)),rep("IBD",nrow(neer_results_pancreas)),rep("Rheumatoid_arthritis",nrow(neer_results_pancreas))
              ,rep("Celiac",nrow(neer_results_pancreas)),
              rep("Hypothyroidism",nrow(pcawg_prs_eur_pancreas))
              ,rep("T1D",nrow(pcawg_prs_eur_pancreas)),rep("Multiple_sclerosis",nrow(pcawg_prs_eur_pancreas))
              ,rep("Psoriasis",nrow(pcawg_prs_eur_pancreas)),rep("IBD",nrow(pcawg_prs_eur_pancreas)),rep("Rheumatoid_arthritis",nrow(pcawg_prs_eur_pancreas)),rep("Celiac",nrow(pcawg_prs_eur_pancreas))
),Cohorts=c(rep("NEER",nrow(neer_results_pancreas)*7),rep("PCAWG",nrow(pcawg_prs_eur_pancreas)*7)))

library(ggplot2)
PRSjoyplot_pancreas <- ggplot(df_joy, aes(x=Z_score, y=Variables,height=..density..,fill=Cohorts,col=Cohorts)) +
  geom_density_ridges(alpha=0.2,scale=1.3)+
  ylab('')+
  theme(text = element_text(size=20))+
  theme_classic()+theme(text = element_text(size = 9),
                        axis.text = element_text(size = 9),
                        axis.title = element_text(size = 9))+ggtitle("Pancreatic cancer")


PRSjoyplot_pancreas


# lung cancer
neer_results_lung <- neer_results[neer_results$ID%in%neer_clinical[neer_clinical$histology_tier2=="Lung",]$ID_match_gen,]
pcawg_prs_eur_lung <- pcawg_prs_eur[pcawg_prs_eur$ID%in%pcawg_clinical[pcawg_clinical$histology_tier2=="Lung",]$normal_specimen_aliquot_id,]


# draw a joy plot
library(ggridges)

df_joy <- data.frame(Z_score = c(neer_results_lung$Hypothyroidism_norm,neer_results_lung$T1D_norm,
                                 neer_results_lung$Multiple_sclerosis_norm,neer_results_lung$Psoriasis_norm,neer_results_lung$IBD_norm,
                                 neer_results_lung$Rheumatoid_arthritis_norm,neer_results_lung$Celiac_norm,pcawg_prs_eur_lung$Hypothyroidism_norm,
                                 pcawg_prs_eur_lung$T1D_norm,pcawg_prs_eur_lung$Multiple_sclerosis_norm,pcawg_prs_eur_lung$Psoriasis_norm,
                                 pcawg_prs_eur_lung$IBD_norm,pcawg_prs_eur_lung$Rheumatoid_arthritis_norm,pcawg_prs_eur_lung$Celiac_norm
), 

Variables = c(rep("Hypothyroidism",nrow(neer_results_lung))
              ,rep("T1D",nrow(neer_results_lung)),rep("Multiple_sclerosis",nrow(neer_results_lung))
              ,rep("Psoriasis",nrow(neer_results_lung)),rep("IBD",nrow(neer_results_lung)),rep("Rheumatoid_arthritis",nrow(neer_results_lung))
              ,rep("Celiac",nrow(neer_results_lung)),
              rep("Hypothyroidism",nrow(pcawg_prs_eur_lung))
              ,rep("T1D",nrow(pcawg_prs_eur_lung)),rep("Multiple_sclerosis",nrow(pcawg_prs_eur_lung))
              ,rep("Psoriasis",nrow(pcawg_prs_eur_lung)),rep("IBD",nrow(pcawg_prs_eur_lung)),rep("Rheumatoid_arthritis",nrow(pcawg_prs_eur_lung)),rep("Celiac",nrow(pcawg_prs_eur_lung))
),Cohorts=c(rep("NEER",nrow(neer_results_lung)*7),rep("PCAWG",nrow(pcawg_prs_eur_lung)*7)))

library(ggplot2)
PRSjoyplot_lung <- ggplot(df_joy, aes(x=Z_score, y=Variables,height=..density..,fill=Cohorts,col=Cohorts)) +
  geom_density_ridges(alpha=0.2,scale=1.3)+
  ylab('')+
  theme(text = element_text(size=20))+
  theme_classic()+theme(text = element_text(size = 9),
                        axis.text = element_text(size = 9),
                        axis.title = element_text(size = 9))+ggtitle("Lung cancer")


PRSjoyplot_lung

library(ggpubr)
joyplot_B_P_L <- ggarrange(PRSjoyplot_breast+theme(legend.position = 'none'),PRSjoyplot_pancreas+theme(legend.position = 'none'),PRSjoyplot_lung+theme(legend.position = 'none'),ncol=3,nrow=1)
joyplot_B_P_L

ggsave("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Fig2_PRS_joyplot_cancertypes.pdf",joyplot_B_P_L,dpi=600,width = 20,height = 10,units = "cm",device = "pdf")




# visualization in separate plots
chr <- c("A","B","C","D","E","F","G")
draw_dist_plot <- function(prs_scores,xlab){
  neer <- prs_scores[prs_scores$Cohorts=="NEER",]
  
  prs_scores <- rbind(prs_scores,neer)
  
  
  
  
  ggplot(prs_scores,aes(x=Z_score,fill=Cohorts))+
    geom_histogram(bins = 30,alpha = 0.5, position = "identity")+
    theme_classic()+
    xlab(xlab)+
    xlim(c(-5,5))+
    scale_y_continuous("PCAWG counts",
                       sec.axis = sec_axis(~ . /2, name = "NEER counts",)) +theme(text = element_text(size = 9),
                                                                                  axis.text = element_text(size = 9),
                                                                                  axis.title = element_text(size = 9))
}
phenos <- unique(df_joy$Variables)
for (i in 1:7){
  assign(paste("fig",chr[i],sep="_"),draw_dist_plot(df_joy[df_joy$Variables==phenos[i],],xlab=paste(phenos[i],"PRS")))
}
library(ggpubr)
sup_fig_hist <- ggarrange(fig_A,fig_B,fig_C,fig_D,fig_E,fig_F,fig_G,ncol=2,nrow=4,labels=chr)
sup_fig_hist
ggsave("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Sup_fig_hist.pdf",sup_fig_hist,dpi=600,width = 18,height = 26,units = "cm",device = "pdf")


# logistic regression to get P-value
hypoT_df <- df_joy[df_joy$Variables=="Hypothyroidism",]
hypoT_df$Y <- ifelse(hypoT_df$Cohorts=="NEER",1,0)
hypoT_logit <- glm(Y~Z_score,family = "binomial",data = hypoT_df)
summary(hypoT_logit)

T1D_df <- df_joy[df_joy$Variables=="T1D",]
T1D_df$Y <- ifelse(T1D_df$Cohorts=="NEER",1,0)
T1D_logit <- glm(Y~Z_score,family = "binomial",data = T1D_df)
summary(T1D_logit)

Multiple_sclerosis_df <- df_joy[df_joy$Variables=="Multiple_sclerosis",]
Multiple_sclerosis_df$Y <- ifelse(Multiple_sclerosis_df$Cohorts=="NEER",1,0)
Multiple_sclerosis_logit <- glm(Y~Z_score,family = "binomial",data = Multiple_sclerosis_df)
summary(Multiple_sclerosis_logit)

Psoriasis_df <- df_joy[df_joy$Variables=="Psoriasis",]
Psoriasis_df$Y <- ifelse(Psoriasis_df$Cohorts=="NEER",1,0)
Psoriasis_logit <- glm(Y~Z_score,family = "binomial",data = Psoriasis_df)
summary(Psoriasis_logit)

Rheumatoid_arthritis_df <- df_joy[df_joy$Variables=="Rheumatoid_arthritis",]
Rheumatoid_arthritis_df$Y <- ifelse(Rheumatoid_arthritis_df$Cohorts=="NEER",1,0)
Rheumatoid_arthritis_logit <- glm(Y~Z_score,family = "binomial",data = Rheumatoid_arthritis_df)
summary(Rheumatoid_arthritis_logit)

IBD_df <- df_joy[df_joy$Variables=="IBD",]
IBD_df$Y <- ifelse(IBD_df$Cohorts=="NEER",1,0)
IBD_logit <- glm(Y~Z_score,family = "binomial",data = IBD_df)
summary(IBD_logit)

Celiac_df <- df_joy[df_joy$Variables=="Celiac",]
Celiac_df$Y <- ifelse(Celiac_df$Cohorts=="NEER",1,0)
Celiac_logit <- glm(Y~Z_score,family = "binomial",data = Celiac_df)
summary(Celiac_logit)




# HOW TO UNDERSTAND THE PART OF T1D PATIENTS WITH HIGH PRS RISK?
row.names(neer_results) <- 1:51
T1D_high <- row.names(neer_results[neer_results$T1D_norm>0,])
T1D_low <- row.names(neer_results[neer_results$T1D_norm<0,])

neer_clinical_H <- neer_clinical[neer_clinical$ID_match_gen%in%T1D_high,]
neer_clinical_L <- neer_clinical[neer_clinical$ID_match_gen%in%T1D_low,]

# Nothing special in histology 
drug_neer <- read.table("clipboard",sep="\t",header=T)
drug_neer_H <- drug_neer[T1D_high,]
drug_neer_low <- drug_neer[T1D_low,]

write.table(drug_neer_H,"T1D_HIGH_NEER_clinical.txt",sep="\t",quote=F)
write.table(drug_neer_low,"T1D_LOW_NEER_clinical.txt",sep="\t",quote=F)

table(drug_neer_H$immunotherapy)
table(drug_neer_low$immunotherapy)

fisher.test(matrix(c(33,6,6,1),nrow=2))







# visualize by every mutation in T1D distribution.
T1D_snps <- read.csv("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Maria_harmonized\\snps_t1d_genes.csv")
head(T1D_snps)
T1D_snps[T1D_snps=="0|0"]<- as.numeric(0)
T1D_snps[T1D_snps=="0|1"]<- as.numeric(1)
T1D_snps[T1D_snps=="1|1"]<- as.numeric(2)
T1D_snps[T1D_snps=="1|0"]<- as.numeric(1)




head(T1D_snps)
T1D_high <- neer_results[neer_results$T1D_norm>0,]$ID
T1D_low <- neer_results[neer_results$T1D_norm<0,]$ID
IBD_low <- neer_results[neer_results$IBD_norm<0,]$ID
length(intersect(T1D_high,IBD_low))

genes <- unique(T1D_snps$Gene)
for(i in 16:71){
  T1D_snps[,i] <- as.numeric(T1D_snps[,i])
}


setwd("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Gene_level_labeled_T1D_PRS")
for(gene_i in genes){
  gene_subset <- as.data.frame(colSums(T1D_snps[T1D_snps$Gene==gene_i,neer_results$ID]))
  with_snp_ER <- row.names(gene_subset)[gene_subset[,1]>0]
  fig <- ggplot(T1D_df,aes(x=Z_score,group=Cohorts,color=Cohorts))+
    geom_density()+
    theme_classic()+
    geom_vline(xintercept = neer_results[neer_results$ID%in%with_snp_ER,]$T1D_norm,color="black",linetype="dashed")+
    xlab("T1D PRS")+
    ggtitle(paste(gene_i,", N with snp =",length(with_snp_ER)))
  png(paste(gene_i,"T1D.png",sep="_"),res=300,units = "cm",width = 18,height=12)
  print(fig)
  dev.off()
  
}


setwd("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Variant_level_labeled_T1D_PRS")

variants <- T1D_snps$snp_id
for(v_i in variants){
  gene_i <- T1D_snps[T1D_snps$snp_id==v_i,]$Gene
  if(T1D_snps[T1D_snps$snp_id==v_i,]$effect_allele == T1D_snps[T1D_snps$snp_id==v_i,]$ALT & T1D_snps[T1D_snps$snp_id==v_i,]$effect_weight>0){
    weight="POS"
  }else if((T1D_snps[T1D_snps$snp_id==v_i,]$effect_allele == T1D_snps[T1D_snps$snp_id==v_i,]$ALT & T1D_snps[T1D_snps$snp_id==v_i,]$effect_weight<0)){
    weight="NEG"
  }else if(T1D_snps[T1D_snps$snp_id==v_i,]$effect_allele == T1D_snps[T1D_snps$snp_id==v_i,]$REF & T1D_snps[T1D_snps$snp_id==v_i,]$effect_weight>0){
    weight="NEG"
  }else if(T1D_snps[T1D_snps$snp_id==v_i,]$effect_allele == T1D_snps[T1D_snps$snp_id==v_i,]$REF & T1D_snps[T1D_snps$snp_id==v_i,]$effect_weight<0){
    weight="POS"
  }
  
  
  variant_subset <- as.data.frame(colSums(T1D_snps[T1D_snps$snp_id==v_i,neer_results$ID]))
  with_snp_ER <- row.names(variant_subset)[variant_subset[,1]>0]
  fig <- ggplot(T1D_df,aes(x=Z_score,group=Cohorts,color=Cohorts))+
    geom_density()+
    theme_classic()+
    geom_vline(xintercept = neer_results[neer_results$ID%in%with_snp_ER,]$T1D_norm,color="black",linetype="dashed")+
    xlab("T1D PRS")+
    ggtitle(paste(v_i,gene_i,weight,", N with snp =",length(with_snp_ER)))
  png(paste(v_i,"T1D.png",sep="_"),res=300,units = "cm",width = 18,height=12)
  print(fig)
  dev.off()
  
}

#T1D_LOW = 9
#TID_HIGH = 42

#run test for each variant

variants <- T1D_snps$snp_id
res <- data.frame(variant_id =variants,weight=NA,OR=NA,p_value=NA,gene=NA)
res$functional_consequence <- T1D_snps$Functional_Consequence
for(v_i in variants){
  gene_i <- T1D_snps[T1D_snps$snp_id==v_i,]$Gene
  if(T1D_snps[T1D_snps$snp_id==v_i,]$effect_allele == T1D_snps[T1D_snps$snp_id==v_i,]$ALT & T1D_snps[T1D_snps$snp_id==v_i,]$effect_weight>0){
    weight="Susceptible"
  }else if((T1D_snps[T1D_snps$snp_id==v_i,]$effect_allele == T1D_snps[T1D_snps$snp_id==v_i,]$ALT & T1D_snps[T1D_snps$snp_id==v_i,]$effect_weight<0)){
    weight="Protective"
  }else if(T1D_snps[T1D_snps$snp_id==v_i,]$effect_allele == T1D_snps[T1D_snps$snp_id==v_i,]$REF & T1D_snps[T1D_snps$snp_id==v_i,]$effect_weight>0){
    weight="Protective"
  }else if(T1D_snps[T1D_snps$snp_id==v_i,]$effect_allele == T1D_snps[T1D_snps$snp_id==v_i,]$REF & T1D_snps[T1D_snps$snp_id==v_i,]$effect_weight<0){
    weight="Susceptible"
  }
  
  
  variant_subset <- as.data.frame(colSums(T1D_snps[T1D_snps$snp_id==v_i,neer_results$ID]))
  T1D_H_alleles <- sum(variant_subset[T1D_high,])
  T1D_L_alleles <- sum(variant_subset[T1D_low,])
  fisher_res <- fisher.test(matrix(c(T1D_H_alleles,length(T1D_high)*2-T1D_H_alleles,T1D_L_alleles,length(T1D_low)*2-T1D_L_alleles),nrow = 2))
  res[res$variant_id==v_i,]$weight <- weight
  res[res$variant_id==v_i,]$OR <- fisher_res$estimate
  res[res$variant_id==v_i,]$p_value<-fisher_res$p.value
  res[res$variant_id==v_i,]$gene <- gene_i

}
res

res$p.adjust <- p.adjust(res$p_value,method = "BH")
res
res[res$p.adjust<0.1,]



###
#skewness test
library(normtest)
set.seed(3398)
skewness.norm.test(neer_results$Hypothyroidism)
set.seed(3398)
skewness.norm.test(neer_results$T1D)
set.seed(3398)
skewness.norm.test(neer_results$Multiple_sclerosis)
set.seed(3398)
skewness.norm.test(neer_results$Psoriasis)
set.seed(3398)
skewness.norm.test(neer_results$IBD)
set.seed(3398)
skewness.norm.test(neer_results$Rheumatoid_arthritis)
set.seed(3398)
skewness.norm.test(neer_results$Celiac)
set.seed(3398)
skewness.norm.test(pcawg_prs_eur$Hypothyroidism)
set.seed(3398)
skewness.norm.test(pcawg_prs_eur$T1D)
set.seed(3398)
skewness.norm.test(pcawg_prs_eur$Multiple_sclerosis)
set.seed(3398)
skewness.norm.test(pcawg_prs_eur$Psoriasis)
set.seed(3398)
skewness.norm.test(pcawg_prs_eur$IBD)
set.seed(3398)
skewness.norm.test(pcawg_prs_eur$Rheumatoid_arthritis)
set.seed(3398)
skewness.norm.test(pcawg_prs_eur$Celiac)











# sup fig
# survival in PCAWG with top 10% and bottom 10% PRS score
sup_Table <- read.csv("E:\\harvard\\Lab\\paper\\PCAWG\\Supplementary Table 1.csv",header=T,stringsAsFactors = F)
clinical <- read.delim("E:\\harvard\\Lab\\PCAWG\\donor.tsv",header=T)
clinical <- merge(clinical,sup_Table,by="icgc_donor_id")



# USE ALL PCAWG PATIENTS???!!!
pcawg_results <- read.delim("./results/PCAWG_all_results.txt",header=T)
library(stringr)

pcawg_results$ID <- str_replace(pcawg_results$ID,pattern = "X","")
pcawg_results$ID <- str_replace_all(pcawg_results$ID,pattern = "[.]","-")
pcawg_prs_eur <- pcawg_results
for(phen in c(colnames(neer_results)[2:8])){
  mean_pcawg <- mean(pcawg_prs_eur[,phen])
  std_pcawg <- sd(pcawg_prs_eur[,phen])
  pcawg_prs_eur[,paste(phen,"norm",sep="_")] <- (pcawg_prs_eur[,phen]-mean_pcawg)/std_pcawg
  neer_results[,paste(phen,"norm",sep="_")] <- (neer_results[,phen]-mean_pcawg)/std_pcawg
}


pcawg_prs_eur_clinical <- merge(pcawg_prs_eur,clinical,by.x="ID",by.y="normal_specimen_aliquot_id")

# IF ONLY INCLUDE one cancer type
pcawg_prs_eur_clinical <- pcawg_prs_eur_clinical[pcawg_prs_eur_clinical$histology_tier2=="Lymphoid",]


# T1D
t1d_lower_cutoff <- quantile(pcawg_prs_eur$T1D_norm,probs=0.5)
t1d_upper_cutoff <- quantile(pcawg_prs_eur$T1D_norm,probs=0.5)
pcawg_prs_eur_clinical_t1d <- pcawg_prs_eur_clinical[pcawg_prs_eur_clinical$T1D_norm<t1d_lower_cutoff|pcawg_prs_eur_clinical$T1D_norm>t1d_upper_cutoff,]
survival_Df <- data.frame(patient_id = pcawg_prs_eur_clinical_t1d$ID,
                          PRS = pcawg_prs_eur_clinical_t1d$T1D_norm, 
                          days_to_death = pcawg_prs_eur_clinical_t1d$donor_survival_time,
                          vital_status=pcawg_prs_eur_clinical_t1d$donor_vital_status,
                          days_to_last_follo_wup = pcawg_prs_eur_clinical_t1d$donor_interval_of_last_followup)

survival_Df$time <- 0
for(i in 1:nrow(survival_Df)){
  if(survival_Df[i,]$vital_status=="deceased"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_death
  }else if(survival_Df[i,]$vital_status=="alive"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_last_follo_wup
  }
}

survival_Df$time <- as.numeric(survival_Df$time)


#install.packages('survival')
library(survival)
#install.packages('survminer')
library(survminer)

survival_Df<- survival_Df[survival_Df$vital_status %in% c("deceased","alive"),]
survival_Df$vital_status <- ifelse(survival_Df$vital_status=="alive",F,T)

survival_Df$T1D_risk <- ifelse(survival_Df$PRS>0,"High","Low")
survival_Df <- survival_Df[is.na(survival_Df$vital_status)==F,]

fit <- survfit(Surv(time,vital_status) ~ T1D_risk,  
               data = survival_Df)

t1d_survival_plot <- ggsurvplot(fit,conf.int = TRUE,pval=T,risk.table = T,risk.table.height = 0.4)+
  xlab("Time (days)")

# IBD
IBD_lower_cutoff <- quantile(pcawg_prs_eur$IBD_norm,probs=0.5)
IBD_upper_cutoff <- quantile(pcawg_prs_eur$IBD_norm,probs=0.5)
pcawg_prs_eur_clinical_IBD <- pcawg_prs_eur_clinical[pcawg_prs_eur_clinical$IBD_norm<IBD_lower_cutoff|pcawg_prs_eur_clinical$IBD_norm>IBD_upper_cutoff,]

survival_Df <- data.frame(patient_id = pcawg_prs_eur_clinical_IBD$ID,
                          PRS = pcawg_prs_eur_clinical_IBD$IBD_norm, 
                          days_to_death = pcawg_prs_eur_clinical_IBD$donor_survival_time,
                          vital_status=pcawg_prs_eur_clinical_IBD$donor_vital_status,
                          days_to_last_follo_wup = pcawg_prs_eur_clinical_IBD$donor_interval_of_last_followup)

survival_Df$time <- 0
for(i in 1:nrow(survival_Df)){
  if(survival_Df[i,]$vital_status=="deceased"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_death
  }else if(survival_Df[i,]$vital_status=="alive"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_last_follo_wup
  }
}

survival_Df$time <- as.numeric(survival_Df$time)


#install.packages('survival')
library(survival)
#install.packages('survminer')
library(survminer)

survival_Df<- survival_Df[survival_Df$vital_status %in% c("deceased","alive"),]
survival_Df$vital_status <- ifelse(survival_Df$vital_status=="alive",F,T)

survival_Df$IBD_risk <- ifelse(survival_Df$PRS>0,"High","Low")
survival_Df <- survival_Df[is.na(survival_Df$vital_status)==F,]

fit <- survfit(Surv(time,vital_status) ~ IBD_risk,  
               data = survival_Df)

ibd_survival_plot <- ggsurvplot(fit,conf.int = TRUE,pval=T,risk.table = T,risk.table.height = 0.4)+
  xlab("Time (days)")
ibd_survival_plot

# Hypothyroidism
Hypothyroidism_lower_cutoff <- quantile(pcawg_prs_eur$Hypothyroidism_norm,probs=0.5)
Hypothyroidism_upper_cutoff <- quantile(pcawg_prs_eur$Hypothyroidism_norm,probs=0.5)
pcawg_prs_eur_clinical_Hypothyroidism <- pcawg_prs_eur_clinical[pcawg_prs_eur_clinical$Hypothyroidism_norm<Hypothyroidism_lower_cutoff|pcawg_prs_eur_clinical$Hypothyroidism_norm>Hypothyroidism_upper_cutoff,]

survival_Df <- data.frame(patient_id = pcawg_prs_eur_clinical_Hypothyroidism$ID,
                          PRS = pcawg_prs_eur_clinical_Hypothyroidism$Hypothyroidism_norm, 
                          days_to_death = pcawg_prs_eur_clinical_Hypothyroidism$donor_survival_time,
                          vital_status=pcawg_prs_eur_clinical_Hypothyroidism$donor_vital_status,
                          days_to_last_follo_wup = pcawg_prs_eur_clinical_Hypothyroidism$donor_interval_of_last_followup)

survival_Df$time <- 0
for(i in 1:nrow(survival_Df)){
  if(survival_Df[i,]$vital_status=="deceased"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_death
  }else if(survival_Df[i,]$vital_status=="alive"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_last_follo_wup
  }
}

survival_Df$time <- as.numeric(survival_Df$time)


#install.packages('survival')
library(survival)
#install.packages('survminer')
library(survminer)

survival_Df<- survival_Df[survival_Df$vital_status %in% c("deceased","alive"),]
survival_Df$vital_status <- ifelse(survival_Df$vital_status=="alive",F,T)

survival_Df$Hypothyroidism_risk <- ifelse(survival_Df$PRS>0,"High","Low")
survival_Df <- survival_Df[is.na(survival_Df$vital_status)==F,]

fit <- survfit(Surv(time,vital_status) ~ Hypothyroidism_risk,  
               data = survival_Df)

hypoT_survival_plot <- ggsurvplot(fit,conf.int = TRUE,pval=T,risk.table = T,risk.table.height = 0.4)+
  xlab("Time (days)")
hypoT_survival_plot

# MS
Multiple_sclerosis_lower_cutoff <- quantile(pcawg_prs_eur$Multiple_sclerosis_norm,probs=0.5)
Multiple_sclerosis_upper_cutoff <- quantile(pcawg_prs_eur$Multiple_sclerosis_norm,probs=0.5)
pcawg_prs_eur_clinical_Multiple_sclerosis <- pcawg_prs_eur_clinical[pcawg_prs_eur_clinical$Multiple_sclerosis_norm<Multiple_sclerosis_lower_cutoff|pcawg_prs_eur_clinical$Multiple_sclerosis_norm>Multiple_sclerosis_upper_cutoff,]

survival_Df <- data.frame(patient_id = pcawg_prs_eur_clinical_Multiple_sclerosis$ID,
                          PRS = pcawg_prs_eur_clinical_Multiple_sclerosis$Multiple_sclerosis_norm, 
                          days_to_death = pcawg_prs_eur_clinical_Multiple_sclerosis$donor_survival_time,
                          vital_status=pcawg_prs_eur_clinical_Multiple_sclerosis$donor_vital_status,
                          days_to_last_follo_wup = pcawg_prs_eur_clinical_Multiple_sclerosis$donor_interval_of_last_followup)

survival_Df$time <- 0
for(i in 1:nrow(survival_Df)){
  if(survival_Df[i,]$vital_status=="deceased"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_death
  }else if(survival_Df[i,]$vital_status=="alive"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_last_follo_wup
  }
}

survival_Df$time <- as.numeric(survival_Df$time)


#install.packages('survival')
library(survival)
#install.packages('survminer')
library(survminer)

survival_Df<- survival_Df[survival_Df$vital_status %in% c("deceased","alive"),]
survival_Df$vital_status <- ifelse(survival_Df$vital_status=="alive",F,T)

survival_Df$Multiple_sclerosis_risk <- ifelse(survival_Df$PRS>0,"High","Low")
survival_Df <- survival_Df[is.na(survival_Df$vital_status)==F,]

fit <- survfit(Surv(time,vital_status) ~ Multiple_sclerosis_risk,  
               data = survival_Df)

MS_survival_plot <- ggsurvplot(fit,conf.int = TRUE,pval=T,risk.table = T,risk.table.height = 0.4)+
  xlab("Time (days)")
MS_survival_plot





# Psoriasis
Psoriasis_lower_cutoff <- quantile(pcawg_prs_eur$Psoriasis_norm,probs=0.5)
Psoriasis_upper_cutoff <- quantile(pcawg_prs_eur$Psoriasis_norm,probs=0.5)
pcawg_prs_eur_clinical_Psoriasis <- pcawg_prs_eur_clinical[pcawg_prs_eur_clinical$Psoriasis_norm<Psoriasis_lower_cutoff|pcawg_prs_eur_clinical$Psoriasis_norm>Psoriasis_upper_cutoff,]

survival_Df <- data.frame(patient_id = pcawg_prs_eur_clinical_Psoriasis$ID,
                          PRS = pcawg_prs_eur_clinical_Psoriasis$Psoriasis_norm, 
                          days_to_death = pcawg_prs_eur_clinical_Psoriasis$donor_survival_time,
                          vital_status=pcawg_prs_eur_clinical_Psoriasis$donor_vital_status,
                          days_to_last_follo_wup = pcawg_prs_eur_clinical_Psoriasis$donor_interval_of_last_followup)

survival_Df$time <- 0
for(i in 1:nrow(survival_Df)){
  if(survival_Df[i,]$vital_status=="deceased"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_death
  }else if(survival_Df[i,]$vital_status=="alive"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_last_follo_wup
  }
}

survival_Df$time <- as.numeric(survival_Df$time)


#install.packages('survival')
library(survival)
#install.packages('survminer')
library(survminer)

survival_Df<- survival_Df[survival_Df$vital_status %in% c("deceased","alive"),]
survival_Df$vital_status <- ifelse(survival_Df$vital_status=="alive",F,T)

survival_Df$Psoriasis_risk <- ifelse(survival_Df$PRS>0,"High","Low")
survival_Df <- survival_Df[is.na(survival_Df$vital_status)==F,]

fit <- survfit(Surv(time,vital_status) ~ Psoriasis_risk,  
               data = survival_Df)

Psoriasis_survival_plot <- ggsurvplot(fit,conf.int = TRUE,pval=T,risk.table = T,risk.table.height = 0.4)+
  xlab("Time (days)")
Psoriasis_survival_plot


# Rheumatoid_arthritis
Rheumatoid_arthritis_lower_cutoff <- quantile(pcawg_prs_eur$Rheumatoid_arthritis_norm,probs=0.5)
Rheumatoid_arthritis_upper_cutoff <- quantile(pcawg_prs_eur$Rheumatoid_arthritis_norm,probs=0.5)
pcawg_prs_eur_clinical_Rheumatoid_arthritis <- pcawg_prs_eur_clinical[pcawg_prs_eur_clinical$Rheumatoid_arthritis_norm<Rheumatoid_arthritis_lower_cutoff|pcawg_prs_eur_clinical$Rheumatoid_arthritis_norm>Rheumatoid_arthritis_upper_cutoff,]

survival_Df <- data.frame(patient_id = pcawg_prs_eur_clinical_Rheumatoid_arthritis$ID,
                          PRS = pcawg_prs_eur_clinical_Rheumatoid_arthritis$Rheumatoid_arthritis_norm, 
                          days_to_death = pcawg_prs_eur_clinical_Rheumatoid_arthritis$donor_survival_time,
                          vital_status=pcawg_prs_eur_clinical_Rheumatoid_arthritis$donor_vital_status,
                          days_to_last_follo_wup = pcawg_prs_eur_clinical_Rheumatoid_arthritis$donor_interval_of_last_followup)

survival_Df$time <- 0
for(i in 1:nrow(survival_Df)){
  if(survival_Df[i,]$vital_status=="deceased"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_death
  }else if(survival_Df[i,]$vital_status=="alive"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_last_follo_wup
  }
}

survival_Df$time <- as.numeric(survival_Df$time)


#install.packages('survival')
library(survival)
#install.packages('survminer')
library(survminer)

survival_Df<- survival_Df[survival_Df$vital_status %in% c("deceased","alive"),]
survival_Df$vital_status <- ifelse(survival_Df$vital_status=="alive",F,T)

survival_Df$Rheumatoid_arthritis_risk <- ifelse(survival_Df$PRS>0,"High","Low")
survival_Df <- survival_Df[is.na(survival_Df$vital_status)==F,]
survival_Df$RA_risk <- survival_Df$Rheumatoid_arthritis_risk
fit <- survfit(Surv(time,vital_status) ~ RA_risk,  
               data = survival_Df)

Rheumatoid_arthritis_survival_plot <- ggsurvplot(fit,conf.int = TRUE,pval=T,risk.table = T,risk.table.height = 0.4)+
  xlab("Time (days)")+
  ylab("")
Rheumatoid_arthritis_survival_plot

# Celiac
Celiac_lower_cutoff <- quantile(pcawg_prs_eur$Celiac_norm,probs=0.5)
Celiac_upper_cutoff <- quantile(pcawg_prs_eur$Celiac_norm,probs=0.5)
pcawg_prs_eur_clinical_Celiac <- pcawg_prs_eur_clinical[pcawg_prs_eur_clinical$Celiac_norm<Celiac_lower_cutoff|pcawg_prs_eur_clinical$Celiac_norm>Celiac_upper_cutoff,]

survival_Df <- data.frame(patient_id = pcawg_prs_eur_clinical_Celiac$ID,
                          PRS = pcawg_prs_eur_clinical_Celiac$Celiac_norm, 
                          days_to_death = pcawg_prs_eur_clinical_Celiac$donor_survival_time,
                          vital_status=pcawg_prs_eur_clinical_Celiac$donor_vital_status,
                          days_to_last_follo_wup = pcawg_prs_eur_clinical_Celiac$donor_interval_of_last_followup)

survival_Df$time <- 0
for(i in 1:nrow(survival_Df)){
  if(survival_Df[i,]$vital_status=="deceased"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_death
  }else if(survival_Df[i,]$vital_status=="alive"){
    survival_Df[i,]$time <- survival_Df[i,]$days_to_last_follo_wup
  }
}

survival_Df$time <- as.numeric(survival_Df$time)


#install.packages('survival')
library(survival)
#install.packages('survminer')
library(survminer)

survival_Df<- survival_Df[survival_Df$vital_status %in% c("deceased","alive"),]
survival_Df$vital_status <- ifelse(survival_Df$vital_status=="alive",F,T)

survival_Df$Celiac_risk <- ifelse(survival_Df$PRS>0,"High","Low")
survival_Df <- survival_Df[is.na(survival_Df$vital_status)==F,]
survival_Df$Celiac_risk <- survival_Df$Celiac_risk
fit <- survfit(Surv(time,vital_status) ~ Celiac_risk,  
               data = survival_Df)

Celiac_survival_plot <- ggsurvplot(fit,conf.int = TRUE,pval=T,risk.table = T,risk.table.height = 0.4)+
  xlab("Time (days)")+
  ylab("")
Celiac_survival_plot




setwd("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission/")

pdf("NO414pcawg_survival_curves.pdf",width = 11,height = 15)
print(ggarrange(hypoT_survival_plot$plot,
                t1d_survival_plot$plot,
                MS_survival_plot$plot,
                Psoriasis_survival_plot$plot,
                ibd_survival_plot$plot,
                Rheumatoid_arthritis_survival_plot$plot,Celiac_survival_plot$plot,ncol=2,nrow=4,labels=c("A","B","C","D","E","F","G")))
dev.off()
setwd("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission/")


for(i in c("hypoT_survival_plot","t1d_survival_plot","MS_survival_plot","Psoriasis_survival_plot","ibd_survival_plot","Rheumatoid_arthritis_survival_plot")){
  assign(i,customize_labels(
    get(i),
    font.title    = c(10),#用长度为3的向量分别指定大小、类型、颜色
    font.subtitle = c(10),
    font.caption  = c(10),
    font.x        = c(10),
    font.y        = c(10),
    font.xtickslab = c(10),
    font.ytickslab = c(10)
    
  ))
  

  png(paste(i,".png",sep=""),res = 300,width = 17,height=8,units="cm")
  print(get(i))
  dev.off()
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
multiplot(hypoT_survival_plot,
           t1d_survival_plot,
           MS_survival_plot,
           Psoriasis_survival_plot,
           ibd_survival_plot,
           Rheumatoid_arthritis_survival_plot,cols=2)





# linkage analysis
# two variants commonly occur in a given patient

T1D_NEER_df <- read.delim("E:\\harvard\\Lab\\PRS_visualization\\T1D\\NEER_PRS_df.txt",row.names = 1)
T1D_NEER_df <- T1D_NEER_df[,neer_results$ID]
T1D_NEER_Res <- as.data.frame(matrix(NA,ncol=66,nrow=66))
row.names(T1D_NEER_Res) <- row.names(T1D_NEER_df)
colnames(T1D_NEER_Res) <- row.names(T1D_NEER_df)

for(i in row.names(T1D_NEER_Res)){
  for(j in colnames(T1D_NEER_Res)){
    double_ <- t(T1D_NEER_df[i,])*t(T1D_NEER_df[j,])
    T1D_NEER_Res[i,j] <- sum(double_ != 0)
  }
}
T1D_NEER_Res <- T1D_NEER_Res/51
T1D_NEER_Res <- round(T1D_NEER_Res,2)
library(pheatmap)
pheatmap(T1D_NEER_Res,cluster_rows = T,cluster_cols = T,display_numbers = T,scale = "none")


T1D_PCAWG_df <- read.delim("E:\\harvard\\Lab\\PRS_visualization\\T1D\\PCAWG_PRS_df.txt",check.names = F)
T1D_PCAWG_df_new <- as.data.frame(matrix(NA,nrow=66,ncol=length(eur_pcawg_id)))
row.names(T1D_PCAWG_df_new) <- unique(T1D_PCAWG_df$Row.names)
colnames(T1D_PCAWG_df_new) <- eur_pcawg_id
for(i in row.names(T1D_PCAWG_df_new)){
  T1D_PCAWG_subset <- T1D_PCAWG_df[T1D_PCAWG_df$Row.names==i,]
  if(nrow(T1D_PCAWG_subset)==1){
    T1D_PCAWG_df_new[i,] <- T1D_PCAWG_subset[,paste("pcawg",eur_pcawg_id,sep="-")]
  }else{
    T1D_PCAWG_df_new[i,] <- colSums(T1D_PCAWG_subset[,paste("pcawg",eur_pcawg_id,sep="-")])
  }
}


T1D_PCAWG_Res <- as.data.frame(matrix(NA,ncol=66,nrow=66))
row.names(T1D_PCAWG_Res) <- row.names(T1D_PCAWG_df_new)
colnames(T1D_PCAWG_Res) <- row.names(T1D_PCAWG_df_new)

for(i in row.names(T1D_PCAWG_Res)){
  for(j in colnames(T1D_PCAWG_Res)){
    double_ <- t(T1D_PCAWG_df_new[i,])*t(T1D_PCAWG_df_new[j,])
    T1D_PCAWG_Res[i,j] <- sum(double_ != 0)
  }
}
T1D_PCAWG_Res <- T1D_PCAWG_Res/414
T1D_PCAWG_Res <- round(T1D_PCAWG_Res,2)
library(pheatmap)
pheatmap(T1D_PCAWG_Res,cluster_rows = T,cluster_cols = T,display_numbers = T,scale = "none")


heatmap_res <- pheatmap(T1D_NEER_Res-T1D_PCAWG_Res,cluster_rows = T,cluster_cols = T,display_numbers = F,scale = "none")

# change row.names and colnames to rs ID
T1D_snps <- read.csv("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\Maria_harmonized\\snps_t1d_genes.csv")






for(i in 1:nrow(T1D_NEER_Res)){
  row.names(T1D_NEER_Res)[i] <- T1D_snps[T1D_snps$Row.names==row.names(T1D_NEER_Res)[i],]$snp_id
}
for(i in 1:ncol(T1D_NEER_Res)){
  colnames(T1D_NEER_Res)[i] <- T1D_snps[T1D_snps$Row.names==colnames(T1D_NEER_Res)[i],]$snp_id
}

for(i in 1:nrow(T1D_PCAWG_Res)){
  row.names(T1D_PCAWG_Res)[i] <- T1D_snps[T1D_snps$Row.names==row.names(T1D_PCAWG_Res)[i],]$snp_id
}
for(i in 1:ncol(T1D_PCAWG_Res)){
  colnames(T1D_PCAWG_Res)[i] <- T1D_snps[T1D_snps$Row.names==colnames(T1D_PCAWG_Res)[i],]$snp_id
}

diff_M <- T1D_NEER_Res-T1D_PCAWG_Res
heatmap_df <- pheatmap(diff_M,cluster_rows = T,cluster_cols = T,display_numbers = F,scale = "none")

od =  hclust(dist(diff_M))$order
diff_M_new <- diff_M[heatmap_df$tree_row$order,heatmap_df$tree_col$order]
diff_M_new[lower.tri(diff_M_new)] <- NA
heatmap_res <- pheatmap(diff_M_new,cluster_rows = F,cluster_cols = F,display_numbers = F,scale = "none", type = "lower",border_color = "transparent",na_col = "transparent")
ggsave("E:\\harvard\\Lab\\Manuscript_figures\\GB_submission\\linkage.pdf",sup_fig_hist,dpi=600,width = 18,height = 26,units = "cm",device = "pdf")



