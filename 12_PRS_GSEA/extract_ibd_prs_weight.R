setwd("/n/data1/hms/dbmi/zaklab/csy/PRS/ibd.npsdat")
for(i in 1:22){
  assign(paste("ibd_snps_chr_",i,sep="_"),read.delim(paste("chrom",i,".maf_0.05.ibdtr.QC3.snpinfo",sep="")))
  
}
head(ibd_snps_chr__1)

ibd_snps <- rbind(ibd_snps_chr__1,ibd_snps_chr__2,ibd_snps_chr__3,
                  ibd_snps_chr__4,ibd_snps_chr__5,ibd_snps_chr__6,
                  ibd_snps_chr__7,ibd_snps_chr__8,ibd_snps_chr__9,
                  ibd_snps_chr__10,ibd_snps_chr__11,ibd_snps_chr__12,
                  ibd_snps_chr__13,ibd_snps_chr__14,ibd_snps_chr__15,
                  ibd_snps_chr__16,ibd_snps_chr__17,ibd_snps_chr__18,
                  ibd_snps_chr__19,ibd_snps_chr__20,ibd_snps_chr__21,
                  ibd_snps_chr__22)

ibd_snps_bed <- data.frame(ibd_snps$chromosome,ibd_snps$position-1,ibd_snps$position)
write.table(ibd_snps_bed,"../IBD_sung_neer/IBD_PRS_weight.bed",quote=F,row.names=F,col.names=F,sep="\t")


# Get weight of IBD model

for(chr in 1:22){
  for(window in c(0,1000,2000,3000)){
    assign(paste("ibd_snpns_chr",chr,"window",window,sep="_"),read.table(paste("maf_0.05.ibdtr.QC3.win_",window,".adjbetahat_pg_qctool.chrom",chr,".txt",sep=""),header = T,stringsAsFactors = F))
  }
}

ibd_snps_0 <- rbind(ibd_snpns_chr_1_window_0,ibd_snpns_chr_2_window_0,ibd_snpns_chr_3_window_0,
                  ibd_snpns_chr_4_window_0,ibd_snpns_chr_5_window_0,ibd_snpns_chr_6_window_0,
                  ibd_snpns_chr_7_window_0,ibd_snpns_chr_8_window_0,ibd_snpns_chr_9_window_0,
                  ibd_snpns_chr_10_window_0,ibd_snpns_chr_11_window_0,ibd_snpns_chr_12_window_0,
                  ibd_snpns_chr_13_window_0,ibd_snpns_chr_14_window_0,ibd_snpns_chr_15_window_0,
                  ibd_snpns_chr_16_window_0,ibd_snpns_chr_17_window_0,ibd_snpns_chr_18_window_0,
                  ibd_snpns_chr_19_window_0,ibd_snpns_chr_20_window_0,ibd_snpns_chr_21_window_0,
                  ibd_snpns_chr_22_window_0)

ibd_snps_1000 <- rbind(ibd_snpns_chr_1_window_1000,ibd_snpns_chr_2_window_1000,ibd_snpns_chr_3_window_1000,
                    ibd_snpns_chr_4_window_1000,ibd_snpns_chr_5_window_1000,ibd_snpns_chr_6_window_1000,
                    ibd_snpns_chr_7_window_1000,ibd_snpns_chr_8_window_1000,ibd_snpns_chr_9_window_1000,
                    ibd_snpns_chr_10_window_1000,ibd_snpns_chr_11_window_1000,ibd_snpns_chr_12_window_1000,
                    ibd_snpns_chr_13_window_1000,ibd_snpns_chr_14_window_1000,ibd_snpns_chr_15_window_1000,
                    ibd_snpns_chr_16_window_1000,ibd_snpns_chr_17_window_1000,ibd_snpns_chr_18_window_1000,
                    ibd_snpns_chr_19_window_1000,ibd_snpns_chr_20_window_1000,ibd_snpns_chr_21_window_1000,
                    ibd_snpns_chr_22_window_1000)

ibd_snps_2000 <- rbind(ibd_snpns_chr_1_window_2000,ibd_snpns_chr_2_window_2000,ibd_snpns_chr_3_window_2000,
                    ibd_snpns_chr_4_window_2000,ibd_snpns_chr_5_window_2000,ibd_snpns_chr_6_window_2000,
                    ibd_snpns_chr_7_window_2000,ibd_snpns_chr_8_window_2000,ibd_snpns_chr_9_window_2000,
                    ibd_snpns_chr_10_window_2000,ibd_snpns_chr_11_window_2000,ibd_snpns_chr_12_window_2000,
                    ibd_snpns_chr_13_window_2000,ibd_snpns_chr_14_window_2000,ibd_snpns_chr_15_window_2000,
                    ibd_snpns_chr_16_window_2000,ibd_snpns_chr_17_window_2000,ibd_snpns_chr_18_window_2000,
                    ibd_snpns_chr_19_window_2000,ibd_snpns_chr_20_window_2000,ibd_snpns_chr_21_window_2000,
                    ibd_snpns_chr_22_window_2000)

ibd_snps_3000 <- rbind(ibd_snpns_chr_1_window_3000,ibd_snpns_chr_2_window_3000,ibd_snpns_chr_3_window_3000,
                    ibd_snpns_chr_4_window_3000,ibd_snpns_chr_5_window_3000,ibd_snpns_chr_6_window_3000,
                    ibd_snpns_chr_7_window_3000,ibd_snpns_chr_8_window_3000,ibd_snpns_chr_9_window_3000,
                    ibd_snpns_chr_10_window_3000,ibd_snpns_chr_11_window_3000,ibd_snpns_chr_12_window_3000,
                    ibd_snpns_chr_13_window_3000,ibd_snpns_chr_14_window_3000,ibd_snpns_chr_15_window_3000,
                    ibd_snpns_chr_16_window_3000,ibd_snpns_chr_17_window_3000,ibd_snpns_chr_18_window_3000,
                    ibd_snpns_chr_19_window_3000,ibd_snpns_chr_20_window_3000,ibd_snpns_chr_21_window_3000,
                    ibd_snpns_chr_22_window_3000)


ibd_snps_weight_matrix <- ibd_snps_0
ibd_snps_weight_matrix$additive_beta <- rowMeans(data.frame(ibd_snps_0$additive_beta,ibd_snps_1000$additive_beta,ibd_snps_2000$additive_beta,ibd_snps_3000$additive_beta))


write.table(ibd_snps_weight_matrix,"../IBD_sung_neer/IBD_snps_weight_sung.txt",sep="\t",row.names=F,col.names=T,quote=F)

