# s<-1000000
# 
# #1-correlation pop. rec. rate by prop of methylated genotype (ALL methy context together)
# methy_prop_index<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.2_unified_methy_index_",s,"_windows.csv", sep = ""))
# 
# 
# Populations<-colnames(methy_prop[[1]][[1]])
# 
# ###### per window #######
# 
# rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))
# 
# genome_table<-matrix(nrow = nrow(rec_list[[2]]), ncol = 7)
# colnames(genome_table)<-paste(1:7, "H", sep = "")
# row.names(genome_table)<-row.names(rec_list[[2]])
# 
# for (c in 1:7){ 
#   windows<-row.names(rec_list[[c]])
#   chr_corr_table<-matrix(ncol = 1, nrow = length(windows))
#   row.names(chr_corr_table)<-windows
#   chr_corr_table[]<-"-"
#   for (w in windows){
#     if (length(unique(as.numeric(rec_list[[c]][w,])))!=1){  
#       COR<-cor.test(as.numeric(rec_list[[c]][w,]), as.numeric(methy_prop_index[[c]][w,]), alternative = "two.sided", method = "pearson")  
#       if (is.na(COR$estimate)==FALSE){ if ((COR$p.value)<0.05){chr_corr_table[w,1]<-round(COR$estimate, digits = 3)}}
#     }
#   }#w
#   genome_table[windows,c]<-chr_corr_table
# }#c
# 
# write.csv(genome_table, paste("/scratch/Federico/3_RILs/5_correlation/Results/E.1.1_unified_methy_corr_",s,"_windows.csv", sep = ""))
# 
# #########################
# 
# ######## per chr ########
# 
# methy_prop_table<-matrix(ncol = 7, nrow = 45)
# colnames(methy_prop_table)<-paste(1:7, "H", sep = "")
# row.names(methy_prop_table)<-Populations
# for (c in 1:7){ for (p in Populations){methy_prop_table[p,c]<-mean(as.numeric(methy_prop_index[[c]][,p]))  }}  
# 
# rec_table<-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/A.5.1_rec_rates_V3_genome_and_chrs.csv")
# row.names(rec_table)<-rec_table[,1]; rec_table<-rec_table[,-1]
# 
# genome_table<-matrix(nrow = 1, ncol = 7)
# colnames(genome_table)<-paste(1:7, "H", sep = "")
# for (c in 1:7){
#   COR<-cor.test(as.numeric(rec_table[,c]), as.numeric(methy_prop_table[,c]), alternative = "two.sided", method = "pearson")  
#   if ((COR$p.value)<0.05){genome_table[,c]<-round(COR$estimate, digits = 3)}
# }#c    
# 
# write.csv(genome_table, paste("/scratch/Federico/3_RILs/5_correlation/Results/E.1.2_unified_methy_corr_",s,"_chrs.csv", sep = ""), row.names = FALSE)

#########################

#2-correlation pop. rec. rate by prop of methy genotype (differentiated by methy context)

library(xlsx)

s<-1000000
methy_prop<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population_win=",s,".RDS", sep = ""))
contexts<-names(methy_prop[[1]])

Populations<-colnames(methy_prop[[1]][[1]])

###### per window #######

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))

genome_table<-matrix(nrow = nrow(rec_list[[2]]), ncol = 7)
colnames(genome_table)<-paste(1:7, "H", sep = "")
row.names(genome_table)<-row.names(rec_list[[2]])

file.remove(paste("/scratch/Federico/3_RILs/5_correlation/Results/E.2.1_methy_corr_",s,"_windows.xlsx", sep = ""))
for (v in contexts){
  for (c in 1:7){ 
    windows<-row.names(rec_list[[c]])
    chr_corr_table<-matrix(ncol = 1, nrow = length(windows))
    row.names(chr_corr_table)<-windows
    chr_corr_table[]<-"-"
    row.names(methy_prop[[c]][[v]])<-windows
    for (w in windows){
      if (length(unique(as.numeric(rec_list[[c]][w,])))!=1){  
        COR<-cor.test(as.numeric(rec_list[[c]][w,]), as.numeric(methy_prop[[c]][[v]][w,]), alternative = "two.sided", method = "pearson")  
        if (is.na(COR$estimate)==FALSE){ if ((COR$p.value)<0.05){chr_corr_table[w,1]<-round(COR$estimate, digits = 3)}}
      }
    }#w
    genome_table[windows,c]<-chr_corr_table
  }#c
  write.xlsx(genome_table, paste("/scratch/Federico/3_RILs/5_correlation/Results/E.2.1_methy_corr_",s,"_windows.xlsx", sep = ""), append = TRUE, sheetName = v)
}#v

#########################

######## per chr ########

rec_table<-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/A.5.1_rec_rates_V3_genome_and_chrs.csv")
row.names(rec_table)<-rec_table[,1]; rec_table<-rec_table[,-1]

genome_table<-matrix(nrow = length(contexts), ncol = 7)
colnames(genome_table)<-paste(1:7, "H", sep = "")
row.names(genome_table)<-contexts

for (v in contexts){
  methy_prop_table<-matrix(ncol = 7, nrow = 45)
  colnames(methy_prop_table)<-paste(1:7, "H", sep = "")
  row.names(methy_prop_table)<-Populations
  for (c in 1:7){ for (p in Populations){methy_prop_table[p,c]<-mean(as.numeric(methy_prop[[c]][[v]][,p]))}}  
  
  for (c in 1:7){
    COR<-cor.test(as.numeric(rec_table[,c]), as.numeric(methy_prop_table[,c]), alternative = "two.sided", method = "pearson")  
    if ((COR$p.value)<0.05){genome_table[v,c]<-round(COR$estimate, digits = 3)}
  }#c 
  
}#v
write.csv(genome_table, paste("/scratch/Federico/3_RILs/5_correlation/Results/E.2.2_methy_corr_",s,"_chrs.csv", sep = ""), row.names = FALSE)

#########################
