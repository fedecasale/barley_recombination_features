
#continue correlation with SV: genotypes are more similar or different to each other based on SVs, 
#and how that impact recombination -> how different are the parents 
#-> how they impact in the thecombination rate of the offspring populations. 
#simple correlation. use the 23 parents. chrs, windows. maybe, divide structural variants in categories. 

#1-correlation pop. rec. rate by prop of sv genotype (ALL SVs together)

s<-1000000
SVs_prop<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))

Populations<-colnames(SVs_prop[[1]][[1]])

###### per window #######

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))

genome_table<-matrix(nrow = nrow(rec_list[[2]]), ncol = 7)
colnames(genome_table)<-paste(1:7, "H", sep = "")
row.names(genome_table)<-row.names(rec_list[[2]])

for (c in 1:7){ 
  windows<-row.names(rec_list[[c]])
  chr_corr_table<-matrix(ncol = 1, nrow = length(windows))
  row.names(chr_corr_table)<-windows
  chr_corr_table[]<-"-"
  for (w in windows){
  if (length(unique(as.numeric(rec_list[[c]][w,])))!=1){  
  COR<-cor.test(as.numeric(rec_list[[c]][w,]), as.numeric(SVs_prop[[c]][w,]), alternative = "two.sided", method = "pearson")  
  if (is.na(COR$estimate)==FALSE){ if ((COR$p.value)<0.05){chr_corr_table[w,1]<-round(COR$estimate, digits = 3)}}
  }
  }#w
  genome_table[windows,c]<-chr_corr_table
}#c

write.csv(genome_table, paste("/scratch/Federico/3_RILs/5_correlation/Results/D.1.1_all_SVs_corr_",s,"_windows.csv", sep = ""))

#########################

######## per chr ########

sv_prop_table<-matrix(ncol = 7, nrow = 45)
colnames(sv_prop_table)<-paste(1:7, "H", sep = "")
row.names(sv_prop_table)<-Populations
for (c in 1:7){ for (p in Populations){sv_prop_table[p,c]<-mean(as.numeric(SVs_prop[[c]][,p]))  }}  

rec_table<-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/A.5.1_rec_rates_V3_genome_and_chrs.csv")
row.names(rec_table)<-rec_table[,1]; rec_table<-rec_table[,-1]
  
genome_table<-matrix(nrow = 1, ncol = 7)
colnames(genome_table)<-paste(1:7, "H", sep = "")
for (c in 1:7){
  COR<-cor.test(as.numeric(rec_table[,c]), as.numeric(sv_prop_table[,c]), alternative = "two.sided", method = "pearson")  
    if ((COR$p.value)<0.05){genome_table[,c]<-round(COR$estimate, digits = 3)}
}#c    

write.csv(genome_table, paste("/scratch/Federico/3_RILs/5_correlation/Results/D.1.2_all_SVs_corr_",s,"_chrs.csv", sep = ""), row.names = FALSE)

#########################



#2-correlation pop. rec. rate by prop of sv genotype (differentiated by SV type)

library(xlsx)

s<-1000000
SVs_prop<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))
SVs<-names(SVs_prop)

Populations<-colnames(SVs_prop[[1]][[1]])

###### per window #######

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))

genome_table<-matrix(nrow = nrow(rec_list[[2]]), ncol = 7)
colnames(genome_table)<-paste(1:7, "H", sep = "")
row.names(genome_table)<-row.names(rec_list[[2]])

file.remove(paste("/scratch/Federico/3_RILs/5_correlation/Results/D.2.1_SVs_corr_",s,"_windows.xlsx", sep = ""))
for (v in SVs){
for (c in 1:7){ 
  windows<-row.names(rec_list[[c]])
  chr_corr_table<-matrix(ncol = 1, nrow = length(windows))
  row.names(chr_corr_table)<-windows
  chr_corr_table[]<-"-"
  for (w in windows){
    if (length(unique(as.numeric(rec_list[[c]][w,])))!=1){  
      COR<-cor.test(as.numeric(rec_list[[c]][w,]), as.numeric(SVs_prop[[v]][[c]][w,]), alternative = "two.sided", method = "pearson")  
      if (is.na(COR$estimate)==FALSE){ if ((COR$p.value)<0.05){chr_corr_table[w,1]<-round(COR$estimate, digits = 3)}}
    }
  }#w
  genome_table[windows,c]<-chr_corr_table
}#c
write.xlsx(genome_table, paste("/scratch/Federico/3_RILs/5_correlation/Results/D.2.1_SVs_corr_",s,"_windows.xlsx", sep = ""), append = TRUE, sheetName = v)
}#v


#########################

######## per chr ########

rec_table<-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/A.5.1_rec_rates_V3_genome_and_chrs.csv")
row.names(rec_table)<-rec_table[,1]; rec_table<-rec_table[,-1]

genome_table<-matrix(nrow = length(SVs), ncol = 7)
colnames(genome_table)<-paste(1:7, "H", sep = "")
row.names(genome_table)<-SVs

for (v in SVs){
sv_prop_table<-matrix(ncol = 7, nrow = 45)
colnames(sv_prop_table)<-paste(1:7, "H", sep = "")
row.names(sv_prop_table)<-Populations
for (c in 1:7){ for (p in Populations){sv_prop_table[p,c]<-mean(as.numeric(SVs_prop[[v]][[c]][,p]))}}  

for (c in 1:7){
  COR<-cor.test(as.numeric(rec_table[,c]), as.numeric(sv_prop_table[,c]), alternative = "two.sided", method = "pearson")  
  if ((COR$p.value)<0.05){genome_table[v,c]<-round(COR$estimate, digits = 3)}
}#c 

}#v
write.csv(genome_table, paste("/scratch/Federico/3_RILs/5_correlation/Results/D.2.2_SVs_corr_",s,"_chrs.csv", sep = ""), row.names = FALSE)

#########################
