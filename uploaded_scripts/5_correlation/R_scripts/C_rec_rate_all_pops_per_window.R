
s<-1000000
rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
Populations<-names(rec_list)  

rec_list2<-list()
for (c in 1:7){ 
  windows<-rec_list[[1]][[c]][,2]
  tabla_chr<-matrix(nrow = length(windows), ncol = 45)
  row.names(tabla_chr)<-windows; colnames(tabla_chr)<-Populations
  for (w in windows){
    for (p in names(rec_list)){
    tabla_chr[w,p]<-rec_list[[p]][[c]][which(rec_list[[p]][[c]][,2]==w),3]
    }#p
  }#w
rec_list2[[c]]<-tabla_chr  
}#c

saveRDS(rec_list2, paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))
