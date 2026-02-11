met_context<-c("cpg", "chg", "chh")

methy_list<-list()
#Methylation variables
METHY_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.1_methylation_per_population_win=1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(METHY_list[[c]][[v]])<-as.numeric(row.names(METHY_list[[c]][[v]]))}}
methy_list[["methylation"]]<-list()
for (c in 1:7){
  tabla<-matrix(nrow = nrow(METHY_list[[c]][[1]]), ncol = 3)
  windows<-as.numeric(row.names(METHY_list[[c]][[1]]))
  row.names(tabla)<-windows
  colnames(tabla)<-met_context
   for (v in met_context){
   for (w in row.names(tabla)){tabla[w,v]<-mean(METHY_list[[c]][[v]][w,], na.rm = TRUE)}    
   }#v
  methy_list[["methylation"]][[c]]<-tabla
}#c

pdif_METHY_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.1_methy_parents_differential_per_pop_win=1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(pdif_METHY_list[[c]][[v]])<-as.numeric(row.names(pdif_METHY_list[[c]][[v]]))}}
methy_list[["methylation_par_dif"]]<-list()
for (c in 1:7){
  tabla<- methy_list[["methylation"]][[c]]; tabla[]<-NA
  for (v in met_context){
    for (w in row.names(tabla)){tabla[w,v]<-mean(pdif_METHY_list[[c]][[v]][w,], na.rm = TRUE)}    
  }
  methy_list[["methylation_par_dif"]][[c]]<-tabla
}#c


DMRs_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.3_DMRs_per_pop_win=1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(DMRs_list[[c]][[v]])<-as.numeric(row.names(DMRs_list[[c]][[v]]))}}
methy_list[["DMRs"]]<-list()
for (c in 1:7){
  tabla<- methy_list[["methylation"]][[c]]; tabla[]<-NA
  for (v in met_context){
    for (w in row.names(tabla)){tabla[w,v]<-mean(DMRs_list[[c]][[v]][w,], na.rm = TRUE)}    
  }
  methy_list[["DMRs"]][[c]]<-tabla
}#c


pdif_DMRs_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.2_DMRs_parents_differential_per_pop_win=1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_list[[c]][[v]])<-as.numeric(row.names(pdif_DMRs_list[[c]][[v]]))}}
methy_list[["DMRs_par_dif"]]<-list()
for (c in 1:7){
  tabla<- methy_list[["methylation"]][[c]]; tabla[]<-NA
  for (v in met_context){
    for (w in row.names(tabla)){tabla[w,v]<-mean(pdif_DMRs_list[[c]][[v]][w,], na.rm = TRUE)}    
  }
  methy_list[["DMRs_par_dif"]][[c]]<-tabla
}#c

DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.3_DMRs_unified_per_pop_wi_1e+06_windows.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(DMRs_UNI_list[[c]])<-as.numeric(row.names(DMRs_UNI_list[[c]]))}}
methy_list[["DMRs_unified"]]<-list()
for (c in 1:7){
  tabla<- methy_list[["methylation"]][[c]]; tabla[]<-NA
    for (w in row.names(tabla)){tabla[w,1]<-mean(DMRs_UNI_list[[c]][w,], na.rm = TRUE)}    
  tabla<-tabla[,1]
  methy_list[["DMRs_unified"]][[c]]<-tabla
}#c

pdif_DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.3_DMRs_parents_differential_unified_per_pop_win=_1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_UNI_list[[c]])<-as.numeric(row.names(pdif_DMRs_UNI_list[[c]]))}}
methy_list[["DMRs_unified_par_dif"]]<-list()
for (c in 1:7){
  tabla<- methy_list[["methylation"]][[c]]; tabla[]<-NA
    for (w in row.names(tabla)){tabla[w,1]<-mean(pdif_DMRs_UNI_list[[c]][w,], na.rm = TRUE)}    
  tabla<-tabla[,1]
  methy_list[["DMRs_unified_par_dif"]][[c]]<-tabla
}#c

saveRDS(methy_list, "/scratch/Federico/3_RILs/4_methylation/Results/D_methylation_pop_average_list.RDS")
