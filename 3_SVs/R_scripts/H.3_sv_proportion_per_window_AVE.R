s<-1000000

SVs_list<-list()
#SVs summed prop
SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
Populations<-colnames(SVs_prop_list[[1]])
SVs_list[["SVs_prop"]]<-list()
for (c in 1:7){
  windows<-row.names(SVs_prop_list[[c]]); names(windows)<-windows; windows[]<-NA
  for (w in names(windows)){windows[w]<-round(mean(as.numeric(SVs_prop_list[[c]][w,]), na.rm = TRUE), digits = 3)}
  SVs_list[["SVs_prop"]][[c]]<-windows
}

#SVs prop by type
TYPES_SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))
SVs_list[["SVs_type_prop"]]<-list()
for (v in names(TYPES_SVs_prop_list)){
  SVs_list[["SVs_type_prop"]][[v]]<-list()
  for (c in 1:7){
    windows<-row.names(SVs_prop_list[[c]]); names(windows)<-windows; windows[]<-NA
    for (w in names(windows)){windows[w]<-round(mean(as.numeric(TYPES_SVs_prop_list[[v]][[c]][w,]), na.rm = TRUE), digits = 3)}
    SVs_list[["SVs_type_prop"]][[v]][[c]]<-windows
  }
}

#Genetic distance among parents
GEN_DIST_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
SVs_list[["gen_dist"]]<-list()
for (c in 1:7){
  windows<-row.names(SVs_prop_list[[c]]); names(windows)<-windows; windows[]<-NA
  for (w in names(windows)){windows[w]<-round(mean(as.numeric(GEN_DIST_list[[c]][w,]), na.rm = TRUE), digits = 3)}
  SVs_list[["gen_dist"]][[c]]<-windows
}

#SVs_list[["gene_prop"]]<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/F.0_gene_proportion_per_window=",s,"_all_pops.RDS", sep = ""))

saveRDS(SVs_list, "/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.3.1_sv_proportion_per_window_AVE.RDS")

################
#with ALTERNATIVE METHOD

s<-1000000

SVs_list<-list()
#SVs summed prop
SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
Populations<-colnames(SVs_prop_list[[1]])
SVs_list[["SVs_prop"]]<-list()
for (c in 1:7){
  windows<-row.names(SVs_prop_list[[c]]); names(windows)<-windows; windows[]<-NA
  for (w in names(windows)){windows[w]<-round(mean(as.numeric(SVs_prop_list[[c]][w,]), na.rm = TRUE), digits = 3)}
  SVs_list[["SVs_prop"]][[c]]<-windows
}

#SVs prop by type
TYPES_SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))
SVs_list[["SVs_type_prop"]]<-list()
for (v in names(TYPES_SVs_prop_list)){
  SVs_list[["SVs_type_prop"]][[v]]<-list()
  for (c in 1:7){
    windows<-row.names(SVs_prop_list[[c]]); names(windows)<-windows; windows[]<-NA
    for (w in names(windows)){windows[w]<-round(mean(as.numeric(TYPES_SVs_prop_list[[v]][[c]][w,]), na.rm = TRUE), digits = 3)}
    SVs_list[["SVs_type_prop"]][[v]][[c]]<-windows
  }
}

#Genetic distance among parents
GEN_DIST_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
SVs_list[["gen_dist"]]<-list()
for (c in 1:7){
  windows<-row.names(SVs_prop_list[[c]]); names(windows)<-windows; windows[]<-NA
  for (w in names(windows)){windows[w]<-round(mean(as.numeric(GEN_DIST_list[[c]][w,]), na.rm = TRUE), digits = 3)}
  SVs_list[["gen_dist"]][[c]]<-windows
}

#SVs_list[["gene_prop"]]<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/F.0_gene_proportion_per_window=",s,"_all_pops.RDS", sep = ""))

saveRDS(SVs_list, "/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.3.2_sv_proportion_per_window_AVE.RDS")