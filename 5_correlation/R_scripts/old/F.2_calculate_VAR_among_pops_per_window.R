#variance calculator
s<-1000000

variables<-list()
variables[["rec_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))
variables[["SVs_prop_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
TYPES_SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))
variables[["TYPES_SVs_prop_list"]]<-list()
for (c in 1:7){
  variables[["TYPES_SVs_prop_list"]][[c]]<-list()
  sub_variables<-names(TYPES_SVs_prop_list)
  for (v in sub_variables){variables[["TYPES_SVs_prop_list"]][[c]][[v]]<-TYPES_SVs_prop_list[[v]][[c]]}
}#c
variables[["GEN_DIST_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
variables[["METHY_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.1_methylation_per_population_win=1e+06.RDS", sep = ""))
variables[["PAR_DIF_METHY_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.1_methy_parents_differential_per_pop_win=1e+06.RDS", sep = ""))
variables[["DMRs_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.3_DMRs_per_pop_win=1e+06.RDS", sep = ""))
variables[["PAR_DIF_DMRs_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.2_DMRs_parents_differential_per_pop_win=1e+06.RDS", sep = ""))
variables[["DMRs_UNI_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.3_DMRs_unified_per_pop_wi_1e+06_windows.RDS", sep = ""))
variables[["PAR_DIF_DMRs_UNI_list"]]<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.3_DMRs_parents_differential_unified_per_pop_win=_",s,".RDS", sep = ""))
names(variables)<-unlist(strsplit(names(variables),"_list"))

var_list<-list()

for (c in 1:7){  
windows<-row.names(variables$rec[[c]])
var_chr<-windows; var_chr[]<-NA; names(var_chr)<-windows
for (i in names(variables)){
if (is.null(names(variables[[i]][[1]]))){ #1 
var_chr_i<-var_chr  
for (w in windows){var_chr_i[[w]]<-round(var(as.numeric(variables[[i]][[c]][w,]), na.rm = TRUE), digits = 5)}#w  
var_list[[i]]<-var_chr_i
} else {
var_list[[i]]<-list()
sub_variables<-names(variables[[i]][[1]])
for (v in sub_variables){
var_chr_i<-var_chr  
row.names(variables[[i]][[c]][[v]])<-as.numeric(row.names(variables[[i]][[c]][[v]])) #to correct a few cases
for (w in windows){var_chr_i[[w]]<-round(var(as.numeric(variables[[i]][[c]][[v]][w,]), na.rm = TRUE), digits = 5)}#w  
var_list[[i]][[v]]<-var_chr_i
}#v
}#if1
}#i
}#c

saveRDS(var_list, paste("/scratch/Federico/3_RILs/5_correlation/Results/F.2_variance_among_pops_per_window=",s,".RDS", sep = ""))