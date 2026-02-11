library("BSDA")
library("xlsx")

###USO LAS VARIABLE PERO NO STANDARIZED
#SVs variables
s<-1000000
rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))
SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
GEN_DIST_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
#Methylation variables
met_context<-c("chg", "chh", "cpg")
DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.3_DMRs_unified_per_pop_wi_1e+06_windows.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(DMRs_UNI_list[[c]])<-as.numeric(row.names(DMRs_UNI_list[[c]]))}}
pdif_DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.3_DMRs_parents_differential_unified_per_pop_win=_",s,".RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_UNI_list[[c]])<-as.numeric(row.names(pdif_DMRs_UNI_list[[c]]))}}

non_standarized_variables<-list()
for (c in 1:7){ cat(c);cat("-")
  non_standarized_variables[[c]]<-list()
  windows<-row.names(rec_list[[c]])
  for (w in windows){   
    rec.rate<-as.numeric(rec_list[[c]][w,])
    variables<-data.frame(
      SVs_prop<-as.numeric(SVs_prop_list[[c]][w,]),
      par.gen.dist<-as.numeric(GEN_DIST_list[[c]][w,]),
      DMRs<-as.numeric(DMRs_UNI_list[[c]][w,]),
      pdif_DMRs<-as.numeric(pdif_DMRs_UNI_list[[c]][w,])
    )
    names(variables)<-unlist(lapply(strsplit(names(variables), "....as."), FUN = function(x){return(x[1])})) 
    non_standarized_variables[[c]][[w]]<-variables
  }#w
}#c

variables_data<-non_standarized_variables
#variables_data<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/I.1.3_REDUCED_standarized_variables_per_window.RDS")
#sw_list<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/I.1.4_REDUCED_step_wise_regression_per_window.RDS")
significant_variables_list<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/I.1.5_REDUCED_sw_reg_per_win_significant_variables.RDS")

variables<-colnames(variables_data[[1]][[1]])


###############################################################

#check if when a variable (v) is not sifnificant when others (o) do,
#the v mean for such window windows is lower than the v mean for the chr.

results_table<-matrix(nrow = 7, ncol = length(variables)); 
row.names(results_table)<-paste(1:7, "H", sep = ""); colnames(results_table)<-variables

for (v in variables){
for (c in 1:7){

  #vector of variable values for all windows
  chr_variable<-c()
  for (w in 1:length(variables_data[[c]])){chr_variable<-c(chr_variable, variables_data[[c]][[w]][[v]])}
  #variables_data[[c]][[1]][[v]]
  
  #vector of variable values for all windows where it doesnt a appear as significant but other variables do 
  tabla<-as.matrix(significant_variables_list[[c]]$coefficients)
  #clean windows with no variables
  tabla<-tabla[names(which(apply(tabla, 1, function(x){return(any(complete.cases(x)==TRUE))})==TRUE)),]
  #clean windows with the variable as significant
  tabla<-tabla[which(is.na(tabla[,v])),]
  #get values from remaining windows for the variable
  sign_wins_variable<-c()
  for (w in row.names(tabla)){ sign_wins_variable<-c(sign_wins_variable, variables_data[[c]][[w]][[v]])}#w
  
  #check if the significant window means are less than the whole set of windows in the chr
  #results_table[c,v]<-t.test(x = sign_wins_variable, y = chr_variable, alternative = "less")$p.value
  #results_table[c,v]<-z.test(x = sign_wins_variable, y = chr_variable, alternative = "less", sigma.x = var(sign_wins_variable), sigma.y = var(chr_variable))$p.value)
  results_table[c,v]<-wilcox.test(x = sign_wins_variable, y = chr_variable, alternative = "less")$p.value

  
}#c
}#v

results_table[]<-round(as.numeric(results_table), digits = 3)
write.csv(results_table, "/scratch/Federico/3_RILs/5_correlation/Results/I.3.1_when_v_doesnt_appear_v_mean_is_low.csv")

##############################

#check if when a variable (v) is sifnificant
#the v mean for such window windows is greater than the v mean for the chr.

results_table<-matrix(nrow = 7, ncol = length(variables)); 
row.names(results_table)<-paste(1:7, "H", sep = ""); colnames(results_table)<-variables

for (v in variables){
  for (c in 1:7){
    
    #vector of variable values for all windows
    chr_variable<-c()
    for (w in 1:length(variables_data[[c]])){chr_variable<-c(chr_variable, variables_data[[c]][[w]][[v]])}
    #variables_data[[c]][[1]][[v]]
    
    #vector of variable values for all windows where it appears as significant 
    tabla<-as.matrix(significant_variables_list[[c]]$coefficients)
    #clean windows with no variables
    tabla<-tabla[names(which(apply(tabla, 1, function(x){return(any(complete.cases(x)==TRUE))})==TRUE)),]
    #clean windows with the variable as significant
    tabla<-tabla[which(is.na(tabla[,v])==FALSE),]
    #get values from remaining windows for the variable
    sign_wins_variable<-c()
    for (w in row.names(tabla)){ sign_wins_variable<-c(sign_wins_variable, variables_data[[c]][[w]][[v]])}#w
    
    #check if the significant window means are greater than the whole set of windows in the chr
    results_table[c,v]<-wilcox.test(x = sign_wins_variable, y = chr_variable, alternative = "greater")$p.value
    
    
  }#c
}#v

results_table[]<-round(as.numeric(results_table), digits = 3)
write.csv(results_table, "/scratch/Federico/3_RILs/5_correlation/Results/I.3.2_when_v_appears_v_mean_is_high.csv")

##################

#check if when a variable (v) is not sifnificant when others (o) do,
#if the v VAR is lower than such of o
#if not, if the v mean for the window is lower than the v mean for the chr.

for (v in variables){
  file.remove(paste("/scratch/Federico/3_RILs/5_correlation/Results/I.3.3_",v,".xlsx", sep = ""))
  for (c in 1:7){
    #vector of variable values for all windows
    chr_variable<-c(); for (w in 1:length(variables_data[[c]])){chr_variable<-c(chr_variable, variables_data[[c]][[w]][[v]])}
    #vector of variable values for all windows where it appears as significant 
    tabla_chr<-as.matrix(significant_variables_list[[c]]$coefficients)
    #clean windows with no variables
    tabla_chr<-tabla_chr[names(which(apply(tabla_chr, 1, function(x){return(any(complete.cases(x)==TRUE))})==TRUE)),]
    #clean windows with the variable as significant
    tabla_chr<-tabla_chr[which(is.na(tabla_chr[,v])),]
    #remaining windows 
    for (w in row.names(tabla_chr)){ 
      other_variables<-variables[-which(is.na(tabla_chr[w,]))]
      for (o in other_variables){ 
      tabla_chr[w,o]<-var.test(x = variables_data[[c]][[w]][[v]], y = variables_data[[c]][[w]][[o]], alternative = "less")$p.value
      if (as.numeric(tabla_chr[w,o])<0.05){tabla_chr[w,o]<-"HIGHER_var"} else {tabla_chr[w,o]<-"EQUAL_var"}
      }#o
    if (any(tabla_chr[w,o]=="EQUAL_var")){
    #check if the significant window is less than the whole set of windows in the chr
    tabla_chr[w,v]<-wilcox.test(x = variables_data[[c]][[w]][[v]], y = chr_variable, alternative = "less")$p.value  
    if (as.numeric(tabla_chr[w,v])<0.05){tabla_chr[w,v]<-"low_MEAN"} else {tabla_chr[w,v]<-"equal_var_equal_mean"}
    } else {tabla_chr[w,v]<-"LOWER_var"}  
    }#w
  write.xlsx(tabla_chr, paste("/scratch/Federico/3_RILs/5_correlation/Results/I.3.3_",v,".xlsx", sep = ""), append = TRUE, sheetName = paste(c, "H", sep = ""))
  }#c
}#v

##############################

