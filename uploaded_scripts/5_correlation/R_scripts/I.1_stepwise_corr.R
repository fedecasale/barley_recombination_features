#STEPWISE REGRESSION

#http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/154-stepwise-regression-essentials-in-r/

#The stepwise regression (or stepwise selection) consists of iteratively adding and removing predictors in the predictive model, in 
#order to find the subset of variables in the data set resulting in the best performing model (a model that lowers prediction error.)
library("MASS")

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

standarize<-function(x){x<-(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE); return(x)}
#normalize<-function(x){x<-(x-mean(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)); return(x)}

s<-1000000

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))

corr_table<-matrix(nrow = nrow(rec_list[[2]]), ncol = 7)
colnames(corr_table)<-paste(1:7, "H", sep = "")
row.names(corr_table)<-row.names(rec_list[[2]])
corr_table[]<-"-"

Populations<-colnames(rec_list[[1]])

#SVs variables
SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
TYPES_SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))
GENE_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/F.0_gene_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
GEN_DIST_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
#Methylation variables
met_context<-c("chg", "chh", "cpg")
METHY_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.1_methylation_per_population_win=1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(METHY_list[[c]][[v]])<-as.numeric(row.names(METHY_list[[c]][[v]]))}}
pdif_METHY_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.1_methy_parents_differential_per_pop_win=1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(pdif_METHY_list[[c]][[v]])<-as.numeric(row.names(pdif_METHY_list[[c]][[v]]))}}
DMRs_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.3_DMRs_per_pop_win=1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(DMRs_list[[c]][[v]])<-as.numeric(row.names(DMRs_list[[c]][[v]]))}}
pdif_DMRs_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.2_DMRs_parents_differential_per_pop_win=1e+06.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_list[[c]][[v]])<-as.numeric(row.names(pdif_DMRs_list[[c]][[v]]))}}
DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.3_DMRs_unified_per_pop_wi_1e+06_windows.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(DMRs_UNI_list[[c]])<-as.numeric(row.names(DMRs_UNI_list[[c]]))}}
pdif_DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.3_DMRs_parents_differential_unified_per_pop_win=_",s,".RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_UNI_list[[c]])<-as.numeric(row.names(pdif_DMRs_UNI_list[[c]]))}}

standarized_variables<-list()

#1-step.wise regression per windows
sw_list<-list()
for (c in 1:7){ cat(c);cat("-")
  standarized_variables[[c]]<-list()
  sw_list[[c]]<-list()
  windows<-row.names(rec_list[[c]])
  for (w in windows){   
    rec.rate<-as.numeric(rec_list[[c]][w,])
    variables<-data.frame(
    #SVs_prop<-as.numeric(SVs_prop_list[[c]][w,]),
    deletions<-as.numeric(TYPES_SVs_prop_list[[1]][[c]][w,]),
    duplications<-as.numeric(TYPES_SVs_prop_list[[2]][[c]][w,]),
    insertions<-as.numeric(TYPES_SVs_prop_list[[3]][[c]][w,]),
    inversions<-as.numeric(TYPES_SVs_prop_list[[4]][[c]][w,]),
    par.gen.dist<-as.numeric(GEN_DIST_list[[c]][w,]),
    cpg_DMRs<-as.numeric(DMRs_list[[c]][["cpg"]][w,]),
    chg_DMRs<-as.numeric(DMRs_list[[c]][["chg"]][w,]),
    chh_DMRs<-as.numeric(DMRs_list[[c]][["chh"]][w,]),
    par.dif.cpg<-as.numeric(pdif_DMRs_list[[c]][["cpg"]][w,]),
    par.dif.chg<-as.numeric(pdif_DMRs_list[[c]][["chg"]][w,]),
    par.dif.chh<-as.numeric(pdif_DMRs_list[[c]][["chh"]][w,])
    )
    names(variables)<-unlist(lapply(strsplit(names(variables), "....as."), FUN = function(x){return(x[1])})) 
    #################
    #normalization
    #if ((length(which(rec.rate==0))==length(rec.rate))==FALSE){rec.rate<-normalize(rec.rate)}
    variables<-apply(variables, 2, FUN = standarize)    
    variables<-as.data.frame(variables)
    #LUEGO INVESTIGAR PORQUE NO HAY GEN.DIST ENTRE ALGUNOS PADRES, BY NOW I ASSUME NO DISTANCE, ALSO FOR SOME SVs
    for (i in 1:ncol(variables)){variables[which(variables[,i]=="NaN"),i]<-0}
    #regression
    full.model<-lm(formula = rec.rate ~. , data = variables)
    if ((full.model$coefficients[[1]]==0)==FALSE){
    step.model<-stepAIC(full.model, direction = "both", trace = FALSE)
    #summary(step.model)
    sw_list[[c]][[w]]<-step.model
    }#if
    standarized_variables[[c]][[w]]<-variables
    }#w
}#c
saveRDS(variables, "/scratch/Federico/3_RILs/5_correlation/Results/I.1.0_standarized_variables_per_window.RDS")
saveRDS(sw_list, "/scratch/Federico/3_RILs/5_correlation/Results/I.1.1_step_wise_regression_per_window.RDS")

######################################

exp_variables<-names(variables)

sw_list<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/I.1.1_step_wise_regression_per_window.RDS")

#I will extract only significant variables and coefficients
alpha<-0.05

significant_variables_list<-list()
for (c in 1:7){ cat(c);cat("-")
  windows<-row.names(rec_list[[c]])
  coeff_table<-matrix(nrow = length(windows), ncol = length(exp_variables))
  row.names(coeff_table)<-windows; colnames(coeff_table)<-exp_variables
  p_value_table<-coeff_table
  significant_variables_list[[c]]<-list()
  chr_table<-matrix(ncol = 4, nrow = 0)
  for (w in windows){ #w
  if (w %in% names(sw_list[[c]])){
  significant_coefficients<-as.data.frame(summary(sw_list[[c]][[w]])$coefficients)[-1,]
  if (any(significant_coefficients[,4]<alpha)){
  significant_coefficients<-significant_coefficients[which(significant_coefficients[,4]<alpha),]
  variables<-row.names(significant_coefficients)
  #save
  for (i in variables){coeff_table[w, i]<-significant_coefficients[i,1]}
  for (i in variables){p_value_table[w, i]<-significant_coefficients[i,4]}
  win_table<-cbind(w, variables, significant_coefficients[,c(1,4)])
  chr_table<-rbind(chr_table, win_table)
  }# if any *
  }#if w 
  }#w
  significant_variables_list[[c]][["coefficients"]]<-coeff_table
  significant_variables_list[[c]][["p_values"]]<-p_value_table
  #CORRECTION FOR REGRESSION COEFFICIENTS OUT OF THE RANGE (-1:1)
  chr_table[which(chr_table[,3]<(-1)),3]<-(-1)
  chr_table[which(chr_table[,3]>1),3]<-1
  chr_table[,3]<-round(chr_table[,3], digits = 2)
  chr_table[,4]<-round(chr_table[,4], digits = 5)
  chr_table<-cbind(paste(c,"H", sep = ""), chr_table)
  colnames(chr_table)<-c("chr", "window", "variable", "reg_coeff", "P_value")
  #save
  significant_variables_list[[c]][["graphic_table"]]<-chr_table
}#c

saveRDS(significant_variables_list, "/scratch/Federico/3_RILs/5_correlation/Results/I.1.2_sw_reg_per_win_significant_variables.RDS")


####################################################

############# REDUCED VARIABLES ####################

library("MASS")

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

standarize<-function(x){x<-(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE); return(x)}
#normalize<-function(x){x<-(x-mean(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)); return(x)}

s<-1000000

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))

corr_table<-matrix(nrow = nrow(rec_list[[2]]), ncol = 7)
colnames(corr_table)<-paste(1:7, "H", sep = "")
row.names(corr_table)<-row.names(rec_list[[2]])
corr_table[]<-"-"

Populations<-colnames(rec_list[[1]])

#SVs variables
SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
GEN_DIST_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
#Methylation variables
met_context<-c("chg", "chh", "cpg")
DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.3_DMRs_unified_per_pop_wi_1e+06_windows.RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(DMRs_UNI_list[[c]])<-as.numeric(row.names(DMRs_UNI_list[[c]]))}}
pdif_DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.3_DMRs_parents_differential_unified_per_pop_win=_",s,".RDS", sep = ""))
for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_UNI_list[[c]])<-as.numeric(row.names(pdif_DMRs_UNI_list[[c]]))}}

standarized_variables<-list()

#1-step.wise regression per windows
sw_list<-list()
for (c in 1:7){ cat(c);cat("-")
  standarized_variables[[c]]<-list()
  sw_list[[c]]<-list()
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
    #################
    #normalization
    #if ((length(which(rec.rate==0))==length(rec.rate))==FALSE){rec.rate<-normalize(rec.rate)}
    variables<-apply(variables, 2, FUN = standarize)    
    variables<-as.data.frame(variables)
    #LUEGO INVESTIGAR PORQUE NO HAY GEN.DIST ENTRE ALGUNOS PADRES, BY NOW I ASSUME NO DISTANCE, ALSO FOR SOME SVs
    for (i in 1:ncol(variables)){variables[which(variables[,i]=="NaN"),i]<-0}
    #regression
    full.model<-lm(formula = rec.rate ~. , data = variables)
    if ((full.model$coefficients[[1]]==0)==FALSE){
      step.model<-stepAIC(full.model, direction = "both", trace = FALSE)
      #summary(step.model)
      sw_list[[c]][[w]]<-step.model
    }#if
    standarized_variables[[c]][[w]]<-variables
  }#w
}#c
saveRDS(standarized_variables, "/scratch/Federico/3_RILs/5_correlation/Results/I.1.3_REDUCED_standarized_variables_per_window.RDS")
saveRDS(sw_list, "/scratch/Federico/3_RILs/5_correlation/Results/I.1.4_REDUCED_step_wise_regression_per_window.RDS")

######################################

exp_variables<-names(variables)

#I will extract only significant variables and coefficients
alpha<-0.05

significant_variables_list<-list()
for (c in 1:7){ cat(c);cat("-")
  windows<-row.names(rec_list[[c]])
  coeff_table<-matrix(nrow = length(windows), ncol = length(exp_variables))
  row.names(coeff_table)<-windows; colnames(coeff_table)<-exp_variables
  p_value_table<-coeff_table
  significant_variables_list[[c]]<-list()
  chr_table<-matrix(ncol = 4, nrow = 0)
  for (w in windows){ #w
    if (w %in% names(sw_list[[c]])){
      significant_coefficients<-as.data.frame(summary(sw_list[[c]][[w]])$coefficients)[-1,]
      if (any(significant_coefficients[,4]<alpha)){
        significant_coefficients<-significant_coefficients[which(significant_coefficients[,4]<alpha),]
        variables<-row.names(significant_coefficients)
        #save
        for (i in variables){coeff_table[w, i]<-significant_coefficients[i,1]}
        for (i in variables){p_value_table[w, i]<-significant_coefficients[i,4]}
        win_table<-cbind(w, variables, significant_coefficients[,c(1,4)])
        chr_table<-rbind(chr_table, win_table)
      }# if any *
    }#if w 
  }#w
  significant_variables_list[[c]][["coefficients"]]<-coeff_table
  significant_variables_list[[c]][["p_values"]]<-p_value_table
  #CORRECTION FOR REGRESSION COEFFICIENTS OUT OF THE RANGE (-1:1)
  chr_table[which(chr_table[,3]<(-1)),3]<-(-1)
  chr_table[which(chr_table[,3]>1),3]<-1
  chr_table[,3]<-round(chr_table[,3], digits = 2)
  chr_table[,4]<-round(chr_table[,4], digits = 5)
  chr_table<-cbind(paste(c,"H", sep = ""), chr_table)
  colnames(chr_table)<-c("chr", "window", "variable", "reg_coeff", "P_value")
  #save
  significant_variables_list[[c]][["graphic_table"]]<-chr_table
}#c

saveRDS(significant_variables_list, "/scratch/Federico/3_RILs/5_correlation/Results/I.1.5_REDUCED_sw_reg_per_win_significant_variables.RDS")

