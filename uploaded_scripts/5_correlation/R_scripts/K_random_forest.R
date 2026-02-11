library(MASS)
library(randomForest)
library(party)
#library(datasets)
#library(caret)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

standarize<-function(x){x<-(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE); return(x)}
#normalize<-function(x){x<-(x-mean(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)); return(x)}

s<-1000000

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))

Populations<-colnames(rec_list[[1]])

#SIMPLIFIED VARIABLES
#dir.create("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest")
# #SVs variables
# SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
# GEN_DIST_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
# #Methylation variables
# met_context<-c("chg", "chh", "cpg")
# DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/viejo/B.4.3_DMRs_unified_per_pop_wi_1e+06_windows.RDS", sep = ""))
# for (c in 1:7){for (v in met_context){row.names(DMRs_UNI_list[[c]])<-as.numeric(row.names(DMRs_UNI_list[[c]]))}}
# pdif_DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.3_DMRs_parents_differential_unified_per_pop_win=_",s,".RDS", sep = ""))
# for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_UNI_list[[c]])<-as.numeric(row.names(pdif_DMRs_UNI_list[[c]]))}}
# 
# standarized_variables<-list()
# #1-variables per windows
# for (c in 1:7){ cat(c);cat("-")
#   standarized_variables[[c]]<-list()
#   windows<-row.names(rec_list[[c]])
#   for (w in windows){
#     rec.rate<-as.numeric(rec_list[[c]][w,])
#     variables<-data.frame(
#       SVs_prop<-as.numeric(SVs_prop_list[[c]][w,]),
#       par.gen.dist<-as.numeric(GEN_DIST_list[[c]][w,]),
#       DMRs<-as.numeric(DMRs_UNI_list[[c]][w,]),
#       pdif_DMRs<-as.numeric(pdif_DMRs_UNI_list[[c]][w,])
#     )
#     names(variables)<-unlist(lapply(strsplit(names(variables), "....as."), FUN = function(x){return(x[1])}))
#     #################
#     #normalization
#     #if ((length(which(rec.rate==0))==length(rec.rate))==FALSE){rec.rate<-normalize(rec.rate)}
#     variables<-apply(variables, 2, FUN = standarize)
#     variables<-as.data.frame(variables)
#     #LUEGO INVESTIGAR PORQUE NO HAY GEN.DIST ENTRE ALGUNOS PADRES, BY NOW I ASSUME NO DISTANCE, ALSO FOR SOME SVs
#     for (i in 1:ncol(variables)){variables[which(variables[,i]=="NaN"),i]<-0}
# 
#     DATA<-cbind(rec.rate, variables)
#     row.names(DATA)<-Populations
# 
#     #training set 70%, #validation set 30%
#     training_pops<-Populations[sample(1:45, 30, replace = FALSE)]
#     validation_pops<-Populations[-which(Populations%in%training_pops)]
# 
#     training_set<-DATA[training_pops,]
#     validation_set<-DATA[validation_pops,]
# 
#     #it is important to generate levels in the data.frame for all the variables
#     #for (i in 1:ncol(training_set)){training_set[[i]]<-as.factor(training_set[[i]])}
#     #for (i in 1:ncol(validation_set)){validation_set[[i]]<-as.factor(validation_set[[i]])}
# 
#     standarized_variables[[c]][[w]]<-list()
#     standarized_variables[[c]][[w]][["full_set"]]<-DATA
#     standarized_variables[[c]][[w]][["training_set"]]<-training_set
#     standarized_variables[[c]][[w]][["validation_set"]]<-validation_set
#   }#w
# }#c
# saveRDS(standarized_variables, "/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.0_standarized_variables_SIMPLIFIED.RDS")
# rm(list = c("DMRs", "DMRs_UNI_list", "par.gen.dist", "GEN_DIST_list","pdif_DMRs", "pdif_DMRs_UNI_list", "SVs_prop", "SVs_prop_list", "rec.rate"))
# rm(list = c("DATA", "training_pops", "training_set", "validation_pops", "validation_set", "variables"))

#FULL VARIABLES
# #SVs variables
# SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
# TYPES_SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))
# GENE_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/F.0_gene_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
# GEN_DIST_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
# #Methylation variables
# met_context<-c("chg", "chh", "cpg")
# METHY_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/viejo/B.1_methylation_per_population_win=1e+06.RDS", sep = ""))
# for (c in 1:7){for (v in met_context){row.names(METHY_list[[c]][[v]])<-as.numeric(row.names(METHY_list[[c]][[v]]))}}
# pdif_METHY_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.1_methy_parents_differential_per_pop_win=1e+06.RDS", sep = ""))
# for (c in 1:7){for (v in met_context){row.names(pdif_METHY_list[[c]][[v]])<-as.numeric(row.names(pdif_METHY_list[[c]][[v]]))}}
# DMRs_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/viejo/B.3_DMRs_per_pop_win=1e+06.RDS", sep = ""))
# for (c in 1:7){for (v in met_context){row.names(DMRs_list[[c]][[v]])<-as.numeric(row.names(DMRs_list[[c]][[v]]))}}
# pdif_DMRs_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.2_DMRs_parents_differential_per_pop_win=1e+06.RDS", sep = ""))
# for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_list[[c]][[v]])<-as.numeric(row.names(pdif_DMRs_list[[c]][[v]]))}}
# DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/viejo/B.4.3_DMRs_unified_per_pop_wi_1e+06_windows.RDS", sep = ""))
# for (c in 1:7){for (v in met_context){row.names(DMRs_UNI_list[[c]])<-as.numeric(row.names(DMRs_UNI_list[[c]]))}}
# pdif_DMRs_UNI_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.3_DMRs_parents_differential_unified_per_pop_win=_",s,".RDS", sep = ""))
# for (c in 1:7){for (v in met_context){row.names(pdif_DMRs_UNI_list[[c]])<-as.numeric(row.names(pdif_DMRs_UNI_list[[c]]))}}
# 
# standarized_variables<-list()
# for (c in 1:7){ cat(c);cat("-")
#   standarized_variables[[c]]<-list()
#   windows<-row.names(rec_list[[c]])
#   for (w in windows){
#     rec.rate<-as.numeric(rec_list[[c]][w,])
#     variables<-data.frame(
#       #SVs_prop<-as.numeric(SVs_prop_list[[c]][w,]),
#       deletions<-as.numeric(TYPES_SVs_prop_list[[1]][[c]][w,]),
#       duplications<-as.numeric(TYPES_SVs_prop_list[[2]][[c]][w,]),
#       insertions<-as.numeric(TYPES_SVs_prop_list[[3]][[c]][w,]),
#       inversions<-as.numeric(TYPES_SVs_prop_list[[4]][[c]][w,]),
#       par.gen.dist<-as.numeric(GEN_DIST_list[[c]][w,]),
#       cpg_DMRs<-as.numeric(DMRs_list[[c]][["cpg"]][w,]),
#       chg_DMRs<-as.numeric(DMRs_list[[c]][["chg"]][w,]),
#       chh_DMRs<-as.numeric(DMRs_list[[c]][["chh"]][w,]),
#       par.dif.cpg<-as.numeric(pdif_DMRs_list[[c]][["cpg"]][w,]),
#       par.dif.chg<-as.numeric(pdif_DMRs_list[[c]][["chg"]][w,]),
#       par.dif.chh<-as.numeric(pdif_DMRs_list[[c]][["chh"]][w,])
#     )
#     names(variables)<-unlist(lapply(strsplit(names(variables), "....as."), FUN = function(x){return(x[1])}))
#     #################
#     #normalization
#     #if ((length(which(rec.rate==0))==length(rec.rate))==FALSE){rec.rate<-normalize(rec.rate)}
#     variables<-apply(variables, 2, FUN = standarize)
#     variables<-as.data.frame(variables)
#     #LUEGO INVESTIGAR PORQUE NO HAY GEN.DIST ENTRE ALGUNOS PADRES, BY NOW I ASSUME NO DISTANCE, ALSO FOR SOME SVs
#     for (i in 1:ncol(variables)){variables[which(variables[,i]=="NaN"),i]<-0}
# 
#     DATA<-cbind(rec.rate, variables)
#     row.names(DATA)<-Populations
# 
#     #training set 70%, #validation set 30%
#     training_pops<-Populations[sample(1:45, 30, replace = FALSE)]
#     validation_pops<-Populations[-which(Populations%in%training_pops)]
# 
#     training_set<-DATA[training_pops,]
#     validation_set<-DATA[validation_pops,]
# 
#     #it is important to generate levels in the data.frame for all the variables
#     #for (i in 1:ncol(training_set)){training_set[[i]]<-as.factor(training_set[[i]])}
#     #for (i in 1:ncol(validation_set)){validation_set[[i]]<-as.factor(validation_set[[i]])}
# 
#     standarized_variables[[c]][[w]]<-list()
#     standarized_variables[[c]][[w]][["full_set"]]<-DATA
#     standarized_variables[[c]][[w]][["training_set"]]<-training_set
#     standarized_variables[[c]][[w]][["validation_set"]]<-validation_set
#   }#w
# }#c
# saveRDS(standarized_variables, "/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.0_standarized_variables_FULL.RDS")

standarized_variables<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.0_standarized_variables_SIMPLIFIED.RDS")
#standarized_variables<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.0_standarized_variables_FULL.RDS")

m<-4 #SIMPLE: 1, 2, 3 ,4  #FULL: 1, 3, 4=default, 6, 9

MAX_NODES<-c(4,8,12,16)

for (n in MAX_NODES){ #n

max_terminal_nodes<-n #4, 8, 12, 16 -> 16 es lo mismo que dejarlo al maximo
j<-paste("mtry=",m,"_maxnodes=",max_terminal_nodes, sep = "")

#train models
models<-list()
for (c in 1:7){ cat(c);cat(": ")   
  windows<-row.names(rec_list[[c]])
  models[[c]]<-list()
  for (w in windows){ #cat(which(windows==w)/length(windows)); cat("-")   
  
    models[[c]][[w]]<-list()
    training_set<-standarized_variables[[c]][[w]][["training_set"]]
    
    full.model<-lm(formula = rec.rate ~., data = training_set)
    if (full.model$coefficients[[1]]!=0){
    
    #STEPWISE regression
    step.model<-stepAIC(full.model, direction = "both", trace = FALSE)
    models[[c]][[w]][["SW"]]<-step.model
    #training_set<-training_set[1:10,]
    #RANDOM FOREST regression
    models[[c]][[w]][["RF"]] <- randomForest(rec.rate~., data=training_set, proximity=TRUE, importance = TRUE, mtry = m, maxnodes = max_terminal_nodes)
    }#if
  }#w
}#c
#saveRDS(models, paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.1_trained_models_",j,".RDS", sep = ""))   

#check predictive accuracy
results_tabla<-matrix(ncol = 3, nrow = 8); colnames(results_tabla)<-c("SW","RF","SW>RF")
MSE_list<-list()
for (c in 1:7){ cat(c);cat("-")   
  windows<-row.names(rec_list[[c]])
  tabla<-matrix(ncol = 3, nrow = length(windows)); row.names(tabla)<-windows; colnames(tabla)<-c("SW","RF","SW>RF")
  for (w in windows){ #cat(which(windows==w)/length(windows)); cat("-")
    if (is.null(models[[c]][[w]]$SW)==FALSE){
    validation_set<-standarized_variables[[c]][[w]][["validation_set"]]
    MSE_SW<-mean((validation_set$rec.rate-predict(models[[c]][[w]]$SW, validation_set))^2); tabla[w,"SW"]<-MSE_SW
    MSE_RF<-mean((validation_set$rec.rate-predict(models[[c]][[w]]$RF, validation_set))^2); tabla[w,"RF"]<-MSE_RF
    tabla[w,"SW>RF"]<-MSE_SW>MSE_RF
    } #if null
  }#w
MSE_list[[c]]<-tabla  
results_tabla[c,1]<-mean(tabla[,"SW"], na.rm = TRUE)
results_tabla[c,2]<-mean(tabla[,"RF"], na.rm = TRUE)
results_tabla[c,3]<-length(which(tabla[,"SW>RF"]==0))/length(which(is.na(tabla[,"SW>RF"])==FALSE))  
}#c
for (r in 1:ncol(results_tabla)){results_tabla[8,r]<-mean(results_tabla[1:7,r])}
row.names(results_tabla)<-c(paste("chr ",1:7,"H", sep = ""), "Genome-wide")
#saveRDS(MSE_list, paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.1.1_SW_vs_RF_tables_",j,".RDS", sep = ""))
write.csv(results_tabla, paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.1.2_SW_vs_RF_summary_",j,".csv", sep = ""))

################

#full data models
models<-list()
for (c in 1:7){ cat(c);cat(": ")   
  windows<-row.names(rec_list[[c]])
  models[[c]]<-list()
  for (w in windows){ #cat(which(windows==w)/length(windows)); cat("-")   
    
    models[[c]][[w]]<-list()
    
    full_set<-standarized_variables[[c]][[w]][["full_set"]]
    
    full.model<-lm(formula = rec.rate ~., data = full_set)
    if (full.model$coefficients[[1]]!=0){
    #STEPWISE regression
    step.model<-stepAIC(full.model, direction = "both", trace = FALSE)
    models[[c]][[w]][["SW"]]<-step.model
    
    #RANDOM FOREST regression
    models[[c]][[w]][["RF"]] <- randomForest(rec.rate~., data=full_set, proximity=TRUE, importance = TRUE,  mtry = m, maxnodes = max_terminal_nodes)
    }#if
  }#w
}#c
#saveRDS(models, paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.3_full_data_models_",j,".RDS", sep = ""))


#I will extract only significant variables and coefficients
alpha<-0.05
exp_variables<-names(standarized_variables[[1]][[1]][[1]])[-1]

significant_variables_list<-list()
for (c in 1:7){ cat(c);cat("-")
  windows<-row.names(rec_list[[c]])
  coeff_table<-matrix(nrow = length(windows), ncol = length(exp_variables))
  row.names(coeff_table)<-windows; colnames(coeff_table)<-exp_variables
  p_value_table<-coeff_table
  significant_variables_list[[c]]<-list()
  chr_table<-matrix(ncol = 4, nrow = 0)
  for (w in windows){ #w  
    if (w %in% names(models[[c]])){
      if (is.null(models[[c]][[w]]$SW)==FALSE){
      significant_coefficients<-as.data.frame(summary(models[[c]][[w]]$SW)$coefficients)[-1,]
      if (nrow(significant_coefficients)!=0){
      if (any(significant_coefficients[,4]<alpha)){
        significant_coefficients<-significant_coefficients[which(significant_coefficients[,4]<alpha),]
        variables<-row.names(significant_coefficients)
        #save
        for (i in variables){coeff_table[w, i]<-significant_coefficients[i,1]}
        for (i in variables){p_value_table[w, i]<-significant_coefficients[i,4]}
        win_table<-cbind(w, variables, significant_coefficients[,c(1,4)])
        chr_table<-rbind(chr_table, win_table)
      }#if null
      }#if any  
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
#saveRDS(significant_variables_list, paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.4_SW_significant_variables_",j,".RDS", sep = ""))


#significant_variables_list<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/I.1.5_REDUCED_sw_reg_per_win_significant_variables.RDS")

results_tabla<-matrix(ncol = 2, nrow = 8); colnames(results_tabla)<-c("All windows","MSEsw>MSErf")
for (c in 1:7){ cat(c);cat(": ")
  SW_table<-significant_variables_list[[c]]$graphic_table
  windows<-unique(as.character(SW_table[,2]))
  chr_index<-c()
  for (w in windows){  
  if (is.null(models[[c]][[w]]$SW)==FALSE){
  SW_sig_var<-SW_table[which(SW_table[,2]==w),]
  significant_variables<-as.character(SW_sig_var$variable)
  n_var<-nrow(SW_sig_var)  
  #get RF importance ranking
  importance<-abs(models[[c]][[w]][["RF"]]$importance[,1])
  important_variables<-names(importance[order(importance, decreasing = TRUE)])[1:n_var]
  index<-length(which(significant_variables%in%important_variables))/length(significant_variables)
  chr_index<-c(chr_index,index)
  }#if null
  }#w
  results_tabla[c,1]<-mean(chr_index)
  
  #for when "MSEsw>MSErf"
  windows<-row.names(MSE_list[[c]][which(MSE_list[[c]][,3]==1),])
  chr_index<-c()
  for (w in windows){  
    if (is.null(models[[c]][[w]]$SW)==FALSE){
      SW_sig_var<-SW_table[which(SW_table[,2]==w),]
      significant_variables<-as.character(SW_sig_var$variable)
      n_var<-nrow(SW_sig_var)  
      #get RF importance ranking
      importance<-abs(models[[c]][[w]][["RF"]]$importance[,1])
      important_variables<-names(importance[order(importance, decreasing = TRUE)])[1:n_var]
      index<-length(which(significant_variables%in%important_variables))/length(significant_variables)
      chr_index<-c(chr_index,index)
    }#if null
  }#w
  results_tabla[c,2]<-mean(chr_index, na.rm = TRUE)
}#c
for (r in 1:ncol(results_tabla)){results_tabla[8,r]<-mean(results_tabla[1:7,r])}
row.names(results_tabla)<-c(paste("chr ",1:7,"H", sep = ""), "Genome-wide")
write.csv(results_tabla, paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.5_SW_vs_RF_selected_variables_summary_",j,".csv", sep = ""))

}#n


m<-4 #SIMPLE: 1, 2, 3 ,4  #FULL: 1, 3, 4=default, 6, 9


  
#get results
mtyrs<-c(1:4)
MAX_NODES<-c(4,8,12,16)
tabla<-matrix(ncol = 4, nrow = 4)
row.names(tabla)<-paste("n=", mtyrs, sep = ""); colnames(tabla)<-paste("n=", MAX_NODES, sep = "")
for (m in mtyrs){
for (n in MAX_NODES){ 
j<-paste("mtry=",m,"_maxnodes=",n, sep = "")  
results_table<-read.csv(paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.1.2_SW_vs_RF_summary_",j,".csv", sep = ""))
tabla[paste("n=", m, sep = ""),paste("n=", n, sep = "")]<-round(results_table[8,4], digits = 2)    
}#n
}#m
write.csv(tabla, "/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.6_SW_vs_RF_summary.csv")

tabla<-matrix(ncol = 4, nrow = 4)
row.names(tabla)<-paste("n=", mtyrs, sep = ""); colnames(tabla)<-paste("n=", MAX_NODES, sep = "")
for (m in mtyrs){
  for (n in MAX_NODES){ 
    j<-paste("mtry=",m,"_maxnodes=",n, sep = "")  
    results_table<-read.csv(paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.5_SW_vs_RF_selected_variables_summary_",j,".csv", sep = ""))
    tabla[paste("n=", m, sep = ""),paste("n=", n, sep = "")]<-round(results_table[8,2], digits = 2)    
  }#n
}#m
write.csv(tabla, "/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.7_SW_vs_RF_variable_chosen.csv")

tabla<-matrix(ncol = 4, nrow = 4)
row.names(tabla)<-paste("n=", mtyrs, sep = ""); colnames(tabla)<-paste("n=", MAX_NODES, sep = "")
for (m in mtyrs){
  for (n in MAX_NODES){ 
    j<-paste("mtry=",m,"_maxnodes=",n, sep = "")  
    results_table<-read.csv(paste("/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.5_SW_vs_RF_selected_variables_summary_",j,".csv", sep = ""))
    tabla[paste("n=", m, sep = ""),paste("n=", n, sep = "")]<-round(results_table[8,3], digits = 2)    
  }#n
}#m
write.csv(tabla, "/scratch/Federico/3_RILs/5_correlation/Results/K_random_forest/K.8_SW_vs_RF_variable_chosen_when_RF<SW.csv")




# confusionMatrix(p2, test$ Species)
# 
# plot(rf)

#y	-> If a factor, classification is assumed, otherwise regression is assumed. If omitted -> unsupervised mode.
#ntree: number of trees in the forest
#mtry -> Number of variables randomly sampled as candidates at each split.
#maxnodes	<- Maximum number of terminal nodes trees in the forest can have.
#importance	Should importance of predictors be assessed?