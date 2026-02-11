library(MASS)
library(randomForest)
library(party)

#0-get data
s<-1000000

rec_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))
SVs_prop_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
GEN_DIST_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
DMRs_UNI_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/4_methylation/Results/B_methylation_per_population/viejo/B.4.3_DMRs_unified_per_pop_wi_1e+06_windows.RDS", sep = ""))
ADD_gen_eff<-readRDS("/home/fcasale/Desktop/Paper_2/2_B_populations_GP/Results/C.2.1_additive_effects_by_window=1Mbp.RDS")
DOM_gen_eff<-readRDS("/home/fcasale/Desktop/Paper_2/2_B_populations_GP/Results/C.2.2_dominant_effects_by_window=1Mbp.RDS")

Populations<-colnames(rec_list[[1]])

#1-accumulate data in data frame (across genome, all windows together)
variables<-matrix(ncol = 5, nrow = 0)
for (c in 1:7){ cat(c); cat("-")
  for (w in 1:nrow(rec_list[[c]])){ 
  variables<-rbind(variables, cbind(rec_list[[c]][w,], SVs_prop_list[[c]][w,], GEN_DIST_list[[c]][w,], DMRs_UNI_list[[c]][w,], ADD_gen_eff[[c]][w,]))
  }#w
}#c
colnames(variables)<-c("rec.rate", "SVs.prop", "gen.dist", "methy.level", "gen.effect")

#2- standarize explanatory variables across genome

variables<-as.data.frame(variables)
for (c in 1:ncol(variables)){variables[[c]]<-as.numeric(as.character(variables[[c]]))}

standarize<-function(x){x<-(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE); return(x)}
#normalize<-function(x){x<-(x-mean(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)); return(x)}
variables[,2:ncol(variables)]<-apply(variables[,2:ncol(variables)], 2, FUN = standarize) 
saveRDS(variables, "/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/j.1.1_standarized_variables_across_genome.RDS")

#3- linear model

full.model<-lm(formula = rec.rate ~ gen.dist + gen.effect + methy.level + SVs.prop , data = variables)

predicted<-predict(object = full.model, variables[,2:ncol(variables)])
cor.test(y = variables$rec.rate, x = predicted)

#4- RF

MAX_NODES<-c(4,8,12,16)

max_terminal_nodes<-n #4, 8, 12, 16 -> 16 es lo mismo que dejarlo al maximo

randomForest(rec.rate~., data=variables, proximity=TRUE, importance = TRUE, mtry = m, maxnodes = max_terminal_nodes)

#5 step-wise

step.model<-stepAIC(full.model, direction = "both", trace = FALSE)
predicted<-predict(object = full.model, variables[,2:ncol(variables)])
cor.test(y = variables$rec.rate, x = predicted)

