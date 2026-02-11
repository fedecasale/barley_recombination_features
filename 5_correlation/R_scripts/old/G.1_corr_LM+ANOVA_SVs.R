make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

s<-1000000

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))

corr_table<-matrix(nrow = nrow(rec_list[[2]]), ncol = 7)
colnames(corr_table)<-paste(1:7, "H", sep = "")
row.names(corr_table)<-row.names(rec_list[[2]])
corr_table[]<-"-"

Populations<-colnames(rec_list[[1]])

SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
TYPES_SVs_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))
GENE_prop_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/F.0_gene_proportion_per_window=",s,"_all_pops.RDS", sep = ""))
GEN_DIST_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))

#plot(GEN_DIST_list[[3]][,1]~row.names(GEN_DIST_list[[3]]))

corr_list<-list()


#1-correlation rec.rate by prop of sv genotype (ALL SVs together)
tabla<-corr_table
for (c in 1:7){ 
  windows<-row.names(rec_list[[c]])
  for (w in windows){
      rec.rate<-as.numeric(rec_list[[c]][w,])
      SVs_prop<-as.numeric(SVs_prop_list[[c]][w,])
      LM<-lm(rec.rate~SVs_prop)
      ANOVA<-anova(LM)
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[1], digits = 3)
  }#w
}#c
corr_list[["1_rr~tot.SVs"]]<-tabla


#2-correlation rec.rate by prop of sv genotype (ALL SVs together), but including genetic distance among parents
tabla<-corr_table
for (c in 1:7){ 
  windows<-row.names(rec_list[[c]])
  for (w in windows){
    rec.rate<-as.numeric(rec_list[[c]][w,])
    SVs_prop<-as.numeric(SVs_prop_list[[c]][w,])
    GEN_DIST<-as.numeric(GEN_DIST_list[[c]][w,])
    LM<-lm(rec.rate~GEN_DIST+SVs_prop)
    ANOVA<-anova(LM)
    tabla[w,c]<-round(ANOVA$`Pr(>F)`[1], digits = 3)
  }#w
}#c
corr_list[["2_rr~gen.dist+tot.SVs"]]<-tabla


#3-correlation rec.rate by prop of sv TYPE genotype
corr_list[["3_SVs"]]<-list()
for (v in names(TYPES_SVs_prop_list)){
tabla<-corr_table
for (c in 1:7){ 
  windows<-row.names(rec_list[[c]])
  for (w in windows){
    rec.rate<-as.numeric(rec_list[[c]][w,])
    SV<-as.numeric(TYPES_SVs_prop_list[[v]][[c]][w,])
    LM<-lm(rec.rate~SV)
    ANOVA<-anova(LM)
    the_last_variable<-nrow(ANOVA)-1
    if (the_last_variable!=0){tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)}
  }#w
}#c
corr_list[["3_SVs"]][[paste("rr~",v, sep = "")]]<-tabla
}#v


#4-correlation rec.rate by prop of sv TYPE genotype, using the SV as last explanatory variable
corr_list[["4_SVs"]]<-list()
for (v in names(TYPES_SVs_prop_list)){
  tabla<-corr_table
  for (c in 1:7){ 
    windows<-row.names(rec_list[[c]])
    for (w in windows){
      rec.rate<-as.numeric(rec_list[[c]][w,])
      j<-which(names(TYPES_SVs_prop_list)==v)
      other_SVs<-names(TYPES_SVs_prop_list)[-j]
      SVs<-list()
      for (i in other_SVs){SVs[[i]]<-as.numeric(TYPES_SVs_prop_list[[i]][[c]][w,])}#i
      SVs[[4]]<-as.numeric(TYPES_SVs_prop_list[[v]][[c]][w,])
      LM<-lm(rec.rate~SVs[[1]]+SVs[[2]]+SVs[[3]]+SVs[[4]])
      ANOVA<-anova(LM)
      the_last_variable<-nrow(ANOVA)-1
      if (the_last_variable!=0){tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)}    
      }#w
  }#c
  corr_list[["4_SVs"]][[paste("rr~",other_SVs[1],"+",other_SVs[2],"+",other_SVs[3],"+",v, sep = "")]]<-tabla
}#v

#5-correlation rec.rate by prop of sv TYPE genotype, using the SV as last explanatory variable, 
#but including genetic distance among parents
corr_list[["5_gen.dist+SVs"]]<-list()
for (v in names(TYPES_SVs_prop_list)){
  tabla<-corr_table
  for (c in 1:7){ 
    windows<-row.names(rec_list[[c]])
    for (w in windows){
      rec.rate<-as.numeric(rec_list[[c]][w,])
      GEN_DIST<-as.numeric(GEN_DIST_list[[c]][w,])
      j<-which(names(TYPES_SVs_prop_list)==v)
      SVs<-list()
      for (i in names(TYPES_SVs_prop_list)[-j]){SVs[[i]]<-as.numeric(TYPES_SVs_prop_list[[i]][[c]][w,])}#i
      SVs[[4]]<-as.numeric(TYPES_SVs_prop_list[[v]][[c]][w,])
      LM<-lm(rec.rate~GEN_DIST+SVs[[1]]+SVs[[2]]+SVs[[3]]+SVs[[4]])
      ANOVA<-anova(LM)
      the_last_variable<-nrow(ANOVA)-1
      if (the_last_variable!=0){tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)}    
    }#w
  }#c
  corr_list[["5_gen.dist+SVs"]][[paste("rr~gen.dist+",other_SVs[1],"+",other_SVs[2],"+",other_SVs[3],"+",v, sep = "")]]<-tabla
}#v
#########################

saveRDS(corr_list, "/scratch/Federico/3_RILs/5_correlation/Results/G.1_corr_rr_SVs.RDS")
