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

corr_list<-list()

#1-correlation rec.rate by methylation level of each met context per window 
corr_list[["1_met"]]<-list()
for (v in met_context){ cat(v); cat(": ")
tabla<-corr_table
for (c in 1:7){ cat(c); cat("-")
  windows<-row.names(rec_list[[c]])
  for (w in windows){   
    rec.rate<-as.numeric(rec_list[[c]][w,])
    METHY_level<-as.numeric(METHY_list[[c]][[v]][w,])
    LM<-lm(rec.rate~METHY_level)
    ANOVA<-anova(LM)
    tabla[w,c]<-round(ANOVA$`Pr(>F)`[1], digits = 3)
  }#w
}#c
corr_list[["1_met"]][[paste("rr~",v, sep = "")]]<-tabla
}#v


#2-correlation rec.rate by methylation level of each met context per window , but taking into account the other contexts in the model
corr_list[["2_met"]]<-list()
for (v in met_context){ cat(v); cat(": ")
  tabla<-corr_table
  other_v<-met_context[-which(met_context==v)]
  for (c in 1:7){ cat(c); cat("-")
    windows<-row.names(rec_list[[c]])
    for (w in windows){
      rec.rate<-as.numeric(rec_list[[c]][w,])
      METHY_level<-as.numeric(METHY_list[[c]][[v]][w,])
      other_v1_level<-as.numeric(METHY_list[[c]][[other_v[1]]][w,])
      other_v2_level<-as.numeric(METHY_list[[c]][[other_v[2]]][w,])
      LM<-lm(rec.rate~other_v1_level+other_v2_level+METHY_level)
      ANOVA<-anova(LM)
      the_last_variable<-nrow(ANOVA)-1
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)
    }#w
  }#c
  corr_list[["2_met"]][[paste("rr~",other_v[1],"+",other_v[2],"+",v,sep = "")]]<-tabla
}#v


#3-correlation rec.rate by methylation level of each met context per window, but including the parental dif in methylation for that context
corr_list[["3_pdif+met"]]<-list()
for (v in met_context){ cat(v); cat(": ")
  tabla<-corr_table
  for (c in 1:7){ cat(c); cat("-")
    windows<-row.names(rec_list[[c]])
    for (w in windows){   
      rec.rate<-as.numeric(rec_list[[c]][w,])
      METHY_level<-as.numeric(METHY_list[[c]][[v]][w,])
      pdif<-as.numeric(pdif_METHY_list[[c]][[v]][w,])
      LM<-lm(rec.rate~pdif+METHY_level)
      ANOVA<-anova(LM)
      the_last_variable<-nrow(ANOVA)-1
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)
    }#w
  }#c
  corr_list[["3_pdif+met"]][[paste("rr~pdif+",v,sep = "")]]<-tabla
}#v



#4-correlation rec.rate by methylation level of each met context per window 
#and including the parental dif in methylation for that context but taking into account the other contexts in the model
corr_list[["4_pdif+met"]]<-list()
for (v in met_context){ cat(v); cat(": ")
  tabla<-corr_table
  other_v<-met_context[-which(met_context==v)]
  for (c in 1:7){ cat(c); cat("-")
    windows<-row.names(rec_list[[c]])
    for (w in windows){   
      rec.rate<-as.numeric(rec_list[[c]][w,])
      METHY_level<-as.numeric(METHY_list[[c]][[v]][w,])
      other_v1_level<-as.numeric(METHY_list[[c]][[other_v[1]]][w,])
      other_v2_level<-as.numeric(METHY_list[[c]][[other_v[2]]][w,])
      pdif<-as.numeric(pdif_METHY_list[[c]][[v]][w,])
      LM<-lm(rec.rate~other_v1_level+other_v2_level+pdif+METHY_level)
      ANOVA<-anova(LM)
      the_last_variable<-nrow(ANOVA)-1
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)
    }#w
  }#c
  corr_list[["4_pdif+met"]][[paste("rr~",other_v[1],"+",other_v[2],"+pdif+",v,sep = "")]]<-tabla
}#v

#5-correlation rec.rate by DMR average level of both parents of each met context per window 
corr_list[["5_DMRs"]]<-list()
for (v in met_context){ cat(v); cat(": ")
  tabla<-corr_table
  for (c in 1:7){ cat(c); cat("-")
    windows<-row.names(rec_list[[c]])
    for (w in windows){   
      rec.rate<-as.numeric(rec_list[[c]][w,])
      DMRs_level<-as.numeric(DMRs_list[[c]][[v]][w,])
      if (any(is.na(DMRs_level))==FALSE){
      LM<-lm(rec.rate~DMRs_level)
      ANOVA<-anova(LM)
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[1], digits = 3)
      } else {tabla[w,c]<-NA}
    }#w
  }#c
  corr_list[["5_DMRs"]][[paste("rr~",v,sep = "")]]<-tabla
}#v

#6-correlation rec.rate by DMR average level of both parents of each met context per window, 
#but taking into account the other contexts in the model
corr_list[["6_DMRs"]]<-list()
for (v in met_context){ cat(v); cat(": ")
  tabla<-corr_table
  other_v<-met_context[-which(met_context==v)]
  for (c in 1:7){ cat(c); cat("-")
    windows<-row.names(rec_list[[c]])
    for (w in windows){   
      rec.rate<-as.numeric(rec_list[[c]][w,])
      DMRs_level<-as.numeric(DMRs_list[[c]][[v]][w,])
      other_v1_level<-as.numeric(DMRs_list[[c]][[other_v[1]]][w,])
      other_v2_level<-as.numeric(DMRs_list[[c]][[other_v[2]]][w,])
      variable_vector<-c(rec.rate, DMRs_level, other_v1_level, other_v2_level)
      if (any(is.na(variable_vector))==FALSE){
      LM<-lm(rec.rate~other_v1_level+other_v2_level+DMRs_level)
      ANOVA<-anova(LM)
      the_last_variable<-nrow(ANOVA)-1
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)
      } else {tabla[w,c]<-NA}
    }#w
  }#c
  corr_list[["6_DMRs"]][[paste("rr~",other_v[1],"+",other_v[2],"+",v,sep = "")]]<-tabla
}#v



#7-correlation rec.rate by DMR level of each met context per window 
#and including the parental dif in methylation for that context but taking into account the other contexts in the model
corr_list[["7_pdif+DMRs"]]<-list()
for (v in met_context){ cat(v); cat(": ")
  tabla<-corr_table
  other_v<-met_context[-which(met_context==v)]
  for (c in 1:7){ cat(c); cat("-")
    windows<-row.names(rec_list[[c]])
    for (w in windows){   
      rec.rate<-as.numeric(rec_list[[c]][w,])
      DMRs_level<-as.numeric(DMRs_list[[c]][[v]][w,])
      other_v1_level<-as.numeric(DMRs_list[[c]][[other_v[1]]][w,])
      other_v2_level<-as.numeric(DMRs_list[[c]][[other_v[2]]][w,])
      pdif<-as.numeric(pdif_DMRs_list[[c]][[v]][w,])
      variable_vector<-c(rec.rate, DMRs_level, other_v1_level, other_v2_level)
      if (any(is.na(variable_vector))==FALSE){
      LM<-lm(rec.rate~other_v1_level+other_v2_level+pdif+DMRs_level)
      ANOVA<-anova(LM)
      the_last_variable<-nrow(ANOVA)-1
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)
      } else {tabla[w,c]<-NA}
    }#w
  }#c
  corr_list[["7_pdif+DMRs"]][[paste("rr~",other_v[1],"+",other_v[2],"+pdif+",v,sep = "")]]<-tabla
}#v


#8-correlation rec.rate by DMR average level of both parents unified (all met context together, weighted) per window 
tabla<-corr_table
for (c in 1:7){ cat(c); cat("-")
  windows<-row.names(rec_list[[c]])
  for (w in windows){   
    rec.rate<-as.numeric(rec_list[[c]][w,])
    DMRs_level<-as.numeric(DMRs_UNI_list[[c]][w,])
    if (any(is.na(DMRs_level))==FALSE){
      LM<-lm(rec.rate~DMRs_level)
      ANOVA<-anova(LM)
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[1], digits = 3)
    } else {tabla[w,c]<-NA}
  }#w
}#c
corr_list[["8_rr~DMRs(ave)"]]<-tabla


#6-correlation rec.rate by DMR average level of both parents of each met context per window, but including the parental dif in DMR for that context
tabla<-corr_table
for (c in 1:7){ cat(c); cat("-")
  windows<-row.names(rec_list[[c]])
  for (w in windows){   
    rec.rate<-as.numeric(rec_list[[c]][w,])
    DMRs_level<-as.numeric(DMRs_UNI_list[[c]][w,])
    pdif_DMRs_UNI<-as.numeric(pdif_DMRs_UNI_list[[c]][w,])
    if (any(is.na(DMRs_level))==FALSE){
      LM<-lm(rec.rate~pdif_DMRs_UNI+DMRs_level)
      ANOVA<-anova(LM)
      the_last_variable<-nrow(ANOVA)-1
      tabla[w,c]<-round(ANOVA$`Pr(>F)`[the_last_variable], digits = 3)
    } else {tabla[w,c]<-NA}
  }#w
}#c
corr_list[["9_rr~pdif+DMRs(ave)"]]<-tabla

#########################

saveRDS(corr_list, "/scratch/Federico/3_RILs/5_correlation/Results/G.2_corr_rr_methylation.RDS")
