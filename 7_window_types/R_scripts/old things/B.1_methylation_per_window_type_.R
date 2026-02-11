library(data.table)
library(xlsx)
library(agricolae)
library(multcompView)

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

dir.create("/scratch/Federico/3_RILs/7_window_types/Results/B_methylation")
#################################################################################

#get windows types and unify
distal_windows<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.0_accum_rec_prob_DISTAL.RDS")
peri_windows<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.0_accum_rec_prob_PERI.RDS")
coldspots<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.1_coldspots_per_pop_WINDOWS.RDS")
hotspots<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.2_hotspots_per_pop_WINDOWS.RDS")

#make list
window_type_list<-list(peri_windows, distal_windows, coldspots, hotspots)
window_types<-c("peri_windows", "distal_windows", "coldspots", "hotspots")
names(window_type_list)<-window_types
saveRDS(window_type_list, "/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.1_window_types.RDS")

#add methylation data #I will work with c sites instead of DMRs because there are windows without DMRs.
methy_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.1_methylation_per_population_win=",win.size,".RDS", sep = ""))
met_contexts<-names(methy_list[[1]])

methy_list_windows<-list()
for (v in met_contexts){ cat(v); cat(": ")
  methy_list_windows[[v]]<-list()
  for (e in window_types){  cat(e); cat("-")
    methy_list_windows[[v]][[e]]<-list()
    for (p in Populations){
      methy_list_windows[[v]][[e]][[p]]<-list()
      for (c in 1:7){    
        windows<-as.character(window_type_list[[e]][[p]][[c]][,2])
        methylation<-methy_list[[c]][[v]][windows,p]
        methy_list_windows[[v]][[e]][[p]][[c]]<-window_type_list[[e]][[p]][[c]]
        row.names(methy_list_windows[[v]][[e]][[p]][[c]])<-methy_list_windows[[v]][[e]][[p]][[c]][,2]
        methy_list_windows[[v]][[e]][[p]][[c]][windows,3]<-methylation
        colnames(methy_list_windows[[v]][[e]][[p]][[c]])[3]<-"methylation_level"
        # #I will take out 5% of the extremes beacuse I think they are creating wrong comparisons
        # tabla<-methy_list_windows[[v]][[e]][[p]][[c]]
        # #limits<-quantile(tabla[,3], probs = seq(0, 1, by = 0.05), na.rm = TRUE)[c(2,20)]
        # limits<-quantile(tabla[,3], na.rm = TRUE)[c(2,4)]
        # tabla<-tabla[-unique(c(which(tabla[,3]<limits[1]), which(tabla[,3]>limits[2]))),]
        # methy_list_windows[[v]][[e]][[p]][[c]]<-tabla
      }#c
    }#p  
  }#e
  cat(fill = TRUE)
}#v
saveRDS(methy_list_windows, "/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.2_methylation_prop_per_window_type.RDS")

# #delete coldspots and hotspots from the regions windows #I WONT DO THIS
# for (m in met_contexts){cat(m); cat(":")
#   for (e in c("peri_windows", "distal_windows")){ cat(e); cat("-")
#     for (p in Populations){
#       for (c in 1:7){
#         coldspots<-which(methy_list_windows[[m]][[e]][[p]][[c]][,1]%in%methy_list_windows[[m]][["coldspots"]][[p]][[c]][,1])
#         hotspots<-which(methy_list_windows[[m]][[e]][[p]][[c]][,1]%in%methy_list_windows[[m]][["hotspots"]][[p]][[c]][,1])
#         to.delete<-c(coldspots, hotspots)
#         methy_list_windows[[m]][[e]][[p]][[c]]<-methy_list_windows[[m]][[e]][[p]][[c]][-to.delete,]
#       }#c
#     }#p
#   }#e
# }#m

#################### check if differences among chrs within populations ###########################

file.remove("/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.3_all_methy_types_windows_DIF_AMONG_CHR_TUKEY.xlsx")

n_tests<-length(window_types)*(length(window_types)-1)/2

tabla<-matrix(ncol = 7+1, nrow = 0)
colnames(tabla)<-c("", paste(1:7, "H", sep = ""))
for (p in Populations){ cat(p); cat(": ")
  tabla_p<-matrix(ncol = 7, nrow = 0)
  for (v in met_contexts){ cat(v); cat("-")
    tabla_v<-matrix(ncol = 7, nrow = length(window_types)+1); row.names(tabla_v)<-c(v, window_types)
    for (e in window_types){
      data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "window_type")
      for (c in 1:7){
        data_table_c<-cbind(methy_list_windows[[v]][[e]][[p]][[c]][,3], rep(c, length(methy_list_windows[[v]][[e]][[p]][[c]][,3])))
        colnames(data_table_c)<-c("methylation", "window_type")
        data_table<-rbind(data_table, data_table_c)
      }#c
      data_table<-as.data.frame(data_table)
      data_table$methylation<-as.numeric(as.character(data_table$methylation))
      LM<-lm(methylation~window_type, data = as.data.frame(data_table))
      #ANOVA + Tukey
      tukey<-HSD.test(y = LM, trt = "window_type", alpha = 0.05/n_tests, unbalanced = TRUE, group = TRUE)$groups
      #PWT<-pairwise.wilcox.test(data_table$methylation, data_table$window_type, p.adjust.method = "bonferroni", paired = FALSE)
      for (c in 1:7){tabla_v[e,c]<-paste(round(tukey[paste(c),1], digits = 3), tukey[paste(c),2])}#e
    }#e
    tabla_p<-rbind(tabla_p, tabla_v)
  }#v
  tabla_p<-cbind(row.names(tabla_p), tabla_p)
  tabla_p<-rbind(c(p, rep("",7)), tabla_p)
  tabla<-rbind(tabla, tabla_p)
  cat(fill = TRUE)
}#p
write.xlsx(tabla, "/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.3_all_methy_types_windows_DIF_AMONG_CHR_TUKEY.xlsx", row.names = FALSE)
###################################################################################################

#unify chrs
for (v in met_contexts){ cat(v); cat(": ")
  for (e in names(methy_list_windows[[1]])){ cat(e); cat(": ")
    for (p in Populations){ cat(p); cat("-")
      pop_data<-c()
      for (c in 1:7){pop_data<-c(pop_data, methy_list_windows[[v]][[e]][[p]][[c]][,3])}#c
      methy_list_windows[[v]][[e]][[p]]<-pop_data
    }#p
  }#e
  cat(fill = TRUE)
}#v

saveRDS(methy_list_windows, "/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.4_methylation_prop_per_window_type_unified.RDS")
#methy_list_windows<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.4_methylation_prop_per_window_type_unified.RDS")

###################### check if means are different with ANOVA + Tukey #############################
file.remove("/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.5.1_all_methy_types_windows_unified_TUKEY.xlsx")

tabla<-matrix(ncol = length(window_types)+1, nrow = 0)
colnames(tabla)<-c("", window_types)
for (p in Populations){ cat(p); cat(": ")
  tabla_p<-matrix(ncol = length(window_types), nrow = 3) 
  colnames(tabla_p)<-window_types
  row.names(tabla_p)<-met_contexts
  for (v in met_contexts){ cat(v); cat("-")
    #gather data to make linear regression
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "window_type")
    for (e in window_types){
      data_table_e<-cbind(methy_list_windows[[v]][[e]][[p]], rep(e, length(methy_list_windows[[v]][[e]][[p]])))
      colnames(data_table_e)<-c("methylation", "window_type")
      data_table<-rbind(data_table, data_table_e)
    }#e
    data_table<-as.data.frame(data_table)
    data_table$methylation<-as.numeric(as.character(data_table$methylation))
    LM<-lm(methylation~window_type, data = as.data.frame(data_table))
    # #ANOVA + Tukey
    tukey<-HSD.test(LM, "window_type", alpha = 0.05/n_tests)$groups
    # #fill results table
    for (e in window_types){tabla_p[v,e]<-paste(round(tukey[e,1], digits = 3), tukey[e,2])}#e
  }#v
  tabla_p<-cbind(c(p, met_contexts), rbind("",tabla_p))
  tabla<-rbind(tabla, tabla_p)
  cat(fill = TRUE)
}#p  

write.xlsx(tabla, "/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.5.1_all_methy_types_windows_unified_TUKEY.xlsx")
####################################################################################################


########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.5.2.1_all_methy_types_windows_unified_WILCOX.xlsx")
file.remove("/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.5.2.2_all_methy_types_windows_unified_WILCOX_GROUPING.xlsx")

tabla_group<-matrix(ncol = length(window_types)+1, nrow = 0); colnames(tabla_group)<-c("",window_types)
for (p in Populations){ cat(p); cat("-") 
  #table for wilcox results
  tabla_p<-matrix(ncol = length(window_types), nrow = 0)
  #table for wilcox grouping
  tabla_g<-matrix(ncol = length(window_types), nrow = 3); colnames(tabla_g)<-window_types; row.names(tabla_g)<-met_contexts
  for (v in met_contexts){ cat(v); cat("-")
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "window_type")
    for (e in window_types){
      data_table_e<-cbind(methy_list_windows[[v]][[e]][[p]], rep(e, length(methy_list_windows[[v]][[e]][[p]])))
      colnames(data_table_e)<-c("methylation", "window_type")
      data_table<-rbind(data_table, data_table_e)
    }#e
    data_table<-as.data.frame(data_table)
    data_table$methylation<-as.numeric(as.character(data_table$methylation))
    #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
    PWT<-pairwise.wilcox.test(data_table$methylation, data_table$window_type)
    p.values_table<-format(round(PWT$p.value, digits = 3), scientific = FALSE)
    pair
    #group by letters
    p_val_list<-list()
    for (c in window_types){ 
      p.val<-c()
      if (any(row.names(p.values_table)==c)){p.val<-c(p.val,p.values_table[c,])}
      if (any(colnames(p.values_table)==c)){p.val<-c(p.val,p.values_table[,c])}
      if (any(p.val=="   NA")){p.val<-p.val[-which(p.val=="   NA")]}
      if (any(p.val=="NA")){p.val<-p.val[-which(p.val=="NA")]}
      p_val_list[[c]]<-list()
      for (r in window_types[-which(window_types==c)]){p_val_list[[c]][[r]]<-as.numeric(p.val[r])}
    }#c  
    p_val_table<-matrix(ncol = length(window_types), nrow = length(window_types))
    colnames(p_val_table)<-window_types; row.names(p_val_table)<-window_types
    for (c in colnames(p_val_table)){ 
      p_val_table[c,c]<-1
      for (r in names(p_val_list[[c]])){p_val_table[r,c]<-as.numeric(p_val_list[[c]][[r]])}
    }
    letter_table<-multcompLetters(p_val_table, compare = "<", threshold = 0.05/n_tests)
    for (c in window_types){tabla_g[v,c]<-gsub(" ", "", letter_table[[2]][c])}
    #add mean
    window_type_means<-window_types; names(window_type_means)<-window_types
    for (c in window_types){window_type_means[c]<-round(mean(data_table$methylation[which(data_table$window_type==c)], na.rm = TRUE), digits = 3)}
    for (c in window_types){tabla_g[v,c]<-paste(window_type_means[c], tabla_g[v,c])}
    
    #arrange table results
    for (c in colnames(p.values_table)){colnames(p.values_table)[which(colnames(p.values_table)==c)]<-paste(c, "=", window_type_means[c], sep = "")}#c 
    for (c in row.names(p.values_table)){row.names(p.values_table)[which(row.names(p.values_table)==c)]<-paste(c, "=", window_type_means[c], sep = "")}#c
    p.values_table<-cbind(c(v,row.names(p.values_table)), rbind(colnames(p.values_table), p.values_table))
    tabla_p<-rbind(tabla_p, "", p.values_table)
  }#v
  colnames(tabla_p)<-NULL
  write.xlsx(tabla_p, "/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.5.2.1_all_methy_types_windows_unified_WILCOX.xlsx", sheetName = p, append = TRUE)
  #arrange table group
  tabla_g<-cbind(c(p, met_contexts), rbind("",tabla_g))
  tabla_group<-rbind(tabla_group, tabla_g)
}#p
write.xlsx(tabla_group, "/scratch/Federico/3_RILs/7_window_types/Results/B_methylation/B.1.5.2.2_all_methy_types_windows_unified_WILCOX_GROUPING.xlsx", row.names = FALSE)

#####################################################################################################

