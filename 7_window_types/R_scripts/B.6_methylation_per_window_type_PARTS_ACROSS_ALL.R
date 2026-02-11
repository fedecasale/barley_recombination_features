library(data.table)
library(xlsx)
library(agricolae)
library(multcompView)

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

methy_list_windows<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.2.4_methylation_prop_per_window_type_unified.RDS")
methy_list_windows$AVE<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.3_methylation_prop_per_window_type_unified_UNI.RDS")

met_contexts<-names(methy_list_windows)
window_types<-names(methy_list_windows[[1]])

window_types_pops<-c()
for (e in window_types){ cat(e); cat(": ")
  for (p in Populations){ 
  window_types_pops<-c(window_types_pops, paste(e, p, sep = "_"))
  }
}

n_tests<-length(window_types_pops)*(length(window_types_pops)-1)/2

###################### check if means are different with ANOVA + Tukey #############################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.2.5.1_all_methy_types_windows_unified_TUKEY.xlsx")

tabla<-matrix(ncol = length(window_types_pops), nrow = length(met_contexts))
row.names(tabla)<-met_contexts
colnames(tabla)<-window_types_pops
  for (v in met_contexts){ cat(v); cat("-")
    #gather data to make linear regression
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "window_type_pop")
    for (e in window_types){ cat(e); cat(": ")
      for (p in Populations){ 
      data_table_e<-cbind(methy_list_windows[[v]][[e]][[p]], rep(paste(e,p,sep = "_"), length(methy_list_windows[[v]][[e]][[p]])))
      colnames(data_table_e)<-c("methylation", "window_type_pop")
      data_table<-rbind(data_table, data_table_e)
      }#p
    }#e
    data_table<-as.data.frame(data_table)
    data_table$methylation<-as.numeric(as.character(data_table$methylation))
    LM<-lm(methylation~window_type_pop, data = as.data.frame(data_table))
    # #ANOVA + Tukey
    tukey<-HSD.test(LM, "window_type_pop", alpha = 0.05/n_tests)$groups
    #fill results table
    for (e in window_types_pops){tabla[v,e]<-paste(round(tukey[e,1], digits = 3), tukey[e,2])}#e
    cat(fill = TRUE)
  }#v

write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.6.1_all_methy_types_windows_unified_TUKEY.xlsx")
####################################################################################################

########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.6.2_all_methy_types_windows_unified_WILCOX_GROUPING.xlsx")

tabla<-matrix(ncol = length(window_types_pops), nrow = length(met_contexts))
row.names(tabla)<-met_contexts
colnames(tabla)<-window_types_pops
for (v in met_contexts){ cat(v); cat("-")
  #gather data to make linear regression
  data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "window_type_pop")
  for (e in window_types){ cat(e); cat(": ")
    for (p in Populations){ 
      data_table_e<-cbind(methy_list_windows[[v]][[e]][[p]], rep(paste(e,p,sep = "_"), length(methy_list_windows[[v]][[e]][[p]])))
      colnames(data_table_e)<-c("methylation", "window_type_pop")
      data_table<-rbind(data_table, data_table_e)
    }#p
  }#e
  data_table<-as.data.frame(data_table)
  data_table$methylation<-as.numeric(as.character(data_table$methylation))  
  PWT<-pairwise.wilcox.test(data_table$methylation, data_table$window_type_pop)
  p.values_table<-format(round(PWT$p.value, digits = 3), scientific = FALSE)

  #group by letters
  p_val_list<-list()
  for (c in window_types_pops){ 
    p.val<-c()
    if (any(row.names(p.values_table)==c)){p.val<-c(p.val,p.values_table[c,])}
    if (any(colnames(p.values_table)==c)){p.val<-c(p.val,p.values_table[,c])}
    p.val<-strsplit(p.val, " ")
    p.val<-unlist(lapply(p.val, function(x){if(any(x=="")){x<-x[-which(x=="")]};return(x)}))
    if (any(p.val=="NA")){p.val<-p.val[-which(p.val=="NA")]}
    p_val_list[[c]]<-list()
    for (r in window_types_pops[-which(window_types_pops==c)]){p_val_list[[c]][[r]]<-as.numeric(p.val[r])}
  }#c  
  
  p_val_table<-matrix(ncol = length(window_types_pops), nrow = length(window_types_pops))
  colnames(p_val_table)<-window_types_pops; row.names(p_val_table)<-window_types_pops
  for (c in colnames(p_val_table)){ 
    p_val_table[c,c]<-1
    for (r in names(p_val_list[[c]])){p_val_table[r,c]<-as.numeric(p_val_list[[c]][[r]])}
  }
  letter_table<-multcompLetters(p_val_table, compare = "<", threshold = 0.05/n_tests)
  for (c in window_types_pops){tabla[v,c]<-gsub(" ", "", letter_table[[2]][c])}
  #add mean
  window_type_means<-window_types_pops; names(window_type_means)<-window_types_pops
  for (c in window_types_pops){window_type_means[c]<-round(mean(data_table$methylation[which(data_table$window_type_pop==c)], na.rm = TRUE), digits = 3)}
  for (c in window_types_pops){tabla[v,c]<-paste(window_type_means[c], tabla[v,c])}
  
  }#v

#save
write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.6.2_all_methy_types_windows_unified_WILCOX_GROUPING.xlsx", row.names = FALSE)

#####################################################################################################

