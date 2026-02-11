library(data.table)
library(xlsx)
library(agricolae)
library(multcompView)

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

svs_list_windows<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.2.2_svs_prop_per_window_type_unified.RDS")
window_types<-names(svs_list_windows)

window_types_pops<-c()
for (e in window_types){ cat(e); cat(": ")
  for (p in Populations){ 
    window_types_pops<-c(window_types_pops, paste(e, p, sep = "_"))
  }
}

###################### check if means are different with ANOVA + Tukey #############################

n_tests<-length(window_types_pops)*(length(window_types_pops)-1)/2

file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.4.1_SVs_prop_windows_unified_TUKEY.xlsx")

tabla<-matrix(ncol = length(window_types_pops), nrow = 1)
colnames(tabla)<-window_types_pops
data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("svs_prop", "window_type_pop")
for (p in Populations){ cat(p); cat("-")
  for (e in window_types){
    data_table_e<-cbind(svs_list_windows[[e]][[p]], rep(paste(e,p,sep = "_"), length(svs_list_windows[[e]][[p]])))
    colnames(data_table_e)<-c("svs_prop", "window_type_pop")
    data_table<-rbind(data_table, data_table_e)
  }#e
}#p
  data_table<-as.data.frame(data_table)
  data_table$svs_prop<-as.numeric(as.character(data_table$svs_prop))
  LM<-lm(svs_prop~window_type_pop, data = data_table)
  #ANOVA + Tukey
  tukey<-HSD.test(LM, "window_type_pop", alpha = 0.05/n_tests)$groups
  #fill results table
  for (e in window_types_pops){tabla[,e]<-paste(round(tukey[e,1], digits = 3), tukey[e,2])}#e
write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.4.1_SVs_prop_windows_unified_TUKEY.xlsx")
####################################################################################################

########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.4.2_SVs_prop_windows_unified_WILCOX_GROUPING.xlsx")

#table for wilcox grouping
tabla<-matrix(ncol = length(window_types_pops), nrow = 1)
colnames(tabla)<-window_types_pops
data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("svs_prop", "window_type_pop")
for (p in Populations){ cat(p); cat("-") 
  for (e in window_types){
    data_table_e<-cbind(svs_list_windows[[e]][[p]], rep(paste(e,p,sep = "_"), length(svs_list_windows[[e]][[p]])))
    colnames(data_table_e)<-c("svs_prop", "window_type_pop")
    data_table<-rbind(data_table, data_table_e)
  }#e
}#p  
data_table<-as.data.frame(data_table)
data_table$svs_prop<-as.numeric(as.character(data_table$svs_prop))
#check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
PWT<-pairwise.wilcox.test(data_table$svs_prop, data_table$window_type_pop)
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
  for (c in window_types_pops){tabla[,c]<-gsub(" ", "", letter_table[[2]][c])}
  #add mean
  window_type_means<-window_types_pops; names(window_type_means)<-window_types_pops
  for (c in window_types_pops){window_type_means[c]<-round(mean(data_table$svs_prop[which(data_table$window_type_pop==c)], na.rm = TRUE), digits = 3)}
  for (c in window_types_pops){tabla[,c]<-paste(window_type_means[c], tabla[,c])}
  
#save
write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.4.2_SVs_prop_windows_unified_WILCOX_GROUPING.xlsx", row.names = FALSE)
