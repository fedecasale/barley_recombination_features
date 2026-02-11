library(data.table)
library(xlsx)
library(agricolae)
library(multcompView)

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

svs_list_windows<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.2.2_svs_prop_per_window_type_unified.RDS")
window_types<-names(svs_list_windows)

###################### check if means are different with ANOVA + Tukey #############################

n_tests<-length(Populations)*(length(Populations)-1)/2

file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.3.1_SVs_prop_windows_unified_TUKEY_across_pops.xlsx")

tabla<-matrix(ncol = length(Populations), nrow = length(window_types))
colnames(tabla)<-Populations
row.names(tabla)<-window_types
for (p in window_types){ cat(p); cat("-")
  #gather data to make linear regression
  data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("svs_prop", "population")
  for (e in Populations){ 
    data_table_e<-cbind(svs_list_windows[[p]][[e]], rep(e, length(svs_list_windows[[p]][[e]])))
    colnames(data_table_e)<-c("svs_prop", "population")
    data_table<-rbind(data_table, data_table_e)
  }#e
  data_table<-as.data.frame(data_table)
  data_table$svs_prop<-as.numeric(as.character(data_table$svs_prop))
  LM<-lm(svs_prop~population, data = as.data.frame(data_table))
  #ANOVA + Tukey
  tukey<-HSD.test(LM, "population", alpha = 0.05/n_tests)$groups
  #fill results table
  for (e in Populations){tabla[p,e]<-paste(round(tukey[e,1], digits = 3), tukey[e,2])}#e
}#p  

write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.3.1_SVs_prop_windows_unified_TUKEY_across_pops.xlsx")
####################################################################################################

########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.3.2_SVs_prop_windows_unified_WILCOX_across_pops.xlsx")
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.3.3_SVs_prop_windows_unified_WILCOX_GROUPING_across_pops.xlsx")

#table for wilcox grouping
tabla_group<-matrix(ncol = length(Populations), nrow = length(window_types))
colnames(tabla_group)<-Populations; row.names(tabla_group)<-window_types
for (p in window_types){ cat(p); cat("-") 
  data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("svs_prop", "population")
  for (e in Populations){
    data_table_e<-cbind(svs_list_windows[[p]][[e]], rep(e, length(svs_list_windows[[p]][[e]])))
    colnames(data_table_e)<-c("svs_prop", "population")
    data_table<-rbind(data_table, data_table_e)
  }#e
  data_table<-as.data.frame(data_table)
  data_table$svs_prop<-as.numeric(as.character(data_table$svs_prop))
  #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
  PWT<-pairwise.wilcox.test(data_table$svs_prop, data_table$population)
  p.values_table<-format(round(PWT$p.value, digits = 3), scientific = FALSE)
  
  #group by letters
  p_val_list<-list()
  for (c in Populations){ 
    p.val<-c()   
    if (any(row.names(p.values_table)==c)){p.val<-c(p.val,p.values_table[c,])}
    if (any(colnames(p.values_table)==c)){p.val<-c(p.val,p.values_table[,c])}
    p.val<-strsplit(p.val, " ")
    p.val<-unlist(lapply(p.val, function(x){if(any(x=="")){x<-x[-which(x=="")]};return(x)}))
    if (any(p.val=="NA")){p.val<-p.val[-which(p.val=="NA")]}
    p_val_list[[c]]<-list()
    for (r in Populations[-which(Populations==c)]){p_val_list[[c]][[r]]<-as.numeric(p.val[r])}
  }#c  
  
  p_val_table<-matrix(ncol = length(Populations), nrow = length(Populations))
  colnames(p_val_table)<-Populations; row.names(p_val_table)<-Populations
  for (c in colnames(p_val_table)){ 
    p_val_table[c,c]<-1
    for (r in names(p_val_list[[c]])){p_val_table[r,c]<-as.numeric(p_val_list[[c]][[r]])}
  }
  letter_table<-multcompLetters(p_val_table, compare = "<", threshold = 0.05/n_tests)
  for (c in Populations){tabla_group[p,c]<-gsub(" ", "", letter_table$Letters[c])}
  #add mean
  population_means<-Populations; names(population_means)<-Populations
  for (c in Populations){population_means[c]<-round(mean(data_table$svs_prop[which(data_table$population==c)], na.rm = TRUE), digits = 3)}
  for (c in Populations){tabla_group[p,c]<-paste(population_means[c], tabla_group[p,c])}
  
  #arrange table results
  for (c in colnames(p.values_table)){colnames(p.values_table)[which(colnames(p.values_table)==c)]<-paste(c, "=", population_means[c], sep = "")}#c 
  for (c in row.names(p.values_table)){row.names(p.values_table)[which(row.names(p.values_table)==c)]<-paste(c, "=", population_means[c], sep = "")}#c
  
  write.xlsx(p.values_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.3.2_SVs_prop_windows_unified_WILCOX_across_pops.xlsx", sheetName = p, append = TRUE)
  
}#p

#format for paper
tabla_group<-cbind(c("Pericentomeric", "Distal proximal", "Distal telomeric", "Distal proximal", "Distal telomeric", "Pericentromeric", "Distal proximal", "Distal telomeric"), tabla_group)
#first check real order match names en la tabla
tabla_group<-cbind(c("Chromosome regions", "", "", "Coldspots", "", "Hotspots", "", ""), tabla_group)
#save
write.xlsx(tabla_group, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.3.3_SVs_prop_windows_unified_WILCOX_GROUPING_across_pops.xlsx", row.names = FALSE)
####################################################################################################
