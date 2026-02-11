library(data.table)
library(xlsx)
library(agricolae)
library(multcompView)

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

methy_list_windows<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.2.4_methylation_prop_per_window_type_unified.RDS")
window_types<-names(methy_list_windows[[1]])

# add average
methy_list_windows$AVE<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.3_methylation_prop_per_window_type_unified_UNI.RDS")
met_contexts<-names(methy_list_windows)[c(2,1,3,4)]

#################### check if differences among chrs within populations ###########################

#need to recalculate
n_tests<-length(Populations)*(length(Populations)-1)/2

###################### check if means are different with ANOVA + Tukey #############################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.4.1_all_methy_types_windows_unified_TUKEY_across_pops.xlsx")

tabla<-matrix(ncol = length(Populations)+1, nrow = 0)
colnames(tabla)<-c("", Populations)
for (v in met_contexts){ cat(v); cat("-")
tabla_v<-matrix(ncol = length(Populations), nrow = length(window_types)) 
colnames(tabla_v)<-Populations
row.names(tabla_v)<-window_types
for (e in window_types){ cat(e); cat(": ")
    #gather data to make linear regression
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "population")
    for (p in Populations){
      data_table_p<-cbind(methy_list_windows[[v]][[e]][[p]], rep(p, length(methy_list_windows[[v]][[e]][[p]])))
      colnames(data_table_p)<-c("methylation", "population")
      data_table<-rbind(data_table, data_table_p)
    }#p
    data_table<-as.data.frame(data_table)
    data_table$methylation<-as.numeric(as.character(data_table$methylation))
    LM<-lm(methylation~population, data = as.data.frame(data_table))
    # #ANOVA + Tukey
    tukey<-HSD.test(LM, "population", alpha = 0.05/n_tests)$groups
    # #fill results table
    for (p in Populations){tabla_v[e,p]<-paste(round(tukey[p,1], digits = 3), tukey[p,2])}#e
  }#e
  tabla_v<-cbind(c(v, window_types), rbind("",tabla_v))
  tabla<-rbind(tabla, tabla_v)
  cat(fill = TRUE)
}#v  

write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.4.1_all_methy_types_windows_unified_TUKEY_across_pops.xlsx")
####################################################################################################

########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.4.2_all_methy_types_windows_unified_WILCOX_across_pops.xlsx")
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.4.2_all_methy_types_windows_unified_WILCOX_GROUPING_across_pops.xlsx")

tabla_group<-matrix(ncol = length(Populations)+1, nrow = 0); colnames(tabla_group)<-c("",Populations)
for (p in window_types){ cat(p); cat("-") 
  #table for wilcox results 
  tabla_p<-matrix(ncol = length(Populations), nrow = 0)
  #table for wilcox grouping
  tabla_g<-matrix(ncol = length(Populations), nrow = length(met_contexts)); colnames(tabla_g)<-Populations; row.names(tabla_g)<-met_contexts
  for (v in met_contexts){ cat(v); cat("-")
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "population")
    for (e in Populations){
      data_table_e<-cbind(methy_list_windows[[v]][[p]][[e]], rep(e, length(methy_list_windows[[v]][[p]][[e]])))
      colnames(data_table_e)<-c("methylation", "population")
      data_table<-rbind(data_table, data_table_e)
    }#e
    data_table<-as.data.frame(data_table)
    data_table$methylation<-as.numeric(as.character(data_table$methylation))
    #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
    PWT<-pairwise.wilcox.test(data_table$methylation, data_table$population)
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
    for (c in Populations){tabla_g[v,c]<-gsub(" ", "", letter_table$Letters[c])}
    #add mean
    population_means<-Populations; names(population_means)<-Populations
    for (c in Populations){population_means[c]<-round(mean(data_table$methylation[which(data_table$population==c)], na.rm = TRUE), digits = 3)}
    for (c in Populations){tabla_g[v,c]<-paste(population_means[c], tabla_g[v,c])}
    
    #arrange table results
    for (c in colnames(p.values_table)){colnames(p.values_table)[which(colnames(p.values_table)==c)]<-paste(c, "=", population_means[c], sep = "")}#c 
    for (c in row.names(p.values_table)){row.names(p.values_table)[which(row.names(p.values_table)==c)]<-paste(c, "=", population_means[c], sep = "")}#c
    p.values_table<-cbind(c(v,row.names(p.values_table)), rbind(colnames(p.values_table), p.values_table))
    tabla_p<-rbind(tabla_p, "", p.values_table)
  }#v
  colnames(tabla_p)<-NULL
  write.xlsx(tabla_p, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.4.2_all_methy_types_windows_unified_WILCOX_across_pops.xlsx", sheetName = p, append = TRUE)
  #arrange table group
  tabla_g<-cbind(c(p, met_contexts), rbind("",tabla_g))
  tabla_group<-rbind(tabla_group, tabla_g)
}#p

#save
write.xlsx(tabla_group, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.4.2_all_methy_types_windows_unified_WILCOX_GROUPING_across_pops.xlsx", row.names = FALSE)

