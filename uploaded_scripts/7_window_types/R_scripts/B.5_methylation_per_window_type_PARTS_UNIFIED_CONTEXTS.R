library(data.table)
library(xlsx)
library(agricolae)
library(multcompView)

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

window_type_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.2.2_methylation_prop_per_window_type.RDS")[[1]]
window_types<-names(window_type_list)

methy_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/4_methylation/Results/B_methylation_per_population/B.2.1_methy_unified_per_pop_win=",win.size,".RDS", sep = ""))

  methy_list_windows<-list()
  for (e in window_types){  cat(e); cat("-")
    methy_list_windows[[e]]<-list()
    for (p in Populations){
      methy_list_windows[[e]][[p]]<-list()
      for (c in 1:7){    
        windows<-as.character(window_type_list[[e]][[p]][[c]][,2])
        methylation<-methy_list[[c]][windows,p]
        methy_list_windows[[e]][[p]][[c]]<-window_type_list[[e]][[p]][[c]]
        row.names(methy_list_windows[[e]][[p]][[c]])<-methy_list_windows[[e]][[p]][[c]][,2]
        methy_list_windows[[e]][[p]][[c]][windows,3]<-methylation
        colnames(methy_list_windows[[e]][[p]][[c]])[3]<-"methylation_level"
      }#c
    }#p  
  }#e
saveRDS(methy_list_windows, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.1_methylation_prop_per_window_type_UNI.RDS")

#################### check if differences among chrs within populations ###########################

file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.2_all_methy_types_windows_DIF_AMONG_CHR_TUKEY_UNI.xlsx")

n_tests<-length(window_types)*(length(window_types)-1)/2

tabla<-matrix(ncol = 7+1, nrow = 0)
colnames(tabla)<-c("", paste(1:7, "H", sep = ""))
for (p in Populations){ cat(p); cat(" - ")   
    tabla_p<-matrix(ncol = 7, nrow = length(window_types)+1); row.names(tabla_p)<-c(p, window_types)
    for (e in window_types){
      data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "window_type")
      for (c in 1:7){
        data_table_c<-cbind(methy_list_windows[[e]][[p]][[c]][,3], rep(c, length(methy_list_windows[[e]][[p]][[c]][,3])))
        colnames(data_table_c)<-c("methylation", "window_type")
        data_table<-rbind(data_table, data_table_c)
      }#c
      data_table<-as.data.frame(data_table)
      data_table$methylation<-as.numeric(as.character(data_table$methylation))
      LM<-lm(methylation~window_type, data = as.data.frame(data_table))
      #ANOVA + Tukey
      tukey<-HSD.test(y = LM, trt = "window_type", alpha = 0.05/n_tests, unbalanced = TRUE, group = TRUE)$groups
      #PWT<-pairwise.wilcox.test(data_table$methylation, data_table$window_type, p.adjust.method = "bonferroni", paired = FALSE)
      for (c in 1:7){tabla_p[e,c]<-paste(round(tukey[paste(c),1], digits = 3), tukey[paste(c),2])}#e
    }#e
  tabla_p[p,]<-""  
  tabla_p<-cbind(row.names(tabla_p), tabla_p)
  tabla<-rbind(tabla, tabla_p)
}#p

write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.2_all_methy_types_windows_DIF_AMONG_CHR_TUKEY_UNI.xlsx", row.names = FALSE)
###################################################################################################

#unify chrs
  for (e in names(methy_list_windows)){ cat(e); cat(": ")
    for (p in Populations){ cat(p); cat("-")
      pop_data<-c()
      for (c in 1:7){pop_data<-c(pop_data, methy_list_windows[[e]][[p]][[c]][,3])}#c
      methy_list_windows[[e]][[p]]<-pop_data
    }#p
  }#e
  
saveRDS(methy_list_windows, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.3_methylation_prop_per_window_type_unified_UNI.RDS")

###################### check if means are different with ANOVA + Tukey #############################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.4_all_methy_types_windows_unified_TUKEY_UNI.xlsx")

tabla<-matrix(ncol = length(window_types), nrow = 3)
row.names(tabla)<-Populations
colnames(tabla)<-window_types
for (p in Populations){ cat(p); cat(" - ")
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "window_type")
    for (e in window_types){
      data_table_e<-cbind(methy_list_windows[[e]][[p]], rep(e, length(methy_list_windows[[e]][[p]])))
      colnames(data_table_e)<-c("methylation", "window_type")
      data_table<-rbind(data_table, data_table_e)
    }#e
    data_table<-as.data.frame(data_table)
    data_table$methylation<-as.numeric(as.character(data_table$methylation))
    LM<-lm(methylation~window_type, data = as.data.frame(data_table))
    # #ANOVA + Tukey
    tukey<-HSD.test(LM, "window_type", alpha = 0.05/n_tests)$groups
    # #fill results table
    for (e in window_types){tabla[p,e]<-paste(round(tukey[e,1], digits = 3), tukey[e,2])}#e
}#p  

write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.4_all_methy_types_windows_unified_TUKEY_UNI.xlsx")
####################################################################################################

########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.5_all_methy_types_windows_unified_WILCOX_UNI.xlsx")
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.6_all_methy_types_windows_unified_WILCOX_GROUPING_UNI.xlsx")

tabla<-matrix(ncol = length(window_types), nrow = 3)
row.names(tabla)<-Populations
colnames(tabla)<-window_types
for (p in Populations){ cat(p); cat("-") 
  #table for wilcox results
  tabla_p<-matrix(ncol = length(window_types), nrow = 1)
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "window_type")
    for (e in window_types){
      data_table_e<-cbind(methy_list_windows[[e]][[p]], rep(e, length(methy_list_windows[[e]][[p]])))
      colnames(data_table_e)<-c("methylation", "window_type")
      data_table<-rbind(data_table, data_table_e)
    }#e
    data_table<-as.data.frame(data_table)
    data_table$methylation<-as.numeric(as.character(data_table$methylation))
    #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
    PWT<-pairwise.wilcox.test(data_table$methylation, data_table$window_type)
    p.values_table<-format(round(PWT$p.value, digits = 3), scientific = FALSE)
    
    #group by letters
    p_val_list<-list()
    for (c in window_types){ 
      p.val<-c()
      if (any(row.names(p.values_table)==c)){p.val<-c(p.val,p.values_table[c,])}
      if (any(colnames(p.values_table)==c)){p.val<-c(p.val,p.values_table[,c])}
      p.val<-strsplit(p.val, " ")
      p.val<-unlist(lapply(p.val, function(x){if(any(x=="")){x<-x[-which(x=="")]};return(x)}))
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
    for (c in window_types){tabla[p,c]<-gsub(" ", "", letter_table[[2]][c])}
    #add mean
    window_type_means<-window_types; names(window_type_means)<-window_types
    for (c in window_types){window_type_means[c]<-round(mean(data_table$methylation[which(data_table$window_type==c)], na.rm = TRUE), digits = 3)}
    for (c in window_types){tabla[p,c]<-paste(window_type_means[c], tabla[p,c])}
    
    #arrange table results
    for (c in colnames(p.values_table)){colnames(p.values_table)[which(colnames(p.values_table)==c)]<-paste(c, "=", window_type_means[c], sep = "")}#c 
    for (c in row.names(p.values_table)){row.names(p.values_table)[which(row.names(p.values_table)==c)]<-paste(c, "=", window_type_means[c], sep = "")}#c
  write.xlsx(p.values_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.5_all_methy_types_windows_unified_WILCOX_UNI.xlsx", sheetName = p, append = TRUE)
}#p
#format por paper
#first check real order match names
tabla<-rbind(c("Pericentomeric", "Distal proximal", "Distal telomeric", "Distal proximal", "Distal telomeric", "Pericentromeric", "Distal proximal", "Distal telomeric"),tabla)
colnames(tabla)<-c("Chromosome regions", "", "", "Coldspots", "", "Hotspots", "", "")
#save
write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.6_all_methy_types_windows_unified_WILCOX_GROUPING_UNI.xlsx", row.names = FALSE)

#####################################################################################################

