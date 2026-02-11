library(data.table)
library(xlsx)
library(agricolae)
library(multcompView)

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes

pop_info<-read.csv("/home/fcasale/Desktop/Paper_2/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

met_context<-c("chg", "cpg", "chh")

#################################################################################

#get windows types
window_type_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.1.2_methylation_prop_per_window_type.RDS")

window_types<-names(window_type_list[[1]])[3:4]

#subset parents methy data with window type per pop per chr
methy_list_RAW<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/A_methylation_per_parent_level_and_Cs_win=",win.size,".RDS", sep = ""))
window_type_parents<-window_type_list
for (e in window_types){ cat(e); cat(": ")
for (p in Populations){ cat(p); cat("-")
#select parents
P1<-pop_info[which(pop_info[,1]==p),2]
P2<-pop_info[which(pop_info[,1]==p),3]
for (c in 1:7){
windows<-as.character(window_type_list[[1]][[e]][[p]][[c]][,2])    
for (v in met_context){
#add P1
window_type_parents[[v]][[e]][[p]][[c]][windows,3]<-methy_list_RAW[[c]][[v]]$methy_data[windows,P1] 
#add P2
window_type_parents[[v]][[e]][[p]][[c]]<-cbind(window_type_parents[[v]][[e]][[p]][[c]], NA)
window_type_parents[[v]][[e]][[p]][[c]][windows,4]<-methy_list_RAW[[c]][[v]]$methy_data[windows,P2] 
colnames(window_type_parents[[v]][[e]][[p]][[c]])[3:4]<-c(P1,P2)
}#v  
}#c  
}#p
cat(fill = TRUE)  
}#e
saveRDS(window_type_parents, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.1_methylation_prop_per_window_type_parents.RDS")

#unify chrs
for (v in met_context){ cat(v); cat(": ")
  for (e in window_types){ cat(e); cat(": ")
    for (p in Populations){ cat(p); cat("-")
      parent_data<-window_type_parents[[v]][[e]][[p]][[1]][,3:4]
      for (c in 2:7){parent_data<-rbind(parent_data, window_type_parents[[v]][[e]][[p]][[c]][,3:4])}#c
      window_type_parents[[v]][[e]][[p]]<-parent_data
    }#p
  }#e
  cat(fill = TRUE)
}#v
saveRDS(window_type_parents, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.2_methylation_prop_per_window_type_parents_unified.RDS")
window_type_parents<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.2_methylation_prop_per_window_type_parents_unified.RDS")

#add distal windows in pop as another group
window_type_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.1.4_methylation_prop_per_window_type_unified.RDS")

###################### check if means are different with ANOVA + Tukey #############################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.3_all_methy_types_windows_unified_TUKEY.xlsx")

window_types<-c("coldspots","hotspots")
parents<-c("HOR8160", "Unumli-Arpa", "SprattArcher")

treatments<-c(parents, "Distal")
n_tests<-length(treatments)*(length(treatments)-1)/2

tabla<-matrix(ncol = length(treatments)+1, nrow = 0)
colnames(tabla)<-c("", treatments)
for (e in window_types){ cat(e, fill = TRUE)
for (p in Populations){ cat(p); cat(": ")   
  tabla_p<-matrix(ncol = length(treatments), nrow = 3) 
  colnames(tabla_p)<-treatments
  row.names(tabla_p)<-met_context
  #select parents
  P1<-pop_info[which(pop_info[,1]==p),2]
  P2<-pop_info[which(pop_info[,1]==p),3]
  for (v in met_context){ cat(v); cat("-")
    #gather data to make linear regression
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "parent")
    for (q in c(P1,P2)){
      data_table_q<-cbind(window_type_parents[[v]][[e]][[p]][,q], rep(q, nrow(window_type_parents[[v]][[e]][[p]])))
      colnames(data_table_q)<-c("methylation", "parent")
      data_table<-rbind(data_table, data_table_q)
    }#e
    #add distal windows in pop
    data_table<-rbind(data_table,cbind(window_type_list[[v]]$distal_windows[[p]], "Distal"))
    data_table<-as.data.frame(data_table)
    data_table$methylation<-as.numeric(as.character(data_table$methylation))
    LM<-lm(methylation~parent, data = as.data.frame(data_table))
    #ANOVA + Tukey
    tukey<-HSD.test(LM, "parent", alpha = 0.05/n_tests)$groups
    #fill results table
    for (q in levels(data_table$parent)){tabla_p[v,q]<-paste(round(tukey[q,1], digits = 3), tukey[q,2])}#e
  }#v
  tabla_p<-cbind(c(p, met_context), rbind("", tabla_p))
  tabla_p<-rbind(c(e, rep("", length(treatments))), tabla_p)
  tabla<-rbind(tabla, tabla_p)
  cat(fill = TRUE)
}#p  
}#e
write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.3_all_methy_types_windows_unified_TUKEY.xlsx")
####################################################################################################


########### check if means are different with multiple kruskall-wallis (WILCOX) ####################

file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.4.1_all_methy_types_windows_unified_WILCOX.xlsx")
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.4.2_all_methy_types_windows_unified_WILCOX_GROUPING.xlsx")

parents<-c("HOR8160", "Unumli-Arpa", "SprattArcher")

tabla_group<-matrix(ncol = length(parents)+1+1, nrow = 0); colnames(tabla_group)<-c("",parents, "Distal")
for (e in window_types){ cat(e, fill = TRUE)
tabla_group<-rbind(tabla_group, cbind(e, "", "", "", ""))  
for (p in Populations){ cat(p); cat("-") 
  #select parents
  P1<-pop_info[which(pop_info[,1]==p),2]
  P2<-pop_info[which(pop_info[,1]==p),3]
  treatments<-c(P1, P2, "Distal")
  #table for wilcox results
  tabla_p<-matrix(ncol = length(treatments), nrow = 0)
  #table for wilcox grouping
  tabla_g<-matrix(ncol = length(treatments), nrow = 3); colnames(tabla_g)<-treatments; row.names(tabla_g)<-met_context
  for (v in met_context){ cat(v); cat("-")
  #gather data to make linear regression
  data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("methylation", "parent")
  for (q in c(P1,P2)){
    data_table_q<-cbind(window_type_parents[[v]][[e]][[p]][,q], rep(q, nrow(window_type_parents[[v]][[e]][[p]])))
    colnames(data_table_q)<-c("methylation", "parent")
    data_table<-rbind(data_table, data_table_q)
  }#e
  #add distal windows in pop
  data_table<-rbind(data_table,cbind(window_type_list[[v]]$distal_windows[[p]], "Distal"))
  data_table<-as.data.frame(data_table)
  data_table$methylation<-as.numeric(as.character(data_table$methylation))
  #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
  PWT<-pairwise.wilcox.test(data_table$methylation, data_table$parent)
  p.values_table<-format(round(PWT$p.value, digits = 3), scientific = FALSE)
  #group by letters
  p_val_list<-list()
  groups_to_compare<-treatments
  for (c in groups_to_compare){ 
    p.val<-c() 
    if (any(row.names(p.values_table)==c)){p.val<-c(p.val,p.values_table[c,])} 
    if (any(colnames(p.values_table)==c)){p.val<-c(p.val,p.values_table[,c])}
    p.val<-strsplit(p.val, " ")
    p.val<-unlist(lapply(p.val, function(x){if(any(x=="")){x<-x[-which(x=="")]};return(x)}))
    if (any(p.val=="NA")){p.val<-p.val[-which(p.val=="NA")]}
    p_val_list[[c]]<-list()
    for (r in groups_to_compare[-which(c(groups_to_compare)==c)]){p_val_list[[c]][[r]]<-as.numeric(p.val[r])}
  }#c  
    p_val_table<-matrix(ncol = length(treatments), nrow = length(treatments))
    colnames(p_val_table)<-treatments; row.names(p_val_table)<-treatments
    for (c in colnames(p_val_table)){ 
      p_val_table[c,c]<-1
      for (r in names(p_val_list[[c]])){p_val_table[r,c]<-as.numeric(p_val_list[[c]][[r]])}
    }
    #paso los NA a FALSE para que funcione la funcion pero despues poner en "-"
    letter_table<-multcompLetters(p_val_table, compare = "<", threshold = 0.05/n_tests)
    for (c in treatments){tabla_g[v,c]<-gsub(" ", "", letter_table[[2]][c])}
    #add mean
    window_type_means<-treatments; names(window_type_means)<-treatments
    for (c in treatments){window_type_means[c]<-round(mean(data_table$methylation[which(data_table$parent==c)], na.rm = TRUE), digits = 3)}
    for (c in treatments){tabla_g[v,c]<-paste(window_type_means[c], tabla_g[v,c])}
    #arrange table results
    for (c in colnames(p.values_table)){colnames(p.values_table)[which(colnames(p.values_table)==c)]<-paste(c, "=", window_type_means[c], sep = "")}#c 
    for (c in row.names(p.values_table)){row.names(p.values_table)[which(row.names(p.values_table)==c)]<-paste(c, "=", window_type_means[c], sep = "")}#c
    p.values_table<-cbind(c(v,row.names(p.values_table)), rbind(colnames(p.values_table), p.values_table))
    tabla_p<-rbind(tabla_p, "", p.values_table)
  }#v
  #arrange table group
  tabla_g2<-matrix(ncol = length(parents)+1+1, nrow = length(met_context)+1)
  row.names(tabla_g2)<-c("",row.names(tabla_g)); colnames(tabla_g2)<-colnames(tabla_group)
  for (c in colnames(tabla_g)){for (v in met_context){tabla_g2[v,c]<-tabla_g[v,c]}}
  tabla_g2[1,]<-""
  tabla_g2[,1]<-c(p, met_context)
  tabla_g2[is.na(tabla_g2)]<-"-"
  tabla_group<-rbind(tabla_group, tabla_g2)
  #save table results
  colnames(tabla_p)<-NULL
  write.xlsx(tabla_p, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.4.1_all_methy_types_windows_unified_WILCOX.xlsx", sheetName = paste(e,p, sep = "_"), append = TRUE)
}#p
}#e 
#save table groups
write.xlsx(tabla_group, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.3.4.2_all_methy_types_windows_unified_WILCOX_GROUPING.xlsx", row.names = FALSE)

#####################################################################################################

############################# graphic boxplots ######################################################

# graphic_table<-matrix(nncol = 5, nrow = 0)
# colnames(graphic_table)<-c("Population", "window_type", "met_context","met_level")