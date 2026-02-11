library(data.table)
library(xlsx)
library(agricolae)
library(multcompView)

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

dir.create("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs")

#load window types list
window_type_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.1.1_window_types.RDS")
window_types<-names(window_type_list)

#add SVs proportion data
svs_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.2_SVs_proportion_per_window=10000.RDS", sep = ""))
#convert to table
for (p in Populations){
  for (c in 1:7){
  values<-unlist(svs_list[[p]][[c]])
  win_bottom<-as.numeric(names(svs_list[[p]][[c]]))-10000
  win_top<-as.numeric(names(svs_list[[p]][[c]]))-1
  tabla<-cbind(win_bottom, win_top, values)
  row.names(tabla)<-win_top
  svs_list[[p]][[c]]<-tabla
  }
}
#add to window types
svs_list_windows<-list()
  for (e in window_types){  cat(e); cat("-")
    svs_list_windows[[e]]<-list()
    for (p in Populations){
      svs_list_windows[[e]][[p]]<-list()
      for (c in 1:7){  
        windows<-as.character(window_type_list[[e]][[p]][[c]][,2])
        svs_list_windows[[e]][[p]][[c]]<-svs_list[[p]][[c]][windows,]
        colnames(svs_list_windows[[e]][[p]][[c]])[3]<-"svs_prop"
        # #I will take out 5% of the extremes because I think they are creating wrong comparisons
        # tabla<-svs_list_windows[[v]][[e]][[p]][[c]]
        # #limits<-quantile(tabla[,3], probs = seq(0, 1, by = 0.05), na.rm = TRUE)[c(2,20)]
        # limits<-quantile(tabla[,3], na.rm = TRUE)[c(2,4)]
        # tabla<-tabla[-unique(c(which(tabla[,3]<limits[1]), which(tabla[,3]>limits[2]))),]
        # svs_list_windows[[v]][[e]][[p]][[c]]<-tabla
      }#c
    }#p  
  }#e
saveRDS(svs_list_windows, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.1.1_svs_prop_per_window_type.RDS")

#unify chrs
for (e in names(svs_list_windows)){ cat(e); cat(": ")
  for (p in Populations){ cat(p); cat("-")
    pop_data<-c()
    for (c in 1:7){pop_data<-c(pop_data, svs_list_windows[[e]][[p]][[c]][,3])}#c
    svs_list_windows[[e]][[p]]<-pop_data
  }#p
}#e
saveRDS(svs_list_windows, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.1.2_svs_prop_per_window_type_unified.RDS")

###################### check if means are different with ANOVA + Tukey #############################

n_tests<-length(window_types)*(length(window_types)-1)/2

file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.1.3_SVs_prop_windows_unified_TUKEY.xlsx")

tabla<-matrix(ncol = length(window_types), nrow = 3)
colnames(tabla)<-window_types
row.names(tabla)<-Populations
for (p in Populations){ cat(p); cat("-")
    #gather data to make linear regression
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("svs_prop", "window_type")
    for (e in window_types){
      data_table_e<-cbind(svs_list_windows[[e]][[p]], rep(e, length(svs_list_windows[[e]][[p]])))
      colnames(data_table_e)<-c("svs_prop", "window_type")
      data_table<-rbind(data_table, data_table_e)
    }#e
    data_table<-as.data.frame(data_table)
    data_table$svs_prop<-as.numeric(as.character(data_table$svs_prop))
    LM<-lm(svs_prop~window_type, data = as.data.frame(data_table))
    #ANOVA + Tukey
    tukey<-HSD.test(LM, "window_type", alpha = 0.05/n_tests)$groups
    #fill results table
    for (e in window_types){tabla[p,e]<-paste(round(tukey[e,1], digits = 3), tukey[e,2])}#e
}#p  

write.xlsx(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.1.3_SVs_prop_windows_unified_TUKEY.xlsx")
####################################################################################################

########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.1.4.1_SVs_prop_windows_unified_WILCOX.xlsx")
file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.1.4.2_SVs_prop_windows_unified_WILCOX_GROUPING.xlsx")

#table for wilcox grouping
tabla_group<-matrix(ncol = length(window_types), nrow = 3)
colnames(tabla_group)<-window_types; row.names(tabla_group)<-Populations
for (p in Populations){ cat(p); cat("-") 
    
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("svs_prop", "window_type")
    for (e in window_types){
      data_table_e<-cbind(svs_list_windows[[e]][[p]], rep(e, length(svs_list_windows[[e]][[p]])))
      colnames(data_table_e)<-c("svs_prop", "window_type")
      data_table<-rbind(data_table, data_table_e)
    }#e
    data_table<-as.data.frame(data_table)
    data_table$svs_prop<-as.numeric(as.character(data_table$svs_prop))
    #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
    PWT<-pairwise.wilcox.test(data_table$svs_prop, data_table$window_type)
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
    for (c in window_types){tabla_group[p,c]<-gsub(" ", "", letter_table[[2]][c])}
    #add mean
    window_type_means<-window_types; names(window_type_means)<-window_types
    for (c in window_types){window_type_means[c]<-round(mean(data_table$svs_prop[which(data_table$window_type==c)], na.rm = TRUE), digits = 3)}
    for (c in window_types){tabla_group[p,c]<-paste(window_type_means[c], tabla_group[p,c])}
    
    #arrange table results
    for (c in colnames(p.values_table)){colnames(p.values_table)[which(colnames(p.values_table)==c)]<-paste(c, "=", window_type_means[c], sep = "")}#c 
    for (c in row.names(p.values_table)){row.names(p.values_table)[which(row.names(p.values_table)==c)]<-paste(c, "=", window_type_means[c], sep = "")}#c
    
  write.xlsx(p.values_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.1.4.1_SVs_prop_windows_unified_WILCOX.xlsx", sheetName = p, append = TRUE)

}#p
#format por paper
#first check real order match names en la tabla
tabla_group<-rbind(c("Pericentomeric", "Distal",  "Coldspots", "Hotspots"), tabla_group)
colnames(tabla_group)<-c("Chromosome regions", "", "", "")
#save
write.xlsx(tabla_group, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.1.4.2_SVs_prop_windows_unified_WILCOX_GROUPING.xlsx")

#####################################################################################################