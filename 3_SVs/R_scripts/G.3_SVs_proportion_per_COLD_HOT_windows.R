library(data.table)
library(xlsx)
library(agricolae)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[,1]

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

#parents.code<-as.matrix(read.table("/scratch/Federico/Paper_2/samplenames.txt"))

met_context<-c("chg", "cpg", "chh")

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
if (win.size == 10000){Populations<-Populations[which(Populations%in%c("HvDRR13","HvDRR27", "HvDRR28"))]} 

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS")
SVs<-names(sv_list)
sv_sizes<-names(sv_list[[1]][[1]][[1]])

#################################################################################

#change sv distance (to the window) by sv proportion in the window
window_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.1_windows_closest_SV.RDS")
SVs_proportion<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.2_SVs_proportion_per_window=10000.RDS")
window_SV_list2<-list()
for (e in names(window_SV_list)){ cat(e, fill = TRUE)
  window_SV_list2[[e]]<-list()  
  for (p in Populations){ cat(p); cat(": ")   
    window_SV_list2[[e]][[p]]<-list()  
    for (c in 1:7){ cat(c); cat("-")
      names(SVs_proportion[[p]][[c]])<-gsub(" ", "",format(as.numeric(names(SVs_proportion[[p]][[c]])), scientific = FALSE))
      windows_table<-window_SV_list[[e]][[1]][[p]][[c]][[1]]
      windows<-as.numeric(windows_table[,1])+(win.size/2)+1
      windows<-gsub(" ", "", as.character(format(windows, scientific = FALSE)))
      row.names(windows_table)<-windows
      colnames(windows_table)[3]<-"win_proportion"
      for (w in windows){windows_table[w,3]<-SVs_proportion[[p]][[c]][[w]]}#w
      window_SV_list2[[e]][[p]][[c]]<-windows_table
    }#c
  }#p
}#e  

saveRDS(window_SV_list2, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.3.1_windows_closestSV_SV_proportion.RDS")
window_SV_list<-window_SV_list2; rm(window_SV_list2)
#window_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.3.1_windows_closestSV_SV_proportion.RDS")
#################################################################################

#unify chrs
for (e in names(window_SV_list)){ cat(e, fill = TRUE)
  for (p in Populations){ cat(p); cat(": ")
    pop_data<-c()
    for (c in 1:7){pop_data<-c(pop_data, window_SV_list[[e]][[p]][[c]][,3])}#c
    window_SV_list[[e]][[p]]<-pop_data
  }#p
}#e

###################### check if means are different with ANOVA + Tukey #############################
file.remove("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.3.2_windows_unified_TUKEY.xlsx")

tabla<-matrix(ncol = length(names(window_SV_list)), nrow = length(Populations))
colnames(tabla)<-names(window_SV_list); row.names(tabla)<-Populations

for (p in Populations){ cat(p); cat("-")
  #run anova + tukey
  data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("distance", "window_type")
  for (e in names(window_SV_list)){
    data_table_e<-cbind(window_SV_list[[e]][[p]], rep(e, length(window_SV_list[[e]][[p]])))
    colnames(data_table_e)<-c("distance", "window_type")
    data_table<-rbind(data_table, data_table_e)
  }#e
  data_table<-as.data.frame(data_table)
  data_table$distance<-as.numeric(as.character(data_table$distance))
  LM<-lm(distance~window_type, data = as.data.frame(data_table))
  #ANOVA + Tukey
  tukey<-HSD.test(LM, "window_type")$groups
  
  #fill results table
  for (e in names(window_SV_list)){tabla[p,e]<-paste(round(tukey[e,1], digits = 3), tukey[e,2])}#e
}#p

write.xlsx(tabla, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.3.2_windows_unified_TUKEY.xlsx")
####################################################################################################


########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.3.3_windows_unified_WILCOX.xlsx")

for (p in Populations){ cat(p); cat("-")
  data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("distance", "window_type")
  for (e in names(window_SV_list)){
    data_table_e<-cbind(window_SV_list[[e]][[p]], rep(e, length(window_SV_list[[e]][[p]])))
    colnames(data_table_e)<-c("distance", "window_type")
    data_table<-rbind(data_table, data_table_e)
  }#e
  data_table<-as.data.frame(data_table)
  data_table$distance<-as.numeric(as.character(data_table$distance))
  #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
  PWT<-pairwise.wilcox.test(data_table$distance, data_table$window_type, p.adjust.method = "BH")
  p.values_table<-format(round(PWT$p.value, digits = 3), scientific = FALSE)
  for (c in colnames(p.values_table)){colnames(p.values_table)[which(colnames(p.values_table)==c)]<-paste(c, "=", round(mean(data_table$distance[which(data_table$window_type==c)]), digits = 3), sep = "")}#c 
  for (c in row.names(p.values_table)){row.names(p.values_table)[which(row.names(p.values_table)==c)]<-paste(c, "=", round(mean(data_table$distance[which(data_table$window_type==c)]), digits = 3), sep = "")}#c
  write.xlsx(p.values_table, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.3.3_windows_unified_WILCOX.xlsx", sheetName = p, append = TRUE)
}#p

####################################################################################################



