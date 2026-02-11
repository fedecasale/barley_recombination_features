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

window_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.2_window_closest_SV.RDS")


###################### check if means are different with ANOVA + Tukey #############################


categories<-c("with_COs", "without_COs", "big_groups", "coldspots")

file.remove(paste("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.1.1_windows_unified_TUKEY.xlsx", sep = ""))

for (d in categories){

if (d == "with_COs"){window_types<-names(window_SV_list)} 
if (d == "without_COs"){window_types<-names(window_SV_list)[-which(names(window_SV_list)=="COs")]} 
if (d == "big_groups"){window_types<-names(window_SV_list)[which(names(window_SV_list)%in%c("distal_windows","hotspots","hotspots_first","coldspots"))]} 
if (d == "coldspots"){window_types<-names(window_SV_list)[which(names(window_SV_list)%in%c("coldspots", "coldspots_high_methy", "coldspots_low_methy"))]} 
  
    
tabla_0<-matrix(ncol = 2+length(names(window_SV_list)), nrow = 0)
colnames(tabla_0)<-c("sv type", "Populations", names(window_SV_list))

for (s in c("small", "big")){ cat(s); cat(":")

tabla<-matrix(ncol = length(names(window_SV_list)), nrow = length(Populations))
colnames(tabla)<-names(window_SV_list); row.names(tabla)<-Populations

for (p in Populations){ cat(p); cat("-")
  #run anova + tukey
  data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("distance", "window_type")
  for (e in window_types){
    data_table_e<-cbind(window_SV_list[[e]][[p]][[s]], rep(e, length(window_SV_list[[e]][[p]][[s]])))
    colnames(data_table_e)<-c("distance", "window_type")
    data_table<-rbind(data_table, data_table_e)
  }#e   
  data_table<-as.data.frame(data_table)
  data_table$distance<-as.numeric(as.character(data_table$distance))
  LM<-lm(distance~window_type, data = as.data.frame(data_table))
  #ANOVA + Tukey
  tukey<-HSD.test(LM, "window_type")$groups

  #fill results table
  for (e in window_types){tabla[p,e]<-paste(round(tukey[e,1], digits = 0), tukey[e,2])}#e
}#p

#add to sv table
tabla_0<-rbind(tabla_0, cbind(c(paste(s),"",""), Populations, tabla))

}#s

cat(fill = TRUE)
write.xlsx(tabla_0, paste("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.1.1_windows_unified_TUKEY.xlsx", sep = ""), sheetName = d, append = TRUE)

}#d
####################################################################################################


########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.1.2_windows_unified_WILCOX.xlsx")

for (p in Populations){ cat(p); cat("-")
  data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("distance", "window_type")
  for (e in names(window_SV_list)){
    data_table_e<-cbind(window_SV_list[[e]][[p]][["big"]], rep(e, length(window_SV_list[[e]][[p]][["big"]])))
    colnames(data_table_e)<-c("distance", "window_type")
    data_table<-rbind(data_table, data_table_e)
  }#e
  data_table<-as.data.frame(data_table)
  data_table$distance<-as.numeric(as.character(data_table$distance))
  #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
  PWT<-pairwise.wilcox.test(data_table$distance, data_table$window_type, p.adjust.method = "BH")
  p.values_table<-format(round(PWT$p.value, digits = 3), scientific = FALSE)
  for (c in colnames(p.values_table)){colnames(p.values_table)[which(colnames(p.values_table)==c)]<-paste(c, "=", round(mean(data_table$distance[which(data_table$window_type==c)]), digits = 0), sep = "")}#c
  for (c in row.names(p.values_table)){row.names(p.values_table)[which(row.names(p.values_table)==c)]<-paste(c, "=", round(mean(data_table$distance[which(data_table$window_type==c)]), digits = 0), sep = "")}#c
  write.xlsx(p.values_table, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.1.2_windows_unified_WILCOX.xlsx", sheetName = p, append = TRUE)
}#p

####################################################################################################