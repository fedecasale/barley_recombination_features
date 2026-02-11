# library(data.table)
# library(xlsx)
library(agricolae)
library(multcompView)

variables_list<-list()
methy<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.2.0_windows_around_break_list.RDS")
for (f in names(methy[[1]])){
variables_list[[f]]<-list()
for (e in names(methy)){variables_list[[f]][[e]]<-methy[[e]][[f]]}#e
}#f  
variables_list[["SVs"]]<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.1.0_windows_around_break_list.RDS")
variables_list[["gene_density"]]<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.5.0_windows_around_break_list.RDS")

n_tests<-10

variables<-names(variables_list)
window_types<-names(variables_list[[1]])
Populations<-names(variables_list[[1]][[1]])
window_around<-as.character(1:11)
results_table<-matrix(ncol = length(window_around)+3, nrow = 0)
colnames(results_table)<-c("variable", "window_type", "population", window_around)

for (v in variables){

for (w in window_types){

tabla_w<-matrix(ncol = length(window_around), nrow = 3)
colnames(tabla_w)<-window_around
row.names(tabla_w)<-Populations

for (p in Populations){

#table for wilcox grouping
windows_around_break_table<-variables_list[[v]][[w]][[p]][[1]]
for (c in 2:7){windows_around_break_table<-rbind(windows_around_break_table, variables_list[[w]][[p]][[c]])}
colnames(windows_around_break_table)<-window_around

data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("variable_value", "window_around")
for (e in window_around){ cat(e); cat("-") 
    data_table_e<-cbind(windows_around_break_table[,e], rep(paste(e), nrow(windows_around_break_table)))
    colnames(data_table_e)<-c("variable_value", "window_around")
    data_table<-rbind(data_table, data_table_e)
  }#e
data_table<-as.data.frame(data_table)
data_table$variable_value<-as.numeric(as.character(data_table$variable_value))

#check with kruskal-wallis pairwise comparisons between group levels with correction for multiple testing
PWT<-pairwise.wilcox.test(data_table$variable_value, data_table$window_around, exact = FALSE) 
p.values_table<-format(round(PWT$p.value, digits = 4), scientific = FALSE)

#group by letters
p_val_list<-list()
for (e in window_around){ 
  p.val<-c()   
  if (any(row.names(p.values_table)==e)){p.val<-c(p.val,p.values_table[e,])}
  if (any(colnames(p.values_table)==e)){p.val<-c(p.val,p.values_table[,e])}
  p.val<-strsplit(p.val, " ")
  p.val<-unlist(lapply(p.val, function(x){if(any(x=="")){x<-x[-which(x=="")]};return(x)}))
  if (any(p.val=="NaN")){p.val[which(p.val=="NaN")]<-"NA"} #when you get cannot compute exact p-value with ties
  if (any(p.val=="NA")){p.val<-p.val[-which(p.val=="NA")]}
  p_val_list[[paste(e)]]<-list()
  for (r in window_around[-which(window_around==e)]){p_val_list[[paste(e)]][[paste(r)]]<-as.numeric(p.val[r])}
}#c  

p_val_table<-matrix(ncol = length(window_around), nrow = length(window_around))
colnames(p_val_table)<-window_around; row.names(p_val_table)<-window_around
for (c in colnames(p_val_table)){ 
  p_val_table[c,c]<-1
  for (r in names(p_val_list[[c]])){p_val_table[r,c]<-as.numeric(p_val_list[[c]][[r]])}
}
#when you get cannot compute exact p-value with ties, it is because values have no variation, so no dif 
p_val_table[which(is.na(p_val_table))]<-1
letter_table<-multcompLetters(p_val_table, compare = "<", threshold = 0.05/n_tests)
for (c in window_around){tabla_w[p,c]<-gsub(" ", "", letter_table[[1]][c])}
#add mean
window_type_means<-window_around; names(window_type_means)<-window_around
for (c in window_around){window_type_means[c]<-round(mean(data_table$variable_value[which(data_table$window_around==c)], na.rm = TRUE), digits = 3)}
for (c in window_around){tabla_w[p,c]<-paste(window_type_means[c], tabla_w[p,c])}

}#p
  
tabla_w<-cbind(v, w, cbind(Populations, tabla_w))
results_table<-rbind(results_table, tabla_w)

}#w 
  
}#v  

write.csv(results_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.6.0_window_groupping_means_and_letters.csv", row.names = FALSE)
