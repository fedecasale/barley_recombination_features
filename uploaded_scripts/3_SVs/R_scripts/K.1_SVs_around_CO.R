layer<-3
win.size<-10000

breakpoint_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")
breakpoint_list<-breakpoint_list[[layer]]
Populations<-names(breakpoint_list)

SVs_proportion<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.2_SVs_proportion_per_window=10000.RDS")

#normalize SV data based on distal windows only
pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere.RDS")
#Antes normalizaba por todo
data<-c()
for (p in Populations){for (c in 1:7){
data_c<-unlist(SVs_proportion[[p]][[c]])
data_c<-data_c[which(as.numeric(names(data_c))>as.numeric(pericentromere_pops[[p]][c,4]))]
distal_windows_data<-data_c[which(as.numeric(names(data_c))<as.numeric(pericentromere_pops[[p]][c,5]))]
data<-c(data, distal_windows_data)
}}
MAX<-max(data, na.rm = TRUE); MIN<-min(data, na.rm = TRUE); MEAN<-mean(data, na.rm = TRUE); SD<-sd(data, na.rm = TRUE)
#normalize<-function(x){x<-(x-MEAN)/(MAX-MIN); return(x)}
standarize<-function(x){x<-(x-MEAN)/SD; return(x)}

#for (p in Populations){for (c in 1:7){SVs_proportion[[p]][[c]]<-lapply(SVs_proportion[[p]][[c]], normalize)}}
for (p in Populations){for (c in 1:7){SVs_proportion[[p]][[c]]<-lapply(SVs_proportion[[p]][[c]], standarize)}}

SVs_proportion[[p]][[c]]

#Ahora normalizo por cada pop
# normalize_list<-list()
# for (p in Populations){
# for (c in 1:7){
# data_c<-unlist(SVs_proportion[[p]][[c]])  
# data_c<-data_c[which(as.numeric(names(data_c))>as.numeric(pericentromere_pops[[p]][c,4]))]
# distal_windows_data<-data_c[which(as.numeric(names(data_c))<as.numeric(pericentromere_pops[[p]][c,5]))]
# data<-c(data, distal_windows_data)
# }#c
# MAX<-max(data, na.rm = TRUE); MIN<-min(data, na.rm = TRUE); MEAN<-mean(data, na.rm = TRUE)
# normalize_list[[p]]<-list(MAX, MIN, MEAN); names(normalize_list[[p]])<-c("MAX", "MIN", "MEAN")
# }#p
# normalize<-function(x){x<-(x-normalize_list[[p]]$MEAN)/(normalize_list[[p]]$MAX-normalize_list[[p]]$MIN); return(x)}
# for (p in Populations){for (c in 1:7){SVs_proportion[[p]][[c]]<-lapply(SVs_proportion[[p]][[c]], normalize)}}



pop_table<-matrix(ncol = 21, nrow = 3); row.names(pop_table)<-Populations

for (p in Populations){cat(p);cat(":")

svs_prop_around_break<-matrix(ncol = 21, nrow = 7)
  
for (c in 1:7){ cat(c);cat("-")

#get windows of CO breakpoint
chr_windows<-as.numeric(names(SVs_proportion[[p]][[c]]))
chr_breaks<-as.numeric(breakpoint_list[[p]][[c]][,6])
break_windows<-c()
for (b in 1:length(chr_breaks)){break_windows<-c(break_windows, chr_windows[which(chr_windows>chr_breaks[b])][1])}
#any((break_windows-chr_breaks)>win.size)

#get windows around
if (any(break_windows<=100000)){break_windows<-break_windows[-which(break_windows<=100000)]}
windows_around_break_table<-matrix(ncol = 21, nrow = 0)
for (b in 1:length(break_windows)){
windows_around_break<-chr_windows[(which(chr_windows==break_windows[b])-10):(which(chr_windows==break_windows[b])+10)]  
windows_around_break_table<-rbind(windows_around_break_table, windows_around_break)
}#b
  
#convert windows to SVs values
values<-unlist(SVs_proportion[[p]][[c]])
window_retrieve_value<-function(x){return(values[as.character(x)])}
windows_around_break_table<-apply(windows_around_break_table, 1, FUN = window_retrieve_value)
windows_around_break_table<-transpose(as.data.frame(windows_around_break_table))

#average values
for (j in 1:ncol(windows_around_break_table)){
svs_prop_around_break[c,j]<-mean(windows_around_break_table[,j], na.rm = TRUE)
}

}#c  

#average values
for (j in 1:ncol(svs_prop_around_break)){pop_table[p,j]<-mean(svs_prop_around_break[,j], na.rm = TRUE)}

cat(fill = TRUE)

}#p

plot(pop_table[2,]~c(1:ncol(pop_table)))
lines(pop_table[2,]~c(1:ncol(pop_table)))

#take means of pericentromere
distal_means<-list()
for (p in Populations){ 
data_distal<-c()
for (c in 1:7){
data_c<-data_c[which(as.numeric(names(data_c))>as.numeric(pericentromere_pops[[p]][c,4]))]
distal_windows_data<-data_c[which(as.numeric(names(data_c))<as.numeric(pericentromere_pops[[p]][c,5]))]
data_distal<-c(data_distal, unlist(SVs_proportion[[p]][[c]][names(distal_windows_data)]))
}#c
distal_means[[p]]<-mean(data_distal, na.rm = TRUE)
}#p

graphic_table<-matrix(ncol = 3, nrow=0); colnames(graphic_table)<-c("Population", "SVs_proportion", "window")
for (p in Populations){graphic_table<-rbind(graphic_table,cbind(p, pop_table[p,4:18], c(-6:8)-1))}
graphic_table<-as.data.frame(graphic_table)
graphic_table$SVs_proportion<-as.numeric(as.character(graphic_table$SVs_proportion))
graphic_table$window<-as.numeric(as.character(graphic_table$window))

ggplot()+
  geom_line(data = graphic_table, aes(x = window, y = SVs_proportion, color = Population)) +
  #geom_hline(yintercept = distal_means$HvDRR13[[1]])+
  
  scale_y_continuous(name = "Window SVs' relative proportion") +
  scale_x_continuous(name = "10 kb windows")+
  
  theme(
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA),
  strip.background = element_blank(),
  #strip.text.x = element_blank(),
  axis.ticks = element_line(),
  legend.key = element_blank(),
  #axis.title.y.left = element_blank(),
  #axis.title.x.bottom = element_blank()
) #+
dev.off()
