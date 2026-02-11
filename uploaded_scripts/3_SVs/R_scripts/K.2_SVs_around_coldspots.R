layer<-3
win.size<-10000

window_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.0_windows_list.RDS")
Populations<-names(window_list[[1]])

SVs_proportion<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.2_SVs_proportion_per_window=10000.RDS")

# #normalize SV data (ANTES: based on distal windows only)
# #pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere.RDS")
# data<-c()
# for (p in Populations){ for (c in 1:7){
#   data_c<-unlist(SVs_proportion[[p]][[c]])
#   #data_c<-data_c[which(as.numeric(names(data_c))>as.numeric(pericentromere_pops[[p]][c,4]))]
#   #distal_windows_data<-data_c[which(as.numeric(names(data_c))<as.numeric(pericentromere_pops[[p]][c,5]))]
#   #data<-c(data, distal_windows_data)
#   data<-c(data, data_c)
# }}
# MAX<-max(data, na.rm = TRUE); MIN<-min(data, na.rm = TRUE); MEAN<-mean(data, na.rm = TRUE); SD<-sd(data, na.rm = TRUE)
# normalize<-function(x){x<-(x-MEAN)/(MAX-MIN); return(x)}
# #normalize<-function(x){x<-x/(MAX-MIN); return(x)}
# #standarize<-function(x){x<-(x-MEAN)/SD; return(x)}
# 
# for (p in Populations){for (c in 1:7){SVs_proportion[[p]][[c]]<-lapply(SVs_proportion[[p]][[c]], normalize)}}
# #for (p in Populations){for (c in 1:7){SVs_proportion[[p]][[c]]<-lapply(SVs_proportion[[p]][[c]], standarize)}}

#get SVs prop of windows around target_windows
window_types<-c("coldspots_low_methy", "hotspots")
pop_table_list<-list()   
for (t in window_types){ cat(t, fill = TRUE)
pop_table<-matrix(ncol = 21, nrow = 3); row.names(pop_table)<-Populations
for (p in Populations){cat(p);cat(":")
  svs_prop_around_break<-matrix(ncol = 21, nrow = 7)
  for (c in 1:7){ cat(c);cat("-")   
    #get windows of CO breakpoint
    chr_windows<-as.numeric(names(SVs_proportion[[p]][[c]]))
    if (t %in% c("coldspots_low_methy", "coldspots_high_methy")){break_windows<-as.numeric(names(window_list[[t]][[p]][[c]]))+1} else { 
    break_windows<-as.numeric(window_list[[t]][[p]][[c]][,2])+1
    }
    #get windows around
    if (any(break_windows<=100000)){break_windows<-break_windows[-which(break_windows<=100000)]}
    windows_around<-lapply(break_windows, function(x){return(chr_windows[(which(chr_windows==x)-10):(which(chr_windows==x)+10)])})
    windows_around<-as.matrix(transpose(as.data.frame(windows_around)))
    #convert windows to SVs values
    values<-unlist(SVs_proportion[[p]][[c]])
    windows_around<-apply(windows_around, 1, function(x){return(values[as.character(x)])})
    windows_around<-transpose(as.data.frame(windows_around))
    #get average value per window
    for (j in 1:ncol(windows_around)){svs_prop_around_break[c,j]<-mean(windows_around[,j], na.rm = TRUE)}
  }#c
  #average values
  for (j in 1:ncol(svs_prop_around_break)){pop_table[p,j]<-mean(svs_prop_around_break[,j], na.rm = TRUE)}
  cat(fill = TRUE)
}#p
pop_table_list[[t]]<-pop_table
}#t

pop_table<-pop_table_list$coldspots_high_methy
plot(pop_table[2,]~c(1:ncol(pop_table)))
lines(pop_table[2,]~c(1:ncol(pop_table)))

# #take means of pericentromere
# distal_means<-list()
# for (p in Populations){ 
#   data_distal<-c()
#   for (c in 1:7){
#     data_c<-data_c[which(as.numeric(names(data_c))>as.numeric(pericentromere_pops[[p]][c,4]))]
#     distal_windows_data<-data_c[which(as.numeric(names(data_c))<as.numeric(pericentromere_pops[[p]][c,5]))]
#     data_distal<-c(data_distal, unlist(SVs_proportion[[p]][[c]][names(distal_windows_data)]))
#   }#c
#   distal_means[[p]]<-mean(data_distal, na.rm = TRUE)
# }#p

mean_table<-matrix(nrow = length(window_types), ncol = ncol(pop_table))
row.names(mean_table)<-window_types
for (t in window_types){

}  
graphic_table<-matrix(ncol = 3, nrow=0); colnames(graphic_table)<-c("Window_type", "SVs_proportion", "window")
for (t in window_types){

  for (p in Populations){
graphic_table<-rbind(graphic_table,cbind(p, pop_table_list[[t]][p,4:18], c(-6:8)-1))
}#p
}#t
graphic_table<-as.data.frame(graphic_table)
graphic_table$SVs_proportion<-as.numeric(as.character(graphic_table$SVs_proportion))
graphic_table$window<-as.numeric(as.character(graphic_table$window))


ggplot()+
  geom_line(data = graphic_table, aes(x = window, y = SVs_proportion, color = Window_type)) +
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

#make average among populations
graphic_table<-matrix(ncol = 3, nrow=0); colnames(graphic_table)<-c("window_type", "window", "SVs_proportion")
for (t in names(pop_table_list)){
promedios<-c(); for (j in 1:ncol(pop_table_list[[t]])){promedios<-c(promedios, mean(pop_table_list[[t]][,j]))}  
t_table<-cbind(names(pop_table_list), 1:ncol(pop_table_list[[t]]), promedios)  
}
graphic_table<-as.data.frame(graphic_table)
graphic_table$SVs_proportion<-as.numeric(as.character(graphic_table$SVs_proportion))
graphic_table$window<-as.numeric(as.character(graphic_table$window))
coldspots<-pop_table_list[["coldspots"]][,]