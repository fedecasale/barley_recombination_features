#estos son los CO accum per 10K win, como los CO intervals son mas grandes que 10K, 
#un CO esta distribuido en several windows (proportionally al physical length que ocupa)

rec_prob<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.4_breakpoints_per_window/D.4_breakpoints_per_window_10K.RDS")

window_sizes<-names(rec_prob)
#voy a usar las tres layers, so uso la layer 3 (que incluye a las otras dos
layers<-names(rec_prob[[1]])[1]
#poder haber usado las tres layers con weight distinto para cada una, pero asi es mas facil
Populations<-names(rec_prob[[1]][[1]])

dir.create("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots")

#1-distribute CO interval in window proportion in each ril #ESTO YA LO HICE EN D.4

#2-sum al windows proportions of rils per window per population
#ya lo hice para window sizes mas grandes, por eso me copio los algoritms del codigo D.5.1 

accumulated_rec_prob<-list()
for (w in window_sizes){
  accumulated_rec_prob[[w]]<-list()
  for (l in layers){ cat(l); cat(" layer"); cat(" - ")
    accumulated_rec_prob[[w]][[l]]<-list()
    for (p in Populations){ cat(p); cat("-")
      rils<-names(rec_prob[[w]][[l]][[p]])
      accumulated_rec_prob[[w]][[l]][[p]]<-list()
      for (c in 1:7) {
        #accumulate COs
        rec_prob_c<-rec_prob[[w]][[l]][[p]][[1]][[c]]; rec_prob_c[,3]<-0
        for (r in rils){rec_prob_c[,3]<-as.numeric(rec_prob_c[,3])+as.numeric(rec_prob[[w]][[l]][[p]][[r]][[c]][,3])}#r
        #generate probability
        rec_prob_c[,3]<-as.numeric(rec_prob_c[,3])/sum(as.numeric(rec_prob_c[,3])) #sum(rec_prob_c[,3])==1
        #save
        accumulated_rec_prob[[w]][[l]][[p]][[c]]<-rec_prob_c
      }#c
    }#p
  }#l
}#w

#simplify list since there is only one window size and layer
accum_rec_prob<-list()
for (w in window_sizes){
  for (l in layers){ cat(l); cat(" layer"); cat(" - ")
    for (p in Populations){ cat(p); cat("-")
      accum_rec_prob[[p]]<-list()
      for (c in 1:7) {accum_rec_prob[[p]][[c]]<-accumulated_rec_prob[[w]][[l]][[p]][[c]]}#c
    }#p
  }#l
}#w
rm(accumulated_rec_prob)

saveRDS(accum_rec_prob, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.0_accumulated_rec_prob_FIRST_LAYER.RDS")

#3-select windows which rec. rate exceeds the genome-wide rec. rate by 2

#3.0-calculate rec. rate per pop (mean COs per win)
#chequear, maybe Michael method...

#3.1 calculate ave. rec. rate per window per pop
mean_wins<-list()
for (p in Populations){
  sum_win<-c()
  for (c in 1:7){sum_win<-c(sum_win, accum_rec_prob[[p]][[c]][,3])}#c    
  mean_wins[[p]]<-mean(sum_win)
}#p

#3.2-get windows that exceed the mean
high_rate_win<-list()
hotspot_table<-matrix(ncol = 7, nrow = 3); row.names(hotspot_table)<-Populations
for (p in Populations){
  high_rate_win[[p]]<-list()
  for (c in 1:7){
    high_rate_win[[p]][[c]]<-accum_rec_prob[[p]][[c]][which(accum_rec_prob[[p]][[c]][,3]>mean_wins[[p]]),]
    hotspot_table[p,c]<-nrow(high_rate_win[[p]][[c]])/nrow(accum_rec_prob[[p]][[c]])
  }#c
}#p
saveRDS(high_rate_win, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.1_windows_exceeding_mean_FIRST_LAYER.RDS")

min(hotspot_table); max(hotspot_table)
#me da que las que exceden son alrededor del 10-15 %

#3.3-check if any high rate window in pericentromeric region
pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere_FIRST_LAYER.RDS")

for (p in Populations){
  centromere<-pericentromere_pops[[p]]
  for (c in 1:7){
    #get windows in pericentromere
    peri_start<-as.numeric(centromere[c,4]); peri_end<-as.numeric(centromere[c,5])
    win_in_peri<-as.data.frame(high_rate_win[[p]][[c]][which(as.numeric(high_rate_win[[p]][[c]][,1])>peri_start),])
    win_in_peri<-as.data.frame(win_in_peri[[p]][[c]][which(as.numeric(win_in_peri[[p]][[c]][,1])<peri_end),])
    cat(c(nrow(win_in_peri), "-"))
  }
}   
#no hay in pericentromeric region

#4-select the top 5% of the windows that exceeded the mean
top_rate_win<-list()
hotspot_table<-matrix(ncol = 7, nrow = 3); row.names(hotspot_table)<-Populations
for (p in Populations){
  top_rate_win[[p]]<-list()
  for (c in 1:7){
    high_rate_win[[p]][[c]]<-high_rate_win[[p]][[c]][order(high_rate_win[[p]][[c]][,3], decreasing = TRUE),]
    top_rate_win[[p]][[c]]<-high_rate_win[[p]][[c]][1:round((nrow(high_rate_win[[p]][[c]])*0.05)),]
    hotspot_table[p,c]<-nrow(top_rate_win[[p]][[c]])/nrow(accum_rec_prob[[p]][[c]])
  }#c
}#p
min(hotspot_table); max(hotspot_table)

saveRDS(top_rate_win, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.1_hotspots_per_pop_FIRST_LAYER.RDS")
hotspot_table[]<-round(as.numeric(hotspot_table[]), digits = 5)
write.csv(hotspot_table, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.2_hotspots_per_pop_stats_FIRST_LAYER.csv")

# #5-intersect selected windows among populations
# hotspots_win<-list()
# hotspot_table<-matrix(ncol = 7, nrow = 2); row.names(hotspot_table)<-c("N wins", "N wins/total wins")
# for (c in 1:7){
#   hotspots_win[[c]]<-top_rate_win[[Populations[1]]][[c]]
#   for (p in Populations[2:length(Populations)]){
#     hotspots_win[[c]]<-hotspots_win[[c]][which(hotspots_win[[c]][,1]%in%top_rate_win[[p]][[c]][,1]),]
#   }#p
#   hotspot_table[1,c]<-nrow(hotspots_win[[c]])
#   hotspot_table[2,c]<-nrow(hotspots_win[[c]])/nrow(accum_rec_prob[[p]][[c]])
# }#c
# min(hotspot_table[2,]); max(hotspot_table[2,])
# 
# saveRDS(hotspots_win, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.3_hotspots_intersected_FIRST_LAYER.RDS")
# 
# hotspot_table<-format(hotspot_table, scientific = FALSE)
# hotspot_table<-cbind(hotspot_table, rbind(sum(as.numeric(hotspot_table[1,])),mean(as.numeric(hotspot_table[2,]))))
# colnames(hotspot_table)<-c(paste(1:7,"H", sep = ""), "genome-wide")
# hotspot_table[1,]<-round(as.numeric(hotspot_table[1,]))
# hotspot_table[2,]<-round(as.numeric(hotspot_table[2,]), digits = 5)
# write.csv(hotspot_table, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.4_hotspots_intersected_stats_FIRST_LAYER.csv")
# 
# #6-get which hotspots were not intersected
# not_intersected<-list()
# hotspot_table<-matrix(ncol = 7, nrow = 3); row.names(hotspot_table)<-Populations
# for (p in Populations){
#   not_intersected[[p]]<-list()
#   for (c in 1:7){
#     shared<-which(row.names(top_rate_win[[p]][[c]])%in%row.names(hotspots_win[[c]]))
#     not_intersected[[p]][[c]]<-top_rate_win[[p]][[c]][-shared,]
#     hotspot_table[p,c]<-nrow(not_intersected[[p]][[c]])/nrow(accum_rec_prob[[p]][[c]])
#   }#c
# }#p
# 
# #
# saveRDS(not_intersected, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.5_hotspots_per_pop_not_intersected_FIRST_LAYER.RDS")
# hotspot_table[]<-round(as.numeric(hotspot_table[]), digits = 5)
# write.csv(hotspot_table, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.6_hotspots_per_pop_stats_not_intersected_FIRST_LAYER.csv")
# #
# 
# # to check #
# accumulated_rec_prob<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.3_CO_prob_per_window_GENETIC_MAPS.RDS", sep = ""))
# 
# 
# ##### check which COs are hotspots
# breakpoints<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")
# hotspots<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.1_hotspots_per_pop_FIRST_LAYER.RDS")
# Populations<-names(breakpoints[[3]])
# 
# make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
#   new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
#   new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}
# 
# COs_not_hotspot_list<-list()
# COs_hotspot_list<-list()
# for (p in Populations){ cat(p); cat("-")
#   COs_not_hotspot_list[[p]]<-list()
#   COs_hotspot_list[[p]]<-list()
#   for (c in 1:7) { 
#     hotspots_c<-hotspots[[p]][[c]]
#     COs_c<-breakpoints[[3]][[p]][[c]]
#     COs_hotspot<-matrix(ncol = 3, nrow = 0); colnames(COs_hotspot)<-c("starts", "ends", "breakpoint")
#     for (r in 1:nrow(hotspots_c)){
#       match<-make.matrix(COs_c[which(COs_c[,4]<hotspots_c[r,1]),])
#       match<-make.matrix(match[which(match[,5]>hotspots_c[r,2]),])
#       COs_hotspot<-rbind(COs_hotspot,make.matrix(match[,c(4:6)]))
#     }#r
#     COs_not_hotspot<-COs_c[-which(COs_c[,6]%in%COs_hotspot[,3]),4:6]
#     COs_hotspot_list[[p]][[c]]<-COs_hotspot
#     COs_not_hotspot_list[[p]][[c]]<-COs_not_hotspot
#   }#c
# }#p
# 
# saveRDS(COs_hotspot_list, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.7_COs_in_hotspots_FIRST_LAYER.RDS")
# saveRDS(COs_not_hotspot_list, "/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.8_COs_not_in_hotspots_FIRST_LAYER.RDS")
