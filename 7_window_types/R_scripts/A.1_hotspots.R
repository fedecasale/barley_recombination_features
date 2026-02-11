make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

#estos son los CO accum per 10K win, como los CO intervals son mas grandes que 10K, 
#un CO esta distribuido en several windows (proportionally al physical length que ocupa)
rec_prob<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.4_breakpoints_per_window/D.4_breakpoints_per_window_10K.RDS")

window_sizes<-names(rec_prob)
#voy a usar las tres layers, so uso la layer 3 (que incluye a las otras dos)
layers<-names(rec_prob[[1]])[3]
#poder haber usado las tres layers con weight distinto para cada una, pero asi es mas facil
Populations<-names(rec_prob[[1]][[1]])

dir.create("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots")

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

saveRDS(accum_rec_prob, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.0_accumulated_rec_prob.RDS")
#accum_rec_prob<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.0_accumulated_rec_prob.RDS")

######### I will normalize accum rec prob here ##############

#I will normalize for each chr in each pop

for (p in Populations){
  for (c in 1:7){
  values<-accum_rec_prob[[p]][[c]][,3]  
  MEAN<-mean(values, na.rm = TRUE); MAX<-max(values, na.rm = TRUE); MIN<-min(values, na.rm = TRUE)
  normalized_values<-unlist(lapply(values, function(x){x<-(x-MEAN)/(MAX-MIN); return(x)}))
  accum_rec_prob[[p]][[c]][,3]<-round(normalized_values, digits = 3)
  }
}

saveRDS(accum_rec_prob, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.0_accumulated_rec_prob_NORMALIZED.RDS")
#accum_rec_prob<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.0_accumulated_rec_prob_NORMALIZED.RDS")
#############################################################

#3-select windows which rec. rate exceeds the genome-wide rec. rate by 2

#3.0-calculate rec. rate per pop (mean COs per win) #chequear, maybe Michael method...
#NO TIENE SENTIDO CREO

#3.1 calculate ave. rec. rate per window per pop
mean_wins<-list()
for (p in Populations){
mean_wins[[p]]<-list()
for (c in 1:7){mean_wins[[p]][[c]]<-mean(accum_rec_prob[[p]][[c]][,3], na.rm = TRUE)}#c    
}#p


#3.2-get windows that exceed the mean
high_rate_win<-list()
hotspot_table<-matrix(ncol = 7, nrow = 3); row.names(hotspot_table)<-Populations
for (p in Populations){
  high_rate_win[[p]]<-list()
  for (c in 1:7){
  high_rate_win[[p]][[c]]<-accum_rec_prob[[p]][[c]][which(accum_rec_prob[[p]][[c]][,3]>mean_wins[[p]][[c]]),]
  hotspot_table[p,c]<-nrow(high_rate_win[[p]][[c]])/nrow(accum_rec_prob[[p]][[c]])
  }#c
}#p
saveRDS(high_rate_win, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.1_windows_exceeding_mean.RDS")

#min(hotspot_table); max(hotspot_table)
#me da que las que exceden son alrededor del 10-15 %

#3.3-check if any high rate window in pericentromeric region
pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere.RDS")

for (p in Populations){
  centromere<-pericentromere_pops[[p]]
  for (c in 1:7){
    #get windows in pericentromere
    peri_start<-as.numeric(centromere[c,4]); peri_end<-as.numeric(centromere[c,5])
    win_in_peri<-as.data.frame(high_rate_win[[p]][[c]][which(as.numeric(high_rate_win[[p]][[c]][,1])>peri_start),])
    win_in_peri<-as.data.frame(win_in_peri[which(as.numeric(win_in_peri[,1])<peri_end),])
    cat(c(nrow(win_in_peri), "-"))
  }
}   

#4.1-select the top 5% of the windows that exceeded the mean
top_rate_win<-list()
hotspot_table<-matrix(ncol = 7, nrow = 3); row.names(hotspot_table)<-Populations
for (p in Populations){
  top_rate_win[[p]]<-list()
  for (c in 1:7){
    #high_rate_win[[p]][[c]]<-high_rate_win[[p]][[c]][order(high_rate_win[[p]][[c]][,3], decreasing = TRUE),]
    #top_rate_win[[p]][[c]]<-high_rate_win[[p]][[c]][1:round((nrow(high_rate_win[[p]][[c]])*0.05)),]
    the_99_percentile<-quantile(high_rate_win[[p]][[c]][,3], probs = seq(0, 1, by = 0.01), na.rm = TRUE)[100]
    top_rate_win[[p]][[c]]<-high_rate_win[[p]][[c]][which(high_rate_win[[p]][[c]][,3]>=the_99_percentile),]
    hotspot_table[p,c]<-nrow(top_rate_win[[p]][[c]])/nrow(accum_rec_prob[[p]][[c]])
  }#c
}#p
saveRDS(top_rate_win, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.2_hotspots_per_pop_WINDOWS.RDS")
#top_rate_win<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.2_hotspots_per_pop_WINDOWS.RDS")

#4.2-check if any top rate window in pericentromeric region
hotspots_peri<-list()
for (p in Populations){
  hotspots_peri[[p]]<-list()
  centromere<-pericentromere_pops[[p]]
  for (c in 1:7){ 
    #get windows in pericentromere
    peri_start<-as.numeric(centromere[c,4]); peri_end<-as.numeric(centromere[c,5])
    win_in_peri<-as.data.frame(top_rate_win[[p]][[c]][which(as.numeric(top_rate_win[[p]][[c]][,1])>peri_start),])
    win_in_peri<-as.data.frame(win_in_peri[which(as.numeric(win_in_peri[,1])<peri_end),])
    hotspots_peri[[p]][[c]]<-win_in_peri
  }
}   
saveRDS(hotspots_peri, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.2_hotspots_per_pop_WINDOWS_pericentromere.RDS")


#3- check if hotspots are regions  - aglutinar todas las windows que esten juntas

hotspots<-list()
for (p in Populations){ 
  hotspots[[p]]<-list()
  for (c in 1:7){
    win_num<-as.numeric(gsub("w_", "",row.names(top_rate_win[[p]][[c]])))
    top_rate_win[[p]][[c]]<-top_rate_win[[p]][[c]][order(win_num),]
    win_num<-win_num[order(win_num)]
    dif<-win_num[2:length(win_num)]-win_num[1:(length(win_num)-1)]
    breaks<-which(dif>1)
    if (length(breaks)!=0){ #1
      hotspots[[p]][[c]]<-list()  
      #el primero
      hotspots[[p]][[c]][[1]]<-top_rate_win[[p]][[c]][1:breaks[1],]
      #los del medio si hay
      if (length(breaks)!=1){ #2
        for (i in 2:(length(breaks)-1)){
          hotspots[[p]][[c]][[length(hotspots[[p]][[c]])+1]]<-top_rate_win[[p]][[c]][(breaks[i]+1):breaks[i+1],]
        }
      }#if2
      #el ultimo
      hotspots[[p]][[c]][[length(hotspots[[p]][[c]])+1]]<-top_rate_win[[p]][[c]][(breaks[length(breaks)]+1):nrow(top_rate_win[[p]][[c]]),]
    }#if1
  }#c
}#p

#

for (p in Populations){
  for (c in 1:7){
    first<-make.matrix(hotspots[[p]][[c]][[1]])
    tabla<-make.matrix(first[1,]); tabla[1,2]<-first[nrow(first),2]
    for (i in 2:length(hotspots[[p]][[c]])){
      following<-make.matrix(hotspots[[p]][[c]][[i]])  
      tabla_i<-make.matrix(following[1,]); tabla_i[1,2]<-following[nrow(following),2]
      tabla<-rbind(tabla, tabla_i)
    }#i
    hotspots[[p]][[c]]<-tabla
    colnames(hotspots[[p]][[c]])[3]<-"length"
    hotspots[[p]][[c]][,3]<-hotspots[[p]][[c]][,2]-hotspots[[p]][[c]][,1]
  }#c
}#p
saveRDS(hotspots, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.3_hotspots_per_pop_REGIONS.RDS")

#5-intersect selected windows among populations

hotspots_intersected<-list()

#intersected among pops, so sharing one parent
combs<-list(c(1,2), c(1,3), c(2,3))
for (q in 1:length(combs)){
int_hotspots<-top_rate_win[[Populations[combs[[q]][1]]]][[1]][which(top_rate_win[[Populations[combs[[q]][1]]]][[1]][,1]%in%top_rate_win[[Populations[combs[[q]][2]]]][[1]][,1]),]
row.names(int_hotspots)<-paste("chr1H_", 1:nrow(int_hotspots), sep = "")  
for (c in 2:7){
chr_table<-top_rate_win[[Populations[combs[[q]][1]]]][[c]][which(top_rate_win[[Populations[combs[[q]][1]]]][[c]][,1]%in%top_rate_win[[Populations[combs[[q]][2]]]][[c]][,1]),]
if (nrow(chr_table)>0){
row.names(chr_table)<-paste("chr", c, "H_", 1:nrow(chr_table), sep = "")  
int_hotspots<-rbind(int_hotspots, chr_table)
}
}#c
hotspots_intersected[[paste(Populations[combs[[q]]], collapse = "_")]]<-int_hotspots
}#q

#intersected among all
int_hotspots<-hotspots_intersected[[1]][which(row.names(hotspots_intersected[[1]])%in%row.names(hotspots_intersected[[2]])),]
int_hotspots<-int_hotspots[which(row.names(int_hotspots)%in%row.names(hotspots_intersected[[3]])),]
hotspots_intersected[[paste(Populations, collapse = "_")]]<-int_hotspots

#save
saveRDS(hotspots_intersected, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.4_hotspots_intersected_WINDOWS.RDS")


#6-get which hotspots were not intersected

hotspots_not_intersected<-list()
for (p in Populations){
  other_pop<-Populations[-which(Populations==p)]
  not_int_hotspots<-top_rate_win[[p]][[1]][-which(top_rate_win[[p]][[1]][,1]%in%c(top_rate_win[[other_pop[1]]][[1]][,1],top_rate_win[[other_pop[2]]][[1]][,1])),]
  row.names(not_int_hotspots)<-paste("chr1H_", 1:nrow(not_int_hotspots), sep = "")  
  for (c in 2:7){
    chr_table<-top_rate_win[[p]][[c]][-which(top_rate_win[[p]][[c]][,1]%in%c(top_rate_win[[other_pop[1]]][[c]][,1],top_rate_win[[other_pop[2]]][[c]][,1])),]
    row.names(chr_table)<-paste("chr", c, "H_", 1:nrow(chr_table), sep = "")  
    not_int_hotspots<-rbind(not_int_hotspots, chr_table)
  }#c
  hotspots_not_intersected[[paste(p)]]<-not_int_hotspots
}#p

#save
saveRDS(hotspots_not_intersected, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.5_hotspots_not_intersected_WINDOWS.RDS")


#7-check which COs are hotspots
breakpoints<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")
hotspots<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.2_hotspots_per_pop_WINDOWS.RDS")
Populations<-names(breakpoints[[3]])
    
COs_not_hotspot_list<-list()
COs_hotspot_list<-list()
for (p in Populations){ cat(p); cat("-")
  COs_not_hotspot_list[[p]]<-list()
  COs_hotspot_list[[p]]<-list()
  for (c in 1:7) { 
    hotspots_c<-hotspots[[p]][[c]]
    COs_c<-breakpoints[[3]][[p]][[c]]
    COs_hotspot<-matrix(ncol = 3, nrow = 0); colnames(COs_hotspot)<-c("starts", "ends", "breakpoint")
    for (r in 1:nrow(hotspots_c)){
    match<-make.matrix(COs_c[which(COs_c[,4]<hotspots_c[r,1]),])
    match<-make.matrix(match[which(match[,5]>hotspots_c[r,2]),])
    COs_hotspot<-rbind(COs_hotspot,make.matrix(match[,c(4:6)]))
    }#r
    COs_not_hotspot<-COs_c[-which(COs_c[,6]%in%COs_hotspot[,3]),4:6]
    COs_hotspot_list[[p]][[c]]<-COs_hotspot
    COs_not_hotspot_list[[p]][[c]]<-COs_not_hotspot
  }#c
}#p

saveRDS(COs_hotspot_list, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.6_COs_in_hotspots.RDS")
saveRDS(COs_not_hotspot_list, "/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.7_COs_not_in_hotspots.RDS")

######################################


