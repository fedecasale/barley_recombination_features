make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

library(ggplot2)
library(impressionist.colors)
library(xlsx)
#####################################################################################################################
dir.create("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots")

#voy a usar las tres layers, so uso la layer 3 (que incluye a las otras dos)
accum_rec_prob<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.0_accumulated_rec_prob.RDS")

window_sizes<-10000

#1-recalcular la peri para las pops...una por una.
#pericentromere_ave<-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede.csv")
pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere.RDS")
Populations<-names(pericentromere_pops)

#2.A-sacar los de la peri -> get distal region 
#2.B-sacar los de los extremos (pegados al peri o pegados al telomere) -> cut distal region

accum_rec_prob_DISTAL<-list()
accum_rec_prob_PERI<-list()
for (p in Populations){
  accum_rec_prob_DISTAL[[p]]<-list()
  accum_rec_prob_PERI[[p]]<-list()
  centromere<-pericentromere_pops[[p]]
  for (c in 1:7){
    #get distal regions
    peri_start<-as.numeric(centromere[c,4]); peri_end<-as.numeric(centromere[c,5])
    arm_1<-as.data.frame(accum_rec_prob[[p]][[c]][which(as.numeric(accum_rec_prob[[p]][[c]][,2])<peri_start),])
    arm_2<-as.data.frame(accum_rec_prob[[p]][[c]][which(as.numeric(accum_rec_prob[[p]][[c]][,1])>peri_end),])
    #cut distal regions
    arm_1_limits<-quantile(arm_1[,1], prob = seq(0,1,by = 0.025))[c(6,36)]  
    arm_1<-arm_1[which(arm_1[,1]>arm_1_limits[1]),]
    arm_1<-arm_1[which(arm_1[,1]<arm_1_limits[2]),]
    arm_2_limits<-quantile(arm_2[,1], prob = seq(0,1,by = 0.025))[c(6,36)]  
    arm_2<-arm_2[which(arm_2[,1]>arm_2_limits[1]),]
    arm_2<-arm_2[which(arm_2[,1]<arm_2_limits[2]),]
    #save
    accum_rec_prob_DISTAL[[p]][[c]]<-rbind(arm_1, arm_2)
    #get also peri windows
    peri<-as.data.frame(accum_rec_prob[[p]][[c]][which(as.numeric(accum_rec_prob[[p]][[c]][,1])>peri_start),])
    peri<-as.data.frame(peri[which(as.numeric(peri[,2])<peri_end),])
    accum_rec_prob_PERI[[p]][[c]]<-peri
  }#c
}#p

saveRDS(accum_rec_prob_DISTAL, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.0_accum_rec_prob_DISTAL.RDS")
saveRDS(accum_rec_prob_PERI, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.0_accum_rec_prob_PERI.RDS")


#3.A-tengo que quedarme con todos los 0
low_rate_win<-list()
coldspot_table<-matrix(ncol = 7, nrow = 3); row.names(coldspot_table)<-Populations
for (p in Populations){
  low_rate_win[[p]]<-list()
  for (c in 1:7){
    low_rate_win[[p]][[c]]<-accum_rec_prob_DISTAL[[p]][[c]][which(accum_rec_prob_DISTAL[[p]][[c]][,3]==0),]
    coldspot_table[p,c]<-nrow(low_rate_win[[p]][[c]])
  }#c
}#p
saveRDS(low_rate_win, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.1_coldspots_per_pop_WINDOWS.RDS")


#3.B-aglutinar todas las windows que esten juntas
coldspots<-list()
#coldspot_table<-matrix(ncol = 7, nrow = 3); row.names(coldspot_table)<-Populations
for (p in Populations){ 
  coldspots[[p]]<-list()
  for (c in 1:7){
    win_num<-as.numeric(gsub("w_", "",row.names(low_rate_win[[p]][[c]])))
    dif<-win_num[2:length(win_num)]-win_num[1:(length(win_num)-1)]
    breaks<-which(dif>1)
    if (length(breaks)!=0){ #1
    coldspots[[p]][[c]]<-list()  
    #el primero
    coldspots[[p]][[c]][[1]]<-low_rate_win[[p]][[c]][1:breaks[1],]
    #los del medio si hay
    if (length(breaks)!=1){ #2
      for (i in 2:(length(breaks)-1)){
        coldspots[[p]][[c]][[length(coldspots[[p]][[c]])+1]]<-low_rate_win[[p]][[c]][(breaks[i]+1):breaks[i+1],]
      }
    }#if2
    #el ultimo
    coldspots[[p]][[c]][[length(coldspots[[p]][[c]])+1]]<-low_rate_win[[p]][[c]][(breaks[length(breaks)]+1):nrow(low_rate_win[[p]][[c]]),]
    }#if1
    #coldspot_table[p,c]<-length(coldspots[[p]][[c]])
  }#c
}#p
#
for (p in Populations){
  for (c in 1:7){
    tabla<-make.matrix(as.matrix(coldspots[[p]][[c]][[1]][1,])); tabla[1,2]<-coldspots[[p]][[c]][[1]][nrow(coldspots[[p]][[c]][[1]]),2]
    for (i in 2:length(coldspots[[p]][[c]])){
    tabla<-rbind(tabla, as.matrix(coldspots[[p]][[c]][[i]][1,])); tabla[i,2]<-coldspots[[p]][[c]][[i]][nrow(coldspots[[p]][[c]][[i]]),2]
    }#i
  coldspots[[p]][[c]]<-tabla
  colnames(coldspots[[p]][[c]])[3]<-"length"
  coldspots[[p]][[c]][,3]<-coldspots[[p]][[c]][,2]-coldspots[[p]][[c]][,1]
  }#c
}#p

saveRDS(coldspots, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.2_coldspots_per_pop_REGION.RDS")




#4-clasify coldspots in size
file.remove("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.3_coldspots_by_size_stats.xlsx")
chrs<-paste("chr",1:7, "H", sep = "")
classification<-c("<=10 kb", "<=500 kb", "<=1 Mbp", "<=5 Mbp", "<=10 Mbp", ">10 Mbp", "mean", "min", "max")
clasiff_tabla<-matrix(ncol = 7, nrow = 9); colnames(clasiff_tabla)<-chrs; row.names(clasiff_tabla)<-classification
coldspots_size<-list()

for (p in Populations){
  
  coldspots_size[[p]]<-list()
  c_tabla_p<-clasiff_tabla
  
  for (c in 1:7){
  
  cp<-as.data.frame(coldspots[[p]][[c]])
  
  #stats
  c_tabla_p[7,c]<-round(mean(cp[,3]))
  c_tabla_p[8,c]<-round(min(cp[,3]))+1
  c_tabla_p[9,c]<-round(max(cp[,3]))+1
  
  #get coldspots by sizes
  coldspots_size[[p]][[c]]<-list()
  
  coldspots_size[[p]][[c]][[classification[1]]]<-cp[which(cp[,3]<=10000),]
  c_tabla_p[1,c]<-length(which(cp[,3]<=10000))
  cp<-cp[-which(cp[,3]<=10000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
  
  coldspots_size[[p]][[c]][[classification[2]]]<-cp[which(cp[,3]<=500000),]
  c_tabla_p[2,c]<-length(which(cp[,3]<=500000))
  cp<-cp[-which(cp[,3]<=500000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp)
  
  coldspots_size[[p]][[c]][[classification[3]]]<-cp[which(cp[,3]<=1000000),]
  c_tabla_p[3,c]<-length(which(cp[,3]<=1000000))
  cp<-cp[-which(cp[,3]<=1000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
  
  coldspots_size[[p]][[c]][[classification[4]]]<-cp[which(cp[,3]<=5000000),]
  c_tabla_p[4,c]<-length(which(cp[,3]<=5000000))
  cp<-cp[-which(cp[,3]<=5000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
  
  if (nrow(cp)!=0){
  coldspots_size[[p]][[c]][[classification[5]]]<-cp[which(cp[,3]<=10000000),]
  c_tabla_p[5,c]<-length(which(cp[,3]<=10000000))
  cp<-cp[-which(cp[,3]<=10000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
  }
  
  if (nrow(cp)!=0){
  coldspots_size[[p]][[c]][[classification[6]]]<-cp[which(cp[,3]>10000000),]
  c_tabla_p[6,c]<-length(which(cp[,3]>10000000))
  }

  }#c
  c_tabla_p[is.na(c_tabla_p)]<-0
  
#save
write.xlsx(c_tabla_p, 
           "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.3_coldspots_by_size_stats.xlsx",
           append = TRUE, sheetName = p)
}#p

saveRDS(coldspots_size, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.4_coldspots_per_pop_per_size.RDS")

#################################################################################

coldspots_intersected<-list()

#intersected among pops, so sharing one parent
combs<-list(c(1,2), c(1,3), c(2,3))
for (q in 1:length(combs)){
  int_coldspots<-low_rate_win[[Populations[combs[[q]][1]]]][[1]][which(low_rate_win[[Populations[combs[[q]][1]]]][[1]][,1]%in%low_rate_win[[Populations[combs[[q]][2]]]][[1]][,1]),]
  row.names(int_coldspots)<-paste("chr1H_", 1:nrow(int_coldspots), sep = "")  
  for (c in 2:7){
    chr_table<-low_rate_win[[Populations[combs[[q]][1]]]][[c]][which(low_rate_win[[Populations[combs[[q]][1]]]][[c]][,1]%in%low_rate_win[[Populations[combs[[q]][2]]]][[c]][,1]),]
    row.names(chr_table)<-paste("chr", c, "H_", 1:nrow(chr_table), sep = "")  
    int_coldspots<-rbind(int_coldspots, chr_table)
  }#c
  coldspots_intersected[[paste(Populations[combs[[q]]], collapse = "_")]]<-int_coldspots
}#q

#intersected among all
int_coldspots<-coldspots_intersected[[1]][which(row.names(coldspots_intersected[[1]])%in%row.names(coldspots_intersected[[2]])),]
int_coldspots<-int_coldspots[which(row.names(int_coldspots)%in%row.names(coldspots_intersected[[3]])),]
coldspots_intersected[[paste(Populations, collapse = "_")]]<-int_coldspots

#save
saveRDS(coldspots_intersected, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.5_coldspots_intersected_WINDOWS.RDS")


#6-get which coldspots were not intersected

coldspots_not_intersected<-list()
for (p in Populations){
  other_pop<-Populations[-which(Populations==p)]
  not_int_coldspots<-low_rate_win[[p]][[1]][-which(low_rate_win[[p]][[1]][,1]%in%c(low_rate_win[[other_pop[1]]][[1]][,1],low_rate_win[[other_pop[2]]][[1]][,1])),]
  row.names(not_int_coldspots)<-paste("chr1H_", 1:nrow(not_int_coldspots), sep = "")  
  for (c in 2:7){
    chr_table<-low_rate_win[[p]][[c]][-which(low_rate_win[[p]][[c]][,1]%in%c(low_rate_win[[other_pop[1]]][[c]][,1],low_rate_win[[other_pop[2]]][[c]][,1])),]
    row.names(chr_table)<-paste("chr", c, "H_", 1:nrow(chr_table), sep = "")  
    not_int_coldspots<-rbind(not_int_coldspots, chr_table)
  }#c
  coldspots_not_intersected[[paste(p)]]<-not_int_coldspots
}#p

#save
saveRDS(coldspots_not_intersected, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.6_coldspots_not_intersected_WINDOWS.RDS")


#COSAS VIEJAS
# #5-intersect coldspots among populations
# 
# #5.1-intersect windows
# coldspots_win<-list()
# coldspot_table<-matrix(ncol = 7, nrow = 2); row.names(coldspot_table)<-c("N wins", "N wins/total wins")
# for (c in 1:7){
#   coldspots_win[[c]]<-low_rate_win[[Populations[1]]][[c]]
#   for (p in Populations[2:length(Populations)]){
#     coldspots_win[[c]]<-coldspots_win[[c]][which(coldspots_win[[c]][,1]%in%low_rate_win[[p]][[c]][,1]),]
#   }#p
#   coldspot_table[1,c]<-nrow(coldspots_win[[c]])
#   coldspot_table[2,c]<-nrow(coldspots_win[[c]])/nrow(accum_rec_prob[[p]][[c]])
# }#c
# saveRDS(coldspots_win, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.5_coldspots_intersected_windows.RDS")
# #5.2-aglutinate windows
# coldspots_int<-list()
# for (c in 1:7){  
#   win_num<-as.numeric(gsub("w_", "",row.names(coldspots_win[[c]])))
#   dif<-win_num[2:length(win_num)]-win_num[1:(length(win_num)-1)]
#   breaks<-which(dif>1)
#   if (length(breaks)!=0){ #1
#     coldspots_int[[c]]<-list()  
#     #el primero
#     coldspots_int[[c]][[1]]<-coldspots_win[[c]][1:breaks[1],]
#     #los del medio si hay
#     if (length(breaks)!=1){ #2
#       for (i in 2:(length(breaks)-1)){
#         coldspots_int[[c]][[length(coldspots_int[[c]])+1]]<-coldspots_win[[c]][(breaks[i]+1):breaks[i+1],]
#       }
#     }#if2
#     #el ultimo
#     coldspots_int[[c]][[length(coldspots_int[[c]])+1]]<-coldspots_win[[c]][(breaks[length(breaks)]+1):nrow(coldspots_win[[c]]),]
#   }#if1
# }#c
# #
# coldspot_table<-rbind(coldspot_table, NA); row.names(coldspot_table)[3]<-"N coldspots"
# for (c in 1:7){
#     tabla<-make.matrix(as.matrix(coldspots_int[[c]][[1]][1,])); tabla[1,2]<-coldspots_int[[c]][[1]][nrow(coldspots_int[[c]][[1]]),2]
#     for (i in 2:length(coldspots_int[[c]])){
#       tabla<-rbind(tabla, as.matrix(coldspots_int[[c]][[i]][1,])); tabla[i,2]<-coldspots_int[[c]][[i]][nrow(coldspots_int[[c]][[i]]),2]
#     }#i
#     coldspots_int[[c]]<-tabla
#     colnames(coldspots_int[[c]])[3]<-"length"
#     coldspots_int[[c]][,3]<-coldspots_int[[c]][,2]-coldspots_int[[c]][,1]
#     #add stats to table  
#     coldspot_table[3,c]<-nrow(coldspots_int[[c]])
#   }#c
# #
# saveRDS(coldspots_int, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.5_coldspots_intersected.RDS")
# chrs<-paste("chr",1:7, "H", sep = "")
# colnames(coldspot_table)<-chrs
# coldspot_table[]<-round(coldspot_table[], digits = 3)
# write.csv(coldspot_table, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.6_coldspots_table_intersected.csv")
# 
# #5.3-clasify intersected coldspots in size
# classification<-c("<=10 kb", "<=500 kb", "<=1 Mbp", "<=5 Mbp", "<=10 Mbp", ">10 Mbp", "mean", "min", "max")
# c_tabla_p<-matrix(ncol = 7, nrow = 9); colnames(c_tabla_p)<-chrs; row.names(c_tabla_p)<-classification
# coldspots_size<-list()
# 
#   coldspots_size<-list()
#   
#   for (c in 1:7){
#     
#     cp<-as.data.frame(coldspots_int[[c]])
#     
#     #stats
#     c_tabla_p[7,c]<-round(mean(cp[,3]))
#     c_tabla_p[8,c]<-round(min(cp[,3]))+1
#     c_tabla_p[9,c]<-round(max(cp[,3]))+1
#     
#     #get coldspots by sizes
#     coldspots_size[[c]]<-list()
#     
#     coldspots_size[[c]][[classification[1]]]<-cp[which(cp[,3]<=10000),]
#     c_tabla_p[1,c]<-length(which(cp[,3]<=10000))
#     cp<-cp[-which(cp[,3]<=10000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
#     
#     coldspots_size[[c]][[classification[2]]]<-cp[which(cp[,3]<=500000),]
#     c_tabla_p[2,c]<-length(which(cp[,3]<=500000))
#     cp<-cp[-which(cp[,3]<=500000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp)
#     
#     coldspots_size[[c]][[classification[3]]]<-cp[which(cp[,3]<=1000000),]
#     c_tabla_p[3,c]<-length(which(cp[,3]<=1000000))
#     cp<-cp[-which(cp[,3]<=1000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
#     
#     coldspots_size[[c]][[classification[4]]]<-cp[which(cp[,3]<=5000000),]
#     c_tabla_p[4,c]<-length(which(cp[,3]<=5000000))
#     cp<-cp[-which(cp[,3]<=5000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
#     
#     if (nrow(cp)!=0){
#       coldspots_size[[c]][[classification[5]]]<-cp[which(cp[,3]<=10000000),]
#       c_tabla_p[5,c]<-length(which(cp[,3]<=10000000))
#       cp<-cp[-which(cp[,3]<=10000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
#     }
#     
#     if (nrow(cp)!=0){
#       coldspots_size[[c]][[classification[6]]]<-cp[which(cp[,3]>10000000),]
#       c_tabla_p[6,c]<-length(which(cp[,3]>10000000))
#     }
#     
#   }#c
#   
# #save
# c_tabla_p[is.na(c_tabla_p)]<-0
# write.csv(c_tabla_p, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.8_coldspots_by_size_stats_intersected.csv")
# saveRDS(coldspots_size, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.7_coldspots_per_size_intersected.RDS")
# 
# 
# ###############
# 
# #6-take out intersected regions from original coldspots, to get not intersected region coldspots
# 
# #6.1-take out intersected windows from total null rec windows - LO DEBO HACER EN A WINDOW BASIS 
# for (p in Populations){
#   for (c in 1:7){
#   shared<-which(row.names(low_rate_win[[p]][[c]])%in%row.names(coldspots_win[[c]]))
#   low_rate_win[[p]][[c]]<-low_rate_win[[p]][[c]][-shared,]
#   }#c
# }#p
# saveRDS(low_rate_win, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.9_coldspots_per_pop_no_intersected_windows.RDS")
# 
# #6.2-aglutinar todas las windows que esten juntas
# coldspots<-list()
# coldspot_table<-matrix(ncol = 7, nrow = 3); row.names(coldspot_table)<-Populations
# for (p in Populations){ 
#   coldspots[[p]]<-list()
#   for (c in 1:7){
#     win_num<-as.numeric(gsub("w_", "",row.names(low_rate_win[[p]][[c]])))
#     dif<-win_num[2:length(win_num)]-win_num[1:(length(win_num)-1)]
#     breaks<-which(dif>1)
#     if (length(breaks)!=0){ #1
#       coldspots[[p]][[c]]<-list()  
#       #el primero
#       coldspots[[p]][[c]][[1]]<-low_rate_win[[p]][[c]][1:breaks[1],]
#       #los del medio si hay
#       if (length(breaks)!=1){ #2
#         for (i in 2:(length(breaks)-1)){
#           coldspots[[p]][[c]][[length(coldspots[[p]][[c]])+1]]<-low_rate_win[[p]][[c]][(breaks[i]+1):breaks[i+1],]
#         }
#       }#if2
#       #el ultimo
#       coldspots[[p]][[c]][[length(coldspots[[p]][[c]])+1]]<-low_rate_win[[p]][[c]][(breaks[length(breaks)]+1):nrow(low_rate_win[[p]][[c]]),]
#     }#if1
#     coldspot_table[p,c]<-length(coldspots[[p]][[c]])
#   }#c
# }#p
# #
# for (p in Populations){
#   for (c in 1:7){
#     tabla<-make.matrix(as.matrix(coldspots[[p]][[c]][[1]][1,])); tabla[1,2]<-coldspots[[p]][[c]][[1]][nrow(coldspots[[p]][[c]][[1]]),2]
#     for (i in 2:length(coldspots[[p]][[c]])){
#       tabla<-rbind(tabla, as.matrix(coldspots[[p]][[c]][[i]][1,])); tabla[i,2]<-coldspots[[p]][[c]][[i]][nrow(coldspots[[p]][[c]][[i]]),2]
#     }#i
#     coldspots[[p]][[c]]<-tabla
#     colnames(coldspots[[p]][[c]])[3]<-"length"
#     coldspots[[p]][[c]][,3]<-coldspots[[p]][[c]][,2]-coldspots[[p]][[c]][,1]
#   }#c
# }#p
# 
# saveRDS(coldspots, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.10_coldspots_per_pop_no_intersected.RDS")
# 
# 
# #6.3-clasify coldspots in size
# file.remove("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.12_coldspots_by_size_stats_no_intersected.xlsx")
# coldspots_size<-list()
# 
# for (p in Populations){
#   
#   coldspots_size[[p]]<-list()
#   c_tabla_p<-clasiff_tabla
#   
#   for (c in 1:7){
#     
#     cp<-as.data.frame(coldspots[[p]][[c]])
#     
#     #stats
#     c_tabla_p[7,c]<-round(mean(cp[,3]))
#     c_tabla_p[8,c]<-round(min(cp[,3]))+1
#     c_tabla_p[9,c]<-round(max(cp[,3]))+1
#     
#     #get coldspots by sizes
#     coldspots_size[[p]][[c]]<-list()
#     
#     coldspots_size[[p]][[c]][[classification[1]]]<-cp[which(cp[,3]<=10000),]
#     c_tabla_p[1,c]<-length(which(cp[,3]<=10000))
#     cp<-cp[-which(cp[,3]<=10000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
#     
#     coldspots_size[[p]][[c]][[classification[2]]]<-cp[which(cp[,3]<=500000),]
#     c_tabla_p[2,c]<-length(which(cp[,3]<=500000))
#     cp<-cp[-which(cp[,3]<=500000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp)
#     
#     coldspots_size[[p]][[c]][[classification[3]]]<-cp[which(cp[,3]<=1000000),]
#     c_tabla_p[3,c]<-length(which(cp[,3]<=1000000))
#     cp<-cp[-which(cp[,3]<=1000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
#     
#     coldspots_size[[p]][[c]][[classification[4]]]<-cp[which(cp[,3]<=5000000),]
#     c_tabla_p[4,c]<-length(which(cp[,3]<=5000000))
#     cp<-cp[-which(cp[,3]<=5000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
#     
#     if (nrow(cp)!=0){
#       coldspots_size[[p]][[c]][[classification[5]]]<-cp[which(cp[,3]<=10000000),]
#       c_tabla_p[5,c]<-length(which(cp[,3]<=10000000))
#       cp<-cp[-which(cp[,3]<=10000000),]; cp<-make.matrix(cp); cp<-as.data.frame(cp) 
#     }
#     
#     if (nrow(cp)!=0){
#       coldspots_size[[p]][[c]][[classification[6]]]<-cp[which(cp[,3]>10000000),]
#       c_tabla_p[6,c]<-length(which(cp[,3]>10000000))
#     }
#     
#   }#c
#   c_tabla_p[is.na(c_tabla_p)]<-0
#   
#   #save
#   write.xlsx(c_tabla_p, 
#              "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.12_coldspots_by_size_stats_no_intersected.xlsx",
#              append = TRUE, sheetName = p)
# }#p
# 
# #save
# saveRDS(coldspots_size, "/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.11_coldspots_per_pop_per_size_no_intersected.RDS")
# ##########





########################## graphic check ##########################

# rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_1e+06.RDS", sep = ""))
# graphic_table2<-matrix(ncol = 5, nrow = 0)
# colnames(graphic_table2)<-c("Population","chr","Mbp","cM_Mbp","region")
# for (p in Populations){for (c in 1:7){graphic_table2<-rbind(graphic_table2, cbind(p,rec_list[[p]][[c]]))}}
# graphic_table2<-as.data.frame(graphic_table2)
# graphic_table2$Mbp=as.numeric(levels(graphic_table2$Mbp))[graphic_table2$Mbp]
# graphic_table2$Mbp<-graphic_table2$Mbp/1000000
# graphic_table2$cM_Mbp=as.numeric(levels(graphic_table2$cM_Mbp))[graphic_table2$cM_Mbp]
# #graphic_table2$fold_change=as.numeric(levels(graphic_table2$fold_change))[graphic_table2$fold_change]
# 
# COLORES<-as.character(graphic_table2$region)
# COLORES[which(COLORES=="D")]<-"red"
# COLORES[which(COLORES%in%c("P","P2"))]<-"blue"
# COLORES[which(COLORES=="I")]<-"green"
# 
# m.f<-0.1
# ancho<-1300*m.f
# alto<-100*m.f
# s<-1000000
# centromere <-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede.csv")
# color_table<-as.matrix(get.color(8,3,c(5,7,12)))
# centro_col<-get.color(8,3,14); peri_col<-get.color(8,3,13)
# centromere[,4:6]<-centromere[,4:6]/s
# 
# centro_peri<-list()
# centro_peri[[1]]<-geom_vline(data = centromere, aes(xintercept = centromere), color = centro_col, size = 0.75)
# centro_peri[[2]]<-geom_vline(data = centromere, aes(xintercept = peri_start), linetype = 2, color = peri_col, size = 0.5)
# centro_peri[[3]]<-geom_vline(data = centromere, aes(xintercept = peri_end), linetype = 2, color = peri_col, size = 0.5)
# 
# facet<-facet_grid(cols = vars(Population), rows = vars(chr), scales = "free", space = "free_y", switch = "y") 
# 
# thema<-theme(strip.background = element_blank(), 
#              panel.background  = element_rect(fill=NA, color = "black", size = 1, linetype = "solid"),
#              panel.grid = element_blank(),
#              panel.spacing = unit(3*m.f, "cm"),
#              #panel.background  = element_rect(fill=NA, size = 3*m.f, linetype = "solid"),
#              strip.text = element_text(size = 150*m.f, margin = margin(50*m.f, 50*m.f, 50*m.f, 50*m.f))) 
# 
# graphic_table3<-matrix(ncol = 4, nrow = 0); colnames(graphic_table3)<-c("Population","chr","start","end")
# for (p in Populations){
#   for (c in 1:7){
#   chr_table<-matrix(ncol = 4, nrow = nrow(coldspots[[p]][[c]]))
#   chr_table[,1]<-p; chr_table[,2]<-paste(c,"H",sep = "") 
#   chr_table[,3:4]<-coldspots[[p]][[c]][,1:2]
#   graphic_table3<-rbind(graphic_table3, chr_table)
#   }
# }
# graphic_table3<-as.data.frame(graphic_table3)
# graphic_table3$start<-as.numeric(as.character(graphic_table3$start))/s
# graphic_table3$end<-as.numeric(as.character(graphic_table3$end))/s
# 
# 
# #pdf(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.2_rec_rates_windows_regions_",s,".pdf", sep = ""), width=ancho, height=alto)
#   ggplot() + 
#     geom_point(data=graphic_table2 ,stat = "identity", aes(x=Mbp, y=cM_Mbp), color = COLORES) +
#     geom_segment(data=graphic_table3, aes(x = start, xend = end, y=2, yend=2), color = "orange")+
#     scale_x_continuous(name = "physical distance (Mbp)", position = "bottom") +
#     scale_y_continuous(name = "recombination rate (cM/Mbp)", position = "right") +    
#     theme(legend.position = "none",
#           axis.text.x.bottom = element_text(size = 110*m.f),
#           axis.text.y.right  = element_text(size = 110*m.f),
#           axis.ticks.length = unit(1*m.f, "cm")
#     ) + 
#     centro_peri + facet + thema
#     
# #dev.off()
###################################################################
