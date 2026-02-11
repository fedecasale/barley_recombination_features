library(data.table)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[,1]

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

#parents.code<-as.matrix(read.table("/scratch/Federico/Paper_2/samplenames.txt"))

#met_context<-c("chg", "cpg", "chh")

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
if (win.size == 10000){Populations<-Populations[which(Populations%in%c("HvDRR13","HvDRR27", "HvDRR28"))]} 
 
#list with windows types
# dir.create("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV")
window_type_list<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/E_window_types/E.4.0_methylation_per_window_type.RDS")
peri_windows<-window_type_list$peri_windows
distal_windows<-window_type_list$distal_windows
# new_list<-list(); for(p in Populations){for (c in 1:7){new_list[[p]][[c]]<-all_windows[[p]][[1]][[c]]}}; all_windows<-new_list; rm(new_list)
coldspots<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.7_coldspots/D.7.0_coldspots_per_pop_windows.RDS")
coldspots_low_methy<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/E_window_types/E.1.3_low_methy_cold_windows_SUMMED_win=10000.RDS")
#coldspots_low_methy_not_int<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/E_coldspots/E.8_low_methy_cold_NOT_INT_windows_SUMMED_win=10000.RDS")
#coldspots_low_methy_int<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/E_coldspots/E.12_low_methy_cold_INT_windows_SUMMED_win=10000.RDS")
coldspots_high_methy<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/E_window_types/E.1.5_high_methy_cold_windows_SUMMED_win=10000.RDS")
#coldspots_high_methy_not_int<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/E_coldspots/E.10_high_methy_cold_NOT_INT_windows_SUMMED_win=10000.RDS")
#coldspots_high_methy_int<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/E_coldspots/E.14_high_methy_cold_INT_windows_SUMMED_win=10000.RDS")
hotspots<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.1_hotspots_per_pop.RDS")
#hotspots_int<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.3_hotspots_intersected.RDS")
#hotspots_not_int<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.5_hotspots_per_pop_not_intersected.RDS")
# hotspots_FIRST<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.2.1_hotspots_per_pop_FIRST_LAYER.RDS")

rec_events<-list(distal_windows, peri_windows, coldspots, coldspots_low_methy, coldspots_high_methy, hotspots)
names(rec_events)<-c("distal_windows", "peri_windows","coldspots", "coldspots_low_methy", "coldspots_high_methy", "hotspots")
# rm(list = c("all_windows","coldspots", "coldspots_low_methy", "coldspots_high_methy", "hotspots"))
saveRDS(rec_events, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.0_windows_list.RDS")

rec_events<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.0_windows_list.RDS")

###################### get closest SVs #############################

get.closest.SV.per.CO<-function(x=x){
  dif<-abs(x-SV_pos)
  return(make.matrix(c(x,min(dif)[1],SV_lengths[names(SV_pos[which(dif==min(dif))[1]])])))
}

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS")
SVs<-names(sv_list)
sv_sizes<-names(sv_list[[1]][[1]][[1]])

#window_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.1_windows_closest_SV.RDS")

# #el del all windows mejor correrlo por el server
rec_events$peri_windows<-NULL
windows_type<-names(rec_events)
window_SV_list<-list()
for (e in windows_type){ cat(e, fill = TRUE)
window_SV_list[[e]]<-list()
for (v in SVs){ cat(v); cat(": ", fill = TRUE)
  window_SV_list[[e]][[v]]<-list()
    for (p in Populations){ cat(p); cat(": ")
      window_SV_list[[e]][[v]][[p]]<-list()
      for (c in 1:7){ cat(c); cat("-")
        window_SV_list[[e]][[v]][[p]][[c]]<-list()
        for (s in sv_sizes){  #cat(s); cat("-")
          if (nrow(sv_list[[v]][[p]][[c]][[s]])!=0){
            SVs_p_c<-sv_list[[v]][[p]][[c]][[s]]
            row.names(SVs_p_c)<-paste(v,"_",c,"_",1:nrow(SVs_p_c), sep = "")
            #add lengths
            SVs_p_c<-cbind(SVs_p_c, NA)
            SVs_p_c[,ncol(SVs_p_c)]<-as.numeric(SVs_p_c[,3])-as.numeric(SVs_p_c[,2])
            #get closer SVs to random COs
            SV_pos<-as.numeric(SVs_p_c[,2:3]); names(SV_pos)<-c(row.names(SVs_p_c),row.names(SVs_p_c))
            SV_lengths<-SVs_p_c[,ncol(SVs_p_c)]; names(SV_lengths)<-row.names(SVs_p_c)
            if (e%in%c("hotspots_int")){ windows<-rec_events[[e]][[c]]; windows_mid_point<-windows[,2]-(win.size/2)}
            if (e%in%c("coldspots_low_methy", "coldspots_high_methy", "coldspots_low_methy_not_int", "coldspots_high_methy_not_int", "coldspots_high_methy_int", "coldspots_low_methy_int")){ #if2
            windows<-names(rec_events[[e]][[p]][[c]])
            windows_mid_point<-as.numeric(windows)-(win.size/2)
            } else {windows<-rec_events[[e]][[p]][[c]]; windows_mid_point<-windows[,2]-(win.size/2)}
            if (length(windows_mid_point)!=0){
            LIST<-lapply(windows_mid_point, FUN = get.closest.SV.per.CO)
            TABLA<-as.matrix(transpose(as.data.frame(matrix(unlist(LIST), nrow = 3, ncol = length(LIST)))))
            colnames(TABLA)<-c("window_mid_point","distance_to_SV", "SV_length")
            window_SV_list[[e]][[v]][[p]][[c]][[s]]<-TABLA
            }#if
          }#if
        }#s
      }#c
    }#p
  cat(fill = TRUE)
 }#v
#saveRDS(window_SV_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1_windows_closest_SV_",e,".RDS", sep = ""))
}#e
# #re-order names
window_SV_list2<-list()
order<-c("distal_windows","hotspots","coldspots","coldspots_low_methy", "coldspots_high_methy")
for (e in order){window_SV_list2[[e]]<-window_SV_list[[e]]}
saveRDS(window_SV_list2, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.1_windows_closest_SV.RDS")
# 
# ####################################################################

window_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.1_windows_closest_SV.RDS")

#################################################################################
#unify chrs and SV_sizes
window_SV_list2<-list()
for (e in names(window_SV_list)){ cat(e, fill = TRUE)
  window_SV_list2[[e]]<-list()  
  for (p in Populations){ cat(p); cat(": ")   
    big_size<-c()
    small_size<-c()
    for (c in 1:7){ cat(c); cat("-")   
      small_size_c<-matrix(nrow = nrow(window_SV_list[[e]][[1]][[p]][[c]][[1]]), ncol = 1)
      row.names(small_size_c)<-window_SV_list[[e]][[1]][[p]][[c]][[1]][,1]
      big_size_c<-small_size_c
      for (v in SVs){ #cat(v); cat(": ", fill = TRUE)
        sv_sizes<-names(window_SV_list[[e]][[v]][[p]][[c]])
        for (s in sv_sizes){ # cat(s); cat("-")
          if (s%in%c("50—299_bp", "0.3—4.9_kb")){
            small_size_c<-cbind(small_size_c, window_SV_list[[e]][[v]][[p]][[c]][[s]][,2])  
          } else {big_size_c<-cbind(big_size_c, window_SV_list[[e]][[v]][[p]][[c]][[s]][,2])}  
        }#s
      }#v    
      big_size<-c(big_size, unlist(apply(big_size_c, 1, function(x){return(min(as.numeric(x), na.rm = TRUE))})))
      small_size<-c(small_size, unlist(apply(small_size_c, 1, function(x){return(min(as.numeric(x), na.rm = TRUE))})))
    }#c
    window_SV_list2[[e]][[p]]<-list(small_size, big_size); names(window_SV_list2[[e]][[p]])<-c("small", "big")   
  }#p
}#e   
saveRDS(window_SV_list2, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.2_window_closest_SV.RDS")
window_SV_list<-window_SV_list2; rm(window_SV_list2)

#add COs
CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO/B.1.2_SV_per_CO_per_SV_size.RDS")
e<-"COs"
window_SV_list[[e]]<-list()
for (p in Populations){ cat(p); cat(": ")
  big_size<-c(); small_size<-c()
  for (c in 1:7){ cat(c); cat("-")
    big_size_c<-matrix(ncol = 2, nrow = 0)
    small_size_c<-matrix(ncol = 2, nrow = 0) 
    for (v in names(CO_SV_list)[-5]){ #cat(v); cat(": ", fill = TRUE)
      sv_sizes<-names(CO_SV_list[[v]][[3]][[p]][[c]])   
      for (s in sv_sizes){ # cat(s); cat("-")
        if (s%in%c("50—299_bp", "0.3—4.9_kb")){small_size_c<-rbind(small_size_c, CO_SV_list[[v]][[3]][[p]][[c]][[s]][,c("breakpoint","distance")])  
        } else {big_size_c<-rbind(big_size_c, CO_SV_list[[v]][[3]][[p]][[c]][[s]][,c("breakpoint","distance")])}  
      }#s
    }#v  
    breakpoints<-unique(big_size_c[,1])    
    closest_svs<-unlist(lapply(breakpoints, function(x){return(min(as.numeric(big_size_c[which(big_size_c[,1]==x),2])))}))
    names(closest_svs)<-breakpoints
    big_size<-c(big_size, closest_svs); 
    #
    breakpoints<-unique(small_size_c[,1])     
    closest_svs<-unlist(lapply(breakpoints, function(x){return(min(as.numeric(small_size_c[which(small_size_c[,1]==x),2])))}))
    names(closest_svs)<-breakpoints
    small_size<-c(small_size, closest_svs)
  }#c
  window_SV_list[[e]][[p]]<-list(small_size, big_size); names(window_SV_list[[e]][[p]])<-c("small", "big")   
}#p
saveRDS(window_SV_list, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.2_window_closest_SV_unified.RDS")

#add COs in hotspots and not in hotspots

# #window_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.2.1_window_closest_SV.RDS")
# COs_types<-list()
# COs_types[["COs_hotspots"]]<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.7_COs_in_hotspots.RDS")
# COs_types[["COs_not_hotspots"]]<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.8_COs_not_in_hotspots.RDS")
# for (e in names(COs_types)){
# window_SV_list[[e]]<-list()
# for (p in Populations){    
#   window_SV_list[[e]][[p]]<-list()
#   small<-c(); big<-c()
#   for (c in 1:7){
#   small<-c(small, unlist(window_SV_list[["COs"]][[p]][["small"]][which(names(window_SV_list[["COs"]][[p]][["small"]])%in%COs_types[[e]][[p]][[c]][,3])]))
#   big<-c(big, unlist(window_SV_list[["COs"]][[p]][["big"]][which(names(window_SV_list[["COs"]][[p]][["big"]])%in%COs_types[[e]][[p]][[c]][,3])]))
#   }#c
#   window_SV_list[[e]][[p]][["small"]]<-small
#   window_SV_list[[e]][[p]][["big"]]<-big
# }#p
# }#e 
# saveRDS(window_SV_list, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.2.1_window_closest_SV.RDS")

#################################################################################






###BIG : keep big sv category with the biggest possible
# small_sizes<-sv_sizes[1:5]
# 
# #unify chrs and SV_sizes
# window_SV_list2<-list()
# for (e in names(window_SV_list)){ cat(e, fill = TRUE)
#   window_SV_list2[[e]]<-list()  
#   for (p in Populations){ cat(p); cat(": ")   
#     big_size<-c()
#     small_size<-c()
#     for (c in 1:7){ cat(c); cat("-")   
#       small_size_c<-matrix(nrow = nrow(window_SV_list[[e]][[1]][[p]][[c]][[1]]), ncol = 1)
#       row.names(small_size_c)<-window_SV_list[[e]][[1]][[p]][[c]][[1]][,1]
#       big_size_c<-small_size_c
#       for (v in SVs){ #cat(v); cat(": ", fill = TRUE)
#         sv_sizes<-names(window_SV_list[[e]][[v]][[p]][[c]])
#         for (s in sv_sizes){ # cat(s); cat("-")
#           if (s%in%c(small_sizes)){
#             small_size_c<-cbind(small_size_c, window_SV_list[[e]][[v]][[p]][[c]][[s]][,2])  
#           } else {big_size_c<-cbind(big_size_c, window_SV_list[[e]][[v]][[p]][[c]][[s]][,2])}  
#         }#s
#       }#v    
#       big_size<-c(big_size, unlist(apply(big_size_c, 1, function(x){return(min(as.numeric(x), na.rm = TRUE))})))
#       small_size<-c(small_size, unlist(apply(small_size_c, 1, function(x){return(min(as.numeric(x), na.rm = TRUE))})))
#     }#c
#     window_SV_list2[[e]][[p]]<-list(small_size, big_size); names(window_SV_list2[[e]][[p]])<-c("small", "big")   
#   }#p
# }#e   
# saveRDS(window_SV_list2, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.3_window_closest_SV_BIG.RDS")
# window_SV_list<-window_SV_list2; rm(window_SV_list2)
# 
# window_SV_list2<-list()
# order<-c("distal_windows","hotspots","coldspots","coldspots_low_methy", "coldspots_high_methy")
# for (e in order){window_SV_list2[[e]]<-window_SV_list[[e]]}
# saveRDS(window_SV_list2, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.3_windows_closest_SV_BIG.RDS")
# # 
# # ####################################################################
# 
# 
# #add COs
# CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO/B.1.2_SV_per_CO_per_SV_size.RDS")
# e<-"COs"
# window_SV_list[[e]]<-list()
# for (p in Populations){ cat(p); cat(": ")
#   big_size<-c(); small_size<-c()
#   for (c in 1:7){ cat(c); cat("-")
#     big_size_c<-matrix(ncol = 2, nrow = 0)
#     small_size_c<-matrix(ncol = 2, nrow = 0) 
#     for (v in names(CO_SV_list)[-5]){ #cat(v); cat(": ", fill = TRUE)
#       sv_sizes<-names(CO_SV_list[[v]][[3]][[p]][[c]])   
#       for (s in sv_sizes){ # cat(s); cat("-")
#         if (s%in%c(small_sizes)){small_size_c<-rbind(small_size_c, CO_SV_list[[v]][[3]][[p]][[c]][[s]][,c("breakpoint","distance")])  
#         } else {big_size_c<-rbind(big_size_c, CO_SV_list[[v]][[3]][[p]][[c]][[s]][,c("breakpoint","distance")])}  
#       }#s
#     }#v  
#     breakpoints<-unique(big_size_c[,1])    
#     closest_svs<-unlist(lapply(breakpoints, function(x){return(min(as.numeric(big_size_c[which(big_size_c[,1]==x),2])))}))
#     names(closest_svs)<-breakpoints
#     big_size<-c(big_size, closest_svs); 
#     #
#     breakpoints<-unique(small_size_c[,1])     
#     closest_svs<-unlist(lapply(breakpoints, function(x){return(min(as.numeric(small_size_c[which(small_size_c[,1]==x),2])))}))
#     names(closest_svs)<-breakpoints
#     small_size<-c(small_size, closest_svs)
#   }#c
#   window_SV_list[[e]][[p]]<-list(small_size, big_size); names(window_SV_list[[e]][[p]])<-c("small", "big")   
# }#p
# saveRDS(window_SV_list, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.4_window_closest_SV_unified_BIG.RDS")
