make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/B.0_SVs_per_windows")

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

window_sizes<-c(10000, 500000, 1000000)[1] #5000000, 10000000)

#get SV per pop
sv_per_pop<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size_ALL_POPS.RDS"))
Populations<-names(sv_per_pop[[1]]) 
SVs<-names(sv_per_pop)
sv_sizes<-names(sv_per_pop[[1]][[1]][[1]])
#sacar indels
sv_sizes<-sv_sizes[-1]

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window")

#for (s in window_sizes){ cat("window size "); cat(s); cat(": ")
s<-window_sizes

if (s == 10000){
Populations<-Populations[which(Populations%in%c("HvDRR13","HvDRR27", "HvDRR28"))]
rec_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.0_windows_list.RDS")[[1]]
} else {
rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
}  
sv_per_window<-list()

for (p in Populations){  cat(p); cat(": ")

sv_per_window[[p]]<-list()

for (c in 1:7) { cat(c); cat("-")

windows<-rec_list[[1]][[c]][,2]
if (s==10000){windows<-windows+1}

#create a register to check if svs in a window are overlapping
window_register<-lapply(windows, function(x){return(x<-list())}); names(window_register)<-windows

for (v in SVs){ #cat(v, fill = TRUE)

#los que empiezan y terminan en una misma window
for (g in sv_sizes){  
if (nrow(sv_per_pop[[v]][[p]][[c]][[g]])!=0){  

#get SV data for that chr
SV_p_c<-sv_per_pop[[v]][[p]][[c]][[g]]

if (nrow(SV_p_c)!=0){ #0

#add SVs count to windows  
for (w in 1:length(windows)){  

#get svs that start in that window
pos_to_get<-SV_p_c[which(as.numeric(SV_p_c[,2])<=as.numeric(windows[w])),]; pos_to_get<-make.matrix(pos_to_get)
pos_to_get<-pos_to_get[which(as.numeric(pos_to_get[,2])>=(as.numeric(windows[w])-s)),]; pos_to_get<-make.matrix(pos_to_get)

if (nrow(pos_to_get)!=0){ #1
  
#si finaliza en la misma window lo adiero
if (any(as.numeric(pos_to_get[,3])<=as.numeric(windows[w]))) { #2
  pos_to_get2<-pos_to_get[which(as.numeric(pos_to_get[,3])<=as.numeric(windows[w])),]; pos_to_get2<-make.matrix(pos_to_get2)
  window_register[[w]][[length(window_register[[w]])+1]]<-pos_to_get2
  pos_to_get<-pos_to_get[-which(as.numeric(pos_to_get[,3])<=as.numeric(windows[w])),]; pos_to_get<-make.matrix(pos_to_get)
}#if2

#si queda algo es porque termina en otra window
if (nrow(pos_to_get)!=0){ #3
  #adiero lo de la window w
  pos_to_get2<-pos_to_get; pos_to_get2[,3]<-windows[w]; pos_to_get2<-make.matrix(pos_to_get2)
  window_register[[w]][[length(window_register[[w]])+1]]<-pos_to_get2

  #lleno las otras
  for (i in 1:nrow(pos_to_get)){
    #la win donde terminan  
    last_win<-windows[which(as.numeric(windows)>as.numeric(pos_to_get[i,3]))][1]
    pos_to_get2<-pos_to_get[i,]; pos_to_get2<-make.matrix(pos_to_get2)
    pos_to_get2[,2]<-as.numeric(last_win)-s; pos_to_get2<-make.matrix(pos_to_get2)
    window_register[[paste(last_win)]][[length(window_register[[paste(last_win)]])+1]]<-pos_to_get2
    
    all_win<-windows[w:which(windows==last_win)]
    if (length(all_win)>2){ #4
      med_win<-all_win[-c(1,length(all_win))]
      for (m in 1:length(med_win)){
        pos_to_get2<-pos_to_get[i,]; pos_to_get2<-make.matrix(pos_to_get2)
        pos_to_get2[,2]<-as.numeric(med_win[m])-s+1; pos_to_get2<-make.matrix(pos_to_get2)
        pos_to_get2[,3]<-as.numeric(med_win[m]); pos_to_get2<-make.matrix(pos_to_get2)
        window_register[[paste(med_win[m])]][[length(window_register[[paste(med_win[m])]])+1]]<-pos_to_get2
      }#m  
      #to.add<-s-1
      #for (m in 1:length(med_win)){sv_table[med_win[m],g]<-sum(as.numeric(sv_table[med_win[m],g]), to.add)}
    }#if4  
  }#i
}#if3    
}#if1
}#w
}#0
}#if
}#g


}#v



sv_per_window[[p]][[c]]<-window_register

}#c

cat(fill = TRUE)

}#p  

#}#s

###############

saveRDS(sv_per_window, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.1_SVs_per_window=",s,".RDS", sep = ""))

################
#saco el occupied percentage per window

sv_per_window<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.1_SVs_per_window=",s,".RDS", sep = ""))
sv_per_window2<-sv_per_window
for (p in Populations){for (c in 1:7){for (w in names(sv_per_window[[p]][[c]])){sv_per_window2[[p]][[c]][[w]]<-NA}}}

for (p in Populations){ cat(p); cat(":")
  for (c in 1:7){ cat(c); cat("-")
    for (w in names(sv_per_window[[p]][[c]])){ #cat(c(w,"-"))
      
      if (length(sv_per_window[[p]][[c]][[w]])!=0){
      #make one table of all tables in window
      window_table<-make.matrix(do.call("rbind", sv_per_window[[p]][[c]][[w]])[,2:3])
      #take ranges positions
      WIN<-apply(window_table, 1, FUN = function(x){return(as.numeric(x[1]):as.numeric(x[2]))})
      #count positions
      WIN<-length(unique(unlist(WIN)))   #CON EL UNIQUE CUENTA CADA POSITION, UNA VEZ SOLA
      #calculate percentage
      sv_per_window2[[p]][[c]][[w]]<-WIN/s
      } else {sv_per_window2[[p]][[c]][[w]]<-0} #if
    
    }#w
  }#c
cat(fill = TRUE)  
}#p  
saveRDS(sv_per_window2, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.2_SVs_proportion_per_window=",s,".RDS", sep = ""))
################

percentage_list<-list()
for (c in 1:7){
  windows<-names(sv_per_window[[p]][[c]])  
  tabla<-matrix(ncol = length(Populations), nrow = length(windows))
  colnames(tabla)<-Populations; row.names(tabla)<-windows
  for (p in Populations){tabla[,p]<-unlist(sv_per_window2[[p]][[c]])}
  percentage_list[[c]]<-tabla
}#c

saveRDS(percentage_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))

################

#HAGO LO MISMO PERO POR SV TYPE POR SEPARADO

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

#get SV per pop
sv_per_pop<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size_ALL_POPS.RDS"))
Populations<-names(sv_per_pop[[1]]) 
SVs<-names(sv_per_pop)
sv_sizes<-names(sv_per_pop[[1]][[1]][[1]])
#sacar indels
sv_sizes<-sv_sizes[-1]

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

window_sizes<-c(500000, 1000000)[2] #5000000, 10000000)

#for (s in window_sizes){ cat("window size "); cat(s); cat(": ")
s<-window_sizes

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))

sv_per_window<-list()

for (v in SVs){ #cat(v, fill = TRUE)

sv_per_window[[v]]<-list()

for (p in Populations){  cat(p); cat(": ")
  
  sv_per_window[[v]][[p]]<-list()
  
  for (c in 1:7) { cat(c); cat("-")
    
    windows<-rec_list[[1]][[c]][,2]
    
    #create a register to check if svs in a window are overlapping
    window_register<-as.list(windows); names(window_register)<-windows; for (w in windows){window_register[[w]]<-list()}
    
      #los que empiezan y terminan en una misma window
      for (g in sv_sizes){  
        if (nrow(sv_per_pop[[v]][[p]][[c]][[g]])!=0){  
          
          #get SV data for that chr
          SV_p_c<-sv_per_pop[[v]][[p]][[c]][[g]]
          
          if (nrow(SV_p_c)!=0){ #0
            
            #add SVs count to windows  
            for (w in 1:length(windows)){  
              
              #get svs that start in that window
              pos_to_get<-SV_p_c[which(as.numeric(SV_p_c[,2])<=as.numeric(windows[w])),]; pos_to_get<-make.matrix(pos_to_get)
              pos_to_get<-pos_to_get[which(as.numeric(pos_to_get[,2])>=(as.numeric(windows[w])-s)),]; pos_to_get<-make.matrix(pos_to_get)
              
              if (nrow(pos_to_get)!=0){ #1
                
                #si finaliza en la misma window lo adiero
                if (any(as.numeric(pos_to_get[,3])<=as.numeric(windows[w]))) { #2
                  pos_to_get2<-pos_to_get[which(as.numeric(pos_to_get[,3])<=as.numeric(windows[w])),]; pos_to_get2<-make.matrix(pos_to_get2)
                  window_register[[w]][[length(window_register[[w]])+1]]<-pos_to_get2
                  pos_to_get<-pos_to_get[-which(as.numeric(pos_to_get[,3])<=as.numeric(windows[w])),]; pos_to_get<-make.matrix(pos_to_get)
                }#if2
                
                #si queda algo es porque termina en otra window
                if (nrow(pos_to_get)!=0){ #3
                  #adiero lo de la window w
                  pos_to_get2<-pos_to_get; pos_to_get2[,3]<-windows[w]; pos_to_get2<-make.matrix(pos_to_get2)
                  window_register[[w]][[length(window_register[[w]])+1]]<-pos_to_get2
                  
                  #lleno las otras
                  for (i in 1:nrow(pos_to_get)){
                    #la win donde terminan  
                    last_win<-windows[which(as.numeric(windows)>as.numeric(pos_to_get[i,3]))][1]
                    pos_to_get2<-pos_to_get[i,]; pos_to_get2<-make.matrix(pos_to_get2)
                    pos_to_get2[,2]<-as.numeric(last_win)-s; pos_to_get2<-make.matrix(pos_to_get2)
                    window_register[[paste(last_win)]][[length(window_register[[paste(last_win)]])+1]]<-pos_to_get2
                    
                    all_win<-windows[w:which(windows==last_win)]
                    if (length(all_win)>2){ #4
                      med_win<-all_win[-c(1,length(all_win))]
                      for (m in 1:length(med_win)){
                        pos_to_get2<-pos_to_get[i,]; pos_to_get2<-make.matrix(pos_to_get2)
                        pos_to_get2[,2]<-as.numeric(med_win[m])-s+1; pos_to_get2<-make.matrix(pos_to_get2)
                        pos_to_get2[,3]<-as.numeric(med_win[m]); pos_to_get2<-make.matrix(pos_to_get2)
                        window_register[[paste(med_win[m])]][[length(window_register[[paste(med_win[m])]])+1]]<-pos_to_get2
                      }#m  
                      #to.add<-s-1
                      #for (m in 1:length(med_win)){sv_table[med_win[m],g]<-sum(as.numeric(sv_table[med_win[m],g]), to.add)}
                    }#if4  
                  }#i
                }#if3    
              }#if1
            }#w
          }#0
        }#if
      }#g
      
    sv_per_window[[v]][[p]][[c]]<-window_register
    
  }#c
  
  cat(fill = TRUE)
  
}#p  
}#v
#}#s

###############

saveRDS(sv_per_window, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.4_SVs_per_window=",s,"_SV_TYPES.RDS", sep = ""))

################
#saco el occupied percentage per window

sv_per_window<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.4_SVs_per_window=",s,"_SV_TYPES.RDS", sep = ""))
sv_per_window2<-sv_per_window
for (v in SVs){for (p in Populations){for (c in 1:7){for (w in names(sv_per_window[[v]][[p]][[c]])){sv_per_window2[[v]][[p]][[c]][[w]]<-NA}}}}

for (v in SVs){ cat(v, fill = TRUE)
for (p in Populations){ cat(p); cat(":")
  for (c in 1:7){ cat(c); cat("-")
    for (w in names(sv_per_window[[v]][[p]][[c]])){ #cat(c(w,"-"))
      
      if (length(sv_per_window[[v]][[p]][[c]][[w]])!=0){
        #make one table of all tables in window
        window_table<-make.matrix(do.call("rbind", sv_per_window[[v]][[p]][[c]][[w]])[,2:3])
        #take ranges positions
        WIN<-apply(window_table, 1, FUN = function(x){return(as.numeric(x[1]):as.numeric(x[2]))})
        #count positions
        WIN<-length(unique(unlist(WIN)))
        #calculate percentage
        sv_per_window2[[v]][[p]][[c]][[w]]<-WIN/s
      } else {sv_per_window2[[v]][[p]][[c]][[w]]<-0} #if
      
    }#w
  }#c
  cat(fill = TRUE)  
}#p 
}#v
saveRDS(sv_per_window2, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.5_SVs_proportion_per_window=",s,"_SV_TYPES.RDS", sep = ""))
################

percentage_list<-list()
for (v in SVs){
percentage_list[[v]]<-list()
for (c in 1:7){
  windows<-names(sv_per_window[[v]][[p]][[c]])  
  tabla<-matrix(ncol = length(Populations), nrow = length(windows))
  colnames(tabla)<-Populations; row.names(tabla)<-windows
  for (p in Populations){tabla[,p]<-unlist(sv_per_window2[[v]][[p]][[c]])}
  percentage_list[[v]][[c]]<-tabla
}#c
}#v

saveRDS(percentage_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))

################