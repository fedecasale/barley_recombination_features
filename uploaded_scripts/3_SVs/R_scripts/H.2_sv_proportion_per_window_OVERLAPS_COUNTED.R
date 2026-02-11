#because different svs overlap, in the normal procedure, a given physical pos is occupied once (regardless if several SVs match it)   
#here I count the matches, which it is biolgically incorrect, but provides better measurement of populations variations...

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/B.0_SVs_per_windows")

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

window_sizes<-c(500000, 1000000)[2] #5000000, 10000000)

#get SV per pop
sv_per_pop<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size_ALL_POPS.RDS"))
Populations<-names(sv_per_pop[[1]]) 
SVs<-names(sv_per_pop)
sv_sizes<-names(sv_per_pop[[1]][[1]][[1]])
#sacar indels
sv_sizes<-sv_sizes[-1]

s<-window_sizes

#saco el occupied percentage per window (counting each position more than one)
sv_per_window<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.1_SVs_per_window=",s,".RDS", sep = ""))
sv_per_window2<-sv_per_window
for (p in Populations){for (c in 1:7){for (w in names(sv_per_window[[p]][[c]])){sv_per_window2[[p]][[c]][[w]]<-NA}}}

for (p in Populations){ cat(p); cat(":")
  for (c in 1:7){ cat(c); cat("-")
    for (w in names(sv_per_window[[p]][[c]])){ #cat(c(w,"-"))
    
      if (length(sv_per_window[[p]][[c]][[w]])!=0){
        #make one table of all tables in window
        window_table<-make.matrix(do.call("rbind", sv_per_window[[p]][[c]][[w]])[,2:3])
        ranges<-apply(window_table  , 1, FUN = function(x){return(abs(as.numeric(x[2])-as.numeric(x[1])))})
        used_bp<-sum(ranges)
        #calculate percentage
        sv_per_window2[[p]][[c]][[w]]<-used_bp/s
      } else {sv_per_window2[[p]][[c]][[w]]<-0} #if
      
    }#w
  }#c
  cat(fill = TRUE)  
}#p  
saveRDS(sv_per_window2, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.2_SVs_proportion_per_window=",s,".RDS", sep = ""))
################

percentage_list<-list()
for (c in 1:7){
  windows<-names(sv_per_window[[p]][[c]])  
  tabla<-matrix(ncol = length(Populations), nrow = length(windows))
  colnames(tabla)<-Populations; row.names(tabla)<-windows
  for (p in Populations){tabla[,p]<-unlist(sv_per_window2[[p]][[c]])}
  percentage_list[[c]]<-tabla
}#c

saveRDS(percentage_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.3_SVs_proportion_per_window=",s,"_all_pops.RDS", sep = ""))

################

#HAGO LO MISMO PERO POR SV TYPE POR SEPARADO

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
          ranges<-apply(window_table  , 1, FUN = function(x){return(abs(as.numeric(x[2])-as.numeric(x[1])))})
          used_bp<-sum(ranges)
          #calculate percentage
          sv_per_window2[[v]][[p]][[c]][[w]]<-used_bp/s
        } else {sv_per_window2[[v]][[p]][[c]][[w]]<-0} #if
        
      }#w
    }#c
    cat(fill = TRUE)  
  }#p 
}#v
saveRDS(sv_per_window2, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.5_SVs_proportion_per_window=",s,"_SV_TYPES.RDS", sep = ""))
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

saveRDS(percentage_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.2.6_SVs_proportion_per_window=",s,"_all_pops_SV_TYPES.RDS", sep = ""))

################