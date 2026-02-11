sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS")
SVs<-names(sv_list)
sv_sizes<-names(sv_list[[1]][[1]][[1]])
Populations<-names(sv_list$deletions)
chrs<-paste("chr", 1:7, sep = "")  #clave sin la H para correr el permTest

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/E_stats")

tabla_0<-matrix(ncol = 7, nrow = 0)
colnames(tabla_0)<-c("sv_type","sv_size","ave_num","ave_length", "median_length", "mean_inter_distance", "median_inter_distance")    

for (v in SVs){ cat(v); cat(" - ")

  tabla<-matrix(ncol = 7, nrow = 7)
  colnames(tabla)<-c("sv_type","sv_size","ave_num","ave_length", "median_length", "mean_inter_distance", "median_inter_distance")    
  row.names(tabla)<-sv_sizes
  
  tabla[,1]<-v  
  
  for (s in sv_sizes){   
  
  sv<-matrix(ncol = 2, nrow = 0)
  inter_distances<-c()
  
    for (p in Populations){ cat(p);cat(": ")
      for (c in 1:7){ cat(c);cat("-")
        if (nrow(sv_list[[v]][[p]][[c]][[s]])>1){
        sv<-rbind(sv, sv_list[[v]][[p]][[c]][[s]][,2:3]) 
        inter_distances<-c(inter_distances, (as.numeric(sv[2:nrow(sv),1])-as.numeric(sv[1:(nrow(sv)-1),2])))
        }#if
      }#c
    }#p
  
    if (is.null(inter_distances)==FALSE){
    tabla[s,2]<-s
    tabla[s,3]<-round(nrow(sv)/length(Populations)/length(chrs))
    tabla[s,4]<-round(mean(as.numeric(sv[,2])-as.numeric(sv[,1])))
    tabla[s,5]<-round(median(as.numeric(sv[,2])-as.numeric(sv[,1])))
    tabla[s,6]<-round(mean(inter_distances))
    tabla[s,7]<-round(median(inter_distances))
    }#if null
  
  }#s

tabla_0<-rbind(tabla_0,tabla)

}#v

tabla_0<-tabla_0[-which(is.na(tabla_0[,2])),]
write.csv(tabla_0, "/scratch/Federico/3_RILs/3_SVs/Results/E_stats/E.3_sv_stats.csv", row.names = FALSE)
