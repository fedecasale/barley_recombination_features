library(data.table)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

random_COs_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.2_random_COs_with_prob/B.2.1_random_COs_with_PROB_per_window_GENETIC_MAPS.RDS")
window_sizes<-names(random_COs_list)
layers<-names(random_COs_list[[1]])
Populations<-names(random_COs_list[[1]][[1]])
###################### get closest SVs #############################

get.closest.SV.per.CO<-function(x=x){
  dif<-abs(x-SV_pos)
  return(make.matrix(c(x,min(dif)[1],SV_lengths[names(SV_pos[which(dif==min(dif))[1]])])))
}

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS")
SVs<-names(sv_list)
sv_sizes<-names(sv_list[[1]][[1]][[1]])

random_CO_SV_list<-list()
for (w in window_sizes){ cat(w); cat(": ", fill = TRUE)
  random_CO_SV_list[[w]]<-list()
  for (v in SVs){ cat(v); cat(": ", fill = TRUE)
    random_CO_SV_list[[w]][[v]]<-list()
    for (l in layers){ cat(l, fill = TRUE)
      random_CO_SV_list[[w]][[v]][[l]]<-list()
      for (p in Populations){ cat(p); cat(": ")   
        random_CO_SV_list[[w]][[v]][[l]][[p]]<-list()
        for (c in 1:7){ cat(c); cat("-")
          random_CO_SV_list[[w]][[v]][[l]][[p]][[c]]<-list()
          for (s in sv_sizes){
            if (nrow(sv_list[[v]][[p]][[c]][[s]])!=0){
              SVs_p_c<-sv_list[[v]][[p]][[c]][[s]]  
              row.names(SVs_p_c)<-paste(v,"_",c,"_",1:nrow(SVs_p_c), sep = "")
              #add lengths 
              SVs_p_c<-cbind(SVs_p_c, NA)
              SVs_p_c[,ncol(SVs_p_c)]<-as.numeric(SVs_p_c[,3])-as.numeric(SVs_p_c[,2])
              #get closer SVs to random COs
              SV_pos<-as.numeric(SVs_p_c[,2:3]); names(SV_pos)<-c(row.names(SVs_p_c),row.names(SVs_p_c))
              SV_lengths<-SVs_p_c[,ncol(SVs_p_c)]; names(SV_lengths)<-row.names(SVs_p_c)
              random_COs<-random_COs_list[[w]][[l]][[p]][[c]]
              LIST<-lapply(random_COs, FUN = get.closest.SV.per.CO) 
              TABLA<-as.matrix(transpose(as.data.frame(matrix(unlist(LIST), nrow = 3, ncol = length(LIST)))))
              colnames(TABLA)<-c("breakpoint","distance_to_SV", "SV_length")
              random_CO_SV_list[[w]][[v]][[l]][[p]][[c]][[s]]<-TABLA 
            }#if
          }#s
        }#c
      }#p
      cat(fill = TRUE)  
    }#l
  }#v
}#w

saveRDS(random_CO_SV_list, "/scratch/Federico/3_RILs/3_SVs/Results/B.2_random_COs_with_prob/B.2.3_SV_per_size_per_RANDOM_CO_GENETIC_MAPS.RDS")

####################################################################