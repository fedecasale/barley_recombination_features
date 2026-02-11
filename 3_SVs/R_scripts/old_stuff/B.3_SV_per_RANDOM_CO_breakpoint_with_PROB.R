library(data.table)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

############## generate random COs with PROB #######################

accumulated_rec_prob<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1_CO_prob_per_window.RDS", sep = ""))
window_sizes<-names(accumulated_rec_prob)
layers<-c("1_first", "2_second", "3_third")
chrs<-paste("chr", 1:7, "H", sep = "")
chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))

breakpoint_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")
Populations<-names(breakpoint_list[[1]])

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/B.3_random_COs_with_prob")

random_COs_list<-list()
for (w in window_sizes){
random_COs_list[[w]]<-list()
  for (l in layers){ cat(l, fill = TRUE)
    random_COs_list[[w]][[l]]<-list()
    for (p in Populations){ cat(p); cat(": ")   
      random_COs_list[[w]][[l]][[p]]<-list()
      for (c in 1:7){ cat(c); cat("-") 
        #create random breakpoints with probability per window      
        N.COs_to_generate<-nrow(breakpoint_list[[l]][[p]][[c]]) #I create the same N COs per chr
        N.windows<-nrow(accumulated_rec_prob[[w]][[l]][[p]][[c]])
        PROB.windows<-as.numeric(accumulated_rec_prob[[w]][[l]][[p]][[c]][,3]) #I used prob per pop per chr per win
        random_win<-sample(x = 1:N.windows, size = N.COs_to_generate, prob = PROB.windows, replace = TRUE)
        random_win<-accumulated_rec_prob[[w]][[l]][[p]][[c]][random_win,]
        random_COs<-apply(random_win, 1, function(x=x){return(sample(x[1]:x[2], size = 1))})
        random_COs_list[[w]][[l]][[p]][[c]]<-random_COs
      }#c
    }#p
  cat(fill = TRUE)  
  }#l
}#w
saveRDS(random_COs_list, "/scratch/Federico/3_RILs/3_SVs/Results/B.3_random_COs_with_prob/B.3.1_random_COs_with_PROB_per_window.RDS")
####################################################################


###################### get closest SVs #############################

get.closest.SV.per.CO<-function(x=x){
  dif<-abs(x-SV_pos)
  return(make.matrix(c(x,min(dif)[1],SV_lengths[names(SV_pos[which(dif==min(dif))[1]])])))
}

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_per_population.RDS")
SVs<-names(sv_list)

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
          SVs_p_c<-sv_list[[v]][[p]][[c]]  
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
          random_CO_SV_list[[w]][[v]][[l]][[p]][[c]]<-TABLA 
        }#c
      }#p
      cat(fill = TRUE)  
    }#l
  }#v
}#w
saveRDS(random_CO_SV_list, "/scratch/Federico/3_RILs/3_SVs/Results/B.3_random_COs_with_prob/B.3.2_SV_per_RANDOM_CO.RDS")
  
####################################################################