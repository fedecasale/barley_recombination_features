library(data.table)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

############## generate random COs with PROB #######################

accumulated_rec_prob<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.3_CO_prob_per_window_GENETIC_MAPS.RDS", sep = ""))
window_sizes<-names(accumulated_rec_prob)
layers<-c("1_first", "2_second", "3_third")
chrs<-paste("chr", 1:7, "H", sep = "")
chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))

breakpoint_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")
Populations<-names(breakpoint_list[[1]])

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/B.2_random_COs_with_prob")

random_COs_list<-list()
for (w in window_sizes){
random_COs_list[[w]]<-list()
  for (l in layers){ cat(l, fill = TRUE)
    random_COs_list[[w]][[l]]<-list()
    for (p in Populations){ cat(p); cat(": ")   
      random_COs_list[[w]][[l]][[p]]<-list()
      for (c in 1:7){ cat(c); cat("-") 
        #create random breakpoints with probability per window      
        #N.COs_to_generate<-nrow(breakpoint_list[[l]][[p]][[c]]) #I create the same N COs per chr
        real_breakpoints<-breakpoint_list[[l]][[p]][[c]]
        #NO DEBERIA HABER MAS
        if (any(real_breakpoints[,6]=="")){real_breakpoints<-real_breakpoints[-which(real_breakpoints[,6]==""),]}
        N.COs_to_generate<-nrow(real_breakpoints) #I create the same N COs per chr
        N.windows<-nrow(accumulated_rec_prob[[w]][[p]][[c]])
        PROB.windows<-as.numeric(accumulated_rec_prob[[w]][[p]][[c]][,3]) #I used prob per pop per chr per win
        random_win<-sample(x = 1:N.windows, size = N.COs_to_generate, prob = PROB.windows, replace = TRUE)
        random_win<-accumulated_rec_prob[[w]][[p]][[c]][random_win,]
        random_COs<-apply(random_win, 1, function(x=x){return(sample(x[1]:x[2], size = 1))})
        random_COs_list[[w]][[l]][[p]][[c]]<-random_COs
      }#c
    }#p
  cat(fill = TRUE)  
  }#l
}#w
saveRDS(random_COs_list, "/scratch/Federico/3_RILs/3_SVs/Results/B.2_random_COs_with_prob/B.2.1_random_COs_with_PROB_per_window_GENETIC_MAPS.RDS")
####################################################################
