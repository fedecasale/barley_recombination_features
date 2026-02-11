library(data.table)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

get.closest.SV.per.CO<-function(x=x){
  dif<-abs(x-SV_pos)
  return(make.matrix(c(x,min(dif)[1],SV_lengths[names(SV_pos[which(dif==min(dif))[1]])])))
}

layers<-c("1_first", "2_second", "3_third")
SVs<-c("deletions", "duplications", "insertions") #"translocations")

chrs<-paste("chr", 1:7, "H", sep = "")
chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_per_population.RDS")
breakpoint_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")

Populations<-names(breakpoint_list[[1]])

random_CO_SV_list<-list()
for (v in SVs){ cat(v); cat(": ", fill = TRUE)
  random_CO_SV_list[[v]]<-list()
  for (l in layers){ cat(l, fill = TRUE)
    random_CO_SV_list[[v]][[l]]<-list()
    for (p in Populations){ cat(p); cat(": ")   
      random_CO_SV_list[[v]][[l]][[p]]<-list()
      for (c in 1:7){ cat(c); cat("-") 
        SVs_p_c<-sv_list[[v]][[p]][[c]]  
        row.names(SVs_p_c)<-paste(v,"_",c,"_",1:nrow(SVs_p_c), sep = "")
        #add lengths
        SVs_p_c<-cbind(SVs_p_c, NA)
        SVs_p_c[,ncol(SVs_p_c)]<-as.numeric(SVs_p_c[,3])-as.numeric(SVs_p_c[,2])
        #create random breakpoints     
        random_COs<-sample(1:as.numeric(chr_lengths[c, 3]), nrow(breakpoint_list[[l]][[p]][[c]]))
        #get closer SVs
        SV_pos<-as.numeric(SVs_p_c[,2:3]); names(SV_pos)<-c(row.names(SVs_p_c),row.names(SVs_p_c))
        SV_lengths<-SVs_p_c[,ncol(SVs_p_c)]; names(SV_lengths)<-row.names(SVs_p_c)
        LIST<-lapply(random_COs, FUN = get.closest.SV.per.CO) 
        TABLA<-as.matrix(transpose(as.data.frame(matrix(unlist(LIST), nrow = 3, ncol = length(LIST)))))
        colnames(TABLA)<-c("breakpoint","distance_to_SV", "SV_length")
        random_CO_SV_list[[v]][[l]][[p]][[c]]<-TABLA 
      }#c
    }#p
  cat(fill = TRUE)  
  }#l
}#v
saveRDS(random_CO_SV_list, "/scratch/Federico/3_RILs/3_SVs/Results/B.2_SV_per_RANDOM_CO.RDS")
