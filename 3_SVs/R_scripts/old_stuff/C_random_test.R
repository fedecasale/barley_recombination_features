wilcox_test<-function(x=x){
  distances<-c(real_distances, x)
  return(wilcox.test(distances~data_type, alternative = "greater", paired = FALSE)$p.value)      
  }

CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO.RDS")

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

layers<-c("1_first", "2_second", "3_third")
layers<-layers[3]
SVs<-c("deletions", "duplications", "insertions") #"translocations")
Populations<-names(CO_SV_list[[1]][[1]])

runs<-1000

#test_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/C_wilcox_test_pvalues.RDS")
test_list<-list()


#dir.create("/scratch/Federico/3_RILs/3_SVs/Results/C_wilcox_text")

for (v in SVs){ cat(v, fill = TRUE)
  test_list[[v]]<-list()  
  for (l in layers){ cat(l); cat(" - ")
    test_list[[v]][[l]]<-list()
    tabla<-matrix(ncol = 3, nrow = 7); row.names(tabla)<-chrs; colnames(tabla)<-Populations
    for (p in Populations){ cat(p);cat(": ")
      random_CO_SV<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/B.2_SV_per_RANDOM_CO/",v,"/",l,"/",p,".RDS", sep = ""))
      for (c in 1:7){ cat(c);cat("-")
        real_distances<-as.numeric(CO_SV_list[[v]][[l]][[p]][[c]][,11])
        n.breakpoints<-length(real_distances)
        #I want to probe that real distances are higher than random distances  
        data_type<-c(rep("REAL", n.breakpoints), rep("RANDOM", n.breakpoints))
        random_distances<-list(); for (r in 1:runs){random_distances[[r]]<-as.numeric(random_CO_SV[[c]][[r]][,2])}
        #p_values<-lapply(random_distances, wilcox_test)
        #tabla[c,p]<-round(median(unlist(p_values)), digits = 4)
        tabla[c,p]<-round(median(unlist(lapply(random_distances, wilcox_test))), digits = 4)
      }#c
      cat(" ")  
    }#p
    cat(" ")
    test_list[[v]][[l]]<-tabla
  }#l
  cat(fill = TRUE)
}#v

saveRDS(test_list, "/scratch/Federico/3_RILs/3_SVs/Results/C_wilcox_test_pvalues.RDS")
