library("xlsx")

CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO/B.1.1_SV_per_CO.RDS")

random_CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.3_random_COs_with_prob/B.3.2_SV_per_RANDOM_CO_GENETIC_MAPS.RDS")

window_sizes<-names(random_CO_SV_list)[order(as.numeric(names(random_CO_SV_list)))]
Populations<-names(CO_SV_list[[1]][[1]])
chrs<-paste("chr",1:7,"H", sep = "")
layers<-c("1_first", "2_second", "3_third")
SVs<-c("deletions", "duplications", "insertions") #"translocations")

tabla<-matrix(ncol = length(layers), nrow = length(SVs)); row.names(tabla)<-SVs; colnames(tabla)<-layers

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/C_wilcox_test")

tabla<-matrix(ncol = 10, nrow = 1)
colnames(tabla)<-c("window_size", "layer", "SV_type", "p.value","ave_real", "med_real", "var_real", "ave_random", "med_random", "var_random")

tabla_0<-tabla

for (w in window_sizes){ cat(w); cat(": ")

    for (v in SVs){ cat(v, fill = TRUE)
    for (l in layers){ cat(l); cat(" - ")
      
      tabla[,c(1:3)]<-c(w,l,v)  
      
      real_distances<-c()
      random_distances<-c()
      
      for (p in Populations){ cat(p);cat(": ")
        for (c in 1:7){ #cat(c);cat("-")
          real_distances<-c(real_distances, as.numeric(CO_SV_list[[v]][[l]][[p]][[c]][,11]))
          random_distances<-c(random_distances, as.numeric(random_CO_SV_list[[w]][[v]][[l]][[p]][[c]][,2])) #elijo la run 1
        }#c
      }#p

      if (length(real_distances)==length(random_distances)){
      
      tabla[,5:7]<-c(round(mean(real_distances)),median(real_distances), round(sd(real_distances)))  
      tabla[,8:10]<-c(round(mean(random_distances)),median(random_distances), round(sd(random_distances)))  
        
      #I want to probe that real distances are higher than random distances  
      n.breakpoints<-length(real_distances)
      data_type<-c(rep("REAL", n.breakpoints), rep("RANDOM", n.breakpoints))
      distances<-c(real_distances, random_distances)
      tabla[,4]<-round(wilcox.test(distances~data_type, alternative = "two.sided", paired = FALSE)$p.value, digits = 4)
      
      } else {cat("DISTANCES LENGTHS DIFFER", fill = TRUE)}
      
    tabla_0<-rbind(tabla_0, tabla)
    
    }#l
    cat(fill = TRUE)
  }#v
}#w

write.csv(tabla_0[-1,], "/scratch/Federico/3_RILs/3_SVs/Results/C_wilcox_test/C.1.1_wilcox_test.csv", row.names = FALSE)

