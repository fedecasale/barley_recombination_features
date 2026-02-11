library("xlsx")

CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO/B.1.2_SV_per_CO_per_SV_size.RDS")

random_CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.2_random_COs_with_prob/B.2.3_SV_per_size_per_RANDOM_CO_GENETIC_MAPS.RDS")
window_sizes<-names(random_CO_SV_list)[order(as.numeric(names(random_CO_SV_list)))]
SVs<-names(random_CO_SV_list[[1]])
layers<-names(random_CO_SV_list[[1]][[1]])
Populations<-names(random_CO_SV_list[[1]][[1]][[1]])

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS")
#sv_sizes<-names(sv_list[[1]][[1]][[1]]); rm(sv_list)

chrs<-paste("chr",1:7,"H", sep = "")

tabla<-matrix(ncol = length(layers), nrow = length(SVs)); row.names(tabla)<-SVs; colnames(tabla)<-layers

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/C_wilcox_test")

tabla<-matrix(ncol = 11, nrow = 1)
colnames(tabla)<-c("window_size", "layer", "SV_type", "SV_size","p.value","ave_real", "med_real", "var_real", "ave_random", "med_random", "var_random")

tabla_0<-tabla

for (w in window_sizes){ cat(w); cat(": ")
  for (v in SVs){ cat(v, fill = TRUE)
    for (l in layers){ cat(l); cat(" - ")
      sv_sizes<-names(CO_SV_list[[v]][[l]][[1]][[1]])
      for (s in sv_sizes){
      
      tabla[]<-NA
      tabla[,c(1:4)]<-c(w,l,v,s)  
      
      real_distances<-c()
      random_distances<-c()
      
      for (p in Populations){ #cat(p);cat(": ")
        for (c in 1:7){ #cat(c);cat("-")
        
            real_distances<-c(real_distances, as.numeric(CO_SV_list[[v]][[l]][[p]][[c]][[s]][,11]))
            random_distances<-c(random_distances, as.numeric(random_CO_SV_list[[w]][[v]][[l]][[p]][[c]][[s]][,2])) #elijo la run 1
                
        }#c
      }#p
      
      if (length(real_distances)==length(random_distances)){
        
        tabla[,6:8]<-c(round(mean(real_distances)),median(real_distances), round(sd(real_distances)))  
        tabla[,9:11]<-c(round(mean(random_distances)),median(random_distances), round(sd(random_distances)))  
        
        #I want to probe that real distances are higher than random distances  
        n.breakpoints<-length(real_distances)
        data_type<-c(rep("REAL", n.breakpoints), rep("RANDOM", n.breakpoints))
        distances<-c(real_distances, random_distances)
        tabla[,5]<-round(wilcox.test(distances~data_type, alternative = "two.sided", paired = FALSE)$p.value, digits = 4)
        
      } else {cat("DISTANCES LENGTHS DIFFER", fill = TRUE)}
      
      tabla_0<-rbind(tabla_0, tabla)

      }#s
    }#l
    cat(fill = TRUE)
  }#v
}#w

write.csv(tabla_0[-1,], "/scratch/Federico/3_RILs/3_SVs/Results/C_wilcox_test/C.2_wilcox_test_with_SV_sizes.csv", row.names = FALSE)
