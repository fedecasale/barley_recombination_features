library(data.table)
library(xlsx)
library(agricolae)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[,1]

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

#parents.code<-as.matrix(read.table("/scratch/Federico/Paper_2/samplenames.txt"))

met_context<-c("chg", "cpg", "chh")

window_sizes<-c(10000, 500000, 1000000)[1]
win.size<-window_sizes
if (win.size == 10000){Populations<-Populations[which(Populations%in%c("HvDRR13","HvDRR27", "HvDRR28"))]} 

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS")
SVs<-names(sv_list)
sv_sizes<-names(sv_list[[1]][[1]][[1]])

#################################################################################

window_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.1.1_windows_closest_SV.RDS")

#################################################################################
#unify chrs and SV_sizes

for (e in names(window_SV_list)){ cat(e, fill = TRUE)
  for (v in SVs){ cat(v); cat(": ", fill = TRUE)
    for (p in Populations){ cat(p); cat(": ")   
      big_size<-c()
      small_size<-c()
      for (c in 1:7){ cat(c); cat("-")
        sv_sizes<-names(window_SV_list[[e]][[v]][[p]][[c]])
        for (s in sv_sizes){ # cat(s); cat("-")
        if (s%in%c("50—299_bp", "0.3—4.9_kb")){small_size<-c(small_size, window_SV_list[[e]][[v]][[p]][[c]][[s]][,2])  
        } else {big_size<-c(big_size, window_SV_list[[e]][[v]][[p]][[c]][[s]][,2])}  
      }#s
    }#c
    window_SV_list[[e]][[v]][[p]]<-list(small_size, big_size); names(window_SV_list[[e]][[v]][[p]])<-c("small", "big")   
  }#p
 }#v
}#e   

#### add COs ####
CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO/B.1.2_SV_per_CO_per_SV_size.RDS")
window_SV_list[["COs"]]<-list()
for (v in names(CO_SV_list)[-5]){ cat(v); cat(": ", fill = TRUE)
  window_SV_list[["COs"]][[v]]<-list()
  for (p in Populations){ cat(p); cat(": ")
    big_size<-c(); small_size<-c()
    for (c in 1:7){ cat(c); cat("-")
      sv_sizes<-names(CO_SV_list[[v]][[3]][[p]][[c]])
      for (s in sv_sizes){ # cat(s); cat("-")
        if (s%in%c("50—299_bp", "0.3—4.9_kb")){small_size<-c(small_size, CO_SV_list[[v]][[3]][[p]][[c]][[s]][,"distance"])  
        } else {big_size<-c(big_size, CO_SV_list[[v]][[3]][[p]][[c]][[s]][,"distance"])}  
      }#s
    }#c
    window_SV_list[["COs"]][[v]][[p]]<-list(small_size, big_size); names(window_SV_list[["COs"]][[v]][[p]])<-c("small", "big")   
  }#p
}#v

###################### check if means are different with ANOVA + Tukey #############################
file.remove("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.2.1_windows_types_TUKEY.xlsx")

for (v in SVs){ cat(v); cat(": ", fill = TRUE)

tabla_0<-matrix(ncol = 2+length(names(window_SV_list)), nrow = 0)
colnames(tabla_0)<-c("sv type", "Populations", names(window_SV_list))

  for (s in c("small", "big")){ cat(s); cat(":")

    tabla<-matrix(ncol = length(names(window_SV_list)), nrow = length(Populations))
    colnames(tabla)<-names(window_SV_list); row.names(tabla)<-Populations

   for (p in Populations){ cat(p); cat("-")
    #run anova + tukey
    data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("distance", "window_type")
    for (e in names(window_SV_list)){
    data_table_e<-cbind(window_SV_list[[e]][[v]][[p]][[s]], rep(e, length(window_SV_list[[e]][[v]][[p]][[s]])))
    colnames(data_table_e)<-c("distance", "window_type")
    data_table<-rbind(data_table, data_table_e)
    }#e
    data_table<-as.data.frame(data_table)
    data_table$distance<-as.numeric(as.character(data_table$distance))
    LM<-lm(distance~window_type, data = as.data.frame(data_table))
    #ANOVA + Tukey
    tukey<-HSD.test(LM, "window_type")$groups

    #fill results table
    for (e in names(window_SV_list)){tabla[p,e]<-paste(round(tukey[e,1], digits = 0), tukey[e,2])}#e
   }#p

  #add to sv table
  tabla_0<-rbind(tabla_0, cbind(c(paste(v,"_",s, sep = ""),"",""), Populations, tabla))

 }#s

cat(fill = TRUE)
write.xlsx(tabla_0, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.2.1_windows_types_TUKEY.xlsx", sheetName = v, append = TRUE)

}#v
####################################################################################################


########### check if means are different with multiple kruskall-wallis (WILCOX) ####################
file.remove("/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.2.2_windows_types_WILCOX.xlsx")

for (v in SVs){ cat(v); cat(": ", fill = TRUE)

  #for (s in c("small", "big")){ cat(s); cat(":")
  
  s<-"big" 
   
    for (p in Populations){ cat(p); cat("-")   
      data_table<-matrix(ncol = 2, nrow = 0); colnames(data_table)<-c("distance", "window_type")
      for (e in names(window_SV_list)){
        data_table_e<-cbind(window_SV_list[[e]][[v]][[p]][[s]], rep(e, length(window_SV_list[[e]][[v]][[p]][[s]])))
        colnames(data_table_e)<-c("distance", "window_type")
        data_table<-rbind(data_table, data_table_e)
      }#e
      data_table<-as.data.frame(data_table)
      data_table$distance<-as.numeric(as.character(data_table$distance))
      #check with kruskal-wallis pairwise comparisons between group levels with corrections for multiple testing.
      PWT<-pairwise.wilcox.test(data_table$distance, data_table$window_type, p.adjust.method = "BH")
      p.values_table<-format(round(PWT$p.value, digits = 3), scientific = FALSE)
      write.xlsx(p.values_table, "/scratch/Federico/3_RILs/3_SVs/Results/G_windows_closest_SV/G.2.2.2_windows_types_WILCOX.xlsx", sheetName = paste(v,"_",p, sep = ""), append = TRUE)
    }#p
    
  #}#s
  
  cat(fill = TRUE)  
  
}#v
####################################################################################################
