window_sizes<-c(500000, 1000000, 5000000, 10000000)
layers<-c("1_first", "2_second", "3_third")
SVs<-c("deletions", "duplications", "insertions") #, "translocations")
Populations<-c("HvDRR13","HvDRR27","HvDRR28")
chrs<-paste("chr",1:7,"H", sep = "")

rec_per_window<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.4_breakpoints_per_window.RDS", sep = ""))

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/C_cor_per_windows")

for (v in SVs){ cat(v, fill = TRUE)
  
  sv_per_window<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/B_SVs_per_windows/B_",v,"_per_window.RDS", sep = ""))
  
  window_list<-list()
  
  for (s in window_sizes){ cat("window size "); cat(s); cat(": ")
    
    window_list[[paste(s)]]<-list()
    
    for (l in layers){ cat(l); cat(" layer"); cat(" - ")
      
      layer_table<-matrix(ncol = 9, nrow = 0)
        
      for (p in Populations){ cat(p); cat(" - ")
        
        rils<-names(rec_per_window[[paste(s)]][[l]][[p]])  
        
        tabla_pop<-matrix(nrow = length(rils), ncol = 7)
        row.names(tabla_pop)<-rils; colnames(tabla_pop)<-chrs
        
        for (r in rils){ #cat(r); cat(" - ")
          
          for (c in 1:7) {  
          
          sv<-as.numeric(sv_per_window[[paste(s)]][[l]][[p]][[r]][[c]][,3])      
          rec<-as.numeric(rec_per_window[[paste(s)]][[l]][[p]][[r]][[c]][,3]) 
          cor<-cor.test(y = rec, x = sv, method = "pearson")
          if (is.na(cor$p.value)==FALSE){if (cor$p.value<0.05){tabla_pop[r,c]<-round(cor$estimate, digits = 2)}}
          
          }#c
          
        }#r
        
        tabla_pop<-cbind(p, row.names(tabla_pop),tabla_pop)
        layer_table<-rbind(layer_table, tabla_pop)
        
      }#p
      
      row.names(layer_table)<-NULL
      colnames(layer_table)<-c("Populations", "rils", chrs)
      window_list[[paste(s)]][[l]]<-layer_table
    
    }#l
  }#s
  
  saveRDS(window_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/C_cor_per_windows/C.1_",v,"_cor_per_window.RDS", sep = ""))
  
}#v



### create excel ###

for (v in SVs){ cat(v, fill = TRUE)
  
window_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/C_cor_per_windows/C.1_",v,"_cor_per_window.RDS", sep = ""))

if (file.exists(paste("/scratch/Federico/3_RILs/3_SVs/Results/C_cor_per_windows/C.2_",v,"_cor_per_window.xlsx", sep = ""))){
file.remove(paste("/scratch/Federico/3_RILs/3_SVs/Results/C_cor_per_windows/C.2_",v,"_cor_per_window.xlsx", sep = ""))  
}

for (s in window_sizes){cat(s); cat(": ")

final_tabla<-as.data.frame(window_list[[paste(s)]][[1]])
for (c in 3:9){final_tabla[[c]]<-as.numeric(as.character(final_tabla[[c]]))}
for (l in layers[2:3]){#cat(l); cat("-")
tabla<-as.data.frame(window_list[[paste(s)]][[l]])[,3:9]
for (c in 1:7){tabla[[c]]<-as.numeric(as.character(tabla[[c]]))}
final_tabla<-cbind(final_tabla, NA,tabla)
}#l

write.xlsx(final_tabla, file = paste("/scratch/Federico/3_RILs/3_SVs/Results/C_cor_per_windows/C.2_",v,"_cor_per_window.xlsx", sep = ""), append = TRUE, sheetName = paste(s), row.names = FALSE, showNA = FALSE)

}#s
cat(fill = TRUE)
}#v

