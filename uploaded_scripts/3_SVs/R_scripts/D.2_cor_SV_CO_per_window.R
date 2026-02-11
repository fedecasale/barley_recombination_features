make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

rec_prob<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1.2_CO_prob_per_window_DENSITY.RDS", sep = ""))
window_sizes<-names(rec_prob)
layers<-c("1_first", "2_second", "3_third")
chrs<-paste("chr", 1:7, "H", sep = "")
chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
Populations<-names(rec_prob[[1]][[1]])
SVs<-c("deletions", "duplications", "insertions", "inversions") #"translocations")

#get barley genes
genes_per_window<-list()
barley_genes<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Barley_Morex_V2_gene_annotation_PGSB.all.descriptions.csv"))
for (s in window_sizes[4]){ cat("window size "); cat(s); cat(": ")  
genes_per_window[[s]]<-list()  
for (c in 1:7){
  #get windows from rec data
  windows<-rec_prob[[paste(s)]][[1]][[1]][[c]]
  windows[,3]<-NA
  #get N genes in window
  chr.table<-barley_genes[which(barley_genes[,3]==chrs[c]),] 
  for (w in 1:nrow(windows)){ #no cuento los genes que empiezan en una y terminan en otra
  pos_to_get<-chr.table[which(as.numeric(chr.table[,5])>=as.numeric(windows[w,1])),]; pos_to_get<-make.matrix(pos_to_get)
  pos_to_get<-pos_to_get[which(as.numeric(pos_to_get[,6])<=as.numeric(windows[w,2])),]; pos_to_get<-make.matrix(pos_to_get)
  windows[w,3]<-nrow(pos_to_get)
  }
  genes_per_window[[s]][[c]]<-windows
}
}
#################

#SIN DIVIDIR PER SV SIZE

file.remove("/scratch/Federico/3_RILs/3_SVs/Results/D_correlation/D.2.1_cor_CO_SV.xlsx")

SVs_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/D_correlation/D.1_SVs_per_size_per_windows.RDS", sep = ""))

for (v in SVs){ #cat(v); cat(": ")

  sv_per_window<-SVs_list[[v]]
    
  w<-window_sizes[4]

  #sv_per_window<-sv_per_window[[w]]
  for (p in Populations){
    for (c in 1:7){
    a<-sv_per_window[[w]][[p]][[c]][[1]]
    if (length(sv_per_window[[w]][[p]][[c]])>1){
    for (i in 2:length(sv_per_window[[w]][[p]][[c]])){a[,3]<-as.numeric(a[,3])+as.numeric(sv_per_window[[w]][[p]][[c]][[i]][,3])}
    }#if
    sv_per_window[[w]][[p]][[c]]<-a
    }}
  
  tabla<-matrix(nrow = 0, ncol = 25)

  for (l in layers){# cat(l); cat(" - ")

  tabla_pop<-matrix(nrow = 3, ncol = 7); row.names(tabla_pop)<-Populations; colnames(tabla_pop)<-chrs
  tabla_pop1<-tabla_pop
  tabla_pop2<-tabla_pop

    for (p in Populations){ #cat(p);cat("-")

      for (c in 1:7){

      rec<-as.numeric(rec_prob[[w]][[l]][[p]][[c]][,3])
      sv<-as.numeric(sv_per_window[[w]][[p]][[c]][,3])
      gene_den<-as.numeric(genes_per_window[[w]][[c]][,3])

      #cor<-cor.test(y = rec, x = sv, method = "pearson")
      #if (is.na(cor$p.value)==FALSE){if (cor$p.value<0.05){tabla_pop[p,c]<-round(cor$estimate, digits = 2)}}
      tabla_pop[p,c]<-round(anova(lm(rec~gene_den))$`Pr(>F)`[[1]], digits = 5)
      tabla_pop1[p,c]<-round(anova(lm(rec~sv))$`Pr(>F)`[[1]], digits = 5)
      tabla_pop2[p,c]<-round(anova(lm(rec~sv+gene_den))$`Pr(>F)`[[1]], digits = 5)

      }#c
    }#p

  tabla<-rbind(tabla, cbind(l, Populations, tabla_pop," ",tabla_pop1, " ", tabla_pop2))

  }#l
#}#w

write.xlsx(tabla, "/scratch/Federico/3_RILs/3_SVs/Results/D_correlation/D.2.1_cor_CO_SV.xlsx", append = TRUE, sheetName = v)

}#v

#################

##ESTO ES SOLO PARA CHEQUEAR QUE POR SIZE ES LO MISMO

sv_per_window<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/D_correlation/D.1_SVs_per_size_per_windows.RDS", sep = ""))

sv_sizes<-names(sv_per_window[[3]][[1]][[1]][[1]])[2:3]

for (s in sv_sizes){
  
file.remove(paste("/scratch/Federico/3_RILs/3_SVs/Results/D_correlation/D.2.2_cor_CO_SV_per_sv_size_",s,".xlsx", sep = ""))
  
for (v in SVs){ #cat(v); cat(": ")

  #for (w in window_sizes){ cat("window size "); #cat(w); cat(": ")
  w<-window_sizes[4]

  tabla<-matrix(nrow = 0, ncol = 25)

  for (l in layers){ #cat(l); cat(" - ")

    tabla_pop<-matrix(nrow = 3, ncol = 7); row.names(tabla_pop)<-Populations; colnames(tabla_pop)<-chrs
    tabla_pop1<-tabla_pop
    tabla_pop2<-tabla_pop

    for (p in Populations){ #cat(p);cat("-")

      for (c in 1:7){

        rec<-as.numeric(rec_prob[[w]][[l]][[p]][[c]][,3])
        sv<-as.numeric(sv_per_window[[v]][[w]][[p]][[c]][[s]][,3])
        gene_den<-as.numeric(genes_per_window[[w]][[c]][,3])

        #cor<-cor.test(y = rec, x = sv, method = "pearson")
        #if (is.na(cor$p.value)==FALSE){if (cor$p.value<0.05){tabla_pop[p,c]<-round(cor$estimate, digits = 2)}}
        tabla_pop[p,c]<-round(anova(lm(rec~gene_den))$`Pr(>F)`[[1]], digits = 5)
        tabla_pop1[p,c]<-round(anova(lm(rec~sv))$`Pr(>F)`[[1]], digits = 5)
        tabla_pop2[p,c]<-round(anova(lm(rec~sv+gene_den))$`Pr(>F)`[[1]], digits = 5)
        
      }#c
    }#p

    tabla<-rbind(tabla, cbind(l, Populations, tabla_pop," ",tabla_pop1, " ", tabla_pop2))

  }#l
  #}#w

  write.xlsx(tabla, paste("/scratch/Federico/3_RILs/3_SVs/Results/D_correlation/D.2.2_cor_CO_SV_per_sv_size_",s,".xlsx", sep = ""), append = TRUE, sheetName = v)

}#v
}#s

#################################

#TABLA PER PAPER, I WILL ONLY CONSIDER LAYER 1 COs
file.remove( "/scratch/Federico/3_RILs/6_graphics_for_paper/Results/tables/D.2_cor_CO_SV.xlsx")

SVs_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/D_correlation/D.1_SVs_per_size_per_windows.RDS", sep = ""))

w<-window_sizes[4]

Populations<-c("HvDRR13", "HvDRR27", "HvDRR28")

TABLA2<-matrix(nrow = 0, ncol = 10)
TABLA3<-matrix(nrow = 0, ncol = 10)
tabla<-matrix(nrow = 3, ncol = 7); row.names(tabla)<-Populations; colnames(tabla)<-paste(1:7,"H", sep = "")

for (v in SVs){ #cat(v); cat(": ")
    
  sv_per_window<-SVs_list[[v]]
  
  w<-window_sizes[4]
  
  #sv_per_window<-sv_per_window[[w]]
  for (p in Populations){
    for (c in 1:7){
      a<-sv_per_window[[w]][[p]][[c]][[1]]
      if (length(sv_per_window[[w]][[p]][[c]])>1){
        for (i in 2:length(sv_per_window[[w]][[p]][[c]])){a[,3]<-as.numeric(a[,3])+as.numeric(sv_per_window[[w]][[p]][[c]][[i]][,3])}
      }#if
      sv_per_window[[w]][[p]][[c]]<-a
    }}
  
    tabla2<-tabla; tabla3<-tabla
    
    for (p in Populations){ #cat(p);cat("-")
      
      for (c in 1:7){
        
        rec<-as.numeric(rec_prob[[w]][[1]][[p]][[c]][,3]) #uso solo first layer
        sv<-as.numeric(sv_per_window[[w]][[p]][[c]][,3])
        gene_den<-as.numeric(genes_per_window[[w]][[c]][,3])
        
        tabla2[p,c]<-round(anova(lm(rec~sv))$`Pr(>F)`[[1]], digits = 5)
        tabla3[p,c]<-round(anova(lm(rec~sv+gene_den))$`Pr(>F)`[[1]], digits = 5)
        
      }#c
    }#p
    
    tabla2<-cbind("rec~sv",v, Populations, tabla2); 
    tabla3<-cbind("rec~sv+geneden",v, Populations, tabla3)
    
    TABLA2<-rbind(TABLA2, tabla2)
    TABLA3<-rbind(TABLA3, tabla3)
}#v

TABLA<-rbind(TABLA2,TABLA3)
colnames(TABLA)<-c("corr","SV","Population", paste(1:7,"H", sep = ""))
write.xlsx(TABLA, "/scratch/Federico/3_RILs/6_graphics_for_paper/Results/tables/D.2_cor_CO_SV.xlsx", row.names = FALSE)

##########################
