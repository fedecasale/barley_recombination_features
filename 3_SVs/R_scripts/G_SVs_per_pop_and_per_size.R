pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[,1]

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/G_corr_sv_rec")

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

SVs<-c("deletions", "duplications", "insertions", "inversions", "translocations")
SVs<-SVs[1:4]

#1-get parents SVs per pop

sv_list<-list()

for (v in SVs){ cat(v); cat(": ")
  
  SV<-as.matrix(fread(paste("/2data/Barley/Hv-DRR-DNAseq2020/mapping_V3/results/SVs/genotyped_",v,"_genome_final_sorted_header_unsplitChroms.csv", sep = ""))) 
  colnames(SV)<-gsub("Unumli_Arpa","Unumli-Arpa", colnames(SV))
  colnames(SV)<-gsub("W23829_803911","W23829/803911", colnames(SV))
  
  if (v%in%"translocations"){  #les voy a poner start y end en el mismo lugar
    SV_out<-SV[,c(1,3,3,5:ncol(SV))]; row.names(SV_out)<-paste(1:nrow(SV_out),"_OUT", sep = "")  
    SV_in<-SV[,c(2,4,4,5:ncol(SV))] ; row.names(SV_in)<-paste(1:nrow(SV_in),"_IN", sep = "")  
    SV<-rbind(SV_out, SV_in)
    colnames(SV)[1:3]<-c("chr","start", "end")
  }
  
  sv_list[[v]]<-list()
  
  for (p in Populations){ cat(p); cat("-")
    
    sv_list[[v]][[p]]<-list()  
    
    #select parents
    P1<-pop_info[which(pop_info[,1]==p),2]
    P2<-pop_info[which(pop_info[,1]==p),3]
    
    #get parent SV data
    cat(length(which(colnames(SV)%in%c(P1,P2)))==2); cat("- ")
    SV_p<-SV[,c(1,2,3,which(colnames(SV)%in%c(P1,P2)))]
    
    #only keep positions where SVs are present
    SV_p<-SV_p[unique(c(which(SV_p[,4]%in%"1"), which(SV_p[,5]%in%"1"))), ]
    #only keep positions where there is a SV between both parents
    to.delete<-which(SV_p[,4]=="1")[which(SV_p[,4]=="1")%in%which(SV_p[,5]=="1")]
    SV_p<-SV_p[-to.delete,]
    
    for (c in 1:7) {
      sv_list[[v]][[p]][[c]]<-SV_p[which(SV_p[,1]%in%chrs[c]),]
    }#c 
    
  }#p 
  cat(fill = TRUE)
}#v
saveRDS(sv_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/G_corr_sv_rec/G.1.1_SVs_per_population.RDS"))


#2-divide in sizes
#SV sizes categories according to Marius: (A: 50 - 300bp; B: 0.3 - 5kb; C: 5 - 50kb; D: 50 - 250kb; E: 0.25 - 1Mb)

source("/scratch/Federico/3_RILs/2_CO_breakpoints/R_scripts/Z_my_functions.R")

#sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_corr_sv_rec/G.1_SVs_per_population.RDS")
#SVs<-names(sv_list)[-which(names(sv_list)=="translocations")]

Populations<-names(sv_list[[1]])

new_sv_list<-list()
for (v in SVs){ cat(c(v,"-"))
  new_sv_list[[v]]<-list()
  for (p in Populations){ 
    new_sv_list[[v]][[p]]<-list()
    for (c in 1:7){  
      
      new_sv_list[[v]][[p]][[c]]<-list()
      
      tabla<-sv_list[[v]][[p]][[c]]
      
      sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
      
      new_sv_list[[v]][[p]][[c]][["<=49_bp"]]<-make.matrix(tabla[which(sv_lengths<=49),])
      if (any(sv_lengths<=49)){tabla<-tabla[-which(sv_lengths<=49),]}; tabla<-make.matrix(tabla)
      sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
      
      new_sv_list[[v]][[p]][[c]][["50—299_bp"]]<-make.matrix(tabla[which(sv_lengths<=299),])
      if (any(sv_lengths<=299)){tabla<-tabla[-which(sv_lengths<=299),]}; tabla<-make.matrix(tabla)
      sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
      
      new_sv_list[[v]][[p]][[c]][["0.3—4.9_kb"]]<-make.matrix(tabla[which(sv_lengths<=4999),])
      if (any(sv_lengths<=4999)){tabla<-tabla[-which(sv_lengths<=4999),]}; tabla<-make.matrix(tabla)
      sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
      
      new_sv_list[[v]][[p]][[c]][["5—49_kb"]]<-make.matrix(tabla[which(sv_lengths<=49999),])
      if (any(sv_lengths<=49999)){tabla<-tabla[-which(sv_lengths<=49999),]}; tabla<-make.matrix(tabla)
      sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
      
      new_sv_list[[v]][[p]][[c]][["50—249_kb"]]<-make.matrix(tabla[which(sv_lengths<=249999),])
      if (any(sv_lengths<=249999)){tabla<-tabla[-which(sv_lengths<=249999),]}; tabla<-make.matrix(tabla)
      sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
      
      new_sv_list[[v]][[p]][[c]][["0.25—1_Mbp"]]<-make.matrix(tabla[which(sv_lengths<=999999),])
      if (any(sv_lengths<=999999)){tabla<-tabla[-which(sv_lengths<=999999),]}; tabla<-make.matrix(tabla)
      sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
      
      new_sv_list[[v]][[p]][[c]][[">1_Mbp"]]<-make.matrix(tabla[which(sv_lengths>=1000000),])
      
    }
  }
}

saveRDS(new_sv_list, "/scratch/Federico/3_RILs/3_SVs/Results/G_corr_sv_rec/G.1.2_SV_per_size.RDS")

################################################# 

#TABLE FOR PAPER
sv_list<-new_sv_list; rm(new_sv_list)
#sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/G_corr_sv_rec/G.2_SV_per_size.RDS")

tabla_0<-matrix(ncol = 2+length(Populations), nrow = 0)
colnames(tabla_0)<-c("SV", "SV size", Populations)

for (v in SVs){ cat(v, fill = TRUE)
  sv_sizes<-names(sv_list[[v]][[1]][[1]])
  tabla<-matrix(ncol = 2+length(Populations), nrow = length(sv_sizes))
  colnames(tabla)<-c("SV", "SV size", Populations)
  row.names(tabla)<-sv_sizes
  tabla[,1]<-v; tabla[,2]<-sv_sizes
  for (s in sv_sizes){  
    for (p in Populations){ #cat(p);cat(": ")
      sv_size_p_c<-c(); for (c in 1:7){ sv_size_p_c<-c(sv_size_p_c, nrow(sv_list[[v]][[p]][[c]][[s]]))}#c
      tabla[s,p]<-sum(sv_size_p_c)  
    }#p
  }#s
  tabla_0<-rbind(tabla_0, tabla)
  cat(fill = TRUE)
}#v

# #add translocations
# translo<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.1_SVs_per_population.RDS")[["translocations"]]
# tabla<-matrix(ncol = 5, nrow = 1); colnames(tabla)<-c("SV", "SV size", Populations)
# tabla[,1]<-"translocations"; tabla[,2]<-"-"
# for (p in Populations){ 
#   sv_size_p_c<-c(); for (c in 1:7){ sv_size_p_c<-c(sv_size_p_c, nrow(translo[[p]][[c]]))}#c
#   tabla[,p]<-sum(sv_size_p_c)  
# }#p
# tabla_0<-rbind(tabla_0, tabla)

write.csv(tabla_0, "/scratch/Federico/3_RILs/3_SVs/Results/G_corr_sv_rec/G.1.3_SV_per_population.csv", row.names = FALSE)

#################################################
