source("/scratch/Federico/3_RILs/2_CO_breakpoints/R_scripts/Z_my_functions.R")

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.1_SVs_per_population_ALL_POPS.RDS")

#SV sizes categories according to Marius: (A: 50 - 300bp; B: 0.3 - 5kb; C: 5 - 50kb; D: 50 - 250kb; E: 0.25 - 1Mb)

SVs<-names(sv_list)[-which(names(sv_list)=="translocations")]
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
saveRDS(new_sv_list, "/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size_ALL_POPS.RDS")

#make shorther list for the 3 populations
sv_list2<-new_sv_list
for (i in 1:length(sv_list)){
  sv_list2[[i]]<-sv_list[[i]][c(10,24,25)]
}
saveRDS(sv_list2, paste("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS"))

################################################# 

#ACA SIGO SOLO CON LAS 3 POPS
Populations<-Populations[c(10,24,25)]

#TABLE FOR PAPER
sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS")

tabla_0<-matrix(ncol = 5, nrow = 0)
colnames(tabla)<-c("SV", "SV size", Populations)

for (v in SVs){ cat(v, fill = TRUE)
sv_sizes<-names(sv_list[[v]][[1]][[1]])
tabla<-matrix(ncol = 5, nrow = length(sv_sizes))
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

#add translocations
translo<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.1_SVs_per_population.RDS")[["translocations"]]
tabla<-matrix(ncol = 5, nrow = 1); colnames(tabla)<-c("SV", "SV size", Populations)
tabla[,1]<-"translocations"; tabla[,2]<-"-"
for (p in Populations){ 
  sv_size_p_c<-c(); for (c in 1:7){ sv_size_p_c<-c(sv_size_p_c, nrow(translo[[p]][[c]]))}#c
  tabla[,p]<-sum(sv_size_p_c)  
}#p
tabla_0<-rbind(tabla_0, tabla)

write.csv(tabla_0, "/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.3_SV_per_population.csv", row.names = FALSE)

#################################################