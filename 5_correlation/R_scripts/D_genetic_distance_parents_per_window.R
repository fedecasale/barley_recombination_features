source("/scratch/Federico/3_RILs/2_CO_breakpoints/R_scripts/Y_gen_distance_Benjamin.R")

dir.create("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances")

s<-1000000

genetic_distances<-list()

#for (c in 1:7){
c<-7
cat(c); cat(": ")
genetic_distances[[c]]<-list()

chr_win_list<-readRDS(paste("/scratch/Federico/3_RILs/1_Marius_data/Results/D_Marius_data_in_windows/D.2_Marius_per_window/D.2_windows_chr",c,"H.RDS", sep = ""))

windows<-names(chr_win_list)

for (w in windows){ cat(w); cat("-")
  
win_table<-chr_win_list[[w]][,-1]
win_table[which(win_table[]=="")]<-NA

genetic_distances[[c]][[w]]<-genetic.distance(win_table)

}#w

saveRDS(genetic_distances, paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D_parent_genetic_distances",c,".RDS", sep = ""))

#}#c
#saveRDS(genetic_distances, paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.1_parent_genetic_distances.RDS", sep = ""))





#gather all chrs in one list
genetic_distances<-list()
for (c in 1:7){
chr<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D_parent_genetic_distances",c,".RDS", sep = ""))
genetic_distances[[c]]<-chr[[c]]
}
saveRDS(genetic_distances, paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.1_parent_genetic_distances.RDS", sep = ""))





genetic_distances<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.1_parent_genetic_distances.RDS", sep = ""))

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[,1]

pop_list<-list()

for (c in 1:7){ cat(c); cat("-")

windows<-names(genetic_distances[[c]])

pop_table<-matrix(nrow = length(windows), ncol = length(Populations))
row.names(pop_table)<-windows
colnames(pop_table)<-Populations

for (p in Populations){  cat(p); cat(": ")
  
#get "P1" "P2"
P1<-pop_info[which(pop_info[,1]==p),2]
P2<-pop_info[which(pop_info[,1]==p),3]
  
for (w in windows){
  pop_table[w,p]<-as.matrix(genetic_distances[[c]][[w]])[P1,P2]
}#w

}#p

pop_list[[c]]<-pop_table

}#c

saveRDS(pop_list, paste("/scratch/Federico/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
