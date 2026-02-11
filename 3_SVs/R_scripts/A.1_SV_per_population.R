pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[,1]
#Populations<-pop_info[c(10,24,25)]

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided")

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

SVs<-c("deletions", "duplications", "insertions", "inversions", "translocations")

sv_list<-list()

for (v in SVs){ cat(v); cat(": ")

SV<-as.matrix(fread(paste("/2data/Barley/Hv-DRR-DNAseq2020/mapping_V3/results/SVs/genotyped_",v,"_genome_final_sorted_header_unsplitChroms.csv", sep = ""))) 
colnames(SV)[which(colnames(SV)=="Unumli_Arpa")]<-"Unumli-Arpa"
colnames(SV)[which(colnames(SV)=="W23829_803911")]<-"W23829/803911"

if (v%in%"translocations"){  #les voy a poner start y end en el mismo lugar
SV_out<-SV[,c(1,3,3,5:ncol(SV))]; row.names(SV_out)<-paste(1:nrow(SV_out),"_OUT", sep = "")  
SV_in<-SV[,c(2,4,4,5:ncol(SV))] ; row.names(SV_in)<-paste(1:nrow(SV_in),"_IN", sep = "")  
SV<-rbind(SV_out, SV_in)
colnames(SV)[1:3]<-c("chr","start", "end")
}

sv_list[[v]]<-list()

for (p in Populations){ cat(p); cat(" - ")

sv_list[[v]][[p]]<-list()  

#select parents
P1<-pop_info[which(pop_info[,1]==p),2]
P2<-pop_info[which(pop_info[,1]==p),3]

#get parent SV data
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

saveRDS(sv_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.1_SVs_per_population_ALL_POPS.RDS"))

#make shorther list for the 3 populations
sv_list2<-sv_list
for (i in 1:length(sv_list)){
sv_list2[[i]]<-sv_list[[i]][c(10,24,25)]
}
saveRDS(sv_list2, paste("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.1_SVs_per_population.RDS"))
