make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_layers")

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

layers<-c("1_first", "2_second", "3_third")

SVs<-c("deletions", "duplications", "insertions") #"translocations")

for (v in SVs){ cat(v, fill = TRUE)

SV<-as.matrix(fread(paste("/2data/Barley/Hv-DRR-DNAseq2020/mapping_V3/results/SVs/genotyped_",v,"_genome_final_sorted_header_unsplitChroms.csv", sep = ""))) 
colnames(SV)<-gsub("_","-", colnames(SV))

sv_layers<-list()

for (l in layers){ cat(l); cat(" layer"); cat(" - "); 
#l<-layers[1]
rec_layer<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.",l,"_layer.RDS", sep = ""))
Populations<-names(rec_layer)  
sv_layers[[paste(l)]]<-list()  

for (p in Populations){ cat(p); cat(": ")
#p<-Populations[1]
rils<-names(rec_layer[[p]])  
sv_layers[[paste(l)]][[p]]<-list()  

#select parents
P1<-pop_info[which(pop_info[,1]==p),2]
P2<-pop_info[which(pop_info[,1]==p),3]

#get parent SV data
SV_p<-SV[,c(1,2,3,which(colnames(SV)%in%c(P1,P2)))]
#only keep positions where SVs are present
SV_p<-SV_p[unique(c(which(SV_p[,4]%in%"1"), which(SV_p[,5]%in%"1"))), ]
  
for (r in rils){
#r<-rils[1]
sv_layers[[paste(l)]][[p]][[r]]<-list()  

for (c in 1:7) { 
#c<-1

chr.length<-as.numeric(chr_lengths[c,3])  

#get rec data  
RIL_table<-rec_layer[[p]][[r]][[c]]
  
#get SV data for that chr
SV_p_c<-SV_p[which(SV_p[,1]==chrs[c]),]    
row.names(SV_p_c)<-paste("pos_", 1:nrow(SV_p_c), sep = "")

SV_table<-SV_p_c; SV_table[,4:5]<-NA; colnames(SV_table)[4:5]<-c("parent","value")

#determine SV info for parental blocks
for (b in 1:nrow(RIL_table)){
pos_to_get<-row.names(SV_p_c)[which(as.numeric(SV_p_c[,2])>=as.numeric(RIL_table[b,5]))]  
pos_to_get<-pos_to_get[which(as.numeric(SV_p_c[pos_to_get,3])<=as.numeric(RIL_table[b,6]))]  
if (RIL_table[b,4]=="P1"){parent<-P1} else if (RIL_table[b,4]=="P2") {parent<-P2}
SV_table[pos_to_get,4]<-parent
SV_table[pos_to_get,5]<-SV_p_c[pos_to_get, parent]  
}#b

#determine SV info for breakpoint regions
#podria hacer un loop donde defino breakpoint por breakpoint, y algo especial para la primera y ultima region
#pero por ahora digo que toda la region que no se es promedio de los dos parents
pos_to_get<-row.names(SV_table)[which(is.na(SV_table[,5]))]  
SV_table[pos_to_get,4]<-"P1_P2"
SV_table[pos_to_get,5]<-(as.numeric(SV_p_c[pos_to_get,P1])+as.numeric(SV_p_c[pos_to_get,P2]))/2

#only keep positions where SVs are present
SV_table<-SV_table[which(as.numeric(SV_table[,5])>0),]

sv_layers[[paste(l)]][[p]][[r]][[c]]<-SV_table

}#c 
}#r
}#p 
cat(fill = TRUE)
}#l
cat(fill = TRUE)

saveRDS(sv_layers, paste("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_layers/A_",v,"_layers.RDS"))

}#v
