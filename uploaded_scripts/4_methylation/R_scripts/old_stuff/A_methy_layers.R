make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)

dir.create("/scratch/Federico/3_RILs/3_Methys/Results/A_Methys_layers")

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

layers<-c("1_first", "2_second", "3_third")

met_context<-c("chg", "chh", "cpg")

parents.code<-as.matrix(read.table("/scratch/Federico/Paper_2/samplenames.txt"))

dir.create("/scratch/Federico/3_RILs/3_Methys/Results/A_Methys_layers/")
  
for (l in layers){ cat(l); cat(" layer"); cat(" - "); 
#l<-layers[1]
rec_layer<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.",l,"_layer.RDS", sep = ""))
Populations<-names(rec_layer)  

dir.create(paste("/scratch/Federico/3_RILs/3_Methys/Results/A_Methys_layers/",l, sep = ""))

for (p in Populations){ cat(p); cat(": ")
#p<-Populations[1]
rils<-names(rec_layer[[p]])  

#select parents
P1<-pop_info[which(pop_info[,1]==p),2]
P2<-pop_info[which(pop_info[,1]==p),3]

P1c<-parents.code[which(parents.code[,2]==P1),1]
P2c<-parents.code[which(parents.code[,2]==P2),1]

dir.create(paste("/scratch/Federico/3_RILs/3_Methys/Results/A_Methys_layers/",l,"/",p, sep = ""))

for (r in rils){
#r<-rils[1]

dir.create(paste("/scratch/Federico/3_RILs/3_Methys/Results/A_Methys_layers/",l,"/",p,"/",r, sep = ""))

for (c in 1:7) { 
#c<-1

chr.length<-as.numeric(chr_lengths[c,3])  

#get rec data  
RIL_table<-rec_layer[[p]][[r]][[c]]

for (v in met_context){ cat(v, fill = TRUE)
  
#get Methy data for that chr
Methy_P1<-as.matrix(fread(paste("/2data/Barley/Barley_Bisulfitseq_HvDRRparents_2020/experiments_mariusk/snp_corrected/7_home_input/",P1c,"_chr",c,"H_",v,".csv", sep = "")))
Methy_P1[,3]<-as.numeric(Methy_P1[,5])/as.numeric(Methy_P1[,6])
Methy_P1<-Methy_P1[,1:3]

Methy_P2<-as.matrix(fread(paste("/2data/Barley/Barley_Bisulfitseq_HvDRRparents_2020/experiments_mariusk/snp_corrected/7_home_input/",P2c,"_chr",c,"H_",v,".csv", sep = "")))
Methy_P2[,3]<-as.numeric(Methy_P2[,5])/as.numeric(Methy_P2[,6])
Methy_P2<-Methy_P2[,1:3]

pos<-unique(c(as.numeric(Methy_P1[,2]),as.numeric(Methy_P2[,2])))  
pos<-pos[order(pos)]
Methy_p_c<-matrix(ncol = 3, nrow = length(pos))
colnames(Methy_p_c)<-c("pos", P1, P2)
Methy_p_c[,1]<-pos
Methy_p_c[,2:3]<-0
row.names(Methy_p_c)<-paste("pos_",pos, sep = "")
row.names(Methy_P1)<-paste("pos_",as.numeric(Methy_P1[,2]), sep = "")
row.names(Methy_P2)<-paste("pos_",as.numeric(Methy_P2[,2]), sep = "")
Methy_p_c[row.names(Methy_P1),2]<-as.numeric(Methy_P1[,3])
Methy_p_c[row.names(Methy_P2),3]<-as.numeric(Methy_P2[,3])

#Determine RIL methylation
Methy_table<-Methy_p_c; Methy_table[,2:3]<-NA; colnames(Methy_table)[2:3]<-c("parent","value")

#determine Methy info for parental blocks
for (b in 1:nrow(RIL_table)){
pos_to_get<-row.names(Methy_p_c)[which(as.numeric(Methy_p_c[,1])>=as.numeric(RIL_table[b,5]))]  
pos_to_get<-pos_to_get[which(as.numeric(Methy_p_c[pos_to_get,1])<=as.numeric(RIL_table[b,6]))]  
if (RIL_table[b,4]=="P1"){parent<-P1} else if (RIL_table[b,4]=="P2") {parent<-P2}
Methy_table[pos_to_get,2]<-parentcat(fill = TRUE)

Methy_table[pos_to_get,3]<-Methy_p_c[pos_to_get, parent]  
}#b

#determine Methy info for breakpoint regions
#podria hacer un loop donde defino breakpoint por breakpoint, y algo especial para la primera y ultima region
#pero por ahora digo que toda la region que no se es promedio de los dos parents
pos_to_get<-row.names(Methy_table)[which(is.na(Methy_table[,3]))]  
Methy_table[pos_to_get,2]<-"P1_P2"
Methy_table[pos_to_get,3]<-(as.numeric(Methy_p_c[pos_to_get,P1])+as.numeric(Methy_p_c[pos_to_get,P2]))/2

#only keep positions where Methys are present
Methy_table<-Methy_table[which(as.numeric(Methy_table[,3])>0),]

#save
fwrite(Methy_table, paste("/scratch/Federico/3_RILs/3_Methys/Results/A_Methys_layers/",l,"/",p,"/",r,"/","chr",c,"H_",v,"_methy_layers.csv"), row.names = FALSE)
rm(Methy_table)

}#v

}#c 
}#r
cat(fill = TRUE)
}#p 
cat(fill = TRUE)
}#l


