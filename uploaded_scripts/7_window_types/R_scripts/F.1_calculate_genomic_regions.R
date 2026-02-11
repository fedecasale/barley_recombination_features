#chr length
chr_length<-read.csv("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede.csv")

#gene segments
exons<-read.csv("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/sources/Hv_Morex.pgsb.Jul2020_OnlyExonRow.gff3", sep = "\t", header = FALSE)
exons<-exons[,c(1,4:5,7)]
colnames(exons)<-c("chr","start","end","sense")
genes<-read.csv("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/sources/Hv_Morex.pgsb.Jul2020_OnlyGeneRow.gff3", sep = "\t", header = FALSE)
genes<-genes[,c(1,4:5,7)]
colnames(genes)<-c("chr","start","end","sense")

chrs<-unique(as.character(exons[,1]))[-8]


exons_list<-list()
for (c in chrs){
exons_list[[c]]<-exons[which(exons$chr==c),-c(1,4)]
}

# #by now I will use only exons with positive sense
# exons<-exons[which(exons$sense=="+"),]
# 
# #get TSS
# TSS<-list()
# for (c in chrs){
# TSS[[c]]<-exons[which(exons$chr==c),"start"]
# }
# 
# #get medium point of exons
# medium<-list()
# for (c in chrs){
# medium[[c]]<-round((exons[which(exons$chr==c),"start"]+exons[which(exons$chr==c),"end"])/2)
# }

#intergenic space (of all genes, including al senses)
intergenic<-list()
for (c in chrs){
genes_chr<-genes[which(genes$chr==c),]
intergenic_chr<-genes_chr[2:3]; intergenic_chr[]<-NA
intergenic_chr[1,1]<-1; intergenic_chr[1,2]<-genes_chr[1,"start"]
for (g in 1:(nrow(genes_chr)-1)){
intergenic_chr[g+1,1]<-genes_chr[g,3]+1
intergenic_chr[g+1,2]<-genes_chr[g+1,2]-1
}#g
intergenic_chr<-rbind(intergenic_chr,NA)
intergenic_chr[nrow(intergenic_chr),1]<-genes_chr[nrow(genes_chr),3]+1
intergenic_chr[nrow(intergenic_chr),2]<-chr_length$chr_end[which(chrs==c)]
intergenic[[c]]<-intergenic_chr
}#c

# #get proximal promoter (puede variar entre 50 y 1000 bp, voy a usar 500 bp)
# promoters<-list()
# for (c in chrs){
# medium_promoter<-intergenic[[c]][,2]-499
# medium_promoter<-medium_promoter[-length(medium_promoter)]  
# #by now I will use only promoters of genes with positive sense
# medium_promoter<-medium_promoter[which(genes$sense=="+")]
# promoters[[c]]<-medium_promoter
# }

#get promoter region (puede variar entre 50 y 1000 bp, voy a usar 500 bp)
promoters<-list()
for (c in chrs){ 
genes_chr<-genes[which(genes$chr==c),]
#by now I will use only promoters of genes with positive sense
genes_chr<-genes_chr[which(genes_chr$sense=="+"),]
promoter_chr<-cbind(genes_chr$start-500, genes_chr$start-1)
colnames(promoter_chr)<-c("start","end")
promoters[[c]]<-promoter_chr
}

#regions near genes (voy a usar 500 bp)
neighboring<-list()
for (c in chrs){
  promoter_chr<-intergenic[[c]]
  promoter_chr[,1]<-promoter_chr[,2]-500
  promoter_chr<-promoter_chr[-nrow(promoter_chr),]
  termination_chr<-intergenic[[c]]
  termination_chr[,2]<-termination_chr[,1]+499
  termination_chr<-termination_chr[-1,]
  neighboring_chr<-rbind(promoter_chr, termination_chr)
  neighboring[[c]]<-neighboring_chr
}

genomic_regions<-list()
genomic_regions[["exons"]]<-exons_list
genomic_regions[["promoters_plus_500bp"]]<-promoters
genomic_regions[["neighboring_500bp"]]<-neighboring
genomic_regions[["intergenic"]]<-intergenic

saveRDS(genomic_regions, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/F_overlap_genomic_regions/F.1_genomic_regions.RDS")

