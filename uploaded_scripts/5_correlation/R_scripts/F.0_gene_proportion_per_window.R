make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

s<-1000000
rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))

#get barley genes
barley_genes<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Barley_Morex_V2_gene_annotation_PGSB.all.descriptions.csv"))
chrs<-unique(barley_genes[,3])[-8]
genes_per_window<-list()
for (c in 1:7){ cat(c); cat("-")
  windows<-row.names(rec_list[[c]])
  gene_prop<-matrix(nrow = length(windows), ncol = 1); row.names(gene_prop)<-windows
  chr.table<-barley_genes[which(barley_genes[,3]==chrs[c]),5:6]
  for (w in 1:length(windows)){ #FALTA HACER: no cuento los genes que empiezan en una y terminan en otra
    pos_to_get<-chr.table[which(as.numeric(chr.table[,1])>=(as.numeric(windows[w])-s)),]; pos_to_get<-make.matrix(pos_to_get)
    pos_to_get<-pos_to_get[which(as.numeric(pos_to_get[,2])<=as.numeric(windows[w])),]; pos_to_get<-make.matrix(pos_to_get)
    if (nrow(pos_to_get)!=0){
      #take ranges positions
      WIN<-apply(pos_to_get, 1, FUN = function(x){return(as.numeric(x[1]):as.numeric(x[2]))})
      #count positions
      WIN<-length(unique(unlist(WIN)))
      #calculate percentage
      gene_prop[w,]<-WIN/s
    } else {gene_prop[w,]<-0}
  }#w
  genes_per_window[[c]]<-gene_prop
}#c

saveRDS(genes_per_window, paste("/scratch/Federico/3_RILs/5_correlation/Results/F.0_gene_proportion_per_window=",s,"_all_pops.RDS", sep = ""))

