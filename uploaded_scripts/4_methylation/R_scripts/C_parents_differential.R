make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[,1]

dir.create("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential")

window_sizes<-c(500000, 1000000)[2]

met_context<-c("chg", "chh", "cpg")

for (s in window_sizes){ cat(s, fill = TRUE)
  
  rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
  
  methy_list<-list()
  
  for (c in 1:7) {  
    
    methy_list[[c]]<-list()
    
    for (v in met_context){ cat(v, fill = TRUE)
      
      methy_list[[c]][[v]]<-list()
      
      methy_data<-read.csv(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A_methylation_per_parent_per_win=",s,"_",v,"_",c,".csv", sep = ""))
      row.names(methy_data)<-methy_data[,1]; methy_data<-methy_data[,-1]
      colnames(methy_data)[which(colnames(methy_data)=="Unumli.Arpa")]<-"Unumli-Arpa"
      colnames(methy_data)[which(colnames(methy_data)=="W23829.803911")]<-"W23829/803911"
      
      pop_methy_data<-matrix(nrow = nrow(methy_data), ncol = 45)
      colnames(pop_methy_data)<-Populations
      row.names(pop_methy_data)<-row.names(methy_data)
      
      for (p in Populations){ cat(p); cat(": ")
        
        #select parents
        P1<-pop_info[which(pop_info[,1]==p),2]
        P2<-pop_info[which(pop_info[,1]==p),3]
        parents_methy_data<-methy_data[,c(P1,P2)]
        
        parents_differential<-abs(as.numeric(parents_methy_data[,1])-as.numeric(parents_methy_data[,2]))
        
        for (i in 1:nrow(pop_methy_data)){pop_methy_data[i,p]<-parents_differential[i]}
        
      }#p
      
      methy_list[[c]][[v]]<-pop_methy_data
      
    }#v
  }#c
  
  saveRDS(methy_list, paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.1_methy_parents_differential_per_pop_win=",s,".RDS", sep = ""))
  
}#s


#######################

parents.code<-as.matrix(read.table("/scratch/Federico/3_RILs/4_methylation/sources/samplenames.txt"))
parents<-parents.code[,2]

#uno general con las differencias en las dmrs
for (s in window_sizes){ cat(s, fill = TRUE)
  
  rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
  
  DMRs_list<-list()
  
  for (c in 1:7) {  cat(c); cat(": ")
    
    DMRs_list[[c]]<-list()
    
    windows<-rec_list[[1]][[c]][,2]
    
    for (v in met_context){ cat(v); cat("-")
    
    DMRs_list[[c]][[v]]<-list()
      
    DMRs_chr<-as.matrix(read.csv(paste("/scratch/Federico/3_RILs/4_methylation/sources/barley_methylation_data/dmrs/DMRs_chr",c,"H_",v,".csv", sep = ""), sep = "")) 
    #correct colnames  
    for (i in 1:nrow(parents.code)){colnames(DMRs_chr)[which(colnames(DMRs_chr)==parents.code[i,1])]<-parents.code[i,2]}#i
      
    pop_DMRs_data<-matrix(nrow = length(windows), ncol = 45)
    colnames(pop_DMRs_data)<-Populations
    row.names(pop_DMRs_data)<-windows
      
    for (w in windows){ 
      
      DMRs_chr_win<-make.matrix(DMRs_chr[which(as.numeric(DMRs_chr[,2])<=as.numeric(w)),])
      DMRs_chr_win<-make.matrix(DMRs_chr_win[which(as.numeric(DMRs_chr_win[,1])>(as.numeric(w)-s)),])
      
      for (p in Populations){ #cat(p); cat(": ")
      
        #select parents
        P1<-pop_info[which(pop_info[,1]==p),2]
        P2<-pop_info[which(pop_info[,1]==p),3]
        parents_DMRs_data<-make.matrix(DMRs_chr_win[,c(P1,P2)])
        parents_DMRs_data<-make.matrix(parents_DMRs_data[complete.cases(parents_DMRs_data),])
        
        pop_DMRs_data[w,p]<-abs(sum(parents_DMRs_data[,1]-parents_DMRs_data[,2])/nrow(parents_DMRs_data))
        
      }#p
      
     }#w
      
    DMRs_list[[c]][[v]]<-pop_DMRs_data
      
   }#v

cat(fill = TRUE)    

}#c
  
saveRDS(DMRs_list, paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.2_DMRs_parents_differential_per_pop_win=",s,".RDS", sep = ""))
  
}#s

#######################
#hago uno unificado pero weighted por la proporcion de cada dmr type en una window over el total de dmrs en esa window

DMRs_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.2_DMRs_parents_differential_per_pop_win=",s,".RDS", sep = ""))
DMRs_weight<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.2_DMRs_weight_win=",s,".RDS", sep = ""))

DMRs_index<-list()
for (c in 1:7){
  DMRs_index[[c]]<-
    ((DMRs_list[[c]][[1]]*DMRs_weight[[c]][,1])+(DMRs_list[[c]][[2]]*DMRs_weight[[c]][,2])+(DMRs_list[[c]][[3]]*DMRs_weight[[c]][,3]))
}#c

saveRDS(DMRs_index, paste("/scratch/Federico/3_RILs/4_methylation/Results/C_parents_differential/C.3_DMRs_parents_differential_unified_per_pop_win=_",s,".RDS", sep = ""))


