for (c in 1:7){
writeLines(
paste(
"library(\"data.table\")

# make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
#   new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
#   new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

get_average<-function(x){return(round(mean(methy_data[which(positions%in%positions[which(positions>(x-s))][which(positions[which(positions>(x-s))]<x)])], na.rm = TRUE), digits = 3))}

chr_lengths<-as.matrix(read.csv(\"/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv\"))
#chrs<-paste(\"chr\",1:7,\"H\", sep = \"\")
chrs<-",c,"

met_context<-c(\"chg\", \"chh\", \"cpg\")

parents.code<-as.matrix(read.table(\"/scratch/Federico/Paper_2/samplenames.txt\"))
parents<-parents.code[,2]

dir.create(\"/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/\")

window_sizes<-c(10000, 500000, 1000000)[1]
met_context<-c(\"cpg\", \"chg\", \"chh\")[3]


#for (s in window_sizes){ cat(s, fill = TRUE)

s<-window_sizes

if (s == 10000){
  rec_list<-readRDS(\"/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.0_accumulated_rec_prob.RDS\")
  parents<-parents[which(parents%in%c(\"HOR8160\",\"SprattArcher\", \"Unumli-Arpa\"))]
} else {
  #rec_list<-readRDS(paste(\"/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_\",s,\".RDS\", sep = \"\"))
  rec_list<-readRDS(paste(\"/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1.1_CO_prob_per_window.RDS\", sep = \"\"))
  rec_list<-rec_list[[paste(s)]]
}

for (c in chrs) { cat(c); cat(\"-\")
  
  #methy_list[[c]]<-list()
  chr.length<-as.numeric(chr_lengths[c,3])  
  windows<-rec_list[[1]][[c]][,2]
  
  for (v in met_context){ cat(v, fill = TRUE)
    
    methy_table<-matrix(ncol = length(parents), nrow = length(windows))
    row.names(methy_table)<-windows
    colnames(methy_table)<-parents
    
    for (p in parents){ cat(p); cat(\": \")
      
      #get Methy data for that parent and chr            
      P1c<-parents.code[which(parents.code[,2]==p),1]
      Methy_P1<-as.matrix(fread(paste(\"/scratch/Federico/3_RILs/4_methylation/sources/barley_methylation_data/meth_state/\",P1c,\"_chr\",c,\"H_\",v,\".csv\", sep = \"\"))) 
      positions<-as.numeric(Methy_P1[,1])
      positions_first_win<-positions[which(positions<s)]
      positions_last_win<-positions[which(positions>windows[length(windows)-1])]
      #eliminate first window positions from data
      positions<-positions[-which(positions%in%positions_first_win)] #[which(positions>=10000)]
      #add the rest of the windows
      positions<-strsplit(as.character(positions), \"\")
      ceros<-length(which(unlist(strsplit(as.character(s),\"\"))==\"0\"))
      positions<-lapply(positions, FUN = function(x){return(x[-c((length(x)-(ceros-1)):length(x))])})
      positions<-unlist(lapply(positions, FUN = function(x){return(paste(x,collapse = \"\"))}))
      positions<-as.numeric(paste(positions, paste(rep(9, ceros), collapse = \"\"), sep = \"\"))
      methy_table[,p]<-unlist(lapply(windows, FUN = function(x){return(length(which(positions==x)))}))
      #add first window and last window
      if (length(positions_first_win)!=0){methy_table[1,p]<-length(positions_first_win)}
      if (length(positions_last_win)!=0){methy_table[nrow(methy_table),p]<-length(positions_last_win)}
      cat(fill = TRUE)
      
    }#p
    
    write.csv(methy_table, paste(\"/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/\",v,\"/A.2_n_measured_Cs_per_parent_per_win=\",s,\"_\",v,\"_\",c,\".csv\", sep = \"\"))
    
    #methy_list[[c]][[v]]<-methy_table
    
  }#v
}#c
#}#s
", sep = ""), 
paste("/scratch/Federico/3_RILs/4_methylation/R_scripts/A.2_methylation_per_parent_per_window_N_MEASURED_Cs_",c,".R", sep = "")
) 
}
