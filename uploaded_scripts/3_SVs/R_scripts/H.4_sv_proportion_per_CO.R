make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

#################################################################################
win_size<-10000

SVs_proportion<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/H_SVs_proportion_per_window/H.1.2_SVs_proportion_per_window=10000.RDS")

CO_SV_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO/B.1.2_SV_per_CO_per_SV_size.RDS")
Populations<-names(CO_SV_list[[1]][[1]])

CO_SV_prop<-list()
for (p in Populations){ cat(p); cat(": ")
 CO_SV_prop[[p]]<-list()
 for (c in 1:7){ cat(c); cat("-")
  windows_end<-as.numeric(names(SVs_proportion[[p]][[c]]))
  windows_start<-windows_end-win_size+1
  COs_table<-CO_SV_list[[1]][[3]][[p]][[c]][[1]][,c(4:6)]
  for (r in 1:nrow(COs_table)){ #cat(r); cat("-")
   start<-as.numeric(COs_table[r,1]); end<-as.numeric(COs_table[r,2]); CO_length<-end-start
   start_win<-windows_end[which(windows_end>start)[1]]
   start_win_segment<-start_win-start
   end_win<-windows_start[which(windows_start<end)[length(which(windows_start<end))]]
   end_win_segment<-end-end_win
   if (sum(start_win_segment, end_win_segment)<CO_length){
   medium_wins<-windows_end[(which(windows_end==start_win)+1):(which(windows_start==end_win)-1)]   
   } else {medium_wins <- NULL}
   start_prop<-start_win_segment*SVs_proportion[[p]][[c]][[paste(start_win)]]
   end_prop<-end_win_segment*SVs_proportion[[p]][[c]][[paste(windows_end[which(windows_start==end_win)])]]
   medium_prop<-sum(win_size*unlist(SVs_proportion[[p]][[c]][paste(medium_wins)]))
   CO_SV_proportion<-(start_prop+medium_prop+end_prop)/CO_length
   COs_table[r,3]<-round(CO_SV_proportion, digits = 3)
  }#r
  colnames(COs_table)[3]<-"SVs_proportion"
  CO_SV_prop[[p]][[c]]<-COs_table 
 }#c   
 cat(fill = TRUE)
}#p

saveRDS(CO_SV_prop, "/scratch/Federico/3_RILs/3_SVs/Results/H.4_SVs_proportion_per_CO_interval.RDS")
#################################################################################