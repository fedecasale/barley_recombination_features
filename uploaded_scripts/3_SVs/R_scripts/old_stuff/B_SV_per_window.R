make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/B.0_SVs_per_windows")

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

window_sizes<-c(500000, 1000000, 5000000, 10000000)

SVs<-c("deletions", "duplications", "insertions") #, "translocations")
#SVs<-SVs[2]

#get SV per pop
sv_per_pop<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_per_population.RDS"))
Populations<-names(sv_per_pop[[1]]) 

rec_prob<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4_CO_prob_per_window_DENSITY.RDS", sep = ""))

for (v in SVs){ cat(v, fill = TRUE)
  
sv_per_window<-list()

for (s in window_sizes){ cat("window size "); cat(s); cat(": ")
#s<-window_sizes[4]
sv_per_window[[paste(s)]]<-list()  

for (p in Populations){ cat(p); cat(" - ")

sv_per_window[[paste(s)]][[p]]<-list()  

for (c in 1:7) {

#get windows from rec data
windows<-rec_prob[[paste(s)]][[1]][[p]][[c]]
windows[,3]<-NA

#get SV data for that chr
SV_p_c<-sv_per_pop[[v]][[p]][[c]]    

#add SVs count to windows  
for (w in 1:nrow(windows)){  
pos_to_get<-SV_p_c[which(as.numeric(SV_p_c[,2])>=as.numeric(windows[w,1])),]; pos_to_get<-make.matrix(pos_to_get)
pos_to_get<-pos_to_get[which(as.numeric(pos_to_get[,3])<=as.numeric(windows[w,2])),]; pos_to_get<-make.matrix(pos_to_get)
windows[w,3]<-nrow(pos_to_get)
}

#NO ADIERO LOS QUE EMPIEZAN EN UNA WINDOW Y TERMINAN EN OTRA POR AHORA, PERO ACA DEJO CODEE CHUNK QUE PODRIA SERVIR
# start_match_window<-row.names(windows)[which(as.numeric(windows[,1])<=as.numeric(RIL_table[b,5]))]
# start_match_window<-start_match_window[length(start_match_window)]
# end_match_window<-row.names(windows)[which(as.numeric(windows[,1])>=as.numeric(RIL_table[b,6]))]
# end_match_window<-end_match_window[1]
# 
# block_length<-as.numeric(RIL_table[b,8])
# 
# if (start_match_window==end_match_window) { #breakpoint in only one window
#   
# occupancy<-(block_length/s)*S #ESTABA ACA
# windows[start_match_window,3]<-as.numeric(windows[start_match_window,3])+occupancy
#   
# } else { #breakpoint is distributed in many windows
# 
#   bkp_length<-as.numeric(breakpoints[f,3])
#   #cuanto ocupa de la primera
#   occupancy<-as.numeric(windows[start_match_window,2])-as.numeric(breakpoints[f,1])
#   occupancy<-round(occupancy/bkp_length, digits = 2)
#   windows[start_match_window,3]<-as.numeric(windows[start_match_window,3])+occupancy
#   #cuanto ocupa de la ultima
#   occupancy<-as.numeric(breakpoints[f,2])-as.numeric(windows[end_match_window,1])
#   occupancy<-round(occupancy/bkp_length, digits = 2)
#   windows[end_match_window,3]<-as.numeric(windows[end_match_window,3])+occupancy
#   #si hay entire windows entre la primera y la ultima
#   if ((as.numeric(gsub("w_","",end_match_window))-as.numeric(gsub("w_","",start_match_window)))>1){
#     freq_per_window<-round(s/bkp_length, digits = 2)
#     start<-which(row.names(windows)==start_match_window)
#     end<-which(row.names(windows)==end_match_window)
#     windows[(start+1):(end-1),3]<-as.numeric(windows[(start+1):(end-1),3])+freq_per_window
#   }#if #si hay entire windows entre la primera y la ultima
# }#else #breakpoint is distributed in many windows

sv_per_window[[paste(s)]][[p]][[c]]<-windows 

}#c 
}#p  
cat(fill = TRUE)
}#s
cat(fill = TRUE)
saveRDS(sv_per_window, paste("/scratch/Federico/3_RILs/3_SVs/Results/B.0_SVs_per_windows/B_",v,"_per_window.RDS", sep = ""))
}#v
