library("data.table")

# make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
#   new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
#   new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

get_average<-function(x){return(round(mean(methy_data[which(positions%in%positions[which(positions>(x-s))][which(positions[which(positions>(x-s))]<x)])], na.rm = TRUE), digits = 3))}

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

met_context<-c("chg", "chh", "cpg")

parents.code<-as.matrix(read.table("/scratch/Federico/Paper_2/samplenames.txt"))
parents<-parents.code[,2]

dir.create("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/")

window_sizes<-c(10000, 500000, 1000000)[1]
met_context<-c("cpg", "chg", "chh")[3]
chrs<-1:7; chrs<-chrs[4]

for (s in window_sizes){ cat(s, fill = TRUE)

if (s == 10000){
rec_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.0_accumulated_rec_prob.RDS")
parents<-parents[which(parents%in%c("HOR8160","SprattArcher", "Unumli-Arpa"))]
} else {
#rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
rec_list<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1.1_CO_prob_per_window.RDS", sep = ""))
rec_list<-rec_list[[paste(s)]]
}

for (c in chrs) { cat(c); cat("-")

#methy_list[[c]]<-list()
chr.length<-as.numeric(chr_lengths[c,3])  
windows<-rec_list[[1]][[c]][,2]

for (v in met_context){ cat(v, fill = TRUE)

methy_table<-matrix(ncol = length(parents), nrow = length(windows))
row.names(methy_table)<-windows
colnames(methy_table)<-parents
  
for (p in parents){ cat(p); cat(": ")

#get Methy data for that parent and chr            
P1c<-parents.code[which(parents.code[,2]==p),1]
Methy_P1<-as.matrix(fread(paste("/scratch/Federico/3_RILs/4_methylation/sources/barley_methylation_data/meth_state/",P1c,"_chr",c,"H_",v,".csv", sep = ""))) 
Methy_P1[,2]<-round(as.numeric(Methy_P1[,2])/as.numeric(Methy_P1[,3]), digits = 3)
Methy_P1<-Methy_P1[,1:2]
positions<-as.numeric(Methy_P1[,1])
methy_data<-as.numeric(Methy_P1[,2])

#porque las positions no son las mismas per parent, hago un average per window per parent

#method 1 tarda mucho
# get_average<-function(x){return(
#   round(mean(
#     methy_data[which(positions%in%positions[which(positions>(x-s))][which(positions[which(positions>(x-s))]<x)])]
#     , na.rm = TRUE), digits = 3)
#   )}
# 
# methy_table[,p]<-unlist(lapply(as.numeric(windows), FUN = get_average))

#method 2 tarda mucho
# for (i in 1:length(windows)){ 
# positions_in_win<-which(positions<as.numeric(windows[i]))  
# methy_table[paste(windows[i]),p]<-round(mean(methy_data[positions_in_win], na.rm = TRUE), digits = 3)
# methy_data<-methy_data[-positions_in_win]  
# positions<-positions[-positions_in_win]
# 

#method 3 tarda ok
# #esto esta OK para 1 Mbp, pero ojo que para 10 kb, la last position de positions<x, podria estar en windows anteriores
# get_last_position<-function(x){pos<-which(positions<x); return(pos[length(pos)])} 
# last_position_list<-lapply(as.numeric(windows), FUN = get_last_position); names(last_position_list)<-windows
# methy_table[1,p]<-round(mean(methy_data[1:last_position_list[[1]]], na.rm = TRUE), digits = 3)
# 
# for (i in 2:length(windows)){ 
# positions_in_win<-(last_position_list[[i-1]]+1):last_position_list[[i]]  
# methy_table[paste(windows[i]),p]<-round(mean(methy_data[positions_in_win], na.rm = TRUE), digits = 3)
# }#i

#method 4, lo mas repido que se me ocurrio y ademas, mas adecuado para 10 kb
#first window
positions_first_win<-positions[which(positions<s)]
if (length(positions_first_win)!=0){methy_first_win<-round(mean(methy_data[which(positions%in%positions_first_win)], na.rm = TRUE), digits = 3)}
positions_last_win<-positions[which(positions>windows[length(windows)-1])]
if (length(positions_last_win)!=0){methy_last_win<-round(mean(methy_data[which(positions%in%positions_last_win)], na.rm = TRUE), digits = 3)}

#eliminate first window positions from data
methy_data<-methy_data[-which(positions%in%positions_first_win)]
positions<-positions[-which(positions%in%positions_first_win)] #[which(positions>=10000)]
#add the rest of the windows
positions<-strsplit(as.character(positions), "")
ceros<-length(which(unlist(strsplit(as.character(s),""))=="0"))
positions<-lapply(positions, FUN = function(x){return(x[-c((length(x)-(ceros-1)):length(x))])})
positions<-unlist(lapply(positions, FUN = function(x){return(paste(x,collapse = ""))}))
positions<-as.numeric(paste(positions, paste(rep(9, ceros), collapse = ""), sep = ""))
methy_table[,p]<-unlist(lapply(windows, FUN = function(x){return(round(mean(methy_data[which(positions==x)], na.rm = TRUE), digits = 3))}))
#add first window and last window
if (length(positions_first_win)!=0){methy_table[1,p]<-methy_first_win}
if (length(positions_last_win)!=0){methy_table[nrow(methy_table),p]<-methy_last_win}
cat(fill = TRUE)
}#p

write.csv(methy_table, paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A_methylation_per_parent_per_win=",s,"_",v,"_",c,".csv", sep = ""))

#methy_list[[c]][[v]]<-methy_table

}#v
}#c
}#s

# v<-"chh"
# s<-1000000
# c<-2
# for (c in 2:7){
# methy_table<-read.csv(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/A_methylation_per_parent_per_win=",s,"_",v,"_",c,"_8.csv", sep = ""))
# row.names(methy_table)<-methy_table[,1]
# methy_table<-methy_table[,-1]
# for (u in c(16,23)){
# methy_table2<-read.csv(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/A_methylation_per_parent_per_win=",s,"_",v,"_",c,"_",u,".csv", sep = ""))[,-1]
# methy_table<-cbind(methy_table, methy_table2)
# }#u
# write.csv(methy_table, paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/chh/A_methylation_per_parent_per_win=",s,"_",v,"_",c,".csv", sep = ""))
# }#c

# met_context<-met_context[1:2]
# for (s in window_sizes){ cat(s, fill = TRUE)
# 
#   if (s == 10000){
#     rec_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.0_accumulated_rec_prob.RDS")
#     parents<-parents[which(parents%in%c("HOR8160","SprattArcher", "Unumli-Arpa"))]
#   } else {
#     #rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
#     rec_list<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1.1_CO_prob_per_window.RDS", sep = ""))
#     rec_list<-rec_list[[paste(s)]]
#   }
# 
#   for (c in chrs) { cat(c); cat("-")
# 
#     #methy_list[[c]]<-list()
#     chr.length<-as.numeric(chr_lengths[c,3])
#     windows<-rec_list[[1]][[c]][,2]
# 
#     for (v in met_context){ cat(v, fill = TRUE)
# 
#       methy_table<-read.csv(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A_methylation_per_parent_per_win=",s,"_",v,"_",c,".csv", sep = ""))
#       colnames(methy_table)[which(colnames(methy_table)=="Unumli.Arpa")]<-"Unumli-Arpa"
#       for (p in parents){ cat(p); cat(": ")
# 
#         #get Methy data for that parent and chr
#         P1c<-parents.code[which(parents.code[,2]==p),1]
#         Methy_P1<-as.matrix(fread(paste("/scratch/Federico/3_RILs/4_methylation/sources/barley_methylation_data/meth_state/",P1c,"_chr",c,"H_",v,".csv", sep = "")))
#         Methy_P1[,2]<-round(as.numeric(Methy_P1[,2])/as.numeric(Methy_P1[,3]), digits = 3)
#         Methy_P1<-Methy_P1[,1:2]
#         positions<-as.numeric(Methy_P1[,1])
#         methy_data<-as.numeric(Methy_P1[,2])
# 
# 
#         #method 4, lo mas repido que se me ocurrio y ademas, mas adecuado para 10 kb
#         #add first window
#         positions_first_win<-positions[which(positions<s)]
#         if (length(positions_first_win)!=0){methy_first_win<-round(mean(methy_data[which(positions%in%positions_first_win)], na.rm = TRUE), digits = 3)}
#         positions_last_win<-positions[which(positions>windows[length(windows)-1])]
#         if (length(positions_last_win)!=0){methy_last_win<-round(mean(methy_data[which(positions%in%positions_last_win)], na.rm = TRUE), digits = 3)}
#         if (length(positions_first_win)!=0){methy_table[1,p]<-methy_first_win}
#         if (length(positions_last_win)!=0){methy_table[nrow(methy_table),p]<-methy_last_win}
# 
#         cat(fill = TRUE)
#       }#p
# 
#       write.csv(methy_table, paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A_methylation_per_parent_per_win=",s,"_",v,"_",c,".csv", sep = ""), row.names = FALSE)
# 
# 
#     }#v
#   }#c
#   #saveRDS(methy_list, paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/A_methylation_per_parent_per_win=",s,"_",v,".RDS", sep = ""))
# }#s
