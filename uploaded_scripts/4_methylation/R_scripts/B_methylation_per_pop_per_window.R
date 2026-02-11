library("data.table")

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[,1]

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
chrs<-paste("chr",1:7,"H", sep = "")

#parents.code<-as.matrix(read.table("/scratch/Federico/Paper_2/samplenames.txt"))

window_sizes<-c(10000, 500000, 1000000)[1]

met_context<-c("chg", "cpg", "chh")

s<-window_sizes
###############################################

#1- Methylation per pop (average of both parents per window) for each met context

#for (s in window_sizes){ cat(s, fill = TRUE)

if (s == 10000){
  rec_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.0_accumulated_rec_prob.RDS")
  Populations<-Populations[which(Populations%in%c("HvDRR13","HvDRR27", "HvDRR28"))]
  #parents<-parents[which(parents%in%c("HOR8160","SprattArcher", "Unumli-Arpa"))]
} else {
  #rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
  rec_list<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1.1_CO_prob_per_window.RDS", sep = ""))
  rec_list<-rec_list[[paste(s)]]
}

methy_list_RAW<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/A_methylation_per_parent_level_and_Cs_win=",s,".RDS", sep = ""))

methy_list<-list()

for (c in 1:7) { cat(c); cat("-")

methy_list[[c]]<-list()

for (v in met_context){ cat(v); cat("-")

methy_list[[c]][[v]]<-list()

methy_data<-methy_list_RAW[[c]][[v]][["methy_data"]]
Cs_counted<-methy_list_RAW[[c]][[v]][["Cs_counted"]]

pop_methy_data<-matrix(nrow = nrow(methy_data), ncol = length(Populations))
colnames(pop_methy_data)<-Populations
row.names(pop_methy_data)<-row.names(methy_data)

for (p in Populations){ cat(p); cat("-")

#select parents
P1<-pop_info[which(pop_info[,1]==p),2]
P2<-pop_info[which(pop_info[,1]==p),3]
#parents_methy_data<-methy_data[,c(P1,P2)]  #VIEJO -> SACAR
#pop_methy_data[,p]<-apply(parents_methy_data,1,function(x){return(mean(as.numeric(x[1:2])))}) #VIEJO -> SACAR

#I will weight values for the quantity of Cs measured per windows
# windows<-row.names(methy_data)
# for (i in windows){
# P1_methy<-as.numeric(methy_data[i,P1])*(Cs_counted[i,P1])   
# P2_methy<-as.numeric(methy_data[i,P2])*(Cs_counted[i,P2])
# pop_methy_data[i,p]<-sum(P1_methy, P2_methy)/sum(Cs_counted[i,P1], Cs_counted[i,P2])
# }#i
#el loop tarda
P1_methy<-methy_data[,P1]*Cs_counted[,P1]   
P2_methy<-methy_data[,P2]*Cs_counted[,P2]
P1_P2_methy<-P1_methy+P2_methy
win_Cs_sum<-Cs_counted[,P1]+Cs_counted[,P2]
pop_methy_data[,p]<-P1_P2_methy/win_Cs_sum

}#p

methy_list[[c]][[v]]<-pop_methy_data

cat(fill = TRUE) 

}#v

cat(fill = TRUE)

}#c

saveRDS(methy_list, paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.1_methylation_per_population_win=",s,".RDS", sep = ""))

#}#s

#################################################################################

#2 - met context unified in one index

Populations<-pop_info[,1]

#for (s in window_sizes){ cat(s, fill = TRUE)

if (s == 10000){
  Populations<-Populations[which(Populations%in%c("HvDRR13","HvDRR27", "HvDRR28"))]
}

methy_list_uni<-list()

for (c in 1:7){
  
    pop_methy_data<-matrix(nrow = nrow(methy_list_RAW[[c]][[1]][[1]]), ncol = length(Populations))
    colnames(pop_methy_data)<-Populations
    row.names(pop_methy_data)<-row.names(methy_list_RAW[[c]][[1]][[1]])
  
    for (p in Populations){
      
      P1<-pop_info[which(pop_info[,1]==p),2]
      P2<-pop_info[which(pop_info[,1]==p),3]
      
      #cpg
      Cs_cpg_P1<-methy_list_RAW[[c]][["cpg"]][["Cs_counted"]][,P1]
      Cs_cpg_P2<-methy_list_RAW[[c]][["cpg"]][["Cs_counted"]][,P2]
      P1_methy_cpg<-methy_list_RAW[[c]][["cpg"]][["methy_data"]][,P1]*Cs_cpg_P1  
      P2_methy_cpg<-methy_list_RAW[[c]][["cpg"]][["methy_data"]][,P2]*Cs_cpg_P2
      #chg
      Cs_chg_P1<-methy_list_RAW[[c]][["chg"]][["Cs_counted"]][,P1]
      Cs_chg_P2<-methy_list_RAW[[c]][["chg"]][["Cs_counted"]][,P2]
      P1_methy_chg<-methy_list_RAW[[c]][["chg"]][["methy_data"]][,P1]*Cs_chg_P1  
      P2_methy_chg<-methy_list_RAW[[c]][["chg"]][["methy_data"]][,P2]*Cs_chg_P2
      #chh
      Cs_chh_P1<-methy_list_RAW[[c]][["chh"]][["Cs_counted"]][,P1]
      Cs_chh_P2<-methy_list_RAW[[c]][["chh"]][["Cs_counted"]][,P2]
      P1_methy_chh<-methy_list_RAW[[c]][["chh"]][["methy_data"]][,P1]*Cs_chh_P1  
      P2_methy_chh<-methy_list_RAW[[c]][["chh"]][["methy_data"]][,P2]*Cs_chh_P2
      
      P1_P2_methy<-P1_methy_cpg+P2_methy_cpg+P1_methy_chg+P2_methy_chg+P1_methy_chh+P2_methy_chh
      win_Cs_sum<-Cs_cpg_P1+Cs_cpg_P2+Cs_chg_P1+Cs_chg_P2+Cs_chh_P1+Cs_chh_P2
      
      pop_methy_data[,p]<-P1_P2_methy/win_Cs_sum
      
    }#p
  methy_list_uni[[c]]<-pop_methy_data
}#c
saveRDS(methy_list_uni, paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.2.1_methy_unified_per_pop_win=",s,".RDS", sep = ""))
#}#s

###########################################

methy_list_RAW<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/A_methylation_per_parent_level_and_Cs_win=",s,".RDS", sep = ""))

#2 - met context unified in one index

Populations<-pop_info[,1]

#for (s in window_sizes){ cat(s, fill = TRUE)

if (s == 10000){
  Populations<-Populations[which(Populations%in%c("HvDRR13","HvDRR27", "HvDRR28"))]
}

methy_list_uni<-list()

for (c in 1:7){
  
  pop_methy_data<-matrix(nrow = nrow(methy_list_RAW[[c]][[1]][[1]]), ncol = length(Populations))
  colnames(pop_methy_data)<-Populations
  row.names(pop_methy_data)<-row.names(methy_list_RAW[[c]][[1]][[1]])
  
  for (p in Populations){
    
    P1<-pop_info[which(pop_info[,1]==p),2]
    P2<-pop_info[which(pop_info[,1]==p),3]
    
    #cpg
    Cs_cpg_P1<-methy_list_RAW[[c]][["cpg"]][["Cs_counted"]][,P1]
    Cs_cpg_P2<-methy_list_RAW[[c]][["cpg"]][["Cs_counted"]][,P2]
    P1_methy_cpg<-methy_list_RAW[[c]][["cpg"]][["methy_data"]][,P1]*Cs_cpg_P1  
    P2_methy_cpg<-methy_list_RAW[[c]][["cpg"]][["methy_data"]][,P2]*Cs_cpg_P2
    #chg
    Cs_chg_P1<-methy_list_RAW[[c]][["chg"]][["Cs_counted"]][,P1]
    Cs_chg_P2<-methy_list_RAW[[c]][["chg"]][["Cs_counted"]][,P2]
    P1_methy_chg<-methy_list_RAW[[c]][["chg"]][["methy_data"]][,P1]*Cs_chg_P1  
    P2_methy_chg<-methy_list_RAW[[c]][["chg"]][["methy_data"]][,P2]*Cs_chg_P2
    
    P1_P2_methy<-P1_methy_cpg+P2_methy_cpg+P1_methy_chg+P2_methy_chg
    win_Cs_sum<-Cs_cpg_P1+Cs_cpg_P2+Cs_chg_P1+Cs_chg_P2
    
    pop_methy_data[,p]<-P1_P2_methy/win_Cs_sum
    
  }#p
  methy_list_uni[[c]]<-pop_methy_data
}#c
saveRDS(methy_list_uni, paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.2.2_methy_unified_CPG_CHG_ONLY_per_pop_win=",s,".RDS", sep = ""))
#}#s



###################################################################################

#check correlation

methy_data<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.1_methylation_per_population_win=10000.RDS")
methy_uni_data<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.2.1_methy_unified_per_pop_win=10000.RDS")
for (c in 1:7){methy_data[[c]][["unified"]]<-methy_uni_data[[c]]}
methy_uni_data2<-readRDS("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.2.2_methy_unified_CPG_CHG_ONLY_per_pop_win=10000.RDS")
for (c in 1:7){methy_data[[c]][["unified_cpg_chg"]]<-methy_uni_data2[[c]]}

#unify chrs
genome_methy_data<-list()
for (i in names(methy_data[[1]])){ 
genome_methy_data[[i]]<-c()
for (c in 1:7){genome_methy_data[[i]]<-c(genome_methy_data[[i]], methy_data[[c]][[i]])}
}#i
#results table
tabla_cor<-matrix(ncol = 5, nrow = 5); contexts<-c("cpg","chg","chh","unified", "unified_cpg_chg")
row.names(tabla_cor)<-contexts; colnames(tabla_cor)<-contexts
for (i in contexts){ for (j in contexts){
test<-cor.test(x = genome_methy_data[[i]], y = genome_methy_data[[j]], alternative = "two.sided", method = "pearson")
tabla_cor[i,j]<-paste(test$estimate, " P=", round(test$p.value, digits = 3), sep = "")  
}}

write.csv(tabla_cor, "/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.3_correlations_table.csv")

###################################################################################

# #################################################################################

#
# #3-DMRs per pop (average of both parents per window) for each met context
#
# parents.code<-as.matrix(read.table("/scratch/Federico/3_RILs/4_methylation/sources/samplenames.txt"))
# parents<-parents.code[,2]
# 
# chrs<-1:7; #chrs<-chrs[1]
# 
# #for (s in window_sizes){ cat(s, fill = TRUE)
# 
#   if (s == 10000){
#   rec_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.0_accumulated_rec_prob.RDS")
#   Populations<-Populations[which(Populations%in%c("HvDRR13","HvDRR27", "HvDRR28"))]
#   #parents<-parents[which(parents%in%c("HOR8160","SprattArcher", "Unumli-Arpa"))]
#   } else {
#   #rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
#   rec_list<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1.1_CO_prob_per_window.RDS", sep = ""))
#   rec_list<-rec_list[[paste(s)]]
#   }
# 
#   DMRs_list<-list()
# 
#   for (c in chrs) {  cat(c); cat(": ")
# 
#     DMRs_list[[c]]<-list()
# 
#     windows<-rec_list[[1]][[c]][,2]
# 
#     for (v in met_context){ cat(v); cat("-")
# 
#       DMRs_list[[c]][[v]]<-list()
# 
#       DMRs_chr<-as.matrix(read.csv(paste("/scratch/Federico/3_RILs/4_methylation/sources/barley_methylation_data/dmrs/DMRs_chr",c,"H_",v,".csv", sep = ""), sep = ""))
#       #correct colnames
#       for (i in 1:nrow(parents.code)){colnames(DMRs_chr)[which(colnames(DMRs_chr)==parents.code[i,1])]<-parents.code[i,2]}#i
# 
#       pop_DMRs_data<-matrix(nrow = length(windows), ncol = length(Populations))
#       colnames(pop_DMRs_data)<-Populations
#       row.names(pop_DMRs_data)<-windows
# 
#       for (w in windows){ #cat(which(windows==w)); cat("-")
# 
#         DMRs_chr_win<-make.matrix(DMRs_chr[which(as.numeric(DMRs_chr[,2])<=as.numeric(w)),])
#         DMRs_chr_win<-make.matrix(DMRs_chr_win[which(as.numeric(DMRs_chr_win[,1])>(as.numeric(w)-s)),])
#         #con esto dejamos afuera los que empiezan en una window y terminan en otra
# 
#         if (length(DMRs_chr_win)!=0){
# 
#           for (p in Populations){ #cat(p); cat(": ")
#             #select parents
#             P1<-pop_info[which(pop_info[,1]==p),2]
#             P2<-pop_info[which(pop_info[,1]==p),3]
#             parents_DMRs_data<-make.matrix(DMRs_chr_win[,c("n.dmc","len",P1,P2)])
#             parents_DMRs_data<-make.matrix(parents_DMRs_data[complete.cases(parents_DMRs_data),])
#             if (nrow(parents_DMRs_data)!=0){
#             #make average index considering dmcs per DMR and the length of DMR in the window
#             parents_DMRs_data<-cbind(parents_DMRs_data,NA); colnames(parents_DMRs_data)[5]<-"index"
#             for (i in 1:nrow(parents_DMRs_data)){
#             parent_mean<-mean(parents_DMRs_data[i,c(P1,P2)])
#             DMR_over_window<-parents_DMRs_data[i,"len"]/s
#             #Cs_over_DMR<-parents_DMRs_data[i,"len"]/4  #no se las c por DMR, asumo 1/4
#             #DMCs_over_Cs<-parents_DMRs_data[i,"n.dmc"]/Cs_over_DMR
#             #parents_DMRs_data[i,5]<-parent_mean*DMCs_over_Cs*DMR_over_window
#             DMCs_over_DMR<-parents_DMRs_data[i,"n.dmc"]/parents_DMRs_data[i,"len"]
#             parents_DMRs_data[i,5]<-parent_mean*DMCs_over_DMR*DMR_over_window
#             }#i
#             pop_DMRs_data[paste(w),p]<-round(sum(parents_DMRs_data[,5]), digits = 4)
#             }#if nrow!=0
#           }#p
# 
#         }#if !=0
# 
#       }#w
# 
#     DMRs_list[[c]][[v]]<-pop_DMRs_data
# 
#     #cat(fill = TRUE)
# 
#     }#v
# 
#     cat(fill = TRUE)
# 
#     #saveRDS(DMRs_list[[c]], paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.3_DMRs_per_pop_win=",s,"_",c,".RDS", sep = ""))
#     #DMRs_list[[c]]<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.3_DMRs_per_pop_win=",s,"_",c,".RDS", sep = ""))
# 
#   }#c
# 
#   saveRDS(DMRs_list, paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.3_DMRs_per_pop_win=",s,".RDS", sep = ""))
# 
# #}#s
# 
# #################################################################################
# 
# #make unified (ahora no necesito ponderar porque ya viene ponderado)
# #el win 10000, me da muchas windows sin DMR...
# DMRs_list_uni<-list()
#   for (c in 1:7){
#     chr_uni_table<-DMRs_list[[c]][[1]]
#     for (v in met_context[2:3]){
#       for (p in Populations){
#         for (w in 1:nrow(chr_uni_table)){
#           chr_uni_table[w,p]<-sum(c(chr_uni_table[w,p], DMRs_list[[c]][[v]][w,p]), na.rm = TRUE)
#         }#w
#       }#p
#     }#v
#     DMRs_list_uni[[c]]<-chr_uni_table
#   }#c
# saveRDS(DMRs_list_uni, paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4_DMRs_unified_per_pop_win=",s,".RDS", sep = ""))
# 
# #########################################################################################################################

#VIEJO
# #4 - hago uno unified de las DMRs  
# DMRs_list<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.3_DMRs_per_pop_win=",s,".RDS", sep = ""))
# rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
#   
# #1 - primero calculo el peso de cada DMR en la window
# #1.1 - la physical proportion
# DMRs_prop<-list()
# for (c in 1:7) {  cat(c); cat(": ")
#   windows<-rec_list[[1]][[c]][,2]
#   weight_table<-matrix(nrow = length(windows), ncol = 3); colnames(weight_table)<-met_context; row.names(weight_table)<-windows
#   for (v in met_context){ cat(v); cat("-")
#     DMRs_chr<-as.matrix(read.csv(paste("/scratch/Federico/3_RILs/4_methylation/sources/barley_methylation_data/dmrs/DMRs_chr",c,"H_",v,".csv", sep = ""), sep = "")) 
#     for (w in windows){ 
#       DMRs_chr_win<-make.matrix(DMRs_chr[which(as.numeric(DMRs_chr[,2])<=as.numeric(w)),])
#       DMRs_chr_win<-make.matrix(DMRs_chr_win[which(as.numeric(DMRs_chr_win[,1])>(as.numeric(w)-s)),])
#       DMRs_length<-sum(as.numeric(DMRs_chr_win[,2])-as.numeric(DMRs_chr_win[,1]))
#       weight_table[w,v]<-DMRs_length/s
#     }#w
#   }#v
#   DMRs_prop[[c]]<-weight_table
# }#c
# saveRDS(DMRs_prop, paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.1_DMRs_proportion_win=",s,".RDS", sep = ""))
#   
# #1.2 - la physical proportion over el total de physical proportion
# DMRs_weight<-DMRs_prop
# for (c in 1:7) {  cat(c); cat(": ")
#   windows<-rec_list[[1]][[c]][,2]
#   for (v in met_context){ 
#     for (w in windows){ 
#     DMRs_weight[[c]][w,v]<-as.numeric(DMRs_prop[[c]][w,v])/sum(as.numeric(DMRs_prop[[c]][w,]))       
#     }#w
#   }#v
# }#c
# saveRDS(DMRs_weight, paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.2_DMRs_weight_win=",s,".RDS", sep = ""))
# 
# #DMRs_weight<-readRDS(paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.2_DMRs_weight_win=",s,".RDS", sep = ""))
# 
# #2 - luego afecto el % de cada met context por su peso en el total (weighted average)
# DMRs_index<-list()
# 
# for (c in 1:7){
#   DMRs_index[[c]]<-
#   ((DMRs_list[[c]][[1]]*DMRs_weight[[c]][,1])+(DMRs_list[[c]][[2]]*DMRs_weight[[c]][,2])+(DMRs_list[[c]][[3]]*DMRs_weight[[c]][,3]))
#                }#c
# 
# saveRDS(DMRs_index, paste("/scratch/Federico/3_RILs/4_methylation/Results/B_methylation_per_population/B.4.3_DMRs_unified_per_pop_wi_",s,"_windows.RDS", sep = ""))
#   
# #################################################################################