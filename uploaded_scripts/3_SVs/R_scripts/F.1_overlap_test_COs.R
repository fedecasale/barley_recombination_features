source("/scratch/Federico/3_RILs/2_CO_breakpoints/R_scripts/Z_my_functions.R")

library(regioneR, verbose = FALSE)

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.1_SVs_per_population.RDS")
breakpoint_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")
centromere <-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede.csv")

pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere.RDS")

layers<-c("1_first", "2_second", "3_third")
Populations<-names(breakpoint_list[[1]])
chrs<-paste("chr", 1:7, sep = "")  #clave sin la H para correr el permTest

SVs<-c("deletions", "duplications", "insertions", "inversions") #"translocations")
SVs<-SVs[4]

# #randomization function based on rec. rate
# accumulated_rec_prob<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.3_CO_prob_per_window_GENETIC_MAPS.RDS", sep = ""))
# 
# random.start<-function(x){
#   random_win<-sample(x = 1:N.windows, size = 1, prob = PROB.windows, replace = TRUE)
#   random_start<-sample(x = rec_prob[random_win,1]:rec_prob[random_win,2], size = 1)
#   return(random_start)
# } 
# 
# 
# randomize_segments<-function(x){
# rec_prob<-accumulated_rec_prob$`5e+06`[[p]][[c]]
# N.windows<-nrow(rec_prob)
# PROB.windows<-as.numeric(rec_prob[,3])
# #
# tabla<-as.data.frame(x)
# random_start_positions<-unlist(lapply(1:nrow(tabla), FUN = random.start))
# tabla[,2]<-random_start_positions
# tabla[,3]<-(random_start_positions+tabla[,4])-1 
# return(GRanges(tabla))
# }
# 
# my_randomize_function<-function(A,B){return(list(randomize_segments(A), randomize_segments(B)))}


#delete small SVs (indels most likely)
for (v in SVs){
for (p in Populations){ 
  for (c in 1:7){ 
    tabla<-sv_list[[v]][[p]][[c]]  
    sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
    if (any(sv_lengths<=49)){tabla<-tabla[-which(sv_lengths<=49),]}; tabla<-make.matrix(tabla)
    sv_list[[v]][[p]][[c]]<-tabla 
    }#p
}#c
}#v

#dir.create("/scratch/Federico/3_RILs/3_SVs/Results/F.1_overlap_test_COs")
pT_list<-list()

for (v in SVs){ cat(v, fill = TRUE)
    
    l<-layers[1] #I use only layer one COs
  
    tabla<-matrix(nrow = 4, ncol = 3); colnames(tabla)<-Populations
    row.names(tabla)<-c("total CO intervals", "total SVs", "overlaping CO-SVs", "P value_OVERLAP")
    
    pT_list[[v]]<-list()
    overlaps<-c()
    
    for (p in Populations){ cat(p);cat(": ")
      
      centromere<-pericentromere_pops[[p]]
      centromere[,1]<-chrs #needs format for permut function
      genome<-centromere[,c(1:3)]; colnames(genome)<-c("seqnames","start","end"); genome<-GRanges(as.data.frame(genome))
      pericentromere<-centromere[,c(1,4:5)]; colnames(pericentromere)<-c("seqnames","start","end"); pericentromere<-GRanges(as.data.frame(pericentromere))
      
      CO_regions_genome<-data.frame(); SV_regions_genome<-data.frame()
        
        for (c in 1:7){ cat(c);cat("-")

        chr<-chrs[c]
        strand<-"+"  #dont know if OK
        
        CO_regions<-breakpoint_list[[l]][[p]][[c]]
        starts<-as.numeric(CO_regions[,which(colnames(CO_regions)=="starts")])+1
        ends<-as.numeric(CO_regions[,which(colnames(CO_regions)=="ends")])-1
        widths<-(ends-starts)+1
        CO_regions <- data.frame(chr, starts, ends, widths, strand)
        colnames(CO_regions)<-c("seqnames","start","end","width","strand")
        #take out peri
        to.delete<-starts[which(starts>centromere[c,4])]; to.delete<-to.delete[which(to.delete<centromere[c,5])]
        to.delete<-which(starts%in%to.delete)
        to.delete2<-ends[which(ends>centromere[c,4])]; to.delete2<-to.delete2[which(to.delete2<centromere[c,5])]
        to.delete2<-which(ends%in%to.delete2)
        to.delete<-unique(to.delete, to.delete2)
        if (is.na(to.delete)==FALSE){CO_regions<-CO_regions[-to.delete,]}
        #
        #tabla[1,c]<-nrow(CO_regions)
        CO_regions_genome<-rbind(CO_regions_genome, CO_regions)
        #CO_regions<-GRanges(CO_regions)
        
        
        SV_regions<-sv_list[[v]][[p]][[c]][,1:4]
        starts<-as.numeric(SV_regions[,2])   #clave poner el c
        ends<-as.numeric(SV_regions[,3])
        widths<-(ends-starts)+1
        SV_regions <- data.frame(chr, starts, ends, widths, strand)
        colnames(SV_regions)<-c("seqnames","start","end","width","strand")
        #take out peri
        to.delete<-starts[which(starts>centromere[c,4])]; to.delete<-to.delete[which(to.delete<centromere[c,5])]
        to.delete<-which(starts%in%to.delete)
        to.delete2<-ends[which(ends>centromere[c,4])]; to.delete2<-to.delete2[which(to.delete2<centromere[c,5])]
        to.delete2<-which(ends%in%to.delete2)
        to.delete<-unique(to.delete, to.delete2)
        if (length(to.delete)!=0){SV_regions<-SV_regions[-to.delete,]}
        #tabla[2,c]<-nrow(SV_regions)
        SV_regions_genome<-rbind(SV_regions_genome, SV_regions)
        #SV_regions<-GRanges(SV_regions)
        
        #el real overlap
        overlaps<-c(overlaps,numOverlaps(CO_regions, SV_regions))
        #tabla[3,c]<-meanDistance(SV_regions, CO_regions)
        }#c
        
        tabla[1,p]<-nrow(CO_regions_genome)
        tabla[2,p]<-nrow(SV_regions_genome)
        tabla[3,p]<-sum(overlaps)
        
        CO_regions_genome<-GRanges(CO_regions_genome)
        SV_regions_genome<-GRanges(SV_regions_genome)
        
        #pero es este numero grande o pequeno, las SV estan asociadas con los CO, o es just chance -> permutation test
        #alternative greater means que el randomizado overlap mas veces -> si es significant es que no estan asociados.
        
        #pT<-permTest(A = SV_regions, B = CO_regions, ntimes=100, alternative = "greater", 
        #             randomize.function = randomizeRegions, genome = genome, mask = pericentromere, per.chromosome = TRUE,
        #             evaluate.function = meanDistance, non.overlapping = FALSE)
        
        #pT<-permTest(A = SV_regions_genome, B = CO_regions_genome, ntimes=1000, alternative = "greater", 
        #             randomize.function = randomizeRegions, genome = genome, mask = pericentromere, per.chromosome = TRUE,
        #             evaluate.function = meanDistance, non.overlapping = FALSE)
        
        #tabla[4,c]<-pT$meanDistance$pval
        
        #pT<-permTest(A = SV_regions, B = CO_regions, ntimes=100, alternative = "less", 
        #            randomize.function = randomizeRegions, genome = genome, mask = pericentromere, per.chromosome = TRUE,
        #            evaluate.function = numOverlaps, non.overlapping = FALSE)
        
        pT<-permTest(A = SV_regions_genome, B = CO_regions_genome, ntimes=1000, alternative = "less", 
                     randomize.function = randomizeRegions, genome = genome, mask = pericentromere, per.chromosome = TRUE,
                     evaluate.function = numOverlaps, non.overlapping = FALSE)
        
        tabla[4,p]<-pT$numOverlaps$pval
        
        pT_list[[v]][[p]]<-pT
      
        #}#c
    
    #write.csv(tabla, paste("/scratch/Federico/3_RILs/3_SVs/Results/F.1_overlap_test_COs/F.1.1_overlap_test_",p,"_",v,".csv", sep = ""))
        
    }#p
  
  tabla[]<-round(tabla[], digits = 4)  
  write.csv(tabla, paste("/scratch/Federico/3_RILs/3_SVs/Results/F.1_overlap_test_COs/F.1.2_overlap_test_",v,".csv", sep = ""), row.names = FALSE)
  
  saveRDS(pT_list[[v]], paste("/scratch/Federico/3_RILs/3_SVs/Results/F.1_overlap_test_COs/F.1.3_permutation_list_test_",v,".RDS", sep = ""))
  
}#v



########## re-unify ##########

#graphic_list<-list()
tabla<-matrix(ncol = 4, nrow = 0)
colnames(tabla)<-c("SV", Populations)

for (v in SVs){ cat(v, fill = TRUE)
  
  #graphic_list[[v]]<-list()
  pT_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/F.1_overlap_test_COs/F.1.3_permutation_list_test_",v,".RDS", sep = ""))
 
  mean_permuted_overlaps<-c()
  observed_overlaps<-c()
  p_values<-c()
  
  for (p in Populations){
  observed_overlaps<-c(observed_overlaps, pT_list[[p]]$numOverlaps$observed)
  mean_permuted_overlaps<-c(mean_permuted_overlaps, mean(pT_list[[p]]$numOverlaps$permuted))
  p_values<-c(p_values, pT_list[[p]]$numOverlaps$pval)
  #graphic_list[[v]][[p]]<-pT_list[[p]]
  }#p
  
  part1<-cbind(c(v,"","",""), c("Total SVs", "Obs. overlaps", "Perm. overlaps","P value"))
  part2<-rbind(tabla_v[2,],observed_overlaps, mean_permuted_overlaps,p_values)
  part2[]<-round(part2[], digits = 4)
  
  tabla<-rbind(tabla,cbind(part1,part2))

}#v  

colnames(tabla)[1:2]<-""
COs<-cbind("Total CO intervals","", tabla_v[1,]); colnames(COs)<-colnames(tabla)
tabla<-rbind(COs, tabla)

write.csv(tabla, paste("/scratch/Federico/3_RILs/3_SVs/Results/F.1_overlap_test_COs/F.1.4_overlap_test_ALL.csv", sep = ""), row.names = FALSE)
saveRDS(graphic_list, "/scratch/Federico/3_RILs/3_SVs/Results/F.1_overlap_test_COs/F.1.5_permutation_list_test.RDS")
##############################
