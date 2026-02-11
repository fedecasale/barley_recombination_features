source("/scratch/Federico/3_RILs/2_CO_breakpoints/R_scripts/Z_my_functions.R")
library(regioneR)

win.size<-10000

#get windows from broad regions
distal_windows<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.0_accum_rec_prob_DISTAL.RDS")
peri_windows<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.0_accum_rec_prob_PERI.RDS")
pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere.RDS")

chrs<-paste("chr", 1:7, sep = "")  #clave sin la H para correr el permTest

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/F.2_overlap_test_COs_and_windows/")

################### SVs ################################

#I will divide COs intervals and SVs in overlap categories
sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.1_SVs_per_population.RDS")
sv_list[["translocations"]]<-NULL
SVs<-names(sv_list)


#divide svs
new_sv_list<-list()
for (h in categories){
new_sv_list[[h]]<-list()
if (h == "all"){j<-50}
if (h == "10K"){j<-10000}
if (h == "100K"){j<-100000}
if (h == "500K"){j<-500000}
for (v in SVs){
  new_sv_list[[h]][[v]]<-list()
  for (p in Populations){
    new_sv_list[[h]][[v]][[p]]<-list()
    for (c in 1:7){
      tabla<-sv_list[[v]][[p]][[c]]
      sv_lengths<-as.numeric(tabla[,which(colnames(tabla)=="stop")])-as.numeric(tabla[,which(colnames(tabla)=="start")])
      if (any(sv_lengths<=j)){tabla<-tabla[-which(sv_lengths<=j),]}; tabla<-make.matrix(tabla)
      new_sv_list[[h]][[v]][[p]][[c]]<-tabla
    }#p
  }#c
}#v
}#h
sv_list<-new_sv_list; rm(new_sv_list)
#########################################################
#create SVs all together
for (h in categories){
sv_list[[h]][["all"]]<-list()  
for (p in Populations){
  sv_list[[h]][["all"]][[p]]<-list()
  for (c in 1:7){
    SVs_all<-matrix(ncol = 5, nrow = 0); colnames(SVs_all)<-colnames(sv_list[[1]][["deletions"]][[p]][[c]])
    for (v in SVs){SVs_all<-rbind(SVs_all, sv_list[[h]][[v]][[p]][[c]])}#v
    sv_list[[h]][["all"]][[p]][[c]]<-SVs_all
  }#c
}#p
}#h
#########################################################

#selection
SVs<-names(sv_list[[1]])
SVs<-SVs[5]

rec_events<-breakpoint_list
##########################################################################################################################

################# overlap coldspots ################

pT_list<-list()

rec_events<-breakpoint_list
  
for (e in names(rec_events)){ cat(paste(e,": ", sep = ""))
#e<-"500K"
  rec_event<-rec_events[[e]]
  
  pT_list[[e]]<-list()
  
  v<-"all"
  
    for (p in Populations){ cat(p);cat(":")
      
      if (e == "intersected"){
        centromere<-average_table
        centromere[,1]<-chrs #needs format for permut function
        genome<-centromere[,c(1:3)]; colnames(genome)<-c("seqnames","start","end"); genome<-GRanges(as.data.frame(genome))
        pericentromere<-centromere[,c(1,4:5)]; colnames(pericentromere)<-c("seqnames","start","end"); pericentromere<-GRanges(as.data.frame(pericentromere))
      } else {
        centromere<-pericentromere_pops[[p]]
        centromere[,1]<-chrs #needs format for permut function
        genome<-centromere[,c(1:3)]; colnames(genome)<-c("seqnames","start","end"); genome<-GRanges(as.data.frame(genome))
        pericentromere<-centromere[,c(1,4:5)]; colnames(pericentromere)<-c("seqnames","start","end"); pericentromere<-GRanges(as.data.frame(pericentromere))
      }
      
      CO_regions_genome<-data.frame(); SV_regions_genome<-data.frame()
      
      for (c in 1:7){ #cat(c);cat("-")
        
        chr<-chrs[c]
        strand<-"+"  #dont know if OK
        starts<-as.numeric(rec_event[[p]][[c]][,1])
        ends<-as.numeric(rec_event[[p]][[c]][,2])
        widths<-(ends-starts)+1
        CO_regions <- data.frame(chr, starts, ends, widths, strand)
        colnames(CO_regions)<-c("seqnames","start","end","width","strand")
        #take out peri  #NO DEBERIA HABER NADA EN PERI, NO ENTIENDO
        to.delete<-starts[which(starts>centromere[c,4])]; to.delete<-to.delete[which(to.delete<centromere[c,5])]
        to.delete<-which(starts%in%to.delete)
        to.delete2<-ends[which(ends>centromere[c,4])]; to.delete2<-to.delete2[which(to.delete2<centromere[c,5])]
        to.delete2<-which(ends%in%to.delete2)
        to.delete<-unique(to.delete, to.delete2)
        if (length(to.delete)!=0){CO_regions<-CO_regions[-to.delete,]}
        CO_regions_genome<-rbind(CO_regions_genome, CO_regions)
        
        SV_regions<-make.matrix(sv_list[[e]][[v]][[p]][[c]][,1:4])
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
        SV_regions_genome<-rbind(SV_regions_genome, SV_regions)
        
      }#c
      
      CO_regions_genome<-GRanges(CO_regions_genome)
      SV_regions_genome<-GRanges(SV_regions_genome)
      
      #pruebo que los coldspots tienen mas overlap con las SVs que by chance
      pT<-permTest(A = SV_regions_genome, B = CO_regions_genome, ntimes=100, alternative = "greater",
                   randomize.function = randomizeRegions, genome = genome, mask = pericentromere, per.chromosome = TRUE,
                   evaluate.function = numOverlaps, non.overlapping = FALSE)
      plot(pT)
      pT_list[[e]][[p]]<-pT
      
      cat(fill = TRUE)
      
    }#p
    
    cat(fill = TRUE)
  
  saveRDS(pT_list, paste("/scratch/Federico/3_RILs/3_SVs/Results/F.2_overlap_test_COs_and_windows/F.1.1_COs_SVs_overlap_by_size.RDS", sep = ""))
  
}#e

# ######### re-unify ##########

pT_list<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/F.2_overlap_test_COs_and_windows/F.1.1_COs_SVs_overlap_by_size.RDS", sep = ""))

v<-"all"
tabla_e<-matrix(ncol = 5, nrow = 0); colnames(tabla_e)<-c("window_type","variable", Populations)
tabla_0<-tabla_e
for (e in names(pT_list)){
  
  tot_svs<-c(); mean_permuted_overlaps<-c(); observed_overlaps<-c(); p_values<-c()
  for (p in Populations){
    tot_svs_p<-nrow(sv_list[[e]][[v]][[p]][[1]]); for (c in 2:7){tot_svs_p<-sum(tot_svs_p, nrow(sv_list[[e]][[v]][[p]][[c]]))}
    tot_svs<-c(tot_svs, tot_svs_p)
    observed_overlaps<-c(observed_overlaps, pT_list[[e]][[p]]$numOverlaps$observed)
    mean_permuted_overlaps<-c(mean_permuted_overlaps, mean(pT_list[[e]][[p]]$numOverlaps$permuted))
    p_values<-c(p_values, pT_list[[e]][[p]]$numOverlaps$pval)
  }#p
  part1<-cbind(c(e,"","",""),c("Total SVs", "Obs. overlaps", "Perm. overlaps","P value"))
  part2<-rbind(tot_svs, observed_overlaps, round(mean_permuted_overlaps),round(p_values, digits = 3))
  tabla<-rbind(tabla_e,cbind(part1,part2))
  tabla_0<-rbind(tabla_0, tabla)
}#e
write.xlsx(tabla_0, paste("/scratch/Federico/3_RILs/3_SVs/Results/F.2_overlap_test_COs_and_windows/F.1.2_COs_SVs_overlap_by_size.xlsx", sep = ""), row.names = FALSE, sheetName = e, append = TRUE)



# #############################################################################################################
