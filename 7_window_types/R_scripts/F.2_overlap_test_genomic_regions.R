#Rscript /scratch/Federico/3_RILs/7_window_types/R_scripts/F.2_overlap_test_genomic_regions.R

library(regioneR)
source("/scratch/Federico/3_RILs/2_CO_breakpoints/R_scripts/Z_my_functions.R")

#get windows from broad regions
distal_windows<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.0_accum_rec_prob_DISTAL.RDS")
peri_windows<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.0_accum_rec_prob_PERI.RDS")
pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere.RDS")

#windows to compare
genomic_regions<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/F_overlap_genomic_regions/F.1_genomic_regions.RDS")

rec_events<-list()
rec_events[["coldspots"]]<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.1_coldspots_per_pop_WINDOWS.RDS")
rec_events[["hotspots"]]<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.2_hotspots_per_pop_WINDOWS.RDS")
Populations<-names(rec_events[[1]])
chrs<-paste("chr", 1:7, sep = "")  #clave sin la H para correr el permTest

breakpoint_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")
breakpoint_list<-breakpoint_list[[3]] #just the very sure COs by now because are few
for (p in Populations){for (c in 1:7){
breakpoints<-breakpoint_list[[p]][[c]][,"breakpoint"]
breakpoint_list[[p]][[c]]<-cbind(breakpoints ,breakpoints)
colnames(breakpoint_list[[p]][[c]])<-c("start","end")
}}
rec_events[["CO_breakpoints"]]<-breakpoint_list

Populations<-names(rec_events[[1]])
chrs<-paste("chr", 1:7, sep = "")  #clave sin la H para correr el permTest

#pT_list<-list()

for (e in names(rec_events)[3]){ cat(paste(e,": ", sep = ""), fill = TRUE)

  #pT_list[[e]]<-list()
  pT_list<-readRDS(paste("/scratch/Federico/3_RILs/7_window_types/Results/F_overlap_genomic_regions/F.2_overlap_gene_regions_",e,".RDS", sep = ""))
  
  for (g in names(genomic_regions)){  cat(paste(e,": ", sep = ""))
  
    # e<-"hotspots"
    # g<-"exons"
    # p<-Populations[1]
    
    rec_event<-rec_events[[e]]
    
    pT_list[[e]][[g]]<-list()
    
    for (p in Populations){ cat(p);cat("-")
      
      centromere<-pericentromere_pops[[p]]
      centromere[,1]<-chrs #needs format for permut function
      genome<-centromere[,c(1:3)]; colnames(genome)<-c("seqnames","start","end"); genome<-GRanges(as.data.frame(genome))
      pericentromere<-centromere[,c(1,4:5)]; colnames(pericentromere)<-c("seqnames","start","end"); pericentromere<-GRanges(as.data.frame(pericentromere))
      
      window_type_genome<-data.frame(); gene_regions_genome<-data.frame()
      
      for (c in 1:7){ #cat(c);cat("-")
      
        chr<-chrs[c]
        strand<-"+"  #dont know if OK
        starts<-as.numeric(rec_event[[p]][[c]][,1])
        ends<-as.numeric(rec_event[[p]][[c]][,2])
        widths<-(ends-starts)+1
        window_type <- data.frame(chr, starts, ends, widths, strand)
        colnames(window_type)<-c("seqnames","start","end","width","strand")
        #take out peri  #NO DEBERIA HABER NADA EN PERI, NO ENTIENDO
        to.delete<-starts[which(starts>centromere[c,4])]; to.delete<-to.delete[which(to.delete<centromere[c,5])]
        to.delete<-which(starts%in%to.delete)
        to.delete2<-ends[which(ends>centromere[c,4])]; to.delete2<-to.delete2[which(to.delete2<centromere[c,5])]
        to.delete2<-which(ends%in%to.delete2)
        to.delete<-unique(to.delete, to.delete2)
        if (length(to.delete)!=0){window_type<-window_type[-to.delete,]}
        window_type_genome<-rbind(window_type_genome, window_type)
        
        gene_regions<-as.matrix(genomic_regions[[g]][[c]])
        starts<-as.numeric(gene_regions[,1])   #clave poner el c
        ends<-as.numeric(gene_regions[,2])
        widths<-(ends-starts)+1
        gene_regions <- data.frame(chr, starts, ends, widths, strand)
        colnames(gene_regions)<-c("seqnames","start","end","width","strand")
        #take out peri
        to.delete<-starts[which(starts>centromere[c,4])]; to.delete<-to.delete[which(to.delete<centromere[c,5])]
        to.delete<-which(starts%in%to.delete)
        to.delete2<-ends[which(ends>centromere[c,4])]; to.delete2<-to.delete2[which(to.delete2<centromere[c,5])]
        to.delete2<-which(ends%in%to.delete2)
        to.delete<-unique(to.delete, to.delete2)
        if (length(to.delete)!=0){gene_regions<-gene_regions[-to.delete,]}
        gene_regions_genome<-rbind(gene_regions_genome, gene_regions)
        
      }#c
      
      window_type_genome<-GRanges(window_type_genome)
      # esto no deberia ser necesitado (es en intergenic que pasa)
      if (any(gene_regions_genome$end<gene_regions_genome$start)){
      gene_regions_genome<-gene_regions_genome[-which(gene_regions_genome$end<gene_regions_genome$start),]
      }
      gene_regions_genome<-GRanges(gene_regions_genome)
      
      #pruebo que los coldspots tienen mas overlap con las genes que by chance
      pT<-permTest(A = gene_regions_genome, B = window_type_genome, ntimes=100, alternative = "greater",
                   randomize.function = randomizeRegions, genome = genome, mask = pericentromere, per.chromosome = TRUE,
                   evaluate.function = numOverlaps, non.overlapping = FALSE)
      
      #plot(pT)
      pT_list[[e]][[g]][[p]]<-pT
      
    }#p
    
    cat(fill = TRUE)
    
    saveRDS(pT_list, paste("/scratch/Federico/3_RILs/7_window_types/Results/F_overlap_genomic_regions/F.2_overlap_gene_regions_",e,".RDS", sep = ""))
    
  }#g
}#e

# ######### re-unify ##########
# pT_list<-list()
# for (e in names(rec_events)){ cat(paste(e,": ", sep = ""), fill = TRUE)
#   rec_event_overlap<-readRDS(paste("/scratch/Federico/3_RILs/7_window_types/Results/F_overlap_genomic_regions/F.2_overlap_gene_regions_",e,".RDS", sep = ""))[[1]]
#   pT_list[[e]]<-rec_event_overlap
# }
# saveRDS(pT_list, paste("/scratch/Federico/3_RILs/7_window_types/Results/F_overlap_genomic_regions/F.2_overlap_gene_regions.RDS", sep = ""))
#####
# 
# pT_list<-readRDS(paste("/scratch/Federico/3_RILs/7_window_types/Results/F_overlap_genomic_regions/F.2_overlap_gene_regions.RDS", sep = ""))
# 
# tabla_final<-matrix(ncol = 6, nrow = 0)
# colnames(tabla_final)<- c("Altered recombination region", "Genomic region","Population", "Observed overlaps", "Permuted overlaps","P value")
# 
# for (e in names(rec_events)){ cat(paste(e,": ", sep = ""), fill = TRUE)
# 
#   tabla_e<-matrix(ncol = 5, nrow = 0)
#   colnames(tabla_e)<- c("Genomic region","Population", "Observed overlaps", "Permuted overlaps","P value")
#   
# for (g in names(genomic_regions)){ cat(paste(e,": ", sep = ""))
#   
#     tabla_g<-matrix(ncol = 3, nrow = 3)
#     colnames(tabla_g)<- c("Observed overlaps", "Permuted overlaps","P value")
#     row.names(tabla_g)<-Populations
#     
#     for (p in Populations){ cat(p);cat("-")
#       tabla_g[p,1]<-round(pT_list[[e]][[g]][[p]]$numOverlaps$observed) 
#       tabla_g[p,2]<-round(mean(pT_list[[e]][[g]][[p]]$numOverlaps$permuted))
#       tabla_g[p,3]<-round(pT_list[[e]][[g]][[p]]$numOverlaps$pval, digits = 4)
#     }
#     
#     tabla_g<-cbind(c(g,"",""), Populations, tabla_g)  
#     tabla_e<-rbind(tabla_e, tabla_g)
# }#g    
# 
# tabla_e<-cbind(c(e, "", "", "", "", "", "","", "", "", "",""), tabla_e)  
# tabla_final<-rbind(tabla_final, tabla_e)
# }#e
# 
# write.csv(tabla_final, "/scratch/Federico/3_RILs/7_window_types/Results/F_overlap_genomic_regions/F.2_overlap_gene_regions.csv", row.names = FALSE)
# 


# #############################################################################################################
