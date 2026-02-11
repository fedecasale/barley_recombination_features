library(ggplot2)
library(impressionist.colors)

breakpoint_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")

layers<-c("1_first", "2_second", "3_third")
chrs<-paste(1:7, "H", sep = "")
Populations<-names(breakpoint_list[[1]])

COs_list<-list()
for(l in layers){ 
breakpoints<-matrix(ncol = 3, nrow = 0); colnames(breakpoints)<-c("Population", "chr", "breakpoint")
for (p in Populations){ for(c in 1:7){breakpoints<-rbind(breakpoints, breakpoint_list[[l]][[p]][[c]][,c(1,3,6)])}}
breakpoints[,2]<-gsub("chr","",breakpoints[,2])
breakpoints<-as.data.frame(breakpoints)
breakpoints$breakpoint<-as.numeric(as.character(breakpoints$breakpoint))
COs_list[[l]]<-breakpoints
}

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.1_SVs_per_population.RDS")
SVs<-names(sv_list)
SVs_list<-list()
for (v in SVs){ 
    SV<-matrix(ncol = 3, nrow = 0); colnames(SV)<-c("Population", "chr", "SV_pos")
    #Estoy usando solo el end, podria hacer lineas que indiquen  el segment
    for (p in Populations){for (c in 1:7){SV<-rbind(SV, cbind(p,sv_list[[v]][[p]][[c]][,c(1,3)]))}}  
    SV[,2]<-gsub("chr","",SV[,2])
    SV<-as.data.frame(SV)
    SV$SV_pos<-as.numeric(as.character(SV$SV_pos))
    SVs_list[[v]]<-SV
}

chr_lengths<-read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv")
chr_lengths[,1]<-chrs
chr_lengths$chr_start<-as.numeric(as.character(chr_lengths$chr_start))
chr_lengths$chr_end<-as.numeric(as.character(chr_lengths$chr_end))
chr_lengths$peri_start<-as.numeric(as.character(chr_lengths$peri_start))
chr_lengths$peri_end<-as.numeric(as.character(chr_lengths$peri_end))

l<-layers[3]

y_seq<-seq(0.1, 1, by = 0.01)


see.palette(8,1)
get.color(8,1,11)
get.color(8,1,9)
get.color(8,2,15)

pdf("/scratch/Federico/3_RILs/3_SVs/Results/E_stats/E.2_distribution_across_genome.pdf", height = 10, width = 15)

ggplot() +
  #COs
  geom_point(data=COs_list[[l]], aes(x=breakpoint, y = sample(y_seq, nrow(COs_list[[l]]), replace = TRUE)+4), shape=20, size=0.0001, color = "#1e3c5a") + 
  #inversions
  geom_point(data=SVs_list[[4]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[4]]), replace = TRUE)+3), shape=20, size=0.0001, color = "#694b5a") + 
  #insertions
  geom_point(data=SVs_list[[3]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[3]]), replace = TRUE)+2), shape=20, size=0.0001, color = "#4b875a") + 
  #deletions
  geom_point(data=SVs_list[[1]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[1]]), replace = TRUE)+1), shape=20, size=0.0001, color = "#e1d25a") + 
  #duplications
  geom_point(data=SVs_list[[2]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[2]]), replace = TRUE)), shape=20, size=0.0001, color = "#c3693c") + 
  #translocations
  #geom_point(data=SVs_list[[5]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[5]]), replace = TRUE)), shape=20, size=0.0001, color = "#694b5a") + 
  #pericentromere
  #geom_segment(data = chr_lengths, aes(x = peri_start, y = 0, xend = peri_end, yend = 0), color = "red")+
  geom_segment(data = chr_lengths, aes(x = centromere, y = 0, xend = centromere, yend = max(y_seq)+4), color = "black", size = 0.75)+
  geom_segment(data = chr_lengths, aes(x = peri_start, y = 0, xend = peri_start, yend = max(y_seq)+4), color = "black", linetype = 2, size = 0.6)+
  geom_segment(data = chr_lengths, aes(x = peri_end, y = 0, xend = peri_end, yend = max(y_seq)+4), color = "black", linetype = 2, size = 0.6)+

  scale_x_continuous(name = "Physical distance (Mbp)", position = "bottom") +
  scale_y_continuous(expand = c(0,0))+
  
  facet_grid(cols = vars(Population), rows = vars(chr), scales = "free", space = "free", switch = "y") +
  
  theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12.5),
        strip.text = element_text(size = 12.5),
        panel.background = element_blank(),
        #panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(),
        )
dev.off()

file.copy("/scratch/Federico/3_RILs/3_SVs/Results/E_stats/E.2_distribution_across_genome.pdf", 
          "/scratch/Federico/3_RILs/6_graphics_for_paper/Results/E_2_distribution_across_genome.pdf")
