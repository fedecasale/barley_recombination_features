accumulated_rec_prob<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.3_CO_prob_per_window_GENETIC_MAPS.RDS", sep = ""))
window_sizes<-names(accumulated_rec_prob)
layers<-c("1_first", "2_second", "3_third")
chrs<-paste("chr", 1:7, "H", sep = "")
Populations<-names(accumulated_rec_prob[[1]])

COs_list<-list()
for (w in window_sizes){
  COs_list[[w]]<-list()
  for (l in layers){ cat(l, fill = TRUE)
    breakpoints<-matrix(ncol = 5, nrow = 0)
    colnames(breakpoints)<-c("Population", "chr", "start","end","rec_prob")
    for (p in Populations){ cat(p); cat(": ")  
      for (c in 1:7){ cat(c); cat("-")   
        PROB.windows<-cbind(p,c,accumulated_rec_prob[[w]][[p]][[c]]) 
        colnames(PROB.windows)<-c("Population", "chr", "start","end","rec_prob")
        breakpoints<-rbind(breakpoints, PROB.windows)
    breakpoints<-as.data.frame(breakpoints)
    breakpoints$start<-as.numeric(as.character(breakpoints$start))
    breakpoints$rec_prob<-as.numeric(as.character(breakpoints$rec_prob))
  }
    }
COs_list[[w]][[l]]<-breakpoints
  }
}

pdf("/scratch/Federico/3_RILs/3_SVs/Results/B.2_random_COs_with_prob/B.2.2.1_distribution_rec_prob_GENETIC_MAPS.pdf")

ggplot() +
  geom_point(data=COs_list[[1]][[1]], aes(x=start, y = rec_prob), shape=20, size=0.1, color = "black") + 
  # geom_point(data=SVs_list[[1]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[1]]), replace = TRUE)+2), shape=20, size=0.0001, color = "#1e3c5a") + 
  # geom_point(data=SVs_list[[2]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[2]]), replace = TRUE)+1), shape=20, size=0.0001, color = "#4b875a") + 
  # geom_point(data=SVs_list[[3]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[3]]), replace = TRUE)), shape=20, size=0.0001, color = "#e1d25a") + 
  # geom_segment(data = chr_lengths, aes(x = peri_start, y = 0, xend = peri_end, yend = 0), color = "red")+
  scale_x_continuous(name = "physical distance (Mbp)", position = "bottom") +
  facet_grid(cols = vars(Population), rows = vars(chr), scales = "free", space = "free", switch = "y") +
  theme(strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
  )

dev.off()

##############################

accumulated_rec_prob<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1_CO_prob_per_window.RDS", sep = ""))

#dir.create("/scratch/Federico/3_RILs/3_SVs/Results/B.2_random_COs_with_prob")

COs_list<-list()
for (w in window_sizes){
  COs_list[[w]]<-list()
  for (l in layers){ cat(l, fill = TRUE)
    breakpoints<-matrix(ncol = 5, nrow = 0)
    colnames(breakpoints)<-c("Population", "chr", "start","end","rec_prob")
    for (p in Populations){ cat(p); cat(": ")  
      for (c in 1:7){ cat(c); cat("-")   
        PROB.windows<-cbind(p,c,accumulated_rec_prob[[w]][[l]][[p]][[c]]) 
        colnames(PROB.windows)<-c("Population", "chr", "start","end","rec_prob")
        breakpoints<-rbind(breakpoints, PROB.windows)
        breakpoints<-as.data.frame(breakpoints)
        breakpoints$start<-as.numeric(as.character(breakpoints$start))
        breakpoints$rec_prob<-as.numeric(as.character(breakpoints$rec_prob))
      }
    }
    COs_list[[w]][[l]]<-breakpoints
  }
}

pdf("/scratch/Federico/3_RILs/3_SVs/Results/B.2_random_COs_with_prob/B.2.2.2_distribution_rec_prob_OLD_L3.pdf")

ggplot() +
  geom_point(data=COs_list[[1]][[3]], aes(x=start, y = rec_prob), shape=20, size=0.1, color = "black") + 
  # geom_point(data=SVs_list[[1]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[1]]), replace = TRUE)+2), shape=20, size=0.0001, color = "#1e3c5a") + 
  # geom_point(data=SVs_list[[2]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[2]]), replace = TRUE)+1), shape=20, size=0.0001, color = "#4b875a") + 
  # geom_point(data=SVs_list[[3]], aes(x=SV_pos, y = sample(y_seq, nrow(SVs_list[[3]]), replace = TRUE)), shape=20, size=0.0001, color = "#e1d25a") + 
  # geom_segment(data = chr_lengths, aes(x = peri_start, y = 0, xend = peri_end, yend = 0), color = "red")+
  scale_x_continuous(name = "physical distance (Mbp)", position = "bottom") +
  facet_grid(cols = vars(Population), rows = vars(chr), scales = "free", space = "free", switch = "y") +
  theme(strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
  )

dev.off()
