library(ggplot2)

Populations<-c("HvDRR13","HvDRR27", "HvDRR28")
win.size<-10000

#get data
#rec_prob<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.0_accumulated_rec_prob.RDS")
rec_prob<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.0_accumulated_rec_prob_NORMALIZED.RDS")
coldspots<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.2_coldspots/A.2.1_coldspots_per_pop_WINDOWS.RDS")
hotspots<-readRDS("/scratch/Federico/3_RILs/7_window_types/Results/A.1_hotspots/A.1.2_hotspots_per_pop_WINDOWS.RDS")

#graphic table
e_list<-list(rec_prob, coldspots, hotspots); names(e_list)<-c("rec_prob", "coldspots", "hotspots")
graphic_table<-matrix(ncol = 5, nrow = 0)
colnames(graphic_table)<-c("Population","chr","Mbp","cM_Mbp","region")
for (e in names(e_list)){ cat(e); cat(": ")
   for (p in Populations){ cat(p); cat("-")
     for (c in 1:7){
       chr_to_add<-cbind(p, paste(c,"H", sep = ""), e_list[[e]][[p]][[c]][,c(1,3)], e)
       colnames(chr_to_add)<-colnames(graphic_table); row.names(chr_to_add)<-NULL 
       graphic_table<-rbind(graphic_table, chr_to_add)
     }#c
   }#p
 cat(fill = TRUE)  
}#e

graphic_table<-as.data.frame(graphic_table)
graphic_table$Mbp=as.numeric(levels(graphic_table$Mbp))[graphic_table$Mbp]
graphic_table$Mbp<-graphic_table$Mbp/1000000
graphic_table$cM_Mbp=as.numeric(levels(graphic_table$cM_Mbp))[graphic_table$cM_Mbp]

#get pericentromere data
#centromere <-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede.csv")
pericentromere_pops<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.4.2_pericentromere.RDS")
centromere<-matrix(ncol = 5, nrow = 0)
colnames(centromere)<-c("Population", "chr", "peri_start", "peri_end", "centromere")
for (p in Populations){centromere<-rbind(centromere, cbind(p, pericentromere_pops[[p]][,c(1,4:6)]))}
centromere<-as.data.frame(centromere)
for (c in 3:5){centromere[,c]<-as.numeric(as.character(centromere[,c]))}
centromere[,3:5]<-centromere[,3:5]/1000000
color_table<-as.matrix(get.color(8,3,c(5,7,12)))
centro_col<-get.color(8,3,14); peri_col<-get.color(8,3,13)

centro_peri<-list()
centro_peri[[1]]<-geom_vline(data = centromere, aes(xintercept = centromere), color = centro_col, size = 0.75)
centro_peri[[2]]<-geom_vline(data = centromere, aes(xintercept = peri_start), linetype = 2, color = peri_col, size = 0.5)
centro_peri[[3]]<-geom_vline(data = centromere, aes(xintercept = peri_end), linetype = 2, color = peri_col, size = 0.5)

facet<-facet_grid(cols = vars(Population), rows = vars(chr), scales = "free", space = "free") 

m.f<-0.1
ancho<-1300*m.f/45*3
alto<-100*m.f

thema<-theme(strip.background = element_blank(),
              panel.background  = element_rect(fill=NA, color = "black", size = 1, linetype = "solid"),
             panel.grid = element_blank(),
             panel.spacing = unit(3*m.f, "cm"),
              #panel.background  = element_rect(fill=NA, size = 3*m.f, linetype = "solid"),
              strip.text = element_text(size = 100*m.f, margin = margin(50*m.f, 50*m.f, 50*m.f, 50*m.f)))

pdf(paste("/scratch/Federico/3_RILs/7_window_types/Results/A.3_hotspots_coldspots.pdf", sep = ""), width=ancho, height=alto)
ggplot() +
     geom_line(data=graphic_table[which(graphic_table$region=="rec_prob"),], aes(x=Mbp, y=cM_Mbp), color = "#aaaaaa", size = 0.25) +
     geom_point(data=graphic_table[which(graphic_table$region=="hotspots"),], aes(x=Mbp, y=cM_Mbp), color = "#e71414", size = 0.6, shape = 19) +
     geom_point(data=graphic_table[which(graphic_table$region=="coldspots"),], aes(x=Mbp, y=cM_Mbp), color = "#1674b1", size = 0.1, shape = 2) +
     
     scale_x_continuous(name = "Physical distance (Mbp)", position = "bottom", expand = c(0.01,0.01)) +
     scale_y_continuous(name = "Normalized CO probability", position = "left") +
 
     centro_peri + facet + thema
dev.off()

file.copy("/scratch/Federico/3_RILs/7_window_types/Results/A.3_hotspots_coldspots.pdf", 
          "/scratch/Federico/3_RILs/6_graphics_for_paper/Results/A.3_hotspots_coldspots.pdf")


#grafico para las pop por separado

alto<-100*m.f/2

pdf(paste("/scratch/Federico/3_RILs/7_window_types/Results/A.3_hotspots_coldspots_per_pop.pdf", sep = ""), width=ancho, height=alto)
for (p in Populations){

#prepare data
graphic_table_pop<-graphic_table[which(graphic_table$Population==p),]
centromere_pop<-centromere[which(centromere[,1]==p),]
centro_peri_pop<-list()
centro_peri_pop[[1]]<-geom_vline(data = centromere_pop, aes(xintercept = centromere), color = centro_col, size = 0.75)
centro_peri_pop[[2]]<-geom_vline(data = centromere_pop, aes(xintercept = peri_start), linetype = 2, color = peri_col, size = 0.5)
centro_peri_pop[[3]]<-geom_vline(data = centromere_pop, aes(xintercept = peri_end), linetype = 2, color = peri_col, size = 0.5)

#make graphic
print(
ggplot() +
  geom_line(data=graphic_table_pop[which(graphic_table_pop$region=="rec_prob"),], aes(x=Mbp, y=cM_Mbp), color = "#aaaaaa", size = 0.25) +
  geom_point(data=graphic_table_pop[which(graphic_table_pop$region=="hotspots"),], aes(x=Mbp, y=cM_Mbp), color = "#e71414", size = 0.6, shape = 19) +
  geom_point(data=graphic_table_pop[which(graphic_table_pop$region=="coldspots"),], aes(x=Mbp, y=cM_Mbp), color = "#1674b1", size = 0.1, shape = 2) +
  
  scale_x_continuous(name = "Physical distance (Mbp)", position = "bottom", expand = c(0.01,0.01)) +
  scale_y_continuous(name = "Normalized CO probability", position = "left") +
  
  centro_peri_pop + facet + thema
)#print
}#p
dev.off()
