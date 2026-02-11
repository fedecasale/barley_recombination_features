window_sizes<-c(500000, 1000000)[2]

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))

ave_rates<-as.matrix(read.csv("/scratch/Federico/3_RILs/5_correlation/Results/A.5.1_rec_rates_V3_genome_and_chrs.csv"))
row.names(ave_rates)<-ave_rates[,1]; ave_rates<-ave_rates[,-1]

#for (s in window_sizes){ cat(s, fill = TRUE)
s<-window_sizes

rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/A.5.2_rec_rates_V3_windows_",s,".RDS", sep = ""))
Populations<-names(rec_list)
#Populations<-Populations[c(10,24,25)]

#ABAJO CALCULABA UNA PERI PARA CADA CHR DE CADA POP, PERO EN VERDAD DEBO HACER UNA GENERAL PARA TODAS LAS POPS PER CHR

ave_rec_list<-list()
for (c in 1:7){ cat(c,"-", sep = "")
  #saco un promedio per window
  ave_rec_chr<-rec_list[[1]][[c]]; ave_rec_chr[]<-0
  for (p in Populations){ ave_rec_chr<-ave_rec_chr+rec_list[[p]][[c]]}
  ave_rec_chr<-ave_rec_chr/length(Populations)
  ave_rec_list[[c]]<-ave_rec_chr  
  #create table
  tabla<-matrix(ncol = 4, nrow = length(ave_rec_chr))
  colnames(tabla)<-c("chr","Mbp","cM/Mbp","region")
  tabla[,1]<-paste(c,"H",sep = "")
  tabla[,2]<-names(ave_rec_chr)
  tabla[,3]<-ave_rec_chr
  tabla[,4]<-"D"
  #calculate peri region  as Baker 2014:
  #"We define the LR-PC region as the continuous region surrounding the centromere for which
  # recombination rate is 20-fold lower than the average for the barley genome (Choo, 1998)".
  ave_rate<-mean(ave_rec_chr)
  low_rec_windows<-which(ave_rec_chr<(ave_rate/20))
  tabla[low_rec_windows,4]<-"P"  
  #make P a continuous trend (en teoria como es promedio, las islands dberian quedar afuera)
  tabla[low_rec_windows[1]:low_rec_windows[length(low_rec_windows)],4]<-"P" 
  #save
  ave_rec_list[[c]]<-tabla
}  
saveRDS(ave_rec_list, paste("/scratch/Federico/3_RILs/5_correlation/Results/B.0_AVERAGE_rec_rates_windows_regions_",s,".RDS", sep = ""))


# #modify summary table
# ave_rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.0_AVERAGE_rec_rates_windows_regions_",s,".RDS", sep = ""))
# centromere<-read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv")
# for (c in 1:7){
# pericentromere<-which(ave_rec_list[[c]][,4]=="P")
# centromere[c,4]<-ave_rec_list[[c]][pericentromere[1],2]
# centromere[c,5]<-ave_rec_list[[c]][pericentromere[length(pericentromere)],2]
# }
# write.csv(centromere, "/scratch/Federico/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede.csv", row.names = FALSE)

#####################################################

#add to pops and look for islands
for (p in Populations){ cat(c(p,":"), sep = "")
for (c in 1:7){ cat(c,"-", sep = "")

windows<-rec_list[[p]][[c]]
tabla<-matrix(ncol = 4, nrow = length(windows))
colnames(tabla)<-c("chr","Mbp","cM/Mbp","region")#,"fold_change")
tabla[,1]<-paste(c,"H",sep = "")
tabla[,2]<-names(windows)
tabla[,3]<-windows
tabla[,4]<-"D"

#add pericentromere region
tabla[which(ave_rec_list[[c]][,4]=="P"),4]<-"P"

#detect islands
ave_rate<-mean(as.numeric(ave_rec_list[[c]][,3]))

low_rec<-which(as.numeric(tabla[,3])<(ave_rate/20))
if (any(tabla[low_rec,4]=="D")){tabla[low_rec[which(tabla[low_rec,4]=="D")],4]<-"I"}

# #asi hice el primer paper
# low_rec_windows<-which(windows<(ave_rate/20))
# tabla[low_rec_windows,4]<-"P"
# 
# #identifico grupos 
# groups<-list()
# cont<-low_rec_windows[2:length(low_rec_windows)]-low_rec_windows[1:(length(low_rec_windows)-1)]
# if (any(cont>1)){
# separators<-which(cont>1)  
# groups[[1]]<-c(names(cont[1:(separators[1]-1)]))
# for (i in 1:length(separators)){
# if (i == length(separators)){groups[[i+1]]<-c(names(cont[separators[i]:length(cont)]))
# } else {
# groups[[i+1]]<-c(names(cont[separators[i]:separators[(i+1)-1]]))  
# }#if
# }#i
# } else {
# groups[[1]]<-c(names(cont[1:length(cont)]))
# }#if
# 
# #check whether groups are pericentromere or islands
# centromere<-as.numeric(chr_lengths[c,6])
# quartile<-as.numeric(chr_lengths[c,3])/4
# interquantile_range_min<-centromere-quartile
# interquantile_range_max<-centromere+quartile
# 
# for (i in 1:length(groups)){ 
# min<-min(as.numeric(groups[[i]]))
# max<-max(as.numeric(groups[[i]]))
# #si el centromere esta adentro del grupo, el grupo es el pericentromere
# if ((centromere>=min)&(centromere<=max)){
# names(groups)[i]<-"pericentromere"  
# tabla[which(tabla[,2]%in%groups[[i]]),4]<-"P"
# } else {
# #si el centromere no esta adentro del grupo, pero estan cerca, el grupo tb es pericentromere
# #para considerar cerca, genere una especie de interquantile range pero alrededor del centromere
# if ((min>=interquantile_range_min)&(max<=interquantile_range_max)) {
# names(groups)[i]<-"pericentromere_2"
# tabla[which(tabla[,2]%in%groups[[i]]),4]<-"P2" #P2
# } else {
# #si el centromere y el grupo estan lejos, el grupo es island
# names(groups)[i]<-"island"  
# tabla[which(tabla[,2]%in%groups[[i]]),4]<-"I"
# }#if
# }#if
# }#i
# 
# ####el pericentromere va a abarcar desde la primera P o P2 hasta la ultima
# first<-which(tabla[,4]%in%c("P","P2"))[1]
# last<-which(tabla[,4]%in%c("P","P2"))[length(which(tabla[,4]%in%c("P","P2")))]
# tabla[first:last,4]<-"P"

# plot
# COLORES<-c(rep("blue",length(which(tabla[,4]=="D"))), rep("red", length(which(tabla[,4]%in%c("P","P2")))),rep("orange",length(which(tabla[,4]=="I"))))
# names(COLORES)<-c(which(tabla[,4]=="D"),which(tabla[,4]%in%c("P","P2")),which(tabla[,4]=="I"))
# COLORES<-COLORES[order(as.numeric(names(COLORES)))]
# plot(tabla[,3]~tabla[,2], col=COLORES)
# title(main = paste(p,c,sep = " "))
 
# #calculate fold change
# tabla[,5]<-windows/ave_rate

#save
rec_list[[p]][[c]]<-tabla

}#c  
}#p  

saveRDS(rec_list, paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))

##############################################

#prepare table
graphic_table2<-matrix(ncol = 5, nrow = 0)
colnames(graphic_table2)<-c("Population","chr","Mbp","cM_Mbp","region")
for (p in Populations){for (c in 1:7){graphic_table2<-rbind(graphic_table2, cbind(p,rec_list[[p]][[c]]))}}
graphic_table2<-as.data.frame(graphic_table2)
graphic_table2$Mbp=as.numeric(levels(graphic_table2$Mbp))[graphic_table2$Mbp]
graphic_table2$Mbp<-graphic_table2$Mbp/1000000
graphic_table2$cM_Mbp=as.numeric(levels(graphic_table2$cM_Mbp))[graphic_table2$cM_Mbp]
#graphic_table2$fold_change=as.numeric(levels(graphic_table2$fold_change))[graphic_table2$fold_change]

COLORES<-as.character(graphic_table2$region)
COLORES[which(COLORES=="D")]<-"red"
COLORES[which(COLORES%in%c("P","P2"))]<-"blue"
COLORES[which(COLORES=="I")]<-"green"

m.f<-0.1
ancho<-1300*m.f
alto<-100*m.f 

pdf(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.2_rec_rates_windows_regions_",s,".pdf", sep = ""), width=ancho, height=alto)
print(
ggplot() + 
  geom_point(data=graphic_table2 ,stat = "identity", aes(x=Mbp, y=cM_Mbp), color = COLORES) +
  scale_x_continuous(name = "physical distance (Mbp)", position = "bottom") +
  scale_y_continuous(name = "recombination rate (cM/Mbp)", position = "right") +    
  theme(legend.position = "none",
        axis.text.x.bottom = element_text(size = 110*m.f),
        axis.text.y.right  = element_text(size = 110*m.f),
        axis.ticks.length = unit(1*m.f, "cm")
  ) + 
  facet_grid(cols = vars(Population), rows = vars(chr), scales = "free", space = "free_y", switch = "y") +
  
  theme(strip.background = element_blank(), 
        panel.background  = element_rect(fill=NA, color = "black", size = 1, linetype = "solid"),
        panel.grid = element_blank(),
        panel.spacing = unit(3*m.f, "cm"),
        #panel.background  = element_rect(fill=NA, size = 3*m.f, linetype = "solid"),
        strip.text = element_text(size = 150*m.f, margin = margin(50*m.f, 50*m.f, 50*m.f, 50*m.f)))
)
dev.off()

##############################################

#}#s

