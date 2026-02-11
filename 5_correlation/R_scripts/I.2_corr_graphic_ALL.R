setwd("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/")

library(ggplot2)
library(impressionist.colors)

s<-1000000

#centromere <-read.csv("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede.csv")
centromere <-read.csv("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede_5.csv")
color_table<-as.matrix(get.color(8,3,c(5,7,12)))
centro_col<-get.color(8,3,14); peri_col<-get.color(8,3,13)
#centro_col<-"blue"; peri_col<-"blue"
centromere[,4:6]<-centromere[,4:6]/s

centro_peri<-list()
centro_peri[[1]]<-geom_vline(data = centromere, aes(xintercept = centromere), color = centro_col, linewidth = 0.75)
centro_peri[[2]]<-geom_vline(data = centromere, aes(xintercept = peri_start), linetype = 2, color = peri_col, linewidth = 0.5)
centro_peri[[3]]<-geom_vline(data = centromere, aes(xintercept = peri_end), linetype = 2, color = peri_col, linewidth = 0.5)

#Prepare table for graphic
sign_var_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.1.2_sw_reg_per_win_significant_variables.RDS")
graphic_table<-sign_var_list[[1]][["graphic_table"]]  
for (c in 2:7){graphic_table<-rbind(graphic_table,sign_var_list[[c]][["graphic_table"]])}

#change names
graphic_table[,3]<-gsub("par.dif.cpg", "Parental methylation difference: CpG", graphic_table[,3])
graphic_table[,3]<-gsub("par.dif.chg", "Parental methylation difference: CHG", graphic_table[,3])
graphic_table[,3]<-gsub("par.dif.chh", "Parental methylation difference: CHH", graphic_table[,3])
graphic_table$variable<-gsub("METHY_cpg","Methylation level: CpG", graphic_table$variable)
graphic_table$variable<-gsub("METHY_chg","Methylation level: CHG", graphic_table$variable)
graphic_table$variable<-gsub("METHY_chh","Methylation level: CHH", graphic_table$variable)
graphic_table$variable<-gsub("genetic_eff","Genetic effects", graphic_table$variable)
graphic_table[,3]<-gsub("par.gen.dist", "Parental sequence divergence", graphic_table[,3])
graphic_table[,3]<-gsub("insertions", "Insertions", graphic_table[,3])
graphic_table[,3]<-gsub("deletions", "Deletions", graphic_table[,3])
graphic_table[,3]<-gsub("inversions", "Inversions", graphic_table[,3])
graphic_table[,3]<-gsub("duplications", "Duplications", graphic_table[,3])

########
unique(graphic_table$variable)

#make transformation of p values
graphic_table<-cbind(graphic_table, NA); colnames(graphic_table)[6]<-"transformed_Pval"
for (i in 1:nrow(graphic_table)){graphic_table[i,6]<-sqrt(as.numeric(graphic_table[i,5])^-1)}
#graphic_table[,6]<-log10(as.numeric(graphic_table[,5]))

unique(graphic_table$variable)

order_for_graphic<-c("Insertions",  "Deletions", "Inversions", "Duplications",
                     "Parental sequence divergence", 
                     "Genetic effects",
                     "Methylation level: CpG", "Methylation level: CHG", "Methylation level: CHH", 
                     "Parental methylation difference: CpG", "Parental methylation difference: CHG", "Parental methylation difference: CHH")                       

order_for_graphic<-order_for_graphic[length(order_for_graphic):1]

#generate relation
pval_rel<-graphic_table[order(as.numeric(graphic_table[,5])), c(5,6)]
p_table<-data.frame(c(0.0001, 0.001, 0.01, 0.05),NA); colnames(p_table)<-c("p_value", "transformed")
for (i in 1:nrow(p_table)){p_table[i,2]<-pval_rel[which(round(pval_rel[,1], digits = 4)==p_table[i,1]),2][1]}

#make data frame and numeric
graphic_table<-as.matrix(graphic_table)
graphic_table[,2]<-as.numeric(graphic_table[,2])/s
graphic_table<-as.data.frame(graphic_table)
for (i in c(2,4,5,6)){graphic_table[,i]<-as.numeric(as.character(graphic_table[,i]))}

LIGHT<-get.color(8,2,15)
DARK<-get.color(8,3,12)
MID<-"white"

pdf("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.2.1_corr_graphic_all.pdf", width = 15, height = 10)

ggplot()+
  geom_tile(data = graphic_table, aes(x = window, y = factor(variable, levels = order_for_graphic), fill = reg_coeff))+
  #geom_point(data = graphic_table, aes(x = window, y = factor(variable, levels = order_for_graphic), color = reg_coeff, size = 0.01*sqrt(P_value^-1))) + #size = log10(P_value))) +
  #scale_size(name = "P value", breaks = p_table[,2]*0.01, labels = c("<0.0001", "0.001", "0.01", "0.05")) +
  scale_fill_gradient2(high = LIGHT, low = DARK, mid = MID, na.value = NA, midpoint = 0, breaks = c(1,0.5,0,-0.5,-1), #name = "Standardized coefficient", 
                       guide = guide_colourbar(title.vjust = 4, ticks.colour = "black", frame.colour = "black")) +
  scale_x_continuous(name = "Physical position (Mbp)", position = "bottom", expand = c(0.01,0.01)) +
  centro_peri +
  
  facet_grid(rows = vars(chr), scales = "free", space = "free",  labeller = labeller(.rows = )) +
  
  theme(
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
    panel.grid.major.y = element_line(size = 0.25, linetype = 'solid', colour = "lightgrey"), 
    panel.grid.major.x = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
    panel.grid.minor.x = element_line(size = 0.25, linetype = 'solid', colour = "grey"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = 'right',
    legend.title = element_blank(),#element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = NA), 
    legend.key.width = unit(2, 'cm'),
  ) 

dev.off()

file.copy("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.2.1_corr_graphic_all.pdf", 
          "/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/")

file.rename("/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/I.2.1_corr_graphic_all.pdf",
            "/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/Figure_2.pdf")

########################################
#Agregated table
s<-1000000
rec_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/B.0_AVERAGE_rec_rates_windows_regions_",s,"_5.RDS", sep = ""))

tabla<-matrix(ncol = 4, nrow = length(order_for_graphic))
row.names(tabla)<-order_for_graphic
colnames(tabla)<-c("Distal_pos","Distal_neg", "Peri_pos", "Peri_neg") #paste(1:7, "H", sep = "")
tabla[]<-0

count_P<-0; count_D<-0
for (c in 1:7){ 
  graphic_table<-sign_var_list[[c]]$graphic_table
  #change names
  graphic_table[,3]<-gsub("par.dif.cpg", "Parental methylation difference: CpG", graphic_table[,3])
  graphic_table[,3]<-gsub("par.dif.chg", "Parental methylation difference: CHG", graphic_table[,3])
  graphic_table[,3]<-gsub("par.dif.chh", "Parental methylation difference: CHH", graphic_table[,3])
  graphic_table$variable<-gsub("METHY_cpg","Methylation level: CpG", graphic_table$variable)
  graphic_table$variable<-gsub("METHY_chg","Methylation level: CHG", graphic_table$variable)
  graphic_table$variable<-gsub("METHY_chh","Methylation level: CHH", graphic_table$variable)
  graphic_table$variable<-gsub("genetic_eff","Genetic effects", graphic_table$variable)
  graphic_table[,3]<-gsub("par.gen.dist", "Parental sequence divergence", graphic_table[,3])
  graphic_table[,3]<-gsub("insertions", "Insertions", graphic_table[,3])
  graphic_table[,3]<-gsub("deletions", "Deletions", graphic_table[,3])
  graphic_table[,3]<-gsub("inversions", "Inversions", graphic_table[,3])
  graphic_table[,3]<-gsub("duplications", "Duplications", graphic_table[,3])
  
  #peri
  peri<-rec_list[[c]][which(rec_list[[c]][,4]=="P"),2]; count_P<-sum(count_P, length(peri)) 
  peri<-which(as.character(graphic_table[,2])%in%peri)
  #peri neg
  peri_neg<-graphic_table[peri[which(as.numeric(graphic_table[peri,4])<0)],]
  tab<-table(peri_neg[,3])  
  for (i in names(tab)){tabla[i,4]<-sum(as.numeric(tabla[i,4]), tab[i])}
  #peri pos
  peri_pos<-graphic_table[peri[which(as.numeric(graphic_table[peri,4])>0)],]
  tab<-table(peri_pos[,3])  
  for (i in names(tab)){tabla[i,3]<-sum(as.numeric(tabla[i,3]), tab[i])}
  
  #distal
  distal<-rec_list[[c]][which(rec_list[[c]][,4]=="D"),2]; count_D<-sum(count_D, length(distal))
  distal<-which(as.character(graphic_table[,2])%in%distal)
  #distal neg
  distal_neg<-graphic_table[distal[which(as.numeric(graphic_table[distal,4])<0)],]
  tab<-table(distal_neg[,3])  
  for (i in names(tab)){tabla[i,2]<-sum(as.numeric(tabla[i,2]), tab[i])}
  #distal pos
  distal_pos<-graphic_table[distal[which(as.numeric(graphic_table[distal,4])>0)],]
  tab<-table(distal_pos[,3])  
  for (i in names(tab)){tabla[i,1]<-sum(as.numeric(tabla[i,1]), tab[i])}
}#c

#make as proportion of windows
tabla[,1]<-as.numeric(tabla[,1])/count_D
tabla[,2]<-as.numeric(tabla[,2])/count_D
tabla[,3]<-as.numeric(tabla[,3])/count_P
tabla[,4]<-as.numeric(tabla[,4])/count_P

tabla[]<-round(as.numeric(tabla[]), digits = 2)

#order variables
tabla<-tabla[nrow(tabla):1,]
tabla<-rbind(c("+","-","+","-"), tabla)
tabla<-rbind(c("Correlation","","Correlation",""), tabla)
tabla<-cbind(row.names(tabla), tabla)
colnames(tabla)<-c("Variable", "", "Distal", "", "Pericentromeric")

write.csv(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.2.2_win_per_significant_variable.csv", row.names = FALSE)

###############################################################

#REDUCED variables
s<-1000000

library(ggplot2)
library(impressionist.colors)

#Prepare table for graphic
sign_var_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.1.5_REDUCED_sw_reg_per_win_significant_variables.RDS")
graphic_table<-sign_var_list[[1]][["graphic_table"]]  

unique(graphic_table$variable)
for (c in 2:7){graphic_table<-rbind(graphic_table,sign_var_list[[c]][["graphic_table"]])}
#change names
graphic_table[,3]<-gsub("SVs_prop", "Total SV load", graphic_table[,3])
graphic_table[,3]<-gsub("par.gen.dist", "Parental sequence divergence", graphic_table[,3])
graphic_table[,3]<-gsub("genetic_eff","Genetic effects", graphic_table[,3])
graphic_table[,3]<-gsub("pdif_METHY", "Ave. parental methylation difference", graphic_table[,3])
graphic_table[,3]<-gsub("METHY","Ave. methylation level", graphic_table[,3])
unique(graphic_table$variable)

#make transformation of p values
graphic_table<-cbind(graphic_table, NA); colnames(graphic_table)[6]<-"transformed_Pval"
for (i in 1:nrow(graphic_table)){graphic_table[i,6]<-sqrt(as.numeric(graphic_table[i,5])^-1)}

order_for_graphic<-c("Total SV load", 
                     "Parental sequence divergence", 
                     "Genetic effects",
                     "Ave. methylation level", 
                     "Ave. parental methylation difference")                       

order_for_graphic<-order_for_graphic[length(order_for_graphic):1]

#generate relation
pval_rel<-graphic_table[order(as.numeric(graphic_table[,5])), c(5,6)]
p_table<-data.frame(c(0.0001, 0.001, 0.01, 0.05),NA); colnames(p_table)<-c("p_value", "transformed")
for (i in 1:nrow(p_table)){p_table[i,2]<-pval_rel[which(round(pval_rel[,1], digits = 4)==p_table[i,1]),2][1]}

#make data frame and numeric
graphic_table<-as.matrix(graphic_table)
graphic_table[,2]<-as.numeric(graphic_table[,2])/s
graphic_table<-as.data.frame(graphic_table)
for (i in c(2,4,5,6)){graphic_table[,i]<-as.numeric(as.character(graphic_table[,i]))}

LIGHT<-get.color(8,2,15)
DARK<-get.color(8,3,12)
MID<-"white"

pdf("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.2.3_REDUCED_corr_graphic_all.pdf", width = 15, height = 10)

ggplot()+
  geom_tile(data = graphic_table, aes(x = window, y = factor(variable, levels = order_for_graphic), fill = reg_coeff))+
  #geom_point(data = graphic_table, aes(x = window, y = factor(variable, levels = order_for_graphic), color = reg_coeff, size = 0.01*sqrt(P_value^-1))) + #size = log10(P_value))) +
  #scale_size(name = "P value", breaks = p_table[,2]*0.01, labels = c("<0.0001", "0.001", "0.01", "0.05")) +
  scale_fill_gradient2(high = LIGHT, low = DARK, mid = MID, na.value = NA, midpoint = 0, breaks = c(1,0.5,0,-0.5,-1), #name = "Standardized coefficient", 
                       guide = guide_colourbar(title.vjust = 4, ticks.colour = "black", frame.colour = "black")) +
  scale_x_continuous(name = "Physical position (Mbp)", position = "bottom", expand = c(0.01,0.01)) +
  centro_peri +
  
  facet_grid(rows = vars(chr), scales = "free", space = "free",  labeller = labeller(.rows = )) +
  
  theme(
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white", linewidth = 0.5, linetype = "solid"),
    panel.grid.major.y = element_line(linewidth = 0.25, linetype = 'solid', colour = "lightgrey"), 
    panel.grid.major.x = element_line(linewidth = 0.25, linetype = 'solid', colour = "grey"), 
    panel.grid.minor.x = element_line(linewidth = 0.25, linetype = 'solid', colour = "grey"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = 'right',
    legend.title = element_blank(),#element_text(size = 9),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA), 
    legend.key.width = unit(2.5, 'cm'),
    text = element_text(size=18)
  ) 

dev.off()

file.copy("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.2.3_REDUCED_corr_graphic_all.pdf", 
          "/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/")

file.rename("/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/I.2.3_REDUCED_corr_graphic_all.pdf",
            "/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/Figure_2_REDUCED.pdf")

########################################
#Agregated table
s<-1000000
rec_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/B.0_AVERAGE_rec_rates_windows_regions_",s,"_5.RDS", sep = ""))

tabla<-matrix(ncol = 4, nrow = length(order_for_graphic))
row.names(tabla)<-order_for_graphic
colnames(tabla)<-c("Distal_pos","Distal_neg", "Peri_pos", "Peri_neg") #paste(1:7, "H", sep = "")
tabla[]<-0

count_P<-0; count_D<-0
for (c in 1:7){ 
  graphic_table<-sign_var_list[[c]]$graphic_table
  graphic_table[,3]<-gsub("SVs_prop", "Total SVs", graphic_table[,3])
  graphic_table[,3]<-gsub("par.gen.dist", "Parental sequence divergence", graphic_table[,3])
  graphic_table[,3]<-gsub("genetic_eff","Genetic effects", graphic_table[,3])
  graphic_table[,3]<-gsub("pdif_METHY", "Ave. parental methylation difference", graphic_table[,3])
  graphic_table[,3]<-gsub("METHY","Ave. methylation level", graphic_table[,3])
  
  #peri
  peri<-rec_list[[c]][which(rec_list[[c]][,4]=="P"),2]; count_P<-sum(count_P, length(peri)) 
  peri<-which(as.character(graphic_table[,2])%in%peri)
  #peri neg
  peri_neg<-graphic_table[peri[which(as.numeric(graphic_table[peri,4])<0)],]
  tab<-table(peri_neg[,3])  
  for (i in names(tab)){tabla[i,4]<-sum(as.numeric(tabla[i,4]), tab[i])}
  #peri pos
  peri_pos<-graphic_table[peri[which(as.numeric(graphic_table[peri,4])>0)],]
  tab<-table(peri_pos[,3])  
  for (i in names(tab)){tabla[i,3]<-sum(as.numeric(tabla[i,3]), tab[i])}
  
  #distal
  distal<-rec_list[[c]][which(rec_list[[c]][,4]=="D"),2]; count_D<-sum(count_D, length(distal))
  distal<-which(as.character(graphic_table[,2])%in%distal)
  #distal neg
  distal_neg<-graphic_table[distal[which(as.numeric(graphic_table[distal,4])<0)],]
  tab<-table(distal_neg[,3])  
  for (i in names(tab)){tabla[i,2]<-sum(as.numeric(tabla[i,2]), tab[i])}
  #distal pos
  distal_pos<-graphic_table[distal[which(as.numeric(graphic_table[distal,4])>0)],]
  tab<-table(distal_pos[,3])  
  for (i in names(tab)){tabla[i,1]<-sum(as.numeric(tabla[i,1]), tab[i])}
}#c

#make as proportion of windows
tabla[,1]<-as.numeric(tabla[,1])/count_D
tabla[,2]<-as.numeric(tabla[,2])/count_D
tabla[,3]<-as.numeric(tabla[,3])/count_P
tabla[,4]<-as.numeric(tabla[,4])/count_P

tabla[]<-round(as.numeric(tabla[]), digits = 2)

#order variables
tabla<-tabla[nrow(tabla):1,]
tabla<-rbind(c("+","-","+","-"), tabla)
tabla<-rbind(c("Correlation","","Correlation",""), tabla)
tabla<-cbind(row.names(tabla), tabla)
colnames(tabla)<-c("Variable", "", "Distal", "", "Pericentromeric")

write.csv(tabla, "/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.2.4_REDUCED_win_per_significant_variable.csv", row.names = FALSE)
