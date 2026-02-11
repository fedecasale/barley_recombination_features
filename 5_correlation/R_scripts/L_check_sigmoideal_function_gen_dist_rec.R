library(ggplot2)
library(impressionist.colors)

#0-get data
s<-1000000

# rec_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))
# GEN_DIST_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/D_genetic_distances/D.2_populations_genetic_distances.RDS", sep = ""))
# 
# variables<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/j.1.1_standarized_variables_across_genome.RDS")
# rec_chr_pop<-variables$rec.rate
# dist_chr_pop<-variables$gen.dist
# 
# y<-as.numeric(rec_list[[5]][])
# x<-as.numeric(GEN_DIST_list[[5]][])
# plot(y~x)
# 
# 
# fit1 <- lm(y~x)
# fit2 <- lm(y~poly(x,2,raw=TRUE))
# fit3 <- lm(y~poly(x,3,raw=TRUE))
# fit4 <- lm(y~poly(x,4,raw=TRUE))
# fit5 <- lm(y~poly(x,5,raw=TRUE))
# fit5 <- lm(y~poly(x,1000,raw=TRUE))
# 
# 
# xaxis<-seq(0, 0.5, length=20)
# lines(xaxis, predict(fit1, data.frame(x=xaxis)), col='green')
# lines(xaxis, predict(fit2, data.frame(x=xaxis)), col='red')
# lines(xaxis, predict(fit3, data.frame(x=xaxis)), col='blue')
# lines(xaxis, predict(fit4, data.frame(x=xaxis)), col='pink')
# lines(xaxis, predict(fit5, data.frame(x=xaxis)), col='yellow')
# 
# summary(fit1)$adj.r.squared
# summary(fit2)$adj.r.squared
# summary(fit3)$adj.r.squared
# summary(fit4)$adj.r.squared
# summary(fit5)$adj.r.squared
# 
# #order<-order(dist_chr_pop, decreasing = FALSE)
# #plot(rec_chr_pop[order]~dist_chr_pop[order])

rec_list<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/C.1_rec_rates_all_pops_per_window_",s,".RDS", sep = ""))
rec_list_AVE<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/B.0_AVERAGE_rec_rates_windows_regions_",s,".RDS", sep = ""))
standarized_variables<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.1.3_REDUCED_standarized_variables_per_window.RDS")
significant_variables_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.1.5_REDUCED_sw_reg_per_win_significant_variables.RDS")
#standarized_variables<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.1.0_standarized_variables_per_window.RDS")
#significant_variables_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/I.1.2_sw_reg_per_win_significant_variables.RDS")

variables<-c("rec","par.gen.dist", "METHY", "genetic_eff", "SVs_prop")

scale<-function(x){(x-min(win_values, na.rm = TRUE))/(max(win_values, na.rm = TRUE)-min(win_values, na.rm = TRUE))}


###################### par.gen.dist ###################################

variable<-"par.gen.dist"

graphic_table<-matrix(ncol = 7, nrow = 0)

for (c in 1:7){
  
chr_win<-significant_variables_list[[c]][[3]]  
chr_win<-as.character(chr_win[which(chr_win[,3]==variable),"window"])

#1-take out pericentromere
distal_wins<-rec_list_AVE[[c]][which(rec_list_AVE[[c]][,4]=="D"),2]
chr_win<-chr_win[which(chr_win%in%distal_wins)]

graphic_table_chr<-matrix(ncol = 7, nrow = 0)

#2-scale within window values
for (w in chr_win){
win_table<-matrix(ncol = 2, nrow = 45)
win_table[,1]<-c
win_table[,2]<-w
for (f in variables){ 
if (f == "rec"){
win_values<-as.numeric(rec_list[[c]][w,])
}else{
win_values<-standarized_variables[[c]][[w]][,f]
}
if (length(which(win_values==0))!=length(win_values)){
win_values<-unlist(lapply(win_values, scale))
}
#win_values<-round(win_values, digits = 2)
win_table<-cbind(win_table, win_values)
}#f
graphic_table_chr <- rbind(graphic_table_chr, win_table)
}#w
colnames(graphic_table_chr)<-c("chr", "win",variables)

#3-clean outliers
for (f in c("rec",variable)){
limits<-quantile(as.numeric(graphic_table_chr[,f]), probs = seq(0,1, 0.05), na.rm = TRUE)[c(2,20)]
graphic_table_chr<-graphic_table_chr[-which(as.numeric(graphic_table_chr[,f])<=limits[1]),]
graphic_table_chr<-graphic_table_chr[-which(as.numeric(graphic_table_chr[,f])>=limits[2]),]
}#f

#4-re-normalize values
for (f in variables){
win_values<-as.numeric(graphic_table_chr[,f])
graphic_table_chr[,f]<-unlist(lapply(win_values, scale))
}

graphic_table <- rbind(graphic_table, graphic_table_chr)

}#c

colnames(graphic_table)<-c("chr", "win","erec","ddist", "bmet","cgen.eff", "aSV")

graphic_table<-as.data.frame(graphic_table)
for (c in 3:7){graphic_table[,c]<-as.numeric(as.character(graphic_table[,c]))}

#second scale
max_rec<-max(graphic_table$erec, na.rm = TRUE)

#define colors
orange<-get.color(8,1,11)
blue<-"#0073CF"
mustard<-"#efdc75" #"#FFD000" #get.color(8,1,12)
green<-get.color(8,2,5)
purple<-"#702963"

target<-"ddist"
variables_2<-colnames(graphic_table)[-c(1,2,which(colnames(graphic_table)==target))]
graphic_table_2<-matrix(ncol = 3, nrow = 0)
for (v in variables_2){graphic_table_2<-rbind(graphic_table_2, cbind(graphic_table[,target],graphic_table[,v],v))}
colnames(graphic_table_2)<-c(target,"var_value","variable")

graphic_table_2<-as.data.frame(graphic_table_2)
graphic_table_2[[target]]<-as.numeric(as.character(graphic_table_2[[target]]))
graphic_table_2$var_value<-as.numeric(as.character(graphic_table_2$var_value))

pdf("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/L.1_seq_div.pdf")
ggplot() +
  geom_point(data = graphic_table, aes(x = ddist, y = erec), size = 0.01, color = "#393A3F")+
  geom_smooth(data = graphic_table_2, aes(x = ddist, y = var_value, color = variable), method = "gam", linewidth = 1.5, se = FALSE) +
  scale_y_continuous(name = "Relative recombination rate", position = "left", expand = c(0.01,0.01), 
                     sec.axis = sec_axis(trans=~./max_rec, name="Relative total SV load/genetic effects/methylation level")) +
  scale_x_continuous(name = "Relative sequence parental divergence", position = "bottom", expand = c(0.005,0.005)) +
  
  scale_color_manual(values = c(blue, mustard, purple, green), name = NULL,
                     breaks = c("erec", c("aSV", "bmet", "cgen.eff")),
                     labels = c("Recombination rate", "Total SVs load", "Methylation level", "Genetic effects")) +

  #facet_grid(rows = vars(chr), scales = "free", space = "free") +
  
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background  =  element_blank(),
        legend.key = element_rect(fill = NA), #needs se=FALSE
        legend.position = "bottom"
        )

dev.off()

###################### SVs_prop ###################################

variable<-"SVs_prop"

graphic_table<-matrix(ncol = 7, nrow = 0)

for (c in 1:7){
  
  chr_win<-significant_variables_list[[c]][[3]]  
  chr_win<-as.character(chr_win[which(chr_win[,3]==variable),"window"])
  
  #1-take out pericentromere
  distal_wins<-rec_list_AVE[[c]][which(rec_list_AVE[[c]][,4]=="D"),2]
  chr_win<-chr_win[which(chr_win%in%distal_wins)]
  
  graphic_table_chr<-matrix(ncol = 7, nrow = 0)
  
  #2-scale within window values
  for (w in chr_win){
    win_table<-matrix(ncol = 2, nrow = 45)
    win_table[,1]<-c
    win_table[,2]<-w
    for (f in variables){ 
      if (f == "rec"){
        win_values<-as.numeric(rec_list[[c]][w,])
      }else{
        win_values<-standarized_variables[[c]][[w]][,f]
      }
      if (length(which(win_values==0))!=length(win_values)){
        win_values<-unlist(lapply(win_values, scale))
      }
      #win_values<-round(win_values, digits = 2)
      win_table<-cbind(win_table, win_values)
    }#f
    graphic_table_chr <- rbind(graphic_table_chr, win_table)
  }#w
  colnames(graphic_table_chr)<-c("chr", "win",variables)
  
  #3-clean outliers
  for (f in c("rec",variable)){
    limits<-quantile(as.numeric(graphic_table_chr[,f]), probs = seq(0,1, 0.05), na.rm = TRUE)[c(2,20)]
    graphic_table_chr<-graphic_table_chr[-which(as.numeric(graphic_table_chr[,f])<=limits[1]),]
    graphic_table_chr<-graphic_table_chr[-which(as.numeric(graphic_table_chr[,f])>=limits[2]),]
  }#f
  
  #4-re-normalize values
  for (f in variables){
    win_values<-as.numeric(graphic_table_chr[,f])
    graphic_table_chr[,f]<-unlist(lapply(win_values, scale))
  }
  
  graphic_table <- rbind(graphic_table, graphic_table_chr)
  
}#c

colnames(graphic_table)<-c("chr", "win","erec","ddist", "bmet","cgen.eff", "aSV")

graphic_table<-as.data.frame(graphic_table)
for (c in 3:7){graphic_table[,c]<-as.numeric(as.character(graphic_table[,c]))}

#second scale
max_rec<-max(graphic_table$erec, na.rm = TRUE)

#define colors
orange<-get.color(8,1,11)
blue<-"#0073CF"
mustard<-"#efdc75" #"#FFD000" #get.color(8,1,12)
musgo<-get.color(8,1,1)
purple<-"#702963"

target<-"aSV"
variables_2<-colnames(graphic_table)[-c(1,2,which(colnames(graphic_table)==target))]
graphic_table_2<-matrix(ncol = 3, nrow = 0)
for (v in variables_2){graphic_table_2<-rbind(graphic_table_2, cbind(graphic_table[,target],graphic_table[,v],v))}
colnames(graphic_table_2)<-c(target,"var_value","variable")

graphic_table_2<-as.data.frame(graphic_table_2)
graphic_table_2[[target]]<-as.numeric(as.character(graphic_table_2[[target]]))
graphic_table_2$var_value<-as.numeric(as.character(graphic_table_2$var_value))

pdf("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/L.2_SV_rec.pdf")
ggplot() +
  geom_point(data = graphic_table, aes(x = aSV, y = erec), size = 0.01, color = "#393A3F")+
  geom_smooth(data = graphic_table_2, aes(x = aSV, y = var_value, color = variable), method = "gam", linewidth = 1.5, se = FALSE) +
  scale_y_continuous(name = "Relative recombination rate", position = "left", expand = c(0.01,0.01), 
                     sec.axis = sec_axis(trans=~./max_rec, name="Relative sequence parental divergence/genetic effects/methylation level")) +
  scale_x_continuous(name = "Relative total SVs load", position = "bottom", expand = c(0.005,0.005)) +
  
  scale_color_manual(values = c(blue, orange, purple, green), name = NULL,
                     breaks = c("erec", c("ddist", "bmet", "cgen.eff")),
                     labels = c("Recombination rate", "Sequence parental divergence", "Methylation level", "Genetic effects")) +
  
  #facet_grid(rows = vars(chr), scales = "free", space = "free") +
  
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background  =  element_blank(),
        legend.key = element_rect(fill = NA), #needs se=FALSE
        legend.position = "bottom"
        )

dev.off()


###################### METHY ###################################

variable<-"METHY"

graphic_table<-matrix(ncol = 7, nrow = 0)

for (c in 1:7){
  
  chr_win<-significant_variables_list[[c]][[3]]  
  chr_win<-as.character(chr_win[which(chr_win[,3]==variable),"window"])
  
  #1-take out pericentromere
  distal_wins<-rec_list_AVE[[c]][which(rec_list_AVE[[c]][,4]=="D"),2]
  chr_win<-chr_win[which(chr_win%in%distal_wins)]
  
  graphic_table_chr<-matrix(ncol = 7, nrow = 0)
  
  #2-scale within window values
  for (w in chr_win){
    win_table<-matrix(ncol = 2, nrow = 45)
    win_table[,1]<-c
    win_table[,2]<-w
    for (f in variables){ 
      if (f == "rec"){
        win_values<-as.numeric(rec_list[[c]][w,])
      }else{
        win_values<-standarized_variables[[c]][[w]][,f]
      }
      if (length(which(win_values==0))!=length(win_values)){
        win_values<-unlist(lapply(win_values, scale))
      }
      #win_values<-round(win_values, digits = 2)
      win_table<-cbind(win_table, win_values)
    }#f
    graphic_table_chr <- rbind(graphic_table_chr, win_table)
  }#w
  colnames(graphic_table_chr)<-c("chr", "win",variables)
  
  #3-clean outliers
  for (f in c("rec",variable)){
    limits<-quantile(as.numeric(graphic_table_chr[,f]), probs = seq(0,1, 0.05), na.rm = TRUE)[c(2,20)]
    graphic_table_chr<-graphic_table_chr[-which(as.numeric(graphic_table_chr[,f])<=limits[1]),]
    graphic_table_chr<-graphic_table_chr[-which(as.numeric(graphic_table_chr[,f])>=limits[2]),]
  }#f
  
  #4-re-normalize values
  for (f in variables){
    win_values<-as.numeric(graphic_table_chr[,f])
    graphic_table_chr[,f]<-unlist(lapply(win_values, scale))
  }
  
  graphic_table <- rbind(graphic_table, graphic_table_chr)
  
}#c

colnames(graphic_table)<-c("chr", "win","erec","ddist", "bmet","cgen.eff", "aSV")

graphic_table<-as.data.frame(graphic_table)
for (c in 3:7){graphic_table[,c]<-as.numeric(as.character(graphic_table[,c]))}

#second scale
max_rec<-max(graphic_table$erec, na.rm = TRUE)

#define colors
orange<-get.color(8,1,11)
blue<-"#0073CF"
mustard<-"#efdc75" #"#FFD000" #get.color(8,1,12)
musgo<-get.color(8,1,1)
purple<-"#702963"

target<-"bmet"
variables_2<-colnames(graphic_table)[-c(1,2,which(colnames(graphic_table)==target))]
graphic_table_2<-matrix(ncol = 3, nrow = 0)
for (v in variables_2){graphic_table_2<-rbind(graphic_table_2, cbind(graphic_table[,target],graphic_table[,v],v))}
colnames(graphic_table_2)<-c(target,"var_value","variable")

graphic_table_2<-as.data.frame(graphic_table_2)
graphic_table_2[[target]]<-as.numeric(as.character(graphic_table_2[[target]]))
graphic_table_2$var_value<-as.numeric(as.character(graphic_table_2$var_value))

pdf("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/L.3_met_rec.pdf")
ggplot() +
  geom_point(data = graphic_table, aes(x = bmet, y = erec), size = 0.01, color = "#393A3F")+
  geom_smooth(data = graphic_table_2, aes(x = bmet, y = var_value, color = variable), method = "gam", linewidth = 1.5, se = FALSE) +
  scale_y_continuous(name = "Relative recombination rate", position = "left", expand = c(0.01,0.01), 
                     sec.axis = sec_axis(trans=~./max_rec, name="Relative sequence parental divergence/genetic effects/total SVs load")) +
  scale_x_continuous(name = "Relative methylation level", position = "bottom", expand = c(0.005,0.005)) +
  
  scale_color_manual(values = c(blue, orange, mustard, green), name = NULL,
                     breaks = c("erec", c("ddist", "aSV", "cgen.eff")),
                     labels = c("Recombination rate", "Sequence parental divergence", "Total SVs load", "Genetic effects")) +
  
  #facet_grid(rows = vars(chr), scales = "free", space = "free") +
  
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background  =  element_blank(),
        legend.key = element_rect(fill = NA), #needs se=FALSE
        legend.position = "bottom"
        )
dev.off()



###################### genetic_eff ###################################

variable<-"genetic_eff"

graphic_table<-matrix(ncol = 7, nrow = 0)

for (c in 1:7){
  
  chr_win<-significant_variables_list[[c]][[3]]  
  chr_win<-as.character(chr_win[which(chr_win[,3]==variable),"window"])
  
  #1-take out pericentromere
  distal_wins<-rec_list_AVE[[c]][which(rec_list_AVE[[c]][,4]=="D"),2]
  chr_win<-chr_win[which(chr_win%in%distal_wins)]
  
  graphic_table_chr<-matrix(ncol = 7, nrow = 0)
  
  #2-scale within window values
  for (w in chr_win){
    win_table<-matrix(ncol = 2, nrow = 45)
    win_table[,1]<-c
    win_table[,2]<-w
    for (f in variables){ 
      if (f == "rec"){
        win_values<-as.numeric(rec_list[[c]][w,])
      }else{
        win_values<-standarized_variables[[c]][[w]][,f]
      }
      if (length(which(win_values==0))!=length(win_values)){
        win_values<-unlist(lapply(win_values, scale))
      }
      #win_values<-round(win_values, digits = 2)
      win_table<-cbind(win_table, win_values)
    }#f
    graphic_table_chr <- rbind(graphic_table_chr, win_table)
  }#w
  colnames(graphic_table_chr)<-c("chr", "win",variables)
  
  #3-clean outliers
  for (f in c("rec",variable)){
    limits<-quantile(as.numeric(graphic_table_chr[,f]), probs = seq(0,1, 0.05), na.rm = TRUE)[c(2,20)]
    graphic_table_chr<-graphic_table_chr[-which(as.numeric(graphic_table_chr[,f])<=limits[1]),]
    graphic_table_chr<-graphic_table_chr[-which(as.numeric(graphic_table_chr[,f])>=limits[2]),]
  }#f
  
  #4-re-normalize values
  for (f in variables){
    win_values<-as.numeric(graphic_table_chr[,f])
    graphic_table_chr[,f]<-unlist(lapply(win_values, scale))
  }
  
  graphic_table <- rbind(graphic_table, graphic_table_chr)
  
}#c

colnames(graphic_table)<-c("chr", "win","erec","ddist", "bmet","cgen.eff", "aSV")

graphic_table<-as.data.frame(graphic_table)
for (c in 3:7){graphic_table[,c]<-as.numeric(as.character(graphic_table[,c]))}

#second scale
max_rec<-max(graphic_table$erec, na.rm = TRUE)

target<-"cgen.eff"
variables_2<-colnames(graphic_table)[-c(1,2,which(colnames(graphic_table)==target))]
graphic_table_2<-matrix(ncol = 3, nrow = 0)
for (v in variables_2){graphic_table_2<-rbind(graphic_table_2, cbind(graphic_table[,target],graphic_table[,v],v))}
colnames(graphic_table_2)<-c(target,"var_value","variable")

graphic_table_2<-as.data.frame(graphic_table_2)
graphic_table_2[[target]]<-as.numeric(as.character(graphic_table_2[[target]]))
graphic_table_2$var_value<-as.numeric(as.character(graphic_table_2$var_value))

pdf("/home/fcasale/Desktop/Paper_2/3_RILs/5_correlation/Results/L.4_gen_eff_rec.pdf")
ggplot() +
  geom_point(data = graphic_table, aes(x = cgen.eff, y = erec), size = 0.01, color = "#393A3F")+
  geom_smooth(data = graphic_table_2, aes(x = cgen.eff, y = var_value, color = variable), method = "gam", linewidth = 1.5, se = FALSE) +
  scale_y_continuous(name = "Relative recombination rate", position = "left", expand = c(0.01,0.01), 
                     sec.axis = sec_axis(trans=~./max_rec, name="Relative sequence parental divergence/total SVs load/methylation level")) +
  scale_x_continuous(name = "Relative genetic effects", position = "bottom", expand = c(0.005,0.005)) +
  
  scale_color_manual(values = c(blue, orange, mustard, purple), name = NULL,
                     breaks = c("erec", c("ddist", "aSV", "bmet")),
                     labels = c("Recombination rate", "Sequence parental divergence", "Total SVs load", "Methylation level")) +
  
  #facet_grid(rows = vars(chr), scales = "free", space = "free") +
  
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background  =  element_blank(),
        legend.key = element_rect(fill = NA), #needs se=FALSE
        legend.position = "bottom"
        )
dev.off()
