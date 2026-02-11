library(xlsx)
library(ggplot2)
library(ggforce)
library(cowplot)
library(ggh4x)
library(ggrepel)
library(dplyr)

results_table<-read.csv("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.6.0_window_groupping_means_and_letters.csv")

variables<-unique(results_table[,1])
window_types<-unique(results_table[,2])
Populations<-unique(results_table[3])

window_around<-as.character(1:11)
colnames(results_table)[4:ncol(results_table)]<-window_around

graphic_table<-cbind(results_table[,1:3], paste(1), results_table[,4]) 
colnames(graphic_table)<-c("variable","window_type", "population", "window","mean")
for (w in 5:14){
w_table<-cbind(results_table[,1:3], paste(w-3), results_table[,w])
colnames(w_table)<-colnames(graphic_table)
graphic_table<-rbind(graphic_table, w_table)
}
mean_letter<-strsplit(as.character(graphic_table$mean), " ")
graphic_table$mean<-unlist(lapply(mean_letter, function(x){return(x[1])}))
graphic_table$letter<-unlist(lapply(mean_letter, function(x){return(x[2])}))
colnames(graphic_table)[6]<-"letter"

# graphic_table$letter_height<-NA
# for (v in unique(graphic_table$variable)){  
# #PERC<-max(as.numeric(graphic_table$mean[which(graphic_table$variable==v)]))*0.01
# graphic_table$letter_height[which(graphic_table$variable==v)]<-as.numeric(graphic_table$mean[which(graphic_table$variable==v)])+PERC
# }

#solo le pongo repel a las pegadas (si pongo todo en repel, repel te pone lejos las que ya estan bien)
graphic_table$letter_reppel<-NA
graphic_table$letter_no_reppel<-NA
new_graphic_table<-graphic_table[1,]; new_graphic_table[]<-NA
for (v in unique(graphic_table$variable)){
if (v %in% c("AVE", "cpg", "chg", "chh")){threshold<-0.006}  
if (v == "SVs"){threshold<-0.004}
if (v == "gene_density"){threshold<-0.01}
for (p in unique(graphic_table$population)){  
for (w in unique(graphic_table$window)){
tabla<-graphic_table %>% filter(variable == v) %>% filter(population == p)  %>% filter(window == w)
tabla_sub<-sapply(as.numeric(tabla$mean),"-",as.numeric(tabla$mean))
tabla_sub[]<-abs(tabla_sub[])
colnames(tabla_sub)<-tabla$window_type; row.names(tabla_sub)<-tabla$window_type
to.repel<-c()
for (i in colnames(tabla_sub)){
to.repel_i<-row.names(tabla_sub)[which(tabla_sub[,i]<threshold)]  
to.repel_i<-to.repel_i[-which(to.repel_i==i)]
if (length(to.repel_i)>0){to.repel_i<-c(to.repel_i, i)}
to.repel<-c(to.repel, to.repel_i)
}
to.repel<-unique(to.repel)
if (length(to.repel)>0){
lowest_overlapping_window_type<-min(as.numeric(tabla$mean[which(tabla$window_type%in%to.repel)]))
to.not.repel<-tabla$window_type[which(tabla$mean==lowest_overlapping_window_type)][1]
to.repel<-to.repel[-which(to.repel%in%to.not.repel)]
tabla$letter_reppel[which(tabla$window_type%in%to.repel)]<-tabla$mean[which(tabla$window_type%in%to.repel)]
}
tabla$letter_no_reppel[which(is.na(tabla$letter_reppel))]<-tabla$mean[which(is.na(tabla$letter_reppel))]
new_graphic_table<-rbind(new_graphic_table, tabla)
}#w
}#p
}#v
graphic_table<-new_graphic_table[-1,]



#I will take out one window from each extreme
graphic_table<-graphic_table[-which(graphic_table$window%in%c("1","11")),]
#add order for window types in each window
graphic_table$window_order<-NA
order_window_types<-c("coldspots_proximal", "coldspots_telomeric", "hotspots_pericentromeric", "hotspots_proximal", "hotspots_telomeric")       
for (o in 1:length(order_window_types)){
o_windows<-graphic_table$window[which(graphic_table$window_type==order_window_types[o])]
o_windows<-(as.numeric(o_windows)*10)-60
graphic_table$window_order[which(graphic_table$window_type==order_window_types[o])]<-paste(o_windows) #paste(o_windows+(2*o), sep = "")
}

graphic_table<-as.data.frame(graphic_table)
graphic_table$mean<-as.numeric(as.character(graphic_table$mean))
graphic_table$window_order<-as.numeric(as.character(graphic_table$window_order))
#graphic_table$letter_height<-as.numeric(as.character(graphic_table$letter_height))
graphic_table$letter_reppel<-as.numeric(as.character(graphic_table$letter_reppel))
graphic_table$letter_no_reppel<-as.numeric(as.character(graphic_table$letter_no_reppel))

write.csv(graphic_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.6.1_graphic_table.csv", row.names = FALSE)

#order variables and set labels
#graphic_table$variable = factor(graphic_table$variable, levels=c('AVE','SVs','gene_density'))
graphic_table$variable = factor(graphic_table$variable, levels=c('cpg','chg','chh'))
graphic_table<-graphic_table[complete.cases(graphic_table$variable),]
#variable_label<-as_labeller(c(`AVE` = "Average methylation level", `SVs` = "Total SVs span", `gene_density` = "Gene density"))
variable_label<-as_labeller(c(`cpg` = "CpG", `chg` = "CHG", `chh` = "CHH"))

#color_cold<-c("#000080", "#4682B4")
color_cold<-c("#1b8a5a", "#1d4877")
color_hot<-c("#f68838", "#ee3e32", "#fbb021")
dark_grey<-"#555555"
light_grey<-"#55555560"

#pdf("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.6.2_windows_around_combination_NEW.pdf", width = 10, height = 7)

ggplot()+
  
  geom_vline(xintercept = 0, linetype = "dashed", color = light_grey) +
  
  geom_line(data = graphic_table, aes(x = window_order, y = mean, color = window_type), linewidth = 0.4) +
  
  geom_point(data = graphic_table, aes(x = window_order, y = mean), size = 3, color = "white") +
  
  geom_text(data = graphic_table, aes(x = window_order, y = letter_no_reppel, label = letter, color = window_type, angle = 30), size = 2.5) +
  geom_text_repel(data = graphic_table, aes(x = window_order, y = letter_reppel, label = letter, color = window_type, angle = 30), 
                  box.padding = 0.01, size = 2.5) + 

  scale_y_continuous(position = "right") +
  
  scale_x_continuous(name = "Physical distance (kb)",
                     limits = c(-40, 40),
                     breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                     labels = c(-40, -30, -20, -10, 0, 10, 20, 30, 40))+
  
  scale_color_manual(values = c(color_cold, color_hot), name = "Differentiated recombination region",
                     breaks = c("coldspots_proximal" ,  "coldspots_telomeric" , "hotspots_pericentromeric", "hotspots_proximal", "hotspots_telomeric"),
                     labels = c("Coldspot-distal proximal" ,  "Coldspot-distal telomeric" , "Hotspot-pericentromeric", "Hotspot-distal proximal", "Hotspot-distal telomeric")
  ) +
  
  facet_grid(variable~population, scales = "free", space = "fixed", switch = "y",
             labeller = labeller(variable = variable_label)) +
  
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_blank(),
    axis.ticks = element_line(),
    legend.key = element_blank(),
    axis.title.y.right = element_blank()
  )

# p + ggh4x::facetted_pos_scales(y = list(
#   variable == "methy" ~ scale_y_continuous(limits = c(0, NA), position = "right"),
#   variable == "SVs" ~ scale_y_continuous(limits = c(0, NA), position = "right"),
#   variable == "genes" ~ scale_y_continuous(limits = c(0,0.05,0.1,0.15,0.2), position = "right")
# ))

# ggdraw(p) +
#   draw_label("Methylation level", x = 0.01, y = 0.58, angle = 90, size = 9) +
#   draw_label("Chromosome region windows", x = 0.21, y = 0.99, size = 9) +
#   draw_label("Coldspot windows", x = 0.505, y = 0.99, size = 9) +
#   draw_label("Hotspot windows", x = 0.79, y = 0.99, size = 9) #+

#dev.off()

file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/E.6.2_windows_around_combination.pdf")
file.copy("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.6.2_windows_around_combination.pdf",
          "/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/E.6.2_windows_around_combination.pdf")


