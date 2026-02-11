library(ggplot2)
library(xlsx)
library(ggforce)
library(cowplot)
library(ggh4x)

SVs<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.1.1_graphic_list.RDS")
methy<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.2.1_graphic_list.RDS")
for (i in names(methy)){methy[[i]]<-methy[[i]][["AVE"]]}
genes<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.5.1_graphic_list.RDS")

graphic_list<-list()  
graphic_list[["SVs"]]<-SVs
graphic_list[["methy"]]<-methy
graphic_list[["genes"]]<-genes

saveRDS(graphic_list, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.6.1_graphic_list.RDS")

variables<-names(graphic_list)
window_types<-names(graphic_list[[1]])
Populations<-names(graphic_list[[1]][[1]])


graphic_table<-matrix(ncol = 5, nrow=0); colnames(graphic_table)<-c("variable","Population", "window_type", "value", "window")
for (v in variables){
for (w in window_types){
  for (p in Populations){
    #for (c in 1:7){
    for (i in 1:ncol(graphic_list[[v]][[w]][[p]])){
      col2win<-(i-6)*10
      #graphic_table<-rbind(graphic_table,cbind(p, w, graphic_list[[w]][[p]][c,i], c, col2win))
      graphic_table<-rbind(graphic_table,cbind(v, p, w, mean(graphic_list[[v]][[w]][[p]][,i]), col2win))
    }#i
    #}#c
  }#p
}#w
}#v  
graphic_table<-as.data.frame(graphic_table)
graphic_table$value<-as.numeric(as.character(graphic_table$value))
graphic_table$value[which(graphic_table$value<0)]<-0
graphic_table$window<-as.numeric(as.character(graphic_table$window))
#graphic_table$chr<-paste(graphic_table$chr, "H", sep = "")
graphic_table<-graphic_table[-which(abs(graphic_table$window)==50),]

#order variables and set labels
graphic_table$variable = factor(graphic_table$variable, levels=c('methy','SVs','genes'))
variable_label<-as_labeller(c(`methy` = "Average methylation level", `SVs` = "Total SVs span", `genes` = "Gene density"))

#color_cold<-c("#000080", "#4682B4")
color_cold<-c("#1b8a5a", "#1d4877")
color_hot<-c("#f68838", "#ee3e32", "#fbb021")



pdf("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.6.2_windows_around_combination.pdf", width = 10, height = 7)

ggplot()+
  geom_line(data = graphic_table, aes(x = window, y = value, color = window_type)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  scale_y_continuous(position = "right") +
  
  scale_x_continuous(name = "Physical distance (kb)",
                     limits = c(-40, 40),
                     breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                     labels = c(-40, -30, -20, -10, 0, 10, 20, 30, 40))+
  
  scale_color_manual(values = c(color_cold, color_hot), name = "Window type",
                     breaks = c("coldspots_proximal" ,  "coldspots_telomeric" , "hotspots_pericentromeric", "hotspots_proximal", "hotspots_telomeric"),
                     labels = c("Coldspots distal proximal" ,  "Coldspots distal telomeric" , "Hotspots pericentromeric", "Hotspots distal proximal", "Hotspots distal telomeric")
  ) +
  
  facet_grid(variable~Population, scales = "free", space = "fixed", switch = "y",
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

dev.off()

file.remove("/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/E.6.2_windows_around_combination.pdf")
file.copy("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.6.2_windows_around_combination.pdf",
          "/home/fcasale/Desktop/Paper_2/3_RILs/6_graphics_for_paper/Results/E.6.2_windows_around_combination.pdf")


