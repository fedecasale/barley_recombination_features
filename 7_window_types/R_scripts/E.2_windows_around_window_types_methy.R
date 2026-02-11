# library(ggplot2)
# library(ggh4x)

Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

methy_level<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/4_methylation/Results/B_methylation_per_population/B.1_methylation_per_population_win=10000.RDS", sep = ""))
methy_level_AVE<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/4_methylation/Results/B_methylation_per_population/B.2.1_methy_unified_per_pop_win=10000.RDS", sep = ""))
for (c in 1:7){
methy_level[[c]][["AVE"]]<-methy_level_AVE[[c]]    
}
met_context<-names(methy_level[[1]])

methy_list_windows<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.2.2_methylation_prop_per_window_type.RDS")
methy_list_windows$AVE<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.5.1_methylation_prop_per_window_type_UNI.RDS")


coldspots_regions<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/A.2_coldspots/A.2.2_coldspots_per_pop_REGION.RDS")
#get coldspots in regions
coldspots_in_regions_list<-list()
for (p in Populations){
  coldspots_in_regions_list[[p]]<-list()
  for (c in 1:7){
    coldspots_in_regions<-c()  
    for (r in 1:nrow(coldspots_regions[[p]][[c]])){
      coldspots_in_regions<-c(coldspots_in_regions,seq(coldspots_regions[[p]][[c]][r,1], coldspots_regions[[p]][[c]][r,2]+1, by = 10000))
    }#r
    coldspots_in_regions_list[[p]][[c]]<-coldspots_in_regions
  }#c
}#p

#divide regions among proximal and telomeric
for (p in Populations){
  for (c in 1:7){
    regions <-coldspots_regions[[p]][[c]]
    proximal<-list(); telomeric<-list()
    for (r in 1:nrow(regions)){
      start<-regions[r,1]
      if (start%in%methy_list_windows[[1]]$coldspots_proximal[[p]][[c]][,1]){proximal[[paste(regions[r,1])]]<-regions[r,]}
      if (start%in%methy_list_windows[[1]]$coldspots_telomeric[[p]][[c]][,1]){telomeric[[paste(regions[r,1])]]<-regions[r,]}
    }#r
    coldspots_regions[[p]][[c]]<-list(proximal, telomeric); names(coldspots_regions[[p]][[c]])<-c("proximal", "telomeric")
  }#c
}#p

#-----

window_types<-c("coldspots_proximal", "coldspots_telomeric", "hotspots_proximal", "hotspots_telomeric", "hotspots_pericentromeric")

graphic_list<-list()

windows_around_break_list<-list()

for (w in window_types){ cat(w, fill = TRUE)

graphic_list[[w]]<-list()
  
windows_around_break_list[[w]]<-list()
  
for (m in met_context){  
  
  graphic_list[[w]][[m]]<-list()
  
  #pop_table<-matrix(ncol = 21, nrow = 3); row.names(pop_table)<-Populations
  windows_around_break_list[[w]][[m]]<-list()
  
  for (p in Populations){cat(p);cat(":")
    
    svs_prop_around_break<-matrix(ncol = 11, nrow = 7)
    
    windows_around_break_list[[w]][[m]][[p]]<-list()
      
    for (c in 1:7){ cat(c);cat("-")

      #get windows of CO breakpoint
      break_windows<-as.numeric(methy_list_windows[[m]][[w]][[p]][[c]][,1])
      
      #take out coldspots in regions and add only one per region 
      if (w%in%c("coldspots_proximal", "coldspots_telomeric")){
        break_windows<-break_windows[-which(break_windows%in%coldspots_in_regions_list[[p]][[c]])]
        #add only one window per region (mid one)  i<-1
        for (i in 1:length(coldspots_regions[[p]][[c]][[gsub("coldspots_","",w)]])){break_windows<-c(break_windows, coldspots_regions[[p]][[c]][[gsub("coldspots_","",w)]][[i]][1])}#i
      }#w
      
      #get windows around
      if (any(break_windows<=100000)){break_windows<-break_windows[-which(break_windows<=100000)]}
      
      #get windows around m<-"cpg"
      chr_windows<-as.numeric(row.names(methy_level[[c]][[m]]))-9999
      similar_windows<-methy_list_windows[[m]][[w]][[p]][[c]][,1]
      windows_around_break_table<-matrix(ncol = 11, nrow = 0)
      for (b in 1:length(break_windows)){  
        #first, take out similar windows when counting windows around a particular type  #solo dejo la break
        chr_windows_wo_similar<-chr_windows[-which(chr_windows%in%(similar_windows[-which(similar_windows==break_windows[b])]))]
        windows_around_break<-chr_windows_wo_similar[(which(chr_windows_wo_similar==break_windows[b])-5):(which(chr_windows_wo_similar==break_windows[b])+5)]  
        windows_around_break_table<-rbind(windows_around_break_table, windows_around_break)
      }#b
  
      #convert windows to SVs values
      values<-methy_level[[c]][[m]][,p]; names(values)<-as.character(as.numeric(names(values))-9999)
      window_retrieve_value<-function(x){return(values[as.character(x)])}
      windows_around_break_table<-apply(windows_around_break_table, 1, FUN = window_retrieve_value)
      windows_around_break_table<-t(windows_around_break_table)
      
      #in the case of the breaks in regions, get the region mean
      if (w%in%c("coldspots_proximal", "coldspots_telomeric")){
        for (b in 1:length(break_windows)){
          if (break_windows[b]%in%coldspots_in_regions_list[[p]][[c]]){
            break_region<-coldspots_regions[[p]][[c]][[gsub("coldspots_","",w)]][[which(names(coldspots_regions[[p]][[c]][[gsub("coldspots_","",w)]])==break_windows[b])]]  
            break_region<-seq(break_region[1], break_region[2]+1, by = 10000)
            break_region_mean<-mean(unlist(lapply(break_region, FUN = window_retrieve_value)))
            windows_around_break_table[b,6]<-break_region_mean
          }#if
        }#b
      }#w
      
      windows_around_break_list[[w]][[m]][[p]][[c]]<-windows_around_break_table
      
      #average values
      for (j in 1:ncol(windows_around_break_table)){svs_prop_around_break[c,j]<-mean(windows_around_break_table[,j], na.rm = TRUE)}
      
    }#c  
    
    #average values
    #for (j in 1:ncol(svs_prop_around_break)){pop_table[p,j]<-mean(svs_prop_around_break[,j], na.rm = TRUE)}
    
    graphic_list[[w]][[m]][[p]]<-svs_prop_around_break
    
    cat(fill = TRUE)
    
  }#p
  
  #graphic_list[[w]][[p]]<-pop_table
  
}#m
}#w

saveRDS(windows_around_break_list, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.2.0_windows_around_break_list.RDS")

saveRDS(graphic_list, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.2.1_graphic_list.RDS")
graphic_list<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.2.1_graphic_list.RDS")


# plot(pop_table[1,]~c(1:ncol(pop_table)))
# lines(pop_table[1,]~c(1:ncol(pop_table)))

# #take means of pericentromere
# distal_means<-list()
# for (p in Populations){ 
#   data_distal<-c()
#   for (c in 1:7){
#     data_c<-data_c[which(as.numeric(names(data_c))>as.numeric(pericentromere_pops[[p]][c,4]))]
#     distal_windows_data<-data_c[which(as.numeric(names(data_c))<as.numeric(pericentromere_pops[[p]][c,5]))]
#     data_distal<-c(data_distal, unlist(methy_level[[p]][[c]][names(distal_windows_data)]))
#   }#c
#   distal_means[[p]]<-mean(data_distal, na.rm = TRUE)
# }#p

#I will average chrs per pop
#graphic_table<-matrix(ncol = 5, nrow=0); colnames(graphic_table)<-c("Population", "window_type", "methy_level","chr", "window")
graphic_table<-matrix(ncol = 5, nrow=0); colnames(graphic_table)<-c("Population", "window_type", "methy_level", "window", "met_contexts")
new_met_contexts<-c("CpG", "CHG", "CHH"); names(new_met_contexts)<-c("cpg", "chg", "chh")
for (w in window_types){
  for (m in met_context){
  for (p in Populations){
    #for (c in 1:7){
    for (i in 1:ncol(graphic_list[[w]][[m]][[p]])){
      col2win<-(i-6)*10
      #graphic_table<-rbind(graphic_table,cbind(p, w, graphic_list[[w]][[p]][c,i], c, col2win))
      graphic_table<-rbind(graphic_table,cbind(p, w, mean(graphic_list[[w]][[m]][[p]][,i]), col2win, new_met_contexts[m]))
    }#i
    #}#c
  }#p
  }#m  
}#w

graphic_table<-as.data.frame(graphic_table)
graphic_table$methy_level<-as.numeric(as.character(graphic_table$methy_level))
graphic_table$window<-as.numeric(as.character(graphic_table$window))
#graphic_table$chr<-paste(graphic_table$chr, "H", sep = "")


#color_cold<-c("#000080", "#4682B4")
color_cold<-c("#1b8a5a", "#1d4877")
color_hot<-c("#fbb021", "#f68838", "#ee3e32") 

pdf("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.2.2_windows_around_methy.pdf", width = 15, height = 10)

ggplot()+
  geom_line(data = graphic_table, aes(x = window, y = methy_level, color = window_type)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  
  scale_y_continuous(name = "Methylation level", position = "left") + 
                     #limits = c(0, 0.05),
                     #breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_x_continuous(name = "Physical distance (kb)",
                     limits = c(-40, 40),
                     breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40),
                     labels = c(-40, -30, -20, -10, 0, 10, 20, 30, 40))+
  
  scale_color_manual(values = c(color_cold, color_hot), name = "Window type",
  breaks = c("coldspots_proximal" ,  "coldspots_telomeric" , "hotspots_pericentromeric", "hotspots_proximal", "hotspots_telomeric"),
  labels = c("Coldspots distal proximal" ,  "Coldspots distal telomeric" , "Hotspots pericentromeric", "Hotspots distal proximal", "Hotspots distal telomeric")
  #guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  #facet_grid(cols = vars(Population), rows = vars(chr), scales = "free", space = "free_y") +
  facet_grid2(rows = vars(Population), cols = vars(factor(met_contexts, levels=c('CpG','CHG','CHH'))),
              scales = "free_y", independent = "y") +
  theme(  
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_blank(),
    #strip.text.x = element_blank(),
    axis.ticks = element_line(),
    legend.key = element_blank(),
    #axis.title.y.left = element_blank(),
    #axis.title.x.bottom = element_blank()
  ) #+

dev.off()