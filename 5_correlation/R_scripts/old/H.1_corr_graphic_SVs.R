library(ggplot2)
library(impressionist.colors)

SVs_corr<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/G.1_corr_rr_SVs.RDS")

COLNAMES<-c("chr", "physical_window", "corr_variable", "corr_p.val", "corr_direction")

p_val<-0.05

graphic_table_list<-list()

for (i in 1:length(SVs_corr)){
  
  #variable<-substr(names(SVs_corr[i]), 3, length(unlist(strsplit(names(SVs_corr[i]),""))))
  variable<-names(SVs_corr)[i]
  
  if (length(SVs_corr[[i]])==4){ #1
    
    for (v in names(SVs_corr[[i]])){
      variable_v<-paste(variable,": ",v, sep = "")
      graphic_table<-matrix(ncol = 5, nrow = 0); colnames(graphic_table)<-COLNAMES
      for (c in 1:7){  
        #make chr graphic table  
        chr_values<-SVs_corr[[i]][[v]][,c]
        if (any(chr_values=="-")){chr_values<-chr_values[-which(chr_values=="-")]}
        chr_table<-matrix(ncol = 5, nrow = length(chr_values)); colnames(chr_table)<-COLNAMES
        chr_table[,1]<-paste(c,"H",sep = "")
        chr_table[,2]<-names(chr_values)
        chr_table[,3]<-variable_v
        chr_table[,4]<-chr_values
        #clean non_significant
        chr_table[which(as.numeric(chr_values)>p_val),4]<-NA  
        graphic_table<-rbind(graphic_table, chr_table)
      }#c
      graphic_table_list[[variable_v]]<-graphic_table    
    }#v
    
  } else { #if1
    
    graphic_table<-matrix(ncol = 5, nrow = 0); colnames(graphic_table)<-COLNAMES
    for (c in 1:7){  
      #make chr graphic table  
      chr_values<-SVs_corr[[i]][,c]
      if (any(chr_values=="-")){chr_values<-chr_values[-which(chr_values=="-")]}
      chr_table<-matrix(ncol = 5, nrow = length(chr_values)); colnames(chr_table)<-COLNAMES
      chr_table[,1]<-paste(c,"H",sep = "")
      chr_table[,2]<-names(chr_values)
      chr_table[,3]<-variable
      chr_table[,4]<-chr_values
      #clean non_significant
      chr_table[which(as.numeric(chr_values)>p_val),4]<-NA
      graphic_table<-rbind(graphic_table, chr_table)
    }#c
    graphic_table_list[[variable]]<-graphic_table    
    
  }#if1 else
  
}#i

#########################################################

#graphic

centromere <-read.csv("/scratch/Federico/3_RILs/5_correlation/Results/B.2_new_chr_length_V3_byFede.csv")
color_table<-as.matrix(get.color(8,3,c(5,7,12)))
centro_col<-get.color(8,3,14); peri_col<-get.color(8,3,13)

#graphic_table_list<-graphic_table_list[1:6]


graphic_table<-matrix(ncol = 5, nrow = 0); colnames(graphic_table)<-COLNAMES
for (i in 1:length(graphic_table_list)){graphic_table<-rbind(graphic_table, graphic_table_list[[i]])}
graphic_table[,3]<-gsub("deletions", "del", graphic_table[,3])
graphic_table[,3]<-gsub("insertions", "ins", graphic_table[,3])
graphic_table[,3]<-gsub("inversions", "inv", graphic_table[,3])
graphic_table[,3]<-gsub("duplications", "dup", graphic_table[,3])
graphic_table<-as.data.frame(graphic_table)
graphic_table$physical_window<-as.numeric(as.character(graphic_table$physical_window))
graphic_table$corr_p.val<-as.numeric(as.character(graphic_table$corr_p.val))

graphic_table[which(graphic_table[,4]=="NaN"),4]<-NA

pdf("/scratch/Federico/3_RILs/5_correlation/Results/H.1.1_corr_graphic_SVs_all.pdf", width = 15, height = 10)

for (c in 1:7){
  
  graphic_table_chr<-graphic_table[which(graphic_table$chr==paste(c,"H",sep = "")),]
  
  plot<-
    ggplot()+
    geom_point(data = graphic_table_chr, aes(x = physical_window, y = corr_variable, color = corr_p.val)) +
    scale_color_gradient(high = "orange", low = "purple", na.value = NA, name = "P value") +
    scale_x_continuous(name = "Physical distance (Mbp)", position = "bottom", expand = c(0.01,0.01)) +
    
    geom_vline(xintercept = centromere$centromere[c], color = centro_col, size = 0.75)+
    geom_vline(xintercept = centromere$peri_start[c], linetype = 2, color = peri_col, size = 0.5)+
    geom_vline(xintercept = centromere$peri_end[c], linetype = 2, color = peri_col, size = 0.5)+
    
    #facet_grid(cols = vars(Population), rows = vars(chr), scales = "free", space = "free", switch = "y") +
    
    theme(
      axis.title.y = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major.y = element_line(size = 0.25, linetype = 'solid', colour = "lightgrey"), 
      panel.grid.major.x = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
      panel.grid.minor.x = element_line(size = 0.25, linetype = 'solid', colour = "grey"),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = 'right',
      legend.text = element_text(size = 7)
    ) +
    
    guides(color = guide_colorbar(title = "P value",
                                  label.position = "right", label.vjust = c(0,0.5,0.5,0.5,0.5,1),
                                  title.position = "top", title.vjust = 1.5,
                                  frame.colour = "black",
                                  barwidth = 1.5,
                                  barheight = 5))
  
  plot<-plot + ggtitle(paste("Chromosome ",c, "H", sep = ""))
  
  print(plot)
  
}#c

dev.off()


####per met context

#saco las variables del promedio
graphic_table<-graphic_table[-which(graphic_table$corr_variable%in%names(SVs_corr)[1:2]),]

list.v<-strsplit(as.character(graphic_table$corr_variable), "")
list.v<-lapply(list.v, function(x){return(paste(x[(length(x)-2):length(x)], collapse = ""))})
list.v<-unlist(list.v)

SVs_types<-unique(list.v)

pdf(paste("/scratch/Federico/3_RILs/5_correlation/Results/H.1.2_corr_graphic_SVs_per_type.pdf", sep = ""), width = 15, height = 10)

for (v in SVs_types){
  
  graphic_table_v<-graphic_table[which(list.v==v),]
  
  plot<-
    ggplot()+
    geom_point(data = graphic_table_v, aes(x = physical_window, y = corr_variable, color = corr_p.val)) +
    scale_color_gradient(high = "orange", low = "purple", na.value = NA, name = "P value") +
    scale_x_continuous(name = "Physical distance (Mbp)", position = "bottom", expand = c(0.01,0.01)) +
    
    geom_vline(xintercept = centromere$centromere[c], color = centro_col, size = 0.75)+
    geom_vline(xintercept = centromere$peri_start[c], linetype = 2, color = peri_col, size = 0.5)+
    geom_vline(xintercept = centromere$peri_end[c], linetype = 2, color = peri_col, size = 0.5)+
    
    facet_grid(rows = vars(chr), scales = "free", space = "free") +
    
    theme(
      axis.title.y = element_blank(),
      strip.background = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major.y = element_line(size = 0.25, linetype = 'solid', colour = "lightgrey"), 
      panel.grid.major.x = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
      panel.grid.minor.x = element_line(size = 0.25, linetype = 'solid', colour = "grey"),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = 'right',
      legend.text = element_text(size = 7)
    ) +
    
    guides(color = guide_colorbar(title = "P value",
                                  label.position = "right", label.vjust = c(0,0.5,0.5,0.5,0.5,1),
                                  title.position = "top", title.vjust = 1.5,
                                  frame.colour = "black",
                                  barwidth = 1.5,
                                  barheight = 5)) 
  
  plot<-plot + ggtitle(paste("SV type: ",v, sep = "")) 
  
  print(plot)
  
}#v
dev.off()
