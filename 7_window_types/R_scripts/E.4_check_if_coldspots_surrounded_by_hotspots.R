library(ggplot2)
library(ggh4x)

Populations<-c("HvDRR13","HvDRR27", "HvDRR28")

svs_list_windows<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.2.2_methylation_prop_per_window_type.RDS")
SVs_proportion<-readRDS(paste("/home/fcasale/Desktop/Paper_2/3_RILs/4_methylation/Results/B_methylation_per_population/B.1_methylation_per_population_win=10000.RDS", sep = ""))

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
      if (start%in%svs_list_windows[[1]]$coldspots_proximal[[p]][[c]][,1]){proximal[[paste(regions[r,1])]]<-regions[r,]}
      if (start%in%svs_list_windows[[1]]$coldspots_telomeric[[p]][[c]][,1]){telomeric[[paste(regions[r,1])]]<-regions[r,]}
    }#r
    coldspots_regions[[p]][[c]]<-list(proximal, telomeric); names(coldspots_regions[[p]][[c]])<-c("proximal", "telomeric")
  }#c
}#p

#window types list
window_types<-c("coldspots_proximal", "coldspots_telomeric", "hotspots_proximal", "hotspots_telomeric", "hotspots_pericentromeric")

types_list<-list()
  for (p in Populations){
    types_list[[p]]<-list()
    for (c in 1:7){
    windows<-c()
    for (w in window_types){ #cat(w, fill = TRUE)
    type_windows<-svs_list_windows[[1]][[w]][[p]][[c]][,1]
    names(type_windows)<-paste(w,1:length(type_windows), sep = "_")
    windows<-c(windows, type_windows)
    }#w
    types_list[[p]][[c]]<-windows
    }#c
  }#p

####

graphic_list<-list()

for (w in window_types){ cat(w, fill = TRUE)
#w<-window_types[3]
graphic_list[[w]]<-list()
  
  for (p in Populations){cat(p);cat(":")
    
    windows_around_break_table_genome<-matrix(ncol = 5, nrow = 0)
    
    for (c in 1:7){ cat(c);cat("-")
      
      #get windows of CO breakpoint
      break_windows<-as.numeric(svs_list_windows[[1]][[w]][[p]][[c]][,1])
      
      #take out coldspots in regions and add only one per region 
      if (w%in%c("coldspots_proximal", "coldspots_telomeric")){
        break_windows<-break_windows[-which(break_windows%in%coldspots_in_regions_list[[p]][[c]])]
        #add only one window per region (mid one)  i<-1
        for (i in 1:length(coldspots_regions[[p]][[c]][[gsub("coldspots_","",w)]])){break_windows<-c(break_windows, coldspots_regions[[p]][[c]][[gsub("coldspots_","",w)]][[i]][1])}#i
      }#w
      
      #get windows around
      if (any(break_windows<=50000)){break_windows<-break_windows[-which(break_windows<=50000)]}
      chr_windows<-as.numeric(row.names(SVs_proportion[[c]][[1]]))-9999
      similar_windows<-svs_list_windows[[1]][[w]][[p]][[c]][,1]
      windows_around_break_table<-matrix(ncol = 5, nrow = 0)
      for (b in 1:length(break_windows)){  
        #first, take out similar windows when counting windows around a particular type  #solo dejo la break
        chr_windows_wo_similar<-chr_windows[-which(chr_windows%in%(similar_windows[-which(similar_windows==break_windows[b])]))]
        windows_around_break<-chr_windows_wo_similar[(which(chr_windows_wo_similar==break_windows[b])-2):(which(chr_windows_wo_similar==break_windows[b])+2)]  
        windows_around_break_table<-rbind(windows_around_break_table, windows_around_break)
      }#b

      #check if surrounding windows are in other window type
      for (b in 1:length(windows_around_break_table)){ 
      x<-windows_around_break_table[b] 
      if(x%in%types_list[[p]][[c]]){
      type<-names(types_list[[p]][[c]])[which(types_list[[p]][[c]]==x)]
      type<-strsplit(type, "_")[[1]][1]
      x<-type
      }
      windows_around_break_table[b]<-x
      }#x
      
      windows_around_break_table_genome<-rbind(windows_around_break_table_genome, windows_around_break_table)
      
    }#c  
    
    graphic_list[[w]][[p]]<-windows_around_break_table_genome
    
    cat(fill = TRUE)
    
  }#p
  
}#w

saveRDS(graphic_list, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.4.1_window_types_around_window_types.RDS")


'%!in%' <- function(x,y)!('%in%'(x,y))
full_table<-matrix(nrow = 0, ncol = 4)
for (p in Populations){
  pop_table<-matrix(nrow = 4, ncol = 2)
  colnames(pop_table)<-c("hotspots", "coldspots")
  row.names(pop_table)<-window_types[1:4]
  for (w in window_types[1:4]){
  types_around<-graphic_list[[w]][[p]][,c(1:2,4:5)]
  hotspots<-length(which(types_around=="hotspots")) 
  coldspots<-length(which(types_around=="coldspots"))
  other<-length(which(types_around%!in%c("hotspots","coldspots")))
  all<-length(types_around)
  pop_table[w,"hotspots"]<-hotspots/all
  pop_table[w,"coldspots"]<-coldspots/all
  }#w
full_table<-rbind(full_table, cbind(c(p,"","",""),row.names(pop_table), pop_table))  
}#p

write.csv(full_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/E_windows_around/E.4.2_window_types_around_window_types.csv")

