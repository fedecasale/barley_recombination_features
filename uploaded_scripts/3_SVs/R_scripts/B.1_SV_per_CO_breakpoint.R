setwd("/scratch/Federico/3_RILs/")

sv_list<-readRDS("/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_divided/A.2_SV_per_size.RDS")
SVs<-names(sv_list)
sv_sizes<-names(sv_list[[1]][[1]][[1]])

breakpoint_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS")
layers<-names(breakpoint_list)
Populations<-names(breakpoint_list[[1]])

chrs<-paste("chr", 1:7, "H", sep = "")

dir.create("/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO")

CO_SV_list<-list()
for (v in SVs){ cat(v); cat(": "); 
  CO_SV_list[[v]]<-list()
for (l in layers){ cat(l); cat(" - "); 
  CO_SV_list[[v]][[l]]<-list()
  for (p in Populations){ #cat(p);cat("-"); 
    CO_SV_list[[v]][[l]][[p]]<-list()
    for (c in 1:7){
    CO_SV_list[[v]][[l]][[p]][[c]]<-list()
    for (s in sv_sizes){ 
    if (nrow(sv_list[[v]][[p]][[c]][[s]])!=0){
    
    SVs_p_c<-sv_list[[v]][[p]][[c]][[s]]  
    row.names(SVs_p_c)<-paste(v,"_",c,"_",1:nrow(SVs_p_c), sep = "")
    breakpoints<-breakpoint_list[[l]][[p]][[c]]  
    #delete non breakpoint lines NO HAY MAS
    if (any(breakpoints[,6]=="")){breakpoints<-breakpoints[-which(breakpoints[,6]==""),]}
    #add SV data
    breakpoints<-cbind(breakpoints, NA, NA, NA, NA, NA) 
    colnames(breakpoints)[7:11]<-c("SV_code","start","end","closest","distance")
      for (b in 1:nrow(breakpoints)){ #cat(b); cat(" - ")
      #I take ths sv which start or ending point is the closest to the CO breakpoint
      dif<-abs(as.numeric(breakpoints[b,6])-as.numeric(SVs_p_c[,2]))  
      closest_start<-which(dif==min(dif)) 
      dif<-abs(as.numeric(breakpoints[b,6])-as.numeric(SVs_p_c[,3]))  
      closest_end<-which(dif==min(dif)) 
      closest_SV<-c(SVs_p_c[closest_start,2],SVs_p_c[closest_end,3])
      names(closest_SV)<-row.names(SVs_p_c)[c(closest_start,closest_end)]
      dif<-abs(as.numeric(breakpoints[b,6])-as.numeric(closest_SV))  
      closest_SV<-closest_SV[which(dif==min(dif))[1]]
      #add info to CO line
      breakpoints[b,7]<-names(closest_SV)
      breakpoints[b,8:9]<-SVs_p_c[names(closest_SV),2:3]
      breakpoints[b,10]<-closest_SV
      breakpoints[b,11]<-min(dif)
      }#b
    CO_SV_list[[v]][[l]][[p]][[c]][[s]]<-breakpoints
    }#if 
    }#s
      
    }#c
  }#p
cat(" ")  
}#l
cat(fill = TRUE)
}#v

saveRDS(CO_SV_list, "/scratch/Federico/3_RILs/3_SVs/Results/B.1_SV_per_CO/B.1.2_SV_per_CO_per_SV_size.RDS")
