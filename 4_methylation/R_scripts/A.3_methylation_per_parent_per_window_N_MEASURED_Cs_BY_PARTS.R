library(data.table)
library(TeachingDemos)
library(dplyr)

#   make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
#   new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
#   new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

get_average<-function(x){return(round(mean(methy_data[which(positions%in%positions[which(positions>(x-s))][which(positions[which(positions>(x-s))]<x)])], na.rm = TRUE), digits = 3))}

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))
#chrs<-paste("chr",1:7,"H", sep = "")

parents.code<-as.matrix(read.table("/scratch/Federico/Paper_2/samplenames.txt"))
parents<-parents.code[,2]

dir.create("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/")

window_sizes<-c(10000, 500000, 1000000)[1]
met_context<-c("cpg", "chg", "chh")[1]
chrs<-1:4; #chrs<-chrs[]
part_size<-1000000

#for (s in window_sizes){ cat(s, fill = TRUE)

s<-window_sizes

ceros<-length(which(unlist(strsplit(as.character(s),""))=="0"))

if (s == 10000){
  rec_list<-readRDS("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.6_hotspots/D.6.0_accumulated_rec_prob.RDS")
  parents<-parents[which(parents%in%c("HOR8160","SprattArcher", "Unumli-Arpa"))]
} else {
  #rec_list<-readRDS(paste("/scratch/Federico/3_RILs/5_correlation/Results/B.1_rec_rates_windows_regions_",s,".RDS", sep = ""))
  rec_list<-readRDS(paste("/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.5_CO_prob_per_window/D.5.1.1_CO_prob_per_window.RDS", sep = ""))
  rec_list<-rec_list[[paste(s)]]
}

for (c in chrs) { cat(c); cat("-")

  #methy_list[[c]]<-list()
  chr.length<-as.numeric(chr_lengths[c,3])  
  windows<-rec_list[[1]][[c]][,2]
  #windows<-windows[1:100]
  for (v in met_context){ cat(v, fill = TRUE)
  
    methy_table<-matrix(ncol = length(parents), nrow = length(windows))
    row.names(methy_table)<-windows
    colnames(methy_table)<-parents
    methy_table[]<-0
    
    methy_table_values<-methy_table
    methy_table_count<-methy_table
    
    for (p in parents){ cat(p); cat(": ")
    
      #get Methy data for that parent and chr            
      P1c<-parents.code[which(parents.code[,2]==p),1]
      Methy_P1<-as.matrix(fread(paste("/scratch/Federico/3_RILs/4_methylation/sources/barley_methylation_data/meth_state/",P1c,"_chr",c,"H_",v,".csv", sep = ""))) 
      Methy_P1[,2]<-round(as.numeric(Methy_P1[,2])/as.numeric(Methy_P1[,3]), digits = 3)
      Methy_P1<-Methy_P1[,1:2]
      positions<-as.numeric(Methy_P1[,1])
      methy_data<-as.numeric(Methy_P1[,2])
      
      #first and last windows
      positions_first_win<-positions[which(positions<s)]
      values_first_win<-methy_data[which(positions<s)]
      positions_last_win<-positions[which(positions>windows[length(windows)-1])]
      values_last_win<-methy_data[which(positions>windows[length(windows)-1])]
      #eliminate first window positions from data
      positions_middle<-positions[-which(positions%in%c(positions_first_win,positions_last_win))] #[which(positions>=10000)]
      values_middle<-methy_data[-which(positions%in%c(positions_first_win,positions_last_win))] #[which(positions>=10000)]
      
      #divide middles positions in parts
      n_parts<-ceiling(length(as.numeric(positions_middle))/part_size); parts<-list(); cat("PARTS:")
      
      for (g in 1:n_parts){ cat(c(paste(g,"/",n_parts, sep = "")), "- ") 
      part_start<-part_size*g-part_size; part_end<-part_size*g-1
      positions<-positions_middle[part_start:part_end]
      values<-values_middle[part_start:part_end]
      
      #delete NA
      if (any(is.na(positions)|is.na(values))){
        to.delete<-unique(which(is.na(positions)|is.na(values)))
        positions<-positions[-to.delete]; values<-values[-to.delete]
      }
      
      #count positions per window
      positions<-strsplit(as.character(positions), "")
      positions<-lapply(positions, FUN = function(x){return(x[-c((length(x)-(ceros-1)):length(x))])})
      positions<-unlist(lapply(positions, FUN = function(x){return(paste(x,collapse = ""))}))
      positions<-as.numeric(paste(positions, paste(rep(9, ceros), collapse = ""), sep = ""))
      positions<-as.data.frame(positions, stringsAsFactors = T)
      positions<-cbind(positions, values); colnames(positions)<-c("pos", "val")
    
      agregated_Cs<-as.matrix(summarise(.data = group_by(.data = positions, pos), val = length(val)))
      accumulated_values<-as.matrix(summarise(.data = group_by(.data = positions, pos), val = sum(val)))
      tabla_part<-cbind(agregated_Cs, accumulated_values[,2])
      row.names(tabla_part)<-tabla_part[,1]
      colnames(tabla_part)<-c("window","Cs","accum_values")

      windows_part<-as.character(as.numeric(format(tabla_part[,1], scientific = FALSE)))
      if (any(((windows_part%in%windows)==FALSE))){windows_part<-windows_part[-which((windows_part%in%windows)==FALSE)]}
      
      #add to tables
      for (w in windows_part){  
        methy_table_count[w,p]<-sum(as.numeric(methy_table_count[w,p]), as.numeric(tabla_part[w,2]), na.rm = T)
        methy_table_values[w,p]<-sum(as.numeric(methy_table_values[w,p]), as.numeric(tabla_part[w,3]), na.rm = T)
      }#w
      
      }#g
      
      #add first window and last window
      if (length(positions_first_win)!=0){
        methy_table_count[1,p]<-length(positions_first_win)
        methy_table_values[1,p]<-sum(values_first_win)
      }
      
      if (length(positions_last_win)!=0){
        methy_table_count[nrow(methy_table),p]<-length(positions_last_win)
        methy_table_values[nrow(methy_table),p]<-sum(values_last_win)
      }
      
      cat(fill = TRUE)
      
    }#p
    
    #create methy table levels
    methy_table_count<-as.data.frame(methy_table_count)
    methy_table_values<-as.data.frame(methy_table_values)
    for (j in 1:ncol(methy_table_count)){
      methy_table_count[,j]<-as.numeric(as.character(methy_table_count[,j]))
      methy_table_values[,j]<-as.numeric(as.character(methy_table_values[,j]))
    }
    methy_table_level<-as.data.frame(methy_table_values)/as.data.frame(methy_table_count)
    
    write.csv(methy_table_count, paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A.3_number_of_Cs_win=",s,"_",v,"_",c,".csv", sep = ""))
    #write.csv(methy_table_values, paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A.3_accum_methy_win=",s,"_",v,"_",c,".csv", sep = ""))
    write.csv(methy_table_level, paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A.3_methylation_level_win=",s,"_",v,"_",c,".csv", sep = ""))
  }#v
}#c
#}#s


#### male list

#for (s in window_sizes){ cat(s, fill = TRUE)

methy_list<-list()

for (c in 1:7) { cat(c); cat("-")
  
  methy_list[[c]]<-list()

  for (v in met_context){ cat(v); cat("-")
    
    methy_list[[c]][[v]]<-list()

    #methy_data<-read.csv(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A_methylation_per_parent_per_win=",s,"_",v,"_",c,".csv", sep = ""))
    methy_data<-read.csv(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A.3_methylation_level_win=",s,"_",v,"_",c,".csv", sep = ""))
    row.names(methy_data)<-methy_data[,1]; methy_data<-methy_data[,-1]
    colnames(methy_data)[which(colnames(methy_data)=="Unumli.Arpa")]<-"Unumli-Arpa"
    colnames(methy_data)[which(colnames(methy_data)=="W23829.803911")]<-"W23829/803911"
    
    Cs_counted<-read.csv(paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/",v,"/A.3_number_of_Cs_win=",s,"_",v,"_",c,".csv", sep = ""))
    row.names(Cs_counted)<-Cs_counted[,1]; Cs_counted<-Cs_counted[,-1]
    colnames(Cs_counted)[which(colnames(Cs_counted)=="Unumli.Arpa")]<-"Unumli-Arpa"
    colnames(Cs_counted)[which(colnames(Cs_counted)=="W23829.803911")]<-"W23829/803911"
    
    methy_list[[c]][[v]][["methy_data"]]<-methy_data
    methy_list[[c]][[v]][["Cs_counted"]]<-Cs_counted
    

  }#v
  
  cat(fill = TRUE)
  
}#c

saveRDS(methy_list, paste("/scratch/Federico/3_RILs/4_methylation/Results/A_methylation_per_parent_per_window/A_methylation_per_parent_level_and_Cs_win=",s,".RDS", sep = ""))

#}#s
