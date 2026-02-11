Populations<-c(paste("HvDRR0",c(2:4,7:9),sep = ""), paste("HvDRR",10:48, sep = ""))

pop_data<-readRDS("/scratch/Federico/Paper_1/A.5_pop_data_curated.RDS")[["2"]]

# 0 ### hay un error en el Paper 1 que es que en el cleaning muchas veces los primeros markers del genetic map se volaron,
#y por eso el mapa empieza en un cM alto..eso genera rec. rates muy grandes para la primer window muchas veces.

#aca hago que tod empiece de 0
#IMPORTANT: el problema que me va a generar es que las funciones que hice para el Paper 1, no las puedo usar mas *

for (p in Populations){
  for (c in 1:7){
    pop_data[[p]][[c]][,4]<-as.numeric(pop_data[[p]][[c]][,4])-as.numeric(pop_data[[p]][[c]][1,4])
  }#c
}#p
#############################

# 1 ### convert V2 to V3 ################

V3<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/iSelect_Marker_V3.csv"))[,c(1,5,8)]

for( p in Populations){ cat(p); cat(": ")
  for (c in 1:7){cat(c); cat(" - ")
#p<-Populations[45];c<-7
    tabla<-pop_data[[p]][[c]][,1:5]

    #add V3 physical value
    V3_chr<-V3[which(V3[,2]==paste("chr",c, "H", sep = "")),]
    tabla<-tabla[which(tabla[,2]%in%V3_chr[,1]),]
    for (m in 1:nrow(tabla)){tabla[m,5]<-V3_chr[which(V3_chr[,1]==tabla[m,2]),3]}
    tabla[,5]<-gsub(",", "", tabla[,5])

    #oder by genetic value
    tabla<-tabla[order(as.numeric(tabla[,4])),]

    #order co-segregating markers
    for (m in 1:nrow(tabla)){
      if (any(tabla[-m,4]==tabla[m,4])){
        to.order<-which(tabla[,4]==tabla[m,4])
        tabla[to.order,]<-tabla[to.order[order(as.numeric(tabla[to.order,5]))],]
      }
    }

    #eliminate not monotonic
    to.discard<-c()
    phys.values<-as.numeric(tabla[,5])
    for (m in 2:nrow(tabla)){if (any(phys.values[m]<phys.values[(m-1):1])){to.discard<-c(to.discard, m)}}
    tabla<-tabla[-to.discard,]

    #correct for negatives
    if (as.numeric(tabla[1,4])<0){
    tabla[2:nrow(tabla),4]<-as.numeric(tabla[2:nrow(tabla),4])+abs(as.numeric(tabla[1,4]))
    tabla[1,4]<-0
    }
    
    #check if both genetic and physical orders are monotonically increasing
    gen_check<-all(as.numeric(tabla[,4]) == cummax(as.numeric(tabla[,4])))
    phys_check<-all(as.numeric(tabla[,5]) == cummax(as.numeric(tabla[,5])))

    if (any(c(gen_check, phys_check))==FALSE){"SOMETHING GOES WRONG"}

    #save
    pop_data[[p]][[c]]<-tabla[,1:5]

  }#c
}#p

saveRDS(pop_data, "/scratch/Federico/3_RILs/5_correlation/Results/A.1_pop_data_V3.RDS")

#############################################

pop_data<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/A.1_pop_data_V3.RDS")

#plot check
library(ggplot2)

#create graphic table
graphic_table<-matrix(ncol = 5, nrow = 0); colnames(graphic_table)<-colnames(pop_data[[1]][[1]])
for (p in Populations){
  tabla<-pop_data[[p]][[1]]
  for (c in 2:7){tabla<-rbind(tabla, pop_data[[p]][[c]])}
  graphic_table<-rbind(graphic_table, tabla)
}

#convert bp to Mbp
colnames(graphic_table)[5]<-"Mbp"
graphic_table[,5]<-as.numeric(as.character(graphic_table[,5]))/1000000

#prepare table
graphic_table<-as.data.frame(graphic_table)
graphic_table$Mbp=as.numeric(levels(graphic_table$Mbp))[graphic_table$Mbp]
#graphic_table$curve=as.numeric(levels(graphic_table$curve))[graphic_table$curve]
graphic_table$cM=as.numeric(levels(graphic_table$cM))[graphic_table$cM]

m.f<-0.1
ancho<-1300*m.f
alto<-130*m.f

pdf("/scratch/Federico/3_RILs/5_correlation/Results/A.2_V3_marey_maps.pdf",width=ancho,height=alto)

ggplot() +
  geom_point(data=graphic_table ,stat = "identity", aes(x=Mbp, y=cM), color = "#663a82") +
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
        strip.text = element_text(size = 150*m.f, margin = margin(50*m.f,50*m.f,50*m.f,50*m.f)))

dev.off()
#############################################


# 2 ####### CUBIC SPLINES ######################


#this parameters belong to Paper1 version, check if needed/work pop by pop, if not adjust manually again
#adjust.parameters<-readRDS("/scratch/Federico/3_RILs/sources/C.2.2_adjust.parameters.RDS")
#*NO PUEDO USAR MAS PARAMETROS DEL PAPER1

cubic_splines<-list()

for(p in Populations){ #print(p)
  cubic_splines[[p]]<-list()
  for(c in 1:7){ #cat(c); cat(" - ")

    physCoord<-c(as.numeric(pop_data[[p]][[c]][,5]))
    genCoord<-c(as.numeric(pop_data[[p]][[c]][,4]))

    modelCubicSpline<-smooth.spline(x = physCoord, y = genCoord, df = 30)

  #  if (is.null(adjust.parameters[[paste(p)]][[paste(c)]])==FALSE){ cat("customized parameters"); cat(" - ")
  #    modelCubicSpline <- adjust.parameters[[paste(p)]][[paste(c)]]
  #  }

    cubic_splines[[p]][[c]]<-modelCubicSpline
  }
  cat(fill = TRUE)
}

#############################################

chr_lengths<-as.matrix(read.csv("/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv"))

window_size<-1000000

new_map_list<-list()
#p<-Populations[24]; c<-5
for (p in Populations){ cat(p); cat(": ")
  new_map_list[[p]]<-list()
  for (c in 1:7){  cat(c); cat(" - ")
    #p<-Populations[24]; c<-5
    tabla<-matrix(nrow = length(seq(1, 700000000, by = window_size)), ncol = 4)
    colnames(tabla)<-colnames(graphic_table)[c(1,3,5,4)]
    tabla[,1]<-p
    tabla[,2]<-c
    tabla[,3]<-seq(1, 700000000, by = window_size)
    model<-cubic_splines[[p]][[c]]
    for (m in 1:nrow(tabla)){tabla[m,4] <- predict(model, as.numeric(tabla[m,3]))$y}

    physical_positions<-as.numeric(tabla[,3])
    new_map<-as.numeric(tabla[,4])

    #### correction by real data #### CORRIJO A TODOS, NO SOLO A LOS NON MONOTONIC
    for (i in length(new_map):1){ #debe ser de atras para adelante
    if (i==length(new_map)){
      new_map[i]<-as.numeric(pop_data[[p]][[c]][nrow(pop_data[[p]][[c]]),4])
    } else {
      if (any(as.numeric(pop_data[[p]][[c]][,5])>physical_positions[i])){
        new_map[i]<-as.numeric(pop_data[[p]][[c]][which(as.numeric(pop_data[[p]][[c]][,5])>physical_positions[i])[1],4])
      } else {new_map[i]<-new_map[i+1]}
    }
    }#for

    ##el mas cercano
    #new_map[i]<-as.numeric(pop_data[[p]][[c]][which(abs(as.numeric(pop_data[[p]][[c]][,5])-physical_positions[i])==min(abs(as.numeric(pop_data[[p]][[c]][,5])-physical_positions[i])))[1],4])

    #arreglo si le hace pegar un salto muy grande #
    big_jump<-which(abs(new_map[-1]-new_map[-length(new_map)])>5)+1  #es clave el +1
    #for (i in big_jump){new_map[i]<-mean(c(new_map[which(new_map[1:i]<new_map[i])], new_map[which(new_map[i:length(new_map)]>new_map[i])]), na.rm = TRUE)}#i
    while (length(big_jump!=0)){
    for (i in big_jump){new_map[i]<-mean(c(new_map[i-1], new_map[i]), na.rm = TRUE)}#i
    big_jump<-which(abs(new_map[-1]-new_map[-length(new_map)])>5)+1
    }

    #plot(c(new_map, pop_data[[p]][[c]][,4])~c(physical_positions, pop_data[[p]][[c]][,5]), col=c(rep("blue", length(new_map)), rep("red", nrow(pop_data[[p]][[c]]))))

    tabla[,4]<-new_map

    #correct for negatives
    if (as.numeric(tabla[1,4])<0){
      tabla[2:nrow(tabla),4]<-as.numeric(tabla[2:nrow(tabla),4])+abs(as.numeric(tabla[1,4]))
      tabla[1,4]<-0
    }  
    
    #corrijo el largo del chr
    tabla<-tabla[-c(which(as.numeric(tabla[,3])>as.numeric(chr_lengths[c,3]))),]

    #plot(c(tabla[,4], pop_data[[p]][[c]][,4])~c(tabla[,3], pop_data[[p]][[c]][,5]), col=c(rep("blue", nrow(tabla)), rep("red", nrow(pop_data[[p]][[c]]))))

    ################################################

    #hago un cubic spline again to smooth the curve, ME INFLA LOS PRIMEROS
    physCoord<-c(as.numeric(tabla[,3]))
    genCoord<-c(as.numeric(tabla[,4]))

    modelCubicSpline<-smooth.spline(x = physCoord, y = genCoord, df = 30)

    for (m in nrow(tabla):1){
      #uso middle window values para smoothear la curva
      if (m!=1){tabla[m,3] <- ((as.numeric(tabla[m,3])-1)+(as.numeric(tabla[m-1,3])-1))/2}
      tabla[m,4] <- predict(modelCubicSpline, as.numeric(tabla[m,3]))$y
      #el cubic spline me genera non-monotonic, so hacer monotonic
      if (m!=nrow(tabla)){if (as.numeric(tabla[m,4])>as.numeric(tabla[m+1,4])){tabla[m,4]<-tabla[m+1,4]}}
    }

    
    #correct for negatives
    if (as.numeric(tabla[1,4])<0){
      tabla[2:nrow(tabla),4]<-as.numeric(tabla[2:nrow(tabla),4])+abs(as.numeric(tabla[1,4]))
      tabla[1,4]<-0
    } 
    ################################################

    #le agrego una extension hasta el fin del chr
    tabla<-rbind(tabla, tabla[nrow(tabla),])
    tabla[nrow(tabla),3]<-chr_lengths[c,3]
    tabla[nrow(tabla),4]<-NA
    
    #si no hay markers en las ultimas windows la pongo como NA, porque el cubic spline ahi tira cualquiera
    last_marker_pos<-as.numeric(pop_data[[p]][[c]][nrow(pop_data[[p]][[c]]),5])
    if (any((as.numeric(tabla[,3])-(window_size/2))>last_marker_pos)){
    tabla[which((as.numeric(tabla[,3])-(window_size/2))>last_marker_pos),4]<-NA
    }
    
    #y para que no queden NA, le pongo un valor generado por el creciemiento promedio de las 5 windows anteriores
    #ESTO ES POSIBLE QUE ME SOBESTIME PERO CASI NINGUNA POP TIENE VALUES EN LAS NA WINDOWS COMO PARA HACER ALGO MEJOR
    if (any(is.na(tabla[,4]))){
    na_win<-which(is.na(tabla[,4]))
    for (i in na_win){
    #Antes hacia un promedio, pero ahora hago un lm
    # growth<-as.numeric(tabla[c((i-5):(i-1)),4])
    # growth<-mean(growth[2:5]-growth[1:4])
    # tabla[i,4]<-as.numeric(tabla[i-1,4])+growth
    last5_gen<-as.numeric(tabla[c((i-5):(i-1)),4])
    last5_pys<-as.numeric(tabla[c((i-5):(i-1)),3])
    LM<-lm(last5_gen~last5_pys)
    tabla[i,4]<-LM$coefficients[1]+LM$coefficients[2]*as.numeric(tabla[i,3])
    }#for
    }#if

    #######################################################

    #lo llevo a 0, me va a generar desfasaje en los viejos puntos y la nueva curva, 
    #pero los viejos puntos estan atados a un inicio ficticio tb  
    tabla[,4]<-as.numeric(tabla[,4])-as.numeric(tabla[1,4])
    
    #######################################################
    
    plot(c(tabla[,4], pop_data[[p]][[c]][,4])~c(tabla[,3], pop_data[[p]][[c]][,5]), col=c(rep("blue", nrow(tabla)), rep("red", nrow(pop_data[[p]][[c]]))))

    #save
    new_map_list[[p]][[c]]<-tabla

  }#c
  cat(fill = TRUE)
}#p

#para las ultimas windows que no tienen markers
##########################################

saveRDS(new_map_list, "/scratch/Federico/3_RILs/5_correlation/Results/A.3_V3_cubic_splines.RDS")

#new_map_list<-readRDS("/scratch/Federico/3_RILs/5_correlation/Results/A.3_V3_cubic_splines.RDS")

#prepare table
graphic_table2<-matrix(ncol = 4, nrow = 0); colnames(graphic_table2)<-c("Population", "chr", "Mbp", "cM")
for (p in Populations){for (c in 1:7){graphic_table2<-rbind(graphic_table2, new_map_list[[p]][[c]])}}
graphic_table2<-as.data.frame(graphic_table2)
graphic_table2$Mbp=as.numeric(levels(graphic_table2$Mbp))[graphic_table2$Mbp]
graphic_table2$cM=as.numeric(levels(graphic_table2$cM))[graphic_table2$cM]
graphic_table2$Mbp<-graphic_table2$Mbp/1000000

# m.f<-0.1
ancho<-1300*m.f
alto<-130*m.f

pdf("/scratch/Federico/3_RILs/5_correlation/Results/A.4_V3_cubic_splines.pdf",width=ancho,height=alto)

ggplot() +
  geom_line(data=graphic_table2 ,stat = "identity", aes(x=Mbp, y=cM), color = "#663a82") +
  geom_point(data=graphic_table ,stat = "identity", aes(x=Mbp, y=cM), size = 0.1) +
  scale_x_continuous(name = "physical distance (Mbp)", position = "bottom") +
  scale_y_continuous(name = "genetic distance (cM)", position = "right") +
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
        strip.text = element_text(size = 150*m.f, margin = margin(50*m.f,50*m.f,50*m.f,50*m.f)))

dev.off()

##############################################

# 3 calculate rec rate

tabla<-matrix(nrow = 45, ncol = 8)
colnames(tabla)<-c(paste(1:7, "H", sep = ""), "Genome")
row.names(tabla)<-Populations
for (p in Populations){
for (c in 1:7){
#chr
map_length<-as.numeric(new_map_list[[p]][[c]][nrow(new_map_list[[p]][[c]]),4])  
physical_length<-as.numeric(new_map_list[[p]][[c]][nrow(new_map_list[[p]][[c]]),3])
tabla[p,c]<-map_length/(physical_length/1000000) #rec rate in cM/Mbp
}
#genome
tabla[p,8]<-mean(as.numeric(tabla[p,1:7]))
}  
write.csv(round(tabla, digits = 3), "/scratch/Federico/3_RILs/5_correlation/Results/A.5.1_rec_rates_V3_genome_and_chrs.csv")

#####



#windows

#first check which window size make sense to catch changes
tabla2<-matrix(nrow = 45, ncol = 8)
colnames(tabla2)<-c(paste(1:7, "H", sep = ""), "Genome"); row.names(tabla2)<-Populations
for (p in Populations){ 
  for (c in 1:7){ 
    physical_positions<-sort(as.numeric(pop_data[[p]][[c]][,5]))
    int_mrk_dis<-physical_positions[2:length(physical_positions)]-physical_positions[1:(length(physical_positions)-1)]
    tabla2[p,c]<-mean(int_mrk_dis)
  }
  tabla2[p,8]<-mean(tabla2[p,1:7])
}
tabla2<-round(tabla2)
mean(tabla2[,8])

rec_list<-list()
window_sizes<-c(500000, 1000000)[2]

for (s in window_sizes){ cat(s, fill = TRUE)

for (p in Populations){ cat(c(p,":"), sep = "")

rec_list[[p]]<-list()

for (c in 1:7){ cat(c,"-", sep = "")

#la ultima window va a contener al end del chr
windows<-1:ceiling(as.numeric(new_map_list[[p]][[c]][nrow(new_map_list[[p]][[c]]),3])/s)
#uso total value of window
windows<-(windows*s)
names(windows)<-windows

#predigo los window values based on, ESTO POR ALGUNA RAZON ME INFLA LOS VALORES AL PRINCIPIO
physCoord<-as.numeric(new_map_list[[p]][[c]][,3])
genCoord<-as.numeric(new_map_list[[p]][[c]][,4])
# modelCubicSpline<-smooth.spline(x = physCoord, y = genCoord, df = 30) 
#
#windows<-round(predict(modelCubicSpline, windows)$y, digits = 2)

genCoord<-genCoord[-1]
physCoord<-physCoord[-1]
physCoord<-ceiling(physCoord)

#calculo el punto medio entre los puntos physicos  #ESTO ASI SOLO FUNCIONA PARA S=1Mbp
windows[]<-NA
for (i in length(genCoord):2) {    
windows[which(as.numeric(names(windows))==((physCoord[i-1]+physCoord[i])/2))]<-genCoord[i-1]+((genCoord[i]-genCoord[i-1])/2)
}

#para las ultimas dos hago un proyectado con la ultima phys position
dif.gen<-genCoord[length(genCoord)]-genCoord[length(genCoord)-1]
dif.phys<-physCoord[length(physCoord)]-physCoord[length(physCoord)-1]
growth_ultima<-dif.gen*s/dif.phys
#para la anteultima hago miti miti
growth_anteultima<-windows[length(windows)-2]-windows[length(windows)-3]
windows[length(windows)-1]<-windows[length(windows)-2]+(growth_anteultima/2+growth_ultima/2)
windows[length(windows)]<-windows[length(windows)-1]+growth_ultima

#transformo los negativos
if (any(windows<0)){
  negatives<-which(windows<0)
  first_positive<-negatives[length(negatives)]+1
  rec_rate<-windows[first_positive]/as.numeric(names(windows[first_positive])) #cM/bp
  for (i in negatives){windows[i]<-as.numeric(names(windows[i]))*rec_rate}
}

#check if monotonic
#which((windows[2:length(windows)]-windows[1:(length(windows)-1)])<0)

#el cubic spline me genera non-monotonic, so hacer monotonic
for (w in (length(windows)-1):1){if(windows[w]>windows[w+1]){windows[w]<-windows[w+1]}}

#calculate rec rate:

#1-calculate genetic map diferences (lo que crecio el map en cada window)
windows[1:length(windows)]<-windows[1:length(windows)]-c(0,windows[1:(length(windows)-1)])
#2-convert to cM/Mbp  (lo que crecio el map en cada window por bp) * 1 Mbp = cM/Mbp
windows[1:length(windows)]<-(windows[1:length(windows)]/s)*1000000

# #hay un problema, la primera gen position ya viene alta de la data del Paper1, #ESTO LO CORREGI MAS ARRIBA
# #so voy a calcular la rec.rate, para la first window, basado en, por ejemplo, si la s = 1Mbp, 
# #antes calculaba la rr basado en las gen distances 0 - 1Mbp, ahora lo hago en la dif 500 Kb - 1Mbp. 
# 
# #1 - primero las otras como antes
# windows[2:length(windows)]<-windows[2:length(windows)]-c(windows[1:(length(windows)-1)])
# #2-convert to cM/Mbp  (lo que crecio el map en cada window por bp) * 1 Mbp = cM/Mbp
# windows[2:length(windows)]<-(windows[2:length(windows)]/s)*1000000
# 
# #3 - luego calculo la first window
# different_s<-(as.numeric(names(windows)[1])-as.numeric(new_map_list[[p]][[c]][2,3]))
# gen_dif<-as.numeric(windows[1])-as.numeric(new_map_list[[p]][[c]][2,4])
# windows[1]<-(gen_dif/different_s)*1000000

#save
rec_list[[p]][[c]]<-windows
}#c 

cat(fill = TRUE)
}#p

saveRDS(rec_list, paste("/scratch/Federico/3_RILs/5_correlation/Results/A.5.2_rec_rates_V3_windows_",s,".RDS", sep = ""))
}#s

##############################################