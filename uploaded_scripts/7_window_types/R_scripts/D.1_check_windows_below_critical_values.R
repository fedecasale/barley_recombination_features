#http://www.ltcconline.net/greenl/courses/201/hyptest/hypmean.htm

# #this is the t critical calculation for unpaired data, unequal sample size, equal variances
# num<-mean(x1)-mean(x2)
# pooled_stdv_num<-(length(x1)-1)*var(x1)+(length(x2)-1)*var(x2)
# pooled_stdv_den<-length(x1)+length(x2)-2
# pooled_stdv<-sqrt(pooled_stdv_num/pooled_stdv_den)
# # #pooled_stdv<-sd_pooled(x1, x2) #me da igual que mi calculo
# den<-pooled_stdv*sqrt(1/length(x1)+1/length(x2))
# t_critical=num/den
#I will dissipate the mean of x1 -> t_critical*den=num 
#t_critical*den+mean(x2)=mean(x1)

library(xlsx)

#load data
methy_list_windows<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/B_methylation/B.2.4_methylation_prop_per_window_type_unified.RDS")
svs_list_windows<-readRDS("/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/C_SVs/C.2.2_svs_prop_per_window_type_unified.RDS")

met_context<-c("cpg", "chg")
Populations<-names(methy_list_windows[[1]][[1]])

#COLDSPOTS

#TEST #1
### coldspots_distal_telomeric VS distal_proximal ###
#colds tienen menos methy que la distal proximal, pero mas SVs
#agarro los below del critical de SVs, y veo si no tienen mas methy
critical<-c()

results_table<-matrix(ncol = 1, nrow = 3)
#colnames(results_table)<-c("Below t-critical", "Above t-critical", "P value")
row.names(results_table)<-Populations
for (m in met_context){
results_table_m<-matrix(ncol = 3, nrow = 3)
colnames(results_table_m)<-c("Below t-critical", "Above t-critical", "P value")
row.names(results_table_m)<-Populations
  for (p in Populations){
    
methy_coldspots_distal_telomeric<-methy_list_windows[[m]]$coldspots_telomeric[[p]]
methy_distal_proximal<-methy_list_windows[[m]]$distal_proximal[[p]]
sv_coldspots_distal_telomeric<-svs_list_windows$coldspots_telomeric[[p]]

x1<-methy_coldspots_distal_telomeric[which(complete.cases(methy_coldspots_distal_telomeric))]
x2<-methy_distal_proximal[which(complete.cases(methy_distal_proximal))]
x1_svs<-sv_coldspots_distal_telomeric[which(complete.cases(sv_coldspots_distal_telomeric))]

# #run wilcox-test just to check
# wilcox.test(x = x1, y = x2, alternative = "greater", paired = FALSE)
# t.test(x = x1, y = x2, alternative = "greater", paired = FALSE)

#I will get the x1 value needed to be under the t critical
pooled_stdv_num<-(length(x1)-1)*var(x1)+(length(x2)-1)*var(x2)
pooled_stdv_den<-length(x1)+length(x2)-2
pooled_stdv<-sqrt(pooled_stdv_num/pooled_stdv_den)
den<-pooled_stdv*sqrt(1/length(x1)+1/length(x2))

#I will dissipate the mean of x1 -> t_critical*den=num
#choose t_critical -> for infinite values, one-tailed test at 0.05 = 	1.645
MEAN_X1_CRITICAL<-1.645*den+mean(x2)
critical<-c(critical, MEAN_X1_CRITICAL)

#I look for the values below the mean_critical
x1_below<-x1[which(x1<MEAN_X1_CRITICAL)]
#check they have increased SVs than the others in the group
x1_svs_below<-x1_svs[which(names(x1_svs)%in%names(x1_below))]
x1_svs_above<-x1_svs[-which(names(x1_svs)%in%names(x1_below))]
#t.test(x = x1_svs_below, y = x1_svs_above, alternative = "greater")
test<-wilcox.test(x = x1_svs_below, y = x1_svs_above, alternative = "greater")

#add values to table
results_table_m[p,1]<-round(mean(x1_svs_below), digits = 3)
results_table_m[p,2]<-round(mean(x1_svs_above), digits = 3)
results_table_m[p,3]<-round(test$p.value, digits = 4)

}#p
results_table<-cbind(results_table, results_table_m)
}#m
results_table<-rbind(NA, results_table)
results_table[1,]<-colnames(results_table)
results_table[,1]<-row.names(results_table)
colnames(results_table)<-c("Population", "CpG","","", "CHG","","")
write.xlsx(results_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/D.1_check_B_and_C_inconsistencies.xlsx", append = TRUE, sheetName = "TEST #1")


#TEST #2
###coldspots_distal_proximal VS distal_telomeric
#colds tienen menos SVs que la distal telomeric, pero mas methy
#agarro los below del critical de methy, y veo si no tienen mas SVs

results_table<-matrix(ncol = 1, nrow = 3)
#colnames(results_table)<-c("Below t-critical", "Above t-critical", "P value")
row.names(results_table)<-Populations
critical<-c()
for (m in met_context){
  results_table_m<-matrix(ncol = 3, nrow = 3)
  colnames(results_table_m)<-c("Below t-critical", "Above t-critical", "P value")
  row.names(results_table_m)<-Populations
  for (p in Populations){
    
    sv_coldspots_distal_proximal<-svs_list_windows$coldspots_proximal[[p]]
    sv_distal_telomeric<-svs_list_windows$distal_telomeric[[p]]
    methy_coldspots_distal_proximal<-methy_list_windows[[m]]$coldspots_proximal[[p]]

    
    x1<-sv_coldspots_distal_proximal[which(complete.cases(sv_coldspots_distal_proximal))]
    x2<-sv_distal_telomeric[which(complete.cases(sv_distal_telomeric))]
    x1_methy<-methy_coldspots_distal_proximal[which(complete.cases(methy_coldspots_distal_proximal))]
    
    # #run wilcox-test just to check
    # wilcox.test(x = x1, y = x2, alternative = "greater", paired = FALSE)
    # t.test(x = x1, y = x2, alternative = "greater", paired = FALSE)
    
    #I will get the x1 value needed to be under the t critical
    pooled_stdv_num<-(length(x1)-1)*var(x1)+(length(x2)-1)*var(x2)
    pooled_stdv_den<-length(x1)+length(x2)-2
    pooled_stdv<-sqrt(pooled_stdv_num/pooled_stdv_den)
    den<-pooled_stdv*sqrt(1/length(x1)+1/length(x2))
    
    #I will dissipate the mean of x1 -> t_critical*den=num
    #choose t_critical -> for infinite values, one-tailed test at 0.05 = 	1.645
    MEAN_X1_CRITICAL<-1.645*den+mean(x2)
    critical<-c(critical, MEAN_X1_CRITICAL)
    #I look for the values below the mean_critical
    x1_below<-x1[which(x1<MEAN_X1_CRITICAL)]
    #check they have increased SVs than the others in the group
    x1_methy_below<-x1_methy[which(names(x1_methy)%in%names(x1_below))]
    x1_methy_above<-x1_methy[-which(names(x1_methy)%in%names(x1_below))]
    #t.test(x = x1_methy_below, y = x1_methy_above, alternative = "greater")
    test<-wilcox.test(x = x1_methy_below, y = x1_methy_above, alternative = "greater")
    
    #add values to table
    results_table_m[p,1]<-round(mean(x1_methy_below), digits = 3)
    results_table_m[p,2]<-round(mean(x1_methy_above), digits = 3)
    results_table_m[p,3]<-round(test$p.value, digits = 4)
    
  }#p
  results_table<-cbind(results_table, results_table_m)
}#m
results_table<-rbind(NA, results_table)
results_table[1,]<-colnames(results_table)
results_table[,1]<-row.names(results_table)
colnames(results_table)<-c("Population", "CpG","","", "CHG","","")
write.xlsx(results_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/D.1_check_B_and_C_inconsistencies.xlsx", append = TRUE, sheetName = "TEST #2")

#SI NO ME DA VEODE SACARLO LOS COLDSPOTS A LA REGION -> NO FUE NECESARIO
#distal_proximal<-distal_proximal[-which(names(distal_proximal)%in%names(coldspots_distal_proximal))]


#HOTSPOTS

#TEST #3
#hotspots in pericentromeric me dieron que igual methylation que 
#distal telomeric e incluso que coldspots telomeric....
#pruebo si los hotspots mayor critical tienen menos SVs que los above critical

results_table<-matrix(ncol = 1, nrow = 3)
#colnames(results_table)<-c("Below t-critical", "Above t-critical", "P value")
row.names(results_table)<-Populations
for (m in met_context){
  results_table_m<-matrix(ncol = 3, nrow = 3)
  colnames(results_table_m)<-c("Below t-critical", "Above t-critical", "P value")
  row.names(results_table_m)<-Populations
  for (p in Populations){
    
    methy_hotspots_pericentromeric<-methy_list_windows[[m]]$hotspots_pericentromeric[[p]]
    methy_distal_telomeric<-methy_list_windows[[m]]$distal_telomeric[[p]]
    svs_hotspots_pericentromeric<-svs_list_windows$hotspots_pericentromeric[[p]]
    
    x1<-methy_hotspots_pericentromeric[which(complete.cases(methy_hotspots_pericentromeric))]
    x2<-methy_distal_telomeric[which(complete.cases(methy_distal_telomeric))]
    x1_svs<-svs_hotspots_pericentromeric[which(complete.cases(svs_hotspots_pericentromeric))]
    
    # #run wilcox-test just to check
    # wilcox.test(x = x1, y = x2, alternative = "less", paired = FALSE)
    # t.test(x = x1, y = x2, alternative = "less", paired = FALSE)
    
    #I will get the x1 value needed to be under the t critical
    pooled_stdv_num<-(length(x1)-1)*var(x1)+(length(x2)-1)*var(x2)
    pooled_stdv_den<-length(x1)+length(x2)-2
    pooled_stdv<-sqrt(pooled_stdv_num/pooled_stdv_den)
    den<-pooled_stdv*sqrt(1/length(x1)+1/length(x2))
    
    #I will dissipate the mean of x1 -> t_critical*den=num
    #choose t_critical -> for infinite values, one-tailed test at 0.05 = 	1.645
    #EN ESTE CASO YO BUSCO VER un lowe-tail test
    #CUANTO MENOS DE LA MEDIA DEBERIAN SER PARA NO SER IGUAL DE METHYLATED
    MEAN_X1_CRITICAL<-(1.645*den-mean(x2))*-1
    
    #I look for the values below the mean_critical
    x1_below<-x1[which(x1<MEAN_X1_CRITICAL)]
    #check they have increased SVs than the others in the group
    x1_svs_below<-x1_svs[which(names(x1_svs)%in%names(x1_below))]
    x1_svs_above<-x1_svs[-which(names(x1_svs)%in%names(x1_below))]
    #t.test(x = x1_svs_above, y = x1_svs_below, alternative = "less")
    test<-wilcox.test(x = x1_svs_above, y = x1_svs_below, alternative = "greater")
    
    #add values to table
    results_table_m[p,1]<-round(mean(x1_svs_above), digits = 3)
    results_table_m[p,2]<-round(mean(x1_svs_below), digits = 3)
    results_table_m[p,3]<-round(test$p.value, digits = 4)
    
    }#p
    results_table<-cbind(results_table, results_table_m)
    }#m
    results_table<-rbind(NA, results_table)
    results_table[1,]<-colnames(results_table)
    results_table[,1]<-row.names(results_table)
    colnames(results_table)<-c("Population", "CpG","","", "CHG","","")
    write.xlsx(results_table, "/home/fcasale/Desktop/Paper_2/3_RILs/7_window_types/Results/D.1_check_B_and_C_inconsistencies.xlsx", append = TRUE, sheetName = "TEST #3")
    

