SVs<-c("deletions", "duplications", "insertions") #\"translocations\")
pop_info<-read.csv("/scratch/Federico/3_RILs/1_Marius_data/sources/A.0_pop.info.csv")
pop_info<-as.matrix(pop_info)
Populations<-pop_info[c(10,24,25),1]
layers<-c("1_first", "2_second", "3_third")


for (v in SVs){
for (l in layers){
for (p in Populations) {
for (r in 1:100) {

writeLines(paste("  

r<-\"",r,"\"
p<-\"",p,"\"
v<-\"",v,"\"
l<-\"",l,"\"

library(data.table)

make.matrix<-function(x){ if (isFALSE(is.matrix(x))){
  new.matrix<-matrix(nrow = 1, ncol = length(x)); colnames(new.matrix)<-names(x)
  new.matrix[,1:ncol(new.matrix)]<-x; return(new.matrix)}else{return(x)}}

get.closest.SV<-function(x=x){
  dif<-abs(x-as.numeric(SVs_p_c[,2:3]))  
  closest_SV<-SVs_p_c[,2:3][which(dif==min(dif))[1]]
  names(closest_SV)<-row.names(which(SVs_p_c==closest_SV, arr.ind=TRUE))[1]   #puede haber dos (muy raro), cuando ambos padres tienen SVs distintas pero que empiezan o terminan en a misma position
  #return(make.matrix(c(x,names(closest_SV), SVs_p_c[names(closest_SV),2:3], closest_SV, min(dif))))
  return(make.matrix(c(x,min(dif), SVs_p_c[names(closest_SV),ncol(SVs_p_c)])))
}

get.distances<-function(y=y){
  cat(which(random_COs[1,]==y[1])); cat(\"-\")
  CO.and.closest.SVs<-lapply(y, FUN = get.closest.SV)
  CO.and.closest.SVs<-matrix(unlist(CO.and.closest.SVs), nrow = ncol(CO.and.closest.SVs[[1]]), ncol = length(CO.and.closest.SVs))
  CO.and.closest.SVs<-transpose(as.data.frame(CO.and.closest.SVs))
  #colnames(CO.and.closest.SVs)<-c(\"breakpoint\", \"SV_code\",\"start\",\"end\",\"closest\",\"distance\")
  colnames(CO.and.closest.SVs)<-c(\"breakpoint\",\"distance_to_SV\", \"SV_length\")
  return(CO.and.closest.SVs)
}

chrs<-paste(\"chr\", 1:7, \"H\", sep = \"\")
chr_lengths<-as.matrix(read.csv(\"/scratch/Federico/3_RILs/sources/Chr_length_V3_byFede.csv\"))

sv_list<-readRDS(\"/scratch/Federico/3_RILs/3_SVs/Results/A_SVs_per_population.RDS\")
breakpoint_list<-readRDS(\"/scratch/Federico/3_RILs/2_CO_breakpoints/Results/D.3_recombination_layers/D.3.4_breakpoints_list.RDS\")

runs<-10

dir.create(\"/scratch/Federico/3_RILs/3_SVs/Results/B.2_SV_per_RANDOM_CO\")

CO_SV_list<-list()
#for (v in SVs){ cat(v); cat(\": \", fill = TRUE)
  CO_SV_list[[v]]<-list()
  dir.create(paste(\"/scratch/Federico/3_RILs/3_SVs/Results/B.2_SV_per_RANDOM_CO/\",v, sep = \"\"))
#  for (l in layers){ cat(l); cat(\" - \")
    CO_SV_list[[v]][[l]]<-list()
    dir.create(paste(\"/scratch/Federico/3_RILs/3_SVs/Results/B.2_SV_per_RANDOM_CO/\",v,\"/\",l, sep = \"\"))
#    for (p in Populations){ cat(p);cat(\": \", fill = TRUE)   
      CO_SV_list[[v]][[l]][[p]]<-list()
      for (c in 1:7){ cat(\"chr\"); cat(c); cat(\": \")
        SVs_p_c<-sv_list[[v]][[p]][[c]]  
        row.names(SVs_p_c)<-paste(v,\"_\",c,\"_\",1:nrow(SVs_p_c), sep = \"\")
        #add lengths
        SVs_p_c<-cbind(SVs_p_c, NA); SVs_p_c[,ncol(SVs_p_c)]<-as.numeric(SVs_p_c[,3])-as.numeric(SVs_p_c[,2])
        #create random breakpoints    ##### USO LA MISMA CANTIDAD DE COs que YA TENGO. 
        
        #by function
        random_COs<-matrix(nrow = nrow(breakpoint_list[[l]][[p]][[c]]), ncol = runs)
        random_COs[,]<-sample(1:as.numeric(chr_lengths[c, 3]), length(random_COs))
        CO_SV_list[[v]][[l]][[p]][[c]]<-apply(random_COs, 2, FUN = get.distances) 
        
        #by loop
        # CO_SV_list[[v]][[l]][[p]][[c]]<-list()
        # for (r in 1:runs){ cat(r); cat(\"-\")
        # random_COs<-sample(1:as.numeric(chr_lengths[c, 3]), nrow(breakpoint_list[[l]][[p]][[c]]))
        # random_COs<-lapply(random_COs, FUN = get.closest.SV)
        # random_COs<-matrix(unlist(random_COs), nrow = ncol(random_COs[[1]]), ncol = length(random_COs))
        # random_COs<-as.matrix(transpose(as.data.frame(random_COs)))
        # colnames(random_COs)<-c(\"breakpoint\", \"SV_code\",\"start\",\"end\",\"closest\",\"distance\")
        # CO_SV_list[[v]][[l]][[p]][[c]][[r]]<-random_COs  
        # }
        
        cat(\" \")
      }#c
      cat(fill = TRUE)
      saveRDS(CO_SV_list, paste(\"/scratch/Federico/3_RILs/3_SVs/Results/B.2_SV_per_RANDOM_CO/\",v,\"/\",l,\"/B.2_\",v,\"_\",l,\"_\",p,\"_\",r,\".RDS\", sep = \"\"))
 #   }#p
    cat(fill = TRUE) 
#  }#l
#}#v

", sep = ""), 
           
paste("/scratch/Federico/3_RILs/3_SVs/R_scripts/R_SV_per_RANDOM_CO_breakpoint_",v,"_",l,"_",p,"_",r,".R", sep = "")
)
  
# writeLines(
# 
# paste(
# "
# #!/bin/bash                                                                     
# #PBS -l select=1:ncpus=1:mem=40gb                                               
# #PBS -l walltime=24:00:00                                                       
# 
# R CMD BATCH /../../scratch/Federico/3_RILs/3_SVs/R_scripts/B.2_SV_per_RANDOM_CO_breakpoint_",v,"_",p,"_",r,".R
# ", sep=""),
#     
# paste("/scratch/Federico/3_RILs/3_SVs/R_scripts/job_",p,"_",r,".sh", sep = "")
#     
# )  
  
}#r
}#p
}#l
}#v


#loop for linux 
for r in {1..10}; do for p in {13,27,28}; do echo "sh job_HvDRR${p}_${r}.sh"; done | bash -;done
for r in {1..100}; do for p in {13,27,28}; do echo "R CMD BATCH R_SV_per_RANDOM_CO_breakpoint_insertions_3_third_HvDRR${p}_${r}.R"; done | bash -;done

#gather runs
for (v in SVs){
for (p in Populations){
pop_list<-list()
for (r in 1:10){}}}
r<-1
run<-readRDS(paste("/scratch/Federico/3_RILs/3_SVs/Results/B.2_SV_per_RANDOM_CO/",v,"/3_third/B.2_",v,"_3_third_",p,"_",r,".RDS", sep = ""))  
}  
}}}