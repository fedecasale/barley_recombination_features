#https://www.r-bloggers.com/2021/04/deep-neural-network-in-r/#:~:text=Neural%20Network%20in%20R%2C%20Neural,the%20strength%20of%20this%20method.
#https://finnstats.com/index.php/2021/04/08/naive-bayes-classification-in-r/
#https://datascienceplus.com/neuralnet-train-and-test-neural-networks-using-r/ -> clave
#library(keras)
library(mlbench)
library(dplyr)
library(magrittr)
library(neuralnet)
#library(tensorflow)

variables<-readRDS("/scratch/Federico/3_RILs/6_graphics_for_paper/Results/C.2.0_variables_per_chr.RDS")

standarize<-function(x){x<-(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE); return(x)}
#normalize<-function(x){x<-(x-mean(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)); return(x)}

variables_geno<-as.data.frame(variables$genome)
variables_geno<-variables_geno[,4:ncol(variables_geno)]
for (i in 1:ncol(variables_geno)){variables_geno[,i]<-as.numeric(as.character(variables_geno[,i]))}

rows.to.delete<-c()
for (i in 1:nrow(variables_geno)){if (any(is.na(variables_geno[i,]))){rows.to.delete<-c(rows.to.delete, i)}}
variables_geno<-variables_geno[-c(rows.to.delete),]

#The neural network requires normalized values for better prediction, 
#just normalize or scale the predictor variables based on below function
variables_geno<-apply(variables_geno, 2, FUN = standarize)    
variables_geno<-as.data.frame(variables_geno)

colnames(variables_geno)<-gsub(" ", "_", colnames(variables_geno))
colnames(variables_geno)
#variable_selection
variables_geno<-variables_geno[,c(1,2,7,8,21,22)]; colnames(variables_geno)
colnames(variables_geno)[3]<-"Par_seq_diver"
colnames(variables_geno)[6]<-"DMRs_par_dif"

training_set<-
n <- neuralnet(Recombination_rate~Total_SVs+Averaged_DMRs+Par_seq_diver+Averaged_DMRs+DMRs_par_dif, #hay que cambiarle el nombre a DMR par dif.
               data = variables_geno,
               hidden = c(2,3), #Deciding on the number of hidden layers in a neural network is not an exact science. In fact, there are instances where accuracy will likely be higher without any hidden layers. Therefore, trial and error plays a significant role in this process. One possibility is to compare how the accuracy of the predictions change as we modify the number of hidden layers. For instance, using a (2,1) configuration ultimately yielded 92.5% classification accuracy for this example.
               linear.output = FALSE,
               stepmax=1e7, #The linear.output variable is set to FALSE, given the impact of the independent variables on the dependent variable (dividend) is assumed to be non-linear.
               threshold = 0.05, #debo pnerlo en 0.01, The threshold is set to 0.01, meaning that if the change in error during an iteration is less than 1%, then no further optimization will be carried out by the model.
               lifesign = 'full',
               rep=1)

plot(n,col.hidden = 'darkgreen',     
     col.hidden.synapse = 'darkgreen',
     show.weights = TRUE,
     information = FALSE,
     fill = 'lightblue')

n$result.matrix

# ###############
# data <- as.matrix(variables_geno)
# dimnames(data) <- NULL
# 
# #data partition
# set.seed(123) 
# ind <- sample(2, nrow(data), replace = T, prob = c(.7, .3))
# training <- data[ind==1,2:3]
# test <- data[ind==2, 2:3]
# trainingtarget <- data[ind==1, 1]
# testtarget <- data[ind==2, 1]
# 
# #ya esta scaled
# m <- colMeans(training)
# s <- apply(training, 2, sd)
# training <- scale(training, center = m, scale = s)
# test <- scale(test, center = m, scale = s)
# 
# 
# use_condaenv("r-tensorflow")
# 
# #model creation
# model <- keras_model_sequential()
# model %>%
#   layer_dense(units = 5, activation = 'relu', input_shape = c(2)) %>%
#   layer_dense(units = 1)
# install_tensorflow()
# 
# #Model Compilation
# model %>% compile(loss = 'mse',
#                   optimizer = 'rmsprop', 
#                   metrics = 'mae') 
# 
# #Model Fitting
# mymodel <- model %>%          
#   fit(training,trainingtarget,
#       epochs = 100,
#       batch_size = 32,
#       validation_split = 0.2)
# 
# #Prediction
# model %>% evaluate(test, testtarget)
# pred <- model %>% predict(test)
# mean((testtarget-pred)^2) 
# 
# #scatter plot
# plot(testtarget, pred) 
