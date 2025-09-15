setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2025")

rm(list = ls())

{
  indep_1 <- read.csv("bayestraits/output/processed_logs/indep_1.txt", sep = "\t")
  indep_2 <- read.csv("bayestraits/output/processed_logs/indep_2.txt", sep = "\t")
  indep_3 <- read.csv("bayestraits/output/processed_logs/indep_3.txt", sep = "\t")
  
  dep_1 <- read.csv("bayestraits/output/processed_logs/dep_1.txt", sep = "\t")
  dep_2 <- read.csv("bayestraits/output/processed_logs/dep_2.txt", sep = "\t")
  dep_3 <- read.csv("bayestraits/output/processed_logs/dep_3.txt", sep = "\t")
  
  dep_nopondloss_1 <- read.csv("bayestraits/output/processed_logs/dep_nopondloss_1.txt", sep = "\t")
  dep_nopondloss_2 <- read.csv("bayestraits/output/processed_logs/dep_nopondloss_2.txt", sep = "\t")
  dep_nopondloss_3 <- read.csv("bayestraits/output/processed_logs/dep_nopondloss_3.txt", sep = "\t")
  
  dep_nogain_1 <- read.csv("bayestraits/output/processed_logs/dep_nogain_1.txt", sep = "\t")
  dep_nogain_2 <- read.csv("bayestraits/output/processed_logs/dep_nogain_2.txt", sep = "\t")
  dep_nogain_3 <- read.csv("bayestraits/output/processed_logs/dep_nogain_3.txt", sep = "\t")
  
  dep_allzero_1 <- read.csv("bayestraits/output/processed_logs/dep_allzero_1.txt", sep = "\t")
  dep_allzero_2 <- read.csv("bayestraits/output/processed_logs/dep_allzero_2.txt", sep = "\t")
  dep_allzero_3 <- read.csv("bayestraits/output/processed_logs/dep_allzero_3.txt", sep = "\t")
  
  dep_res_1 <- read.csv("bayestraits/output/processed_logs/dep_restricted_1.txt", sep = "\t")
  dep_res_2 <- read.csv("bayestraits/output/processed_logs/dep_restricted_2.txt", sep = "\t")
  dep_res_3 <- read.csv("bayestraits/output/processed_logs/dep_restricted_3.txt", sep = "\t")
}

cols_dep <- match(c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"), colnames(dep_1))
cols_indep <- match(c("alpha1", "beta1", "alpha2", "beta2"), colnames(indep_1))

library(vioplot)

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(indep_1[,cols_indep], colMed = "black", main = "indep")
  vioplot(indep_2[,cols_indep], colMed = "black")
  vioplot(indep_3[,cols_indep], colMed = "black")
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(dep_1[,cols_dep], colMed = "black", main = "Dep")
  vioplot(dep_2[,cols_dep], colMed = "black")
  vioplot(dep_3[,cols_dep], colMed = "black")
  par(mfrow = c(1,1))}


{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(dep_nopondloss_1[,cols_dep], colMed = "black", main = "Dep_no_pond_loss")
  vioplot(dep_nopondloss_2[,cols_dep], colMed = "black")
  vioplot(dep_nopondloss_3[,cols_dep], colMed = "black")
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(dep_nogain_1[,cols_dep], colMed = "black", main = "Dep no gain")
  vioplot(dep_nogain_2[,cols_dep], colMed = "black")
  vioplot(dep_nogain_3[,cols_dep], colMed = "black")
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(dep_allzero_1[,cols_dep], colMed = "black", main = "Dep all zero")
  vioplot(dep_allzero_2[,cols_dep], colMed = "black")
  vioplot(dep_allzero_3[,cols_dep], colMed = "black")
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  vioplot(dep_res_1[,cols_dep], colMed = "black", main = "Dep_res")
  vioplot(dep_res_2[,cols_dep], colMed = "black")
  vioplot(dep_res_3[,cols_dep], colMed = "black")
  par(mfrow = c(1,1))}


################### Compare Models with 2logBFs ###########################################
rm(list = ls())

library(readr) #version ‘2.1.4’ used for function read_file
mods <- c("indep", "dep", "dep_nopondloss", "dep_nogain", "dep_allzero",
          "dep_restricted")

BF_matrix <- as.data.frame(matrix(nrow = length(mods)*3, ncol = 3))
colnames(BF_matrix) <- c("Model_class", "replicate", "LogMargLik")
BF_matrix$Model_class <- rep(mods, each = 3)
BF_matrix$replicate <- c(rep(1:3, length(mods)))

for(i in 1:dim(BF_matrix)[1]){
  A <- read_file(paste("bayestraits/output/processed_stones/", BF_matrix$Model_class[i], 
                       "_", BF_matrix$replicate[i], ".stones.txt", sep = ""))
  BF_matrix$LogMargLik[i] <- round(as.numeric(gsub("Log marginal likelihood:\t", "", 
                                             gsub("\r\n", "", A))),2)
}

BF_matrix


# This matrix is a simple compilation of estimated Log Marginal Likelihoods (LMLs)
# for each individual model run in BayesTraits. The code depends on how you set
# up your directories for where to store individual BayesTraits outputs and requires
# the steppingstones utility to estimate LMLs. The highest log marginal likelihood 
# is the model run that did the best job fitting to the data, and the highest individual 
# run for each model is the thing we are most interested in. 

avgs <- sapply(split(BF_matrix$LogMargLik, BF_matrix$Model_class), mean) 
SDs <- sapply(split(BF_matrix$LogMargLik, BF_matrix$Model_class), sd) 
LMLs <- sapply(split(BF_matrix$LogMargLik, BF_matrix$Model_class), max) 

cbind(round(avgs, 2), round(SDs, 2))


# here we go through and pick the best individual run from each model group to compare

BF_test <- as.data.frame(matrix("-", nrow = length(mods), ncol = length(mods)))
colnames(BF_test) <- gsub("eight_state_", "", names(LMLs))
rownames(BF_test) <- gsub("eight_state_", "", names(LMLs))

for(i in 1:dim(BF_test)[1]){
  for(j in 1:dim(BF_test)[1]){
    BF_test[i,j] <- round(2*(LMLs[i] - LMLs[j]),2)
  }}

BF_test
# this matrix does the explicit comparison, and creates 2*log Bayes Factors (2LogBFs)
# for each comparison between any two models. A 2LogBF is calculated with this 
# equation: 2*(LML[Model1] - LML[Model2]). This setup has several published 
# recommendations for interpretation, with higher numbers providing greater support
# for model 1, and negative numbers implying varying degrees of support for model 2. 


###################### get rate estimates for different models###############################


source("lung_loss_git/scripts/functions/summary_posterior_function.R")

read_in_model = function(model){
  name <- paste("bayestraits/output/processed_logs/", model, ".txt", sep = "")
  data <- read.csv(name, sep = "\t")
  return(data)}

BF_matrix

{
  model <- "dep_allzero_1"
  summary_posterior(read_in_model(model), grep("q", colnames(read_in_model(model))))
}


