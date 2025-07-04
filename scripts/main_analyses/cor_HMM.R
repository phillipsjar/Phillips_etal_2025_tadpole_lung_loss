library(corHMM)
library(ape)
library(phytools)
packageVersion("corHMM")
setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2024")

rm(list = ls())

#BT_tree_set <- read.nexus(file = "bayestraits/trees/tree_set.nex")
ML_tree_aqu <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_aqu.nex")

data_aqu <- read.csv("lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")
data_aqu$ecology[is.na(data_aqu$ecology)] <- "?"

#trim data to taxa present in trees
data_corHMM <- data_aqu[data_aqu$Taxa %in% ML_tree_aqu$tip.label,]

dataset_corHMM <- as.data.frame(cbind(data_corHMM$Taxa, data_corHMM$ecology, data_corHMM$lung))
colnames(dataset_corHMM) <- c("Genus_sp", "eco", "lung")

dataset_corHMM$eco[which(dataset_corHMM$eco == 1)] <- "lentic"
dataset_corHMM$eco[which(dataset_corHMM$eco == 0)] <- "lotic"

corHMMtest <- fitCorrelationTest(ML_tree_aqu, dataset_corHMM, simplified_models = TRUE)

{
AIC_indep_simp <- corHMMtest$simplified_independent_model_fit$AIC
AIC_indep_HMM_simp <- corHMMtest$simplified_hidden_Markov_independent_model_fit$AIC
AIC_dep_simp <- corHMMtest$simplified_correlated_model_fit$AIC
AIC_dep_HMM_simp <- corHMMtest$simplified_hidden_Markov_correlated_model_fit$AIC
AIC_indep <- corHMMtest$independent_model_fit$AIC
AIC_indep_HMM <- corHMMtest$hidden_Markov_independent_model_fit$AIC
AIC_dep <- corHMMtest$correlated_model_fit$AIC
AIC_dep_HMM <- corHMMtest$hidden_Markov_correlated_model_fit$AIC}


AIC_values <- c(AIC_indep_simp, AIC_indep_HMM_simp, AIC_dep_simp, AIC_dep_HMM_simp, 
           AIC_indep, AIC_indep_HMM, AIC_dep, AIC_dep_HMM)
names(AIC_values) <- c("AIC_indep_simp", "AIC_indep_HMM_simp", "AIC_dep_simp", "AIC_dep_HMM_simp", 
                       "AIC_indep", "AIC_indep_HMM", "AIC_dep", "AIC_dep_HMM")

maxLH_deltas <- rep(NA, length(AIC_values))
names(maxLH_deltas) <- names(AIC_values)
for(i in 1:length(AIC_values)){
  maxLH_deltas[i] <- AIC_values[i] - AIC_values[which.min(AIC_values)]
}

maxLH_deltas

corHMMtest$simplified_hidden_Markov_correlated_model_fit

corHMMtest$simplified_hidden_Markov_correlated_model_fit




