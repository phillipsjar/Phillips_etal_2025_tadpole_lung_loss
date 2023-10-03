# Libraries used
library(readr) #version ‘2.1.4’ used for function read_file

#######################################################################################
# 1 - correlated evolution test
#            - use Log Marginal likelihoods to calculate Bayesfactors for comparison of independent and dependent models

# 2 - Dollo's Law test of reversible evolution - how much support is there for some rates to be zero?
#     4a - examine model strings to see which rates are estimated at zero
#     4b - create figure of a representative posterior

# 3 - perform stochastic character mapping using bayestraits rate matrices


#######################################################################################
############################ upload data #################################

# first locate cleaned up and exported BayesTraits posteriors

load(file = "bayestraits_exports/master_dep_export.Rdata")
load(file = "bayestraits_exports/master_indep_export.Rdata")
load(file = "bayestraits_exports/master_six_state_dep_export.Rdata")
load(file = "bayestraits_exports/master_six_state_indep_export.Rdata")

master_dep_export <- master_dep
master_indep_export <- master_indep
rm(list = c("master_dep", "master_indep"))
#######################################################################################
########################### Test of correlated evolution ##############################
library(readr) #version ‘2.1.4’ used for function read_file

# Bayestraits outputs a .stones file upon the stones command with an estimate log marginal likl
# those files have  been edited by a shell script (process_stones.sh in extra_scripts) to pull
# out only the log marginal likelihood as a single number and provide a specific naming convention
# model_replicate_LML.txt -> dep_1_LML.txt for example

# master matrix of marginal likelihoods for Bayesfactor testing

BF_matrix <- as.data.frame(matrix(nrow = 12, ncol = 3))
colnames(BF_matrix) <- c("Model_class", "replicate", "LogMargLik")
BF_matrix$Model_class <- c(rep("dep", 3), rep("indep", 3), rep("six_state_dep", 3), rep("six_state_indep", 3))
BF_matrix$replicate <- c(rep(1:3, 2))

for(i in 1:dim(BF_matrix)[1]){
              A <- read_file(paste("bayestraits_exports/stones/", BF_matrix$Model_class[i], 
                                      "_", BF_matrix$replicate[i], ".stones.txt", sep = ""))
              BF_matrix$LogMargLik[i] <- as.numeric(gsub("Log marginal likelihood:\t", "", 
                                                         gsub("\r\n", "", A)))
      }







# For actual Bayesfactor test, we use averages of each model class (if converged)
dep <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "dep")]))
indep <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "indep")]))
six_dep <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "six_state_dep")]))
six_indep <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "six_state_indep")]))

# Bayesfactor test:
#                     2logBF = 2*(LML_complex - LML_simple)
#                     higher BF's mean more support for more complex models

2*(dep - indep)
2*(six_dep - six_indep)

#just for fun

BF_vis = function(complex, simple){
  BF = (complex - simple)
  BF = 2*BF
  if(BF < 2.2){paste("no support for complex model")}
  if(6 > BF & BF >= 2.2){paste("positive evidence for complex model")}
  if(10 > BF & BF > 6){paste("strong evidence for complex model")}
  if(BF > 10){paste("very strong evidence for complex model")}}

BF_vis(dep, indep)
BF_vis(six_dep, six_indep)


##### we get strong evidence for the dependent model, suggesting there is correlated evolution happening

rm(list = c("BF_matrix", "BF", "dep", "indep", "six_dep", "six_indep", "i",
            "complex", "simple", "BF_vis", "A"))


#####################################################################################################

########## Dollo's Law check

#use reverse jump to test what proportion of the posterior is spent in a true Dollo's Law frame (no regains),
# or else no regains in certain conditions, and if lungs are ever lost outside streams

source("lung_loss_git/scripts/functions/Dollo_check.R")

RJ_model_testing(master_dep_export, model = "dependent")
RJ_model_testing(master_indep_export, model = "independent")
RJ_model_testing(master_six_state_dep_export, model = "six_state")
RJ_model_testing(master_six_state_indep_export, model = "six_state")


