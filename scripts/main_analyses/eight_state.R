setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2024")

rm(list = ls())
par(mar = c(5,3,2,1))


# import processed BayesTraits posteriors
# use provided shell script to process the raw BT output files to remove extra text

{
eight_state_full_1 <- read.csv("bayestraits/output/processed_logs/eight_state_full_1.txt", sep = "\t")
eight_state_full_2 <- read.csv("bayestraits/output/processed_logs/eight_state_full_2.txt", sep = "\t")
eight_state_full_3 <- read.csv("bayestraits/output/processed_logs/eight_state_full_3.txt", sep = "\t")

eight_state_indep_1 <- read.csv("bayestraits/output/processed_logs/eight_state_indep_1.txt", sep = "\t")
eight_state_indep_2 <- read.csv("bayestraits/output/processed_logs/eight_state_indep_2.txt", sep = "\t")
eight_state_indep_3 <- read.csv("bayestraits/output/processed_logs/eight_state_indep_3.txt", sep = "\t")

eight_state_nopondloss_1 <- read.csv("bayestraits/output/processed_logs/eight_state_nopondloss_1.txt", sep = "\t")
eight_state_nopondloss_2 <- read.csv("bayestraits/output/processed_logs/eight_state_nopondloss_2.txt", sep = "\t")
eight_state_nopondloss_3 <- read.csv("bayestraits/output/processed_logs/eight_state_nopondloss_3.txt", sep = "\t")

eight_state_nostreamloss_1 <- read.csv("bayestraits/output/processed_logs/eight_state_nostreamloss_1.txt", sep = "\t")
eight_state_nostreamloss_2 <- read.csv("bayestraits/output/processed_logs/eight_state_nostreamloss_2.txt", sep = "\t")
eight_state_nostreamloss_3 <- read.csv("bayestraits/output/processed_logs/eight_state_nostreamloss_3.txt", sep = "\t")

eight_state_nogains_1 <- read.csv("bayestraits/output/processed_logs/eight_state_nogains_1.txt", sep = "\t")
eight_state_nogains_2 <- read.csv("bayestraits/output/processed_logs/eight_state_nogains_2.txt", sep = "\t")
eight_state_nogains_3 <- read.csv("bayestraits/output/processed_logs/eight_state_nogains_3.txt", sep = "\t")

eight_state_res_1 <- read.csv("bayestraits/output/processed_logs/eight_state_res_1.txt", sep = "\t")
eight_state_res_2 <- read.csv("bayestraits/output/processed_logs/eight_state_res_2.txt", sep = "\t")
eight_state_res_3 <- read.csv("bayestraits/output/processed_logs/eight_state_res_3.txt", sep = "\t")

eight_state_bestmod_1 <- read.csv("bayestraits/output/processed_logs/eight_state_bestmod_1.txt", sep = "\t")
eight_state_bestmod_2 <- read.csv("bayestraits/output/processed_logs/eight_state_bestmod_2.txt", sep = "\t")
eight_state_bestmod_3 <- read.csv("bayestraits/output/processed_logs/eight_state_bestmod_3.txt", sep = "\t")

eight_state_nopondgain_1 <- read.csv("bayestraits/output/processed_logs/eight_state_nopondgain_1.txt", sep = "\t")
eight_state_nopondgain_2 <- read.csv("bayestraits/output/processed_logs/eight_state_nopondgain_2.txt", sep = "\t")
eight_state_nopondgain_3 <- read.csv("bayestraits/output/processed_logs/eight_state_nopondgain_3.txt", sep = "\t")

}

# columns of interest (lung losses and lung gains)
columns <- match(c("q65", "q31", "q42", "q87", "q56", "q13", "q24", "q78"), colnames(eight_state_full_1))


#visualize the posteriors of three independent runs of a given model

library(vioplot)
{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(log(eight_state_full_1[,columns]+1), colMed = "black", main = "Eight full",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_full_2[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_full_3[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(log(eight_state_indep_1[,columns]+1), colMed = "black", main = "Eight indep",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_indep_2[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_indep_3[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(log(eight_state_nopondloss_1[,columns]+1), colMed = "black", main = "no pond loss",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_nopondloss_2[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_nopondloss_3[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  par(mfrow = c(1,1))}


{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(log(eight_state_nostreamloss_1[,columns]+1), colMed = "black", main = "no stream loss",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_nostreamloss_2[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_nostreamloss_3[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(log(eight_state_nogains_1[,columns]+1), colMed = "black", main = "no regains",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_nogains_2[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_nogains_3[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(log(eight_state_bestmod_1[,columns]+1), colMed = "black", main = "Eight best mod",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_bestmod_2[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_bestmod_3[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  par(mfrow = c(1,1))}

{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(log(eight_state_res_1[,columns]+1), colMed = "black", main = "Eight res",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_res_2[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_res_3[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  par(mfrow = c(1,1))}


{par(mfrow = c(3,1))
  par(mar = c(3,3,2,1))
  vioplot(log(eight_state_nopondgain_1[,columns]+1), colMed = "black", main = "Eight no pond gain",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_nopondgain_2[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  vioplot(log(eight_state_nopondgain_3[,columns]+1), colMed = "black",
          yaxt = "n")
  axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  par(mfrow = c(1,1))}

################### Compare Models with 2logBFs ###########################################

# to compare different models, we can look at the estimated log Marginal likelihoods
# this is estimated using the stepping stones method in Bayestraits
# use provided shell script to process .stones files from Bayestraits in bulk


library(readr) #version ‘2.1.4’ used for function read_file
rm(list = ls())
mods <- c("eight_state_full", "eight_state_indep", "eight_state_nopondloss",
          "eight_state_nostreamloss", "eight_state_nogains", "eight_state_nopondgain",
          "eight_state_res", "eight_state_bestmod")
BF_matrix <- as.data.frame(matrix(nrow = length(mods)*3, ncol = 3))
colnames(BF_matrix) <- c("Model_class", "replicate", "LogMargLik")
BF_matrix$Model_class <- rep(mods, each = 3)
BF_matrix$replicate <- c(rep(1:3, length(mods)))

# this loop reads processed shell files. This code is built around using the 
# same naming scheme everywhere, so be sure to save Bayestraits files with a 
# specific naming scheme that can be called later.

for(i in 1:dim(BF_matrix)[1]){
  A <- read_file(paste("bayestraits/output/processed_stones/", BF_matrix$Model_class[i], 
                       "_", BF_matrix$replicate[i], ".stones.txt", sep = ""))
  BF_matrix$LogMargLik[i] <- as.numeric(gsub("Log marginal likelihood:\t", "", 
                                             gsub("\r\n", "", A)))
}

BF_matrix

# This matrix is a simple compilation of estimated Log Marginal Likelihoods (LMLs)
# for each individual model run in BayesTraits. The code depends on how you set
# up your directories for where to store individual BayesTraits outputs and requires
# the steppingstones utility to estimate LMLs. The highest log marginal likelihood 
# is the model run that did the best job fitting to the data, and the highest individual 
# run for each model is the thing we are most interested in. 


avgs <- sapply(split(BF_matrix$LogMargLik, BF_matrix$Model_class), mean) 
# average LML for each model type
SDs <- sapply(split(BF_matrix$LogMargLik, BF_matrix$Model_class), sd) 
# SD of LML for each model type. High SDs imply that different runs are not converging
# on a single solution.

LMLs <- sapply(split(BF_matrix$LogMargLik, BF_matrix$Model_class), max) 
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

paste(BF_matrix$Model_class[which.max(BF_matrix$LogMargLik)], BF_matrix$replicate[which.max(BF_matrix$LogMargLik)], sep = "_")
# here is the best individual model

# and a quick visual look at what the model says about lung evolution
# the first four violin plots show different types of lung loss 
# (lentic, generalized lotic, specialized lotic, and terrestrial), while 
# the final four plots show rates of regain in those contexts.

eight_state_bestmod_1 <- read.csv("bayestraits/output/processed_logs/eight_state_bestmod_1.txt", sep = "\t")
columns <- match(c("q65", "q31", "q42", "q87", "q56", "q13", "q24", "q78"), colnames(eight_state_bestmod_1))
library(vioplot)

vioplot(log(eight_state_bestmod_1[,columns]+1), colMed = "black", main = "Eight best mod",
        yaxt = "n")
axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
     tck = -.01, lwd.ticks = 2, lwd = 1)
abline(h = 0, col = "red", lty = 2)

# this model actually does not converge across runs, with very different LMLs across 
# some of the independent runs. Additionally, trace plots for even the best individual 
# run do not converge on a single answer. 

colfunc <- colorRampPalette(c("#0000FF01", "#AAAAAA01", "#FF000001"))

colors <- colfunc(1000)
library(scales)
lh_scaled <- round(rescale(eight_state_bestmod_1$Lh, to = c(1, 1000), from = range(eight_state_bestmod_1$Lh)))
# added color scheme to show higher likelihood peaks in red and lower lH valleys in blue

plot(eight_state_bestmod_1$Iteration, eight_state_bestmod_1$q12, pch = 21, col = "transparent",cex = .5,
     bg = colors[lh_scaled])

# this makes trusting parameter estimates from this model difficult, even if the model
# provides evidence for specific hypotheses. Instead, we can use the best model that converged
# to get potential parameter ranges. The best model that converged is the third restricted 
# regains model iteration.

eight_state_res_3 <- read.csv("bayestraits/output/processed_logs/eight_state_res_3.txt", sep = "\t")

vioplot(log(eight_state_res_3[,columns]+1), colMed = "black", main = "Eight regains res",
        yaxt = "n")
axis(side = 2, at = log(c(1,2,6,21,51)), labels = c(0,.01,.05,.2,.5),
     tck = -.01, lwd.ticks = 2, lwd = 1)
abline(h = 0, col = "red", lty = 2)

lh_scaled <- round(rescale(eight_state_res_3$Lh, to = c(1, 1000), from = range(eight_state_res_3$Lh)))
plot(eight_state_res_3$Iteration, eight_state_res_3$q21, pch = 21, col = "transparent",cex = .5,
     bg = colors[lh_scaled])


###################### get rate estimates for different models###############################


source("lung_loss_git/scripts/functions/summary_posterior_function.R")

read_in_model = function(model){
  name <- paste("bayestraits/output/processed_logs/", model, ".txt", sep = "")
  data <- read.csv(name, sep = "\t")
  return(data)}

BF_matrix

{
model <- "eight_state_full_1"
summary_posterior(read_in_model(model), grep("q", colnames(read_in_model(model))))
}






