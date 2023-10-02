# Libraries used
library(readr) #version ‘2.1.4’ used for function read_file

#######################################################################################
#                             Code overview

# 1 - upload data 
# perform burnin

# 2 - check for convergence of main models
# visually examine posterior distribution of different model runs to check convergence
# create_single_master posterior for converged models

# 3 - correlated evolution test
#            - use Log Marginal likelihoods to calculate Bayesfactors for comparison of independent and dependent models

# 4 - summarize posterior of the supported model
#     4a - examine model strings to see which rates are estimated at zero
#     4b - create figure of a representative posterior

# 5 - perform stochastic character mapping using bayestraits rate matrices

load("Bayestraits/data/master_export.Rdata")
load("Bayestraits/results/maxLH_dep_master_export.Rdata")
load("Bayestraits/results/tree_set_dep_master_export.Rdata")
load("Bayestraits/results/portik_dep_master_export.Rdata")


#######################################################################################
############################ upload and clean up data #################################
# these files have been edited using shell scripts (process_logs.sh in extra_scripts) to remove other text, 
# allowing direct importing of Bayestraits output files as .csv files
# there are three runs per model and each model is run on three different trees

dep_mcmc_1 <- read.csv("git/bayestraits/output/processed_logs/dep_mcmc_1.txt",sep = "\t")
dep_mcmc_2 <- read.csv("git/bayestraits/output/processed_logs/dep_mcmc_1.txt",sep = "\t")
dep_mcmc_3 <- read.csv("git/bayestraits/output/processed_logs/dep_mcmc_1.txt",sep = "\t")

indep_mcmc_1 <- read.csv("git/bayestraits/output/processed_logs/indep_mcmc_1.txt",sep = "\t")
indep_mcmc_2 <- read.csv("git/bayestraits/output/processed_logs/indep_mcmc_2.txt",sep = "\t")
indep_mcmc_3 <- read.csv("git/bayestraits/output/processed_logs/indep_mcmc_3.txt",sep = "\t")



#perform burnins
{
  
dep_mcmc_1      <- dep_mcmc_1[round(dim(dep_mcmc_1)[1]*.25) : dim(dep_mcmc_1)[1],]
dep_mcmc_2      <- dep_mcmc_2[round(dim(dep_mcmc_2)[1]*.25) : dim(dep_mcmc_2)[1],]
dep_mcmc_3      <- dep_mcmc_3[round(dim(dep_mcmc_3)[1]*.25) : dim(dep_mcmc_3)[1],]

indep_mcmc_1    <- indep_mcmc_1[round(dim(indep_mcmc_1)[1]*.25) : dim(indep_mcmc_1)[1],]
indep_mcmc_2    <- indep_mcmc_2[round(dim(indep_mcmc_2)[1]*.25) : dim(indep_mcmc_2)[1],]
indep_mcmc_3    <- indep_mcmc_3[round(dim(indep_mcmc_3)[1]*.25) : dim(indep_mcmc_3)[1],]
}

  
{
six_state <- six_state[round(dim(six_state)[1]*.25) : dim(six_state)[1],]
  
  
max_LH_dep1      <- max_LH_dep1[round(dim(max_LH_dep1)[1]*.25) : dim(max_LH_dep1)[1],]
max_LH_dep2      <- max_LH_dep2[round(dim(max_LH_dep2)[1]*.25) : dim(max_LH_dep2)[1],]
max_LH_dep3      <- max_LH_dep3[round(dim(max_LH_dep3)[1]*.25) : dim(max_LH_dep3)[1],]
max_LH_indep1    <- max_LH_indep1[round(dim(max_LH_indep1)[1]*.25) : dim(max_LH_indep1)[1],]
max_LH_indep2    <- max_LH_indep2[round(dim(max_LH_indep2)[1]*.25) : dim(max_LH_indep2)[1],]
max_LH_indep3    <- max_LH_indep3[round(dim(max_LH_indep3)[1]*.25) : dim(max_LH_indep3)[1],]
tree_set_dep1    <- tree_set_dep1[round(dim(tree_set_dep1)[1]*.25) : dim(tree_set_dep1)[1],]
tree_set_dep2    <- tree_set_dep2[round(dim(tree_set_dep2)[1]*.25) : dim(tree_set_dep2)[1],]
tree_set_dep3    <- tree_set_dep3[round(dim(tree_set_dep3)[1]*.25) : dim(tree_set_dep3)[1],]
tree_set_indep1  <- tree_set_indep1[round(dim(tree_set_indep1)[1]*.25) : dim(tree_set_indep1)[1],]
tree_set_indep2  <- tree_set_indep2[round(dim(tree_set_indep2)[1]*.25) : dim(tree_set_indep2)[1],]
tree_set_indep3  <- tree_set_indep3[round(dim(tree_set_indep3)[1]*.25) : dim(tree_set_indep3)[1],]
portik_dep1      <- portik_dep1[round(dim(portik_dep1)[1]*.25) : dim(portik_dep1)[1],]
portik_dep2      <- portik_dep2[round(dim(portik_dep2)[1]*.25) : dim(portik_dep2)[1],]
portik_dep3      <- portik_dep3[round(dim(portik_dep3)[1]*.25) : dim(portik_dep3)[1],]
portik_indep1    <- portik_indep1[round(dim(portik_indep1)[1]*.25) : dim(portik_indep1)[1],]
portik_indep2    <- portik_indep2[round(dim(portik_indep2)[1]*.25) : dim(portik_indep2)[1],]
portik_indep3    <- portik_indep3[round(dim(portik_indep3)[1]*.25) : dim(portik_indep3)[1],]}


#######################################################################################
############################# Check within-run model convergence #################################
#######################################################################################

posterior_summary = function(data,column.names){
  if (sum(dim(data)[2]) != 0){
    cols <- match(column.names, colnames(data))
    par(mfrow = c(dim(data[cols])[2]/4, 4))
    for(i in 1:dim(data[,cols])[2]){
      smoothScatter(data$Iteration, data[,cols[i]], ylim = c(0,max(data[,cols])),
                    cex=4,nr=500, xlab = "iter", 
                    ylab = colnames(data)[cols[i]])}} else{
                      par(mfrow = c(1,1))
                      smoothScatter(data$Iteration, data[,cols],cex=4,nr=500, xlab = "iter", ylab = colnames(rows))
                    }
}

posterior_summary(dep_mcmc_1, c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"))
posterior_summary(dep_mcmc_2, c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"))
posterior_summary(dep_mcmc_3, c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"))

posterior_summary(indep_mcmc_1, c("alpha1", "alpha2", "beta1", "beta2"))
posterior_summary(indep_mcmc_2, c("alpha1", "alpha2", "beta1", "beta2"))
posterior_summary(indep_mcmc_3, c("alpha1", "alpha2", "beta1", "beta2"))



#######################################################################################
############################# Check across run model convergence #################################
#######################################################################################

#little function to reduce clutter (found in extra_scripts)

source("git/scripts/extra_scripts/violin_comparison_function.R")

violin_comparison_dep(dep_mcmc_1, dep_mcmc_2, dep_mcmc_3, "dep")
violin_comparison(max_LH_indep1, max_LH_indep2, max_LH_indep3, "Max_LH_indep")



######################################################################################################
# for each converged run, create a master 

master_export = function(data1, data2, data3){
  master <- as.data.frame(matrix(ncol = dim(data1)[2], nrow = (dim(data1)[1]*3)))
  colnames(master) <- colnames(data1)
  master$run_number <- c(rep(1,dim(data1)[1]), rep(2,dim(data2)[1]), rep(3,dim(data3)[1]))
  
  master[which(master$run_number == 1),1:dim(data1)[2]] <- data1[1:dim(data1)[2]]
  master[which(master$run_number == 2),1:dim(data2)[2]] <- data2[1:dim(data2)[2]]
  master[which(master$run_number == 3),1:dim(data3)[2]] <- data3[1:dim(data3)[2]]
  master <- as.data.frame(master)
  return(master)
}

master_dep <- master_export(dep_mcmc_1, dep_mcmc_2, dep_mcmc_3)
master_indep <- master_export(indep_mcmc_1, indep_mcmc_2, indep_mcmc_3)

master_dep_export  <-  master_dep[seq(1,dim(master_dep)[1], length.out = 100000),]
master_indep_export <- master_indep[seq(1,dim(master_indep)[1], length.out = 100000),]


save(master_dep_export, file = "git/data/master_dep_export.Rdata")
save(master_indep_export, file = "git/data/master_indep_export.Rdata")



#remove actual runs to de-clutter

rm(list = c("dep_mcmc_1", "dep_mcmc_2", "dep_mcmc_3", "indep_mcmc_1", "indep_mcmc_2", 
            "indep_mcmc_3", "posterior_summary", "violin_comparison"))

#####################################################################################################

########## Dollo's Law check

#use reverse jump to test what proportion of the posterior is spent in a true Dollo's Law frame (no regains),
# or else no regains in certain conditions, and if lungs are ever lost outside streams

#function meant for cleaned up BT output with model string intact
RJ_model_testing = function(data){
  models <- unique(data$Model.string)
  model_matrix <- matrix(nrow = length(models), ncol = 3)
  colnames(model_matrix) <- c("model_string", "number_observed", "percent_of_posterior")
  rownames(model_matrix) <- 1:length(models)
  model_matrix[,1] <- models
  for(i in 1:length(models)){
    model_matrix[i,2] <- length(which(data$Model.string %in% model_matrix[i,1]))
  }
  runs <- dim(data)[1]
  
  model_matrix[,3] <- (as.numeric(model_matrix[,2])/runs)
  output <- matrix(NA,nrow = 2, ncol = 3)
  colnames(output) <- c("No regains", "No lentic regains", "No lentic losses")
  rownames(output) <- c("posterior odds", "BayesFactor")
  output[1,1] <- sum(as.numeric(model_matrix[grep("Z . . . . Z . .", model_matrix[,1]),3]))
  output[1,2] <- sum(as.numeric(model_matrix[grep(". . . . . Z . .", model_matrix[,1]),3]))
  output[1,3] <- sum(as.numeric(model_matrix[grep(". . . . . . . Z", model_matrix[,1]),3]))
  output[2,1] <- output[1,1]/(1-output[1,1])
  output[2,2] <- output[1,2]/(1-output[1,2])
  output[2,3] <- output[1,3]/(1-output[1,3])
  
  return(output)
}

RJ_model_testing(master_dep)


#######################################################################################
########################### Test of correlated evolution ##############################
library(readr) #version ‘2.1.4’ used for function read_file

# Bayestraits outputs a .stones file upon the stones command with an estimate log marginal likl
# those files have  been edited by a shell script (process_stones.sh in extra_scripts) to pull
# out only the log marginal likelihood as a single number and provide a specific naming convention
# tree_model_replicate_LML.txt -> portik_dep_1_LML.txt for example

# master matrix of marginal likelihoods for Bayesfactor testing

BF_matrix <- as.data.frame(matrix(nrow = 18, ncol = 4))
colnames(BF_matrix) <- c("Model_class", "tree", "replicate", "LogMargLik")
BF_matrix$Model_class <- rep(c(rep("dep", 3), rep("indep", 3)),3)
BF_matrix$tree <- c(rep("maxLH", 6), rep("tree_set", 6), rep("portik", 6))
BF_matrix$replicate <- c(rep(1:3, 2))

for(i in 1:6){
  BF_matrix$LogMargLik[i] <- read_file(paste("stones/", BF_matrix$tree[i], "_", BF_matrix$Model_class[i], 
                                  "_", BF_matrix$replicate[i], "_LML.txt", sep = ""))}

# For actual Bayesfactor test, we use averages of each model class (if converged)
dep <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "dep")]))
indep <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "indep")]))

# Bayesfactor test:
#                     BF = 2*(LML_complex - mod2_simple)
#                     higher BF's mean more support for more complex models

BF = 2*(dep - indep)

#just for fun
BF_vis = function(complex, simple){
  BF = 2*(complex - simple)
  if(BF < 2){paste("no support for complex model")}
  if(5 > BF & BF >= 2){paste("positive evidence for complex model")}
  if(10 > BF & BF > 5){paste("strong evidence for complex model")}
  if(BF > 10){paste("very strong evidence for complex model")}}

BF_vis(dep, indep)
##### we get strong evidence for the dependent model, suggesting there is correlated evolution happening

rm(list = c("BF_matrix", "BF", "dep", "indep", "i"))
########################################################################################################
### Examine the dependent run in more detail, using the master log file of the posterior distribution
########################################################################################################

load("Bayestraits/data/master_export.Rdata")
master_maxLH_dep <- master_export


rm(master_export)



#################################################################################
####                visualize posterior for figure 3A
library(HDInterval)
library(vioplot)
# calculate highest density intervals of the master posterior
#BT_data <- master_maxLH_dep
#states_OI <- c("q12", "q13", "q15", "q21", "q24", "q31", "q34", "q36", "q42", "q43", "q51", "q56", "q63", "q65")
states_OI <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
BT_data <- master_dep_export[,match(states_OI, colnames(master_dep_export))]*.001

posterior_intervals <- hdi(BT_data, credMass = .95)

#lungless

{pdf(file = "figures/rate_violin_plots.pdf", bg = "transparent", width = 5, height = 3.5)
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                
labels <- c("gain lungs (lotic)", "lotic to lentic (no lungs)", "lose lungs (lotic)",
            "lotic to lentic (lunged)", "lentic to lotic (no lungs)", "gain lungs (lentic)",
            "lentic to lotic (lunged)", "lose lungs (lentic)")
par(mar = c(5.1,4.1,4.1,2.1)) #reset default margins
par(mfrow=c(1,1))
  par(xpd = TRUE)
  vioplot(BT_data, lineCol = "transparent", rectCol = "transparent",  xaxt = "n", yaxt = "n",
          colMed = "transparent", col = "darkgrey", main = "",
          ylab = "")
  #segments(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3], col = "white", lwd = 2)  
  text(x = (1:length(labels))+((par("usr")[2]-par("usr")[1])*.05), y = par("usr")[3]-(par("usr")[4]-par("usr")[3])*.05,
       labels = labels, xpd = NA, adj = 1, srt = 25, cex = .65)
  
  axis(2,cex.axis=.75)
  
  #legend(0,((par("usr")[3]) - (par("usr")[3])*.025), legend = labels, ncol = 3, fill = "darkgray", bty = "n", cex = .75)
  for(i in 1:dim(posterior_intervals)[2]){
    #segments(i, 0, i, par("usr")[3]-(par("usr")[4]-par("usr")[3])*.05, lty = "dashed")
    segments(i, posterior_intervals[2,i], i,posterior_intervals[1,i], lwd = 1.75, col = "white")
    segments(i, posterior_intervals[2,i], i,posterior_intervals[1,i], lwd = 1.5, col = "black")
    segments(i - ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[1,i], 
             i + ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[1,i], lwd = 1.75, col = "white")
    segments(i - ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[2,i], 
             i + ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[2,i], lwd = 1.75, col = "white")
    segments(i - ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[1,i], 
             i + ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[1,i], lwd = 1.5, col = "black")
    segments(i - ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[2,i], 
             i + ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[2,i], lwd = 1.5, col = "black")
    points(i, median(BT_data[,i]), pch = 21, bg = "white", cex = 1.25)
  }
dev.off()}
   

avgs_dep <- sapply(master_dep_export[,grep("q12", colnames(master_dep_export)):(grep("q12", colnames(master_dep_export))+7)], "mean");
avgs_indep <- sapply(master_indep_export[,grep("alpha1", colnames(master_indep_export)):(grep("alpha1", colnames(master_indep_export))+3)], "mean");




require(scales)
avgs_scaled <- round(rescale(avgs_maxLH_dep,c(.5,7.5),c(0,5)),2)

round(avgs_scaled,2)

















