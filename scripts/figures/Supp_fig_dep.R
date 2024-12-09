
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



source("lung_loss_git/scripts/functions/combine_posteriors.R")
indep <- combine_posteriors(indep_1, indep_2, indep_3)
dep <- combine_posteriors(dep_1, dep_2, dep_3)
dep_nopondloss <- combine_posteriors(dep_nopondloss_1, dep_nopondloss_2, dep_nopondloss_3)
dep_nogain <- combine_posteriors(dep_nogain_1, dep_nogain_2, dep_nogain_3)
dep_allzero <- combine_posteriors(dep_allzero_1, dep_allzero_2, dep_allzero_3)
dep_res <- combine_posteriors(dep_res_1, dep_res_2, dep_res_3)

library(readr) #version ‘2.1.4’ used for function read_file

BF_matrix <- as.data.frame(matrix(nrow = 6*3, ncol = 3))
colnames(BF_matrix) <- c("Model_class", "replicate", "LogMargLik")
BF_matrix$Model_class <- c(rep("indep", 3), rep("dep", 3), 
                           rep("dep_nopondloss", 3) , rep("dep_nogain", 3), 
                           rep("dep_allzero", 3),  rep("dep_restricted", 3))
BF_matrix$replicate <- c(rep(1:3, 6))

for(i in 1:dim(BF_matrix)[1]){
  A <- read_file(paste("bayestraits_exports/stones/", BF_matrix$Model_class[i], 
                       "_", BF_matrix$replicate[i], ".stones.txt", sep = ""))
  BF_matrix$LogMargLik[i] <- as.numeric(gsub("Log marginal likelihood:\t", "", 
                                             gsub("\r\n", "", A)))
}

indep_LML <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "indep")]))
dep_LML <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "dep")]))
dep_nopondloss_LML <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "dep_nopondloss")]))
dep_nogain_LML <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "dep_nogain")]))
dep_allzero_LML <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "dep_allzero")]))
dep_res_LML <- mean(as.numeric(BF_matrix$LogMargLik[which(BF_matrix$Model_class == "dep_restricted")]))

2*(dep_LML - indep_LML) 
#dep supported

2*(dep_LML - dep_nopondloss_LML) 
# no pond loss supported

2*(dep_LML - dep_nogain_LML) 
# no regains supported

2*(dep_LML - dep_allzero_LML) 
# both no regains and no pond losses

2*(dep_LML - dep_res_LML) 
#dep_res supported (same as all zeroes)

2*(dep_res_LML - dep_allzero_LML) 
# dep_allzero is best model

dep_final <- dep_allzero

save(dep_final, file = "lung_loss_git/processed_data/BT_output/best_four_model.Rdata")


indep_final <- indep

rm(list = c("dep_1", "dep_2", "dep_3", 
            "dep_nogain_1", "dep_nogain_2", "dep_nogain_3", 
            "dep_nopondloss_1", "dep_nopondloss_2", "dep_nopondloss_3", 
            "dep_allzero_1", "dep_allzero_2", "dep_allzero_3", 
            "dep_res_1", "dep_res_2", "dep_res_3", 
            "indep_1", "indep_2", "indep_3"))

source("lung_loss_git/scripts/functions/summary_posterior_function.R")

summary_dep <- summary_posterior(dep_final, grep("q", colnames(dep_final)))

summary_indep <- summary_posterior(indep_final, match(c("alpha1", "beta1", "alpha2", "beta2"), 
                                                      colnames(indep_final)))

ML_tree_4state <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_aqu.nex") 

data_dep <- read.csv(file = "lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")
data_dep <- data_dep[data_dep$Taxa %in% ML_tree_4state$tip.label,]

data_indep <- read.csv(file = "lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")
data_indep <- data_dep[data_dep$Taxa %in% ML_tree_4state$tip.label,]


state_labels <- four_state_state <- c("X_S", "Lu_S", "X_P", "Lu_P")
trait_col_names <- four_state_col <- c("ecology", "lung")
mode = "dependent"
model = "dependent"


source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")
cols <- setNames(c("purple", "blue", "red", "#B9FFB4"), four_state_state)

Q <- custom_Q(rate_sum_dep[,2], model = "dependent", scale = .001)
sim_dat_dep <- make_simmap_data(data_dep, ML_tree_4state, mode = "double binary", 
                            four_state_state, four_state_col)

A <- make.simmap(ML_tree_4state, sim_dat_dep, nsim = 100,
                 Q = Q, state_labels = four_state_state, pi=rep(1/4,4))









plotSimmap(A[[2]], cols, pts=F, fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)



{jpeg(file = "figures/supp_fig_dep_phylo.jpg", bg = "transparent", units = "in", res = 1000, width = 10, height = 10)
    
  plotTree(A[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(A[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) {
    plot(A[[i]], colors=sapply(cols, make.transparent,alpha=0.02),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  
  dev.off()}




mode = "independent"
model = "independent"

source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")
cols <- setNames(c("purple", "blue", "red", "#B9FFB4"), four_state_state)

Q <- custom_Q(rate_sum_indep[,2], model = "independent", scale = .001)
sim_dat_indep <- make_simmap_data(data_indep, ML_tree_4state, mode = "double binary", 
                                four_state_state, four_state_col)

indep_simmap <- make.simmap(ML_tree_4state, sim_dat_indep, nsim = 100,
                 Q = Q, state_labels = four_state_state, pi=rep(1/4,4))

plotSimmap(indep_simmap[[1]], cols, pts=F, fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)



{jpeg(file = "figures/supp_fig_indep_phylo.jpg", bg = "transparent", units = "in", res = 1000, width = 10, height = 10)
  
  plotTree(indep_simmap[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(indep_simmap[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) {
    plot(indep_simmap[[i]], colors=sapply(cols, make.transparent,alpha=0.02),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  
  dev.off()}



                  






