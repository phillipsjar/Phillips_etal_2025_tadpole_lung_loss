setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2024")
rm(list = ls())

indep_1 <- read.csv("bayestraits/output/processed_logs/indep_1.txt", sep = "\t")
dep_1 <- read.csv("bayestraits/output/processed_logs/dep_1.txt", sep = "\t")

source("lung_loss_git/scripts/functions/summary_posterior_function.R")

rate_sum_dep <- summary_posterior(dep_1, grep("q", colnames(dep_1)))
a <- summary_posterior(indep_1, 7:10)
rate_sum_indep <- rate_sum_dep
rate_sum_indep[,] <- NA

{rate_sum_indep[c(1),] <- a[3,]
 rate_sum_indep[c(6),] <- a[3,]
 rate_sum_indep[c(2),] <- a[1,]
 rate_sum_indep[c(4),] <- a[1,]
 rate_sum_indep[c(3),] <- a[4,]
 rate_sum_indep[c(8),] <- a[4,]
 rate_sum_indep[c(5),] <- a[2,]
 rate_sum_indep[c(7),] <- a[2,]}

rm(list = c("dep_1", "indep_1", "a"))
library(ape)
library(phytools)
ML_aqua_tree <- read.tree(file = "lung_loss_git/trees/edited_trees/maxLH_aqu_vis_tree.tre") 
state_labels <- c("X_S", "Lu_S", "X_P", "Lu_P")
trait_col_names <- c("ecology", "lung")
#mode = "eight_state"
#model = "eight_state"

Nruns <- 10
data <- read.csv(file = "lung_loss_git/processed_data/lung_data/full_data.csv")
data <- data[data$Taxa %in% ML_aqua_tree$tip.label,]




source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")

Q_indep <- custom_Q(rate_sum_indep[,2], model = "dependent", scale = .001)
Q_dep <- custom_Q(rate_sum_dep[,2], model = "dependent", scale = .001)

sim_dat <- make_simmap_data(data, ML_aqua_tree, mode = "double binary", 
                            state_labels, trait_col_names)


A <- make.simmap(ML_aqua_tree, sim_dat, nsim = 100,
                 Q = Q_indep, state_labels = state_labels)
B <- make.simmap(ML_aqua_tree, sim_dat, nsim = 100,
                 Q = Q_dep, state_labels = state_labels)

cols1 <- setNames(c("magenta", "#99CCFF", "red", "#B9FFB4"), 
                         c("X_S", "Lu_S", "X_P", "Lu_P"))


tiplabels(tip = tips, pch = 16, cex = .4, col = "black")

{png(file = "figures/supp_figure_2_indep.png", bg = "transparent", units = "in", res = 1200, width = 4, height = 10)
plotTree(A[[1]],fsize=0.1,offset=0.15, #black tree with wide edges to be outlines later
         lwd=1.7)
for(i in 1:4){
  tips <- match(names(which(sim_dat[,i] == 1)), A[[1]]$tip.label)
  tiplabels(tip = tips, pch = 16, cex = .2, offset = 1, col = cols1[i])
}
par(fg="transparent",lend=1)
plotTree(A[[1]],fsize=0.1,offset=0.15,     #white tree to create the outline (smaller lwd)
         lwd=1,color="white",add=TRUE)
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency

for(i in 1:100) {
  plot(A[[i]], colors=sapply(cols1, make.transparent,alpha=0.03),
       add=TRUE, lwd=1, fsize=0.1, offset=0.15)
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}
dev.off()}

{png(file = "figures/supp_figure_2_dep.png", bg = "transparent", units = "in", res = 1000, width = 4, height = 10)
  plotTree(B[[1]],fsize=0.1,offset=0.15, #black tree with wide edges to be outlines later
           lwd=2)
  par(fg="transparent",lend=1)
  plotTree(B[[1]],fsize=0.1,offset=0.15,     #white tree to create the outline (smaller lwd)
           lwd=1,color="white",add=TRUE)
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) {
    plot(B[[i]], colors=sapply(cols1, make.transparent,alpha=0.03),
         add=TRUE, lwd=1, fsize=0.1, offset=0.15)
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  for(i in 1:4){
    tips <- match(names(which(sim_dat[,i] == 1)), A[[1]]$tip.label)
    tiplabels(tip = tips, pch = 16, cex = .2, offset = 1, col = cols1[i])
  }
  dev.off()}

rm(list = ls())


eight_state_res_3 <- read.csv("bayestraits/output/processed_logs/eight_state_res_3.txt", sep = "\t")
source("lung_loss_git/scripts/functions/summary_posterior_function.R")

rate_sum <- summary_posterior(eight_state_res_3, grep("q", colnames(eight_state_res_3)))

ML_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/MaxLh_tree_full.nex") 
state_labels <- eight_state_state <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P", "X_terr", "Lu_terr")
trait_col_names <- eight_state_col <- c("ecology", "lung", "Spec_lotic", "terrestrial")
mode = "eight_state"
model = "eight_state"

data <- read.csv(file = "lung_loss_git/processed_data/lung_data/full_data.csv")
data <- data[data$Taxa %in% ML_tree$tip.label,]

source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")


cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#B9FFB4",
                          "brown", "#CC6600"), 
                          c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                            "X_terr", "Lu_terr"))

Q <- custom_Q(rate_sum[,2], model = "eight_state", scale = .01)
#Q <- custom_Q(rates, model = "eight_state", scale = .01)

sim_dat <- make_simmap_data(data, ML_tree, mode = "eight_state", 
                            eight_state_state, eight_state_col)



tip_1 <- getMRCA(ML_tree, c("Ascaphus_truei", "Leiopelma_hochstetteri"))


ML_tree <- bind.tip(ML_tree, "tip_1" ,edge.length=0,where=tip_1)


sim_dat2 <- rbind(sim_dat, c(0,0,.5,0,0,.5,0,0))
rownames(sim_dat2) <- c(rownames(sim_dat), "tip_1")
sim_dat <- sim_dat2; rm(sim_dat2)

C <- make.simmap(ML_tree, sim_dat, nsim = 100,
                 Q = Q, state_labels = eight_state_state, 
                 pi= c(.25, 0, .25, 0, .25, .25, 0, 0))

C <- lapply(C,drop.tip.simmap,tip="tip_1")
class(C)<-"multiPhylo"


{png(file = "figures/supp_figure_2_eight.png", bg = "transparent", units = "in", res = 1000, width = 4, height = 10)
  plotTree(C[[1]],fsize=0.1,offset=0.15, #black tree with wide edges to be outlines later
           lwd=2)
  par(fg="transparent",lend=1)
  plotTree(C[[1]],fsize=0.1,offset=0.15,     #white tree to create the outline (smaller lwd)
           lwd=1,color="white",add=TRUE)
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) {
    plot(C[[i]], colors=sapply(cols2, make.transparent,alpha=0.03),
         add=TRUE, lwd=1, fsize=0.1, offset=0.15)
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  
  for(i in 1:8){
  tips <- match(names(which(sim_dat[,i] == 1)), C[[1]]$tip.label)
  tiplabels(tip = tips, pch = 16, cex = .2, offset = 1, col = cols2[i])
  }
  dev.off()}








