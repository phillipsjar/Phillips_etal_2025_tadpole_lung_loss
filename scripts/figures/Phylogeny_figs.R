setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2025")



rm(list = ls())
library(ape)
library(phytools)

# read in bayestraits posteriors for select models

eight_bestmod_1 <- read.csv("bayestraits/output/processed_logs/eight_state_bestmod_1.txt", sep = "\t")
eight_state_res_3 <- read.csv("bayestraits/output/processed_logs/eight_state_res_3.txt", sep = "\t")

#load in function that extracts summary statistics from posterior
source("lung_loss_git/scripts/functions/summary_posterior_function.R")

#extract summary statistics for all markov parameters (eg, q12, q13)
rate_sum <- summary_posterior(eight_bestmod_1, grep("q", colnames(eight_bestmod_1)))
rate_sum_full <- summary_posterior(eight_state_res_3, grep("q", colnames(eight_bestmod_1)))

#load in ML tree
ML_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/MaxLh_tree_full.nex") 

#parameters to be used in later functions
state_labels <- eight_state_state <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P", "X_terr", "Lu_terr")
trait_col_names <- eight_state_col <- c("ecology", "lung", "Spec_lotic", "terrestrial")
mode = "eight_state"
model = "eight_state"
Nruns <- 10

#load in tip data and trim to ML tree
data <- read.csv(file = "lung_loss_git/processed_data/lung_data/full_data.csv")
data <- data[data$Taxa %in% ML_tree$tip.label,]

#load in custom functions to create Q matrices and format tip data for simmaps
source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")

#color scheme
cols1 <- setNames(c("red", "purple", "white", "white", "red", "white",
                         "#8B5E3C", "white"), 
                         c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                           "X_terr", "Lu_terr"))
cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#B9FFB4",
                          "#8B5E3C", "#CC6600"), 
                          c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                            "X_terr", "Lu_terr"))

#make Q matrices for simmap simulations
Q <- custom_Q(rate_sum[,2], model = "eight_state", scale = .01)
Q_full  <- custom_Q(rate_sum_full[,1], model = "eight_state", scale = .01)

#format data for simmaping
sim_dat <- make_simmap_data(data, ML_tree, mode = "eight_state", 
                            eight_state_state, eight_state_col)

#define some nodes to be fossilized with taxa that define a certain node
node_1 <- c("Ascaphus_truei", "Rana_sauteri")
node_2 <- c("Heleophryne_regis", "Rana_sauteri")
node_3 <- c("Oophaga_pumilio", "Allobates_nidicola")
node_4 <- c("Leptodactylus_fuscus", "Adenomera_marmorata")
node_5 <- c("Thoropa_miliaris", "Alsodes_monticola")
node_6 <- c("Bufo_bufo", "Melanophryniscus_pachyrhynus")

nodes <- list(node_1, node_2, node_3, node_4, node_5, node_6)
nodes_n <- c("tip_1", "tip_2", "tip_3", "tip_4", "tip_5", "tip_6")

# add these tip data to the simmap data and ML tree
for(i in 1:(length(nodes)-1)){
  n <- getMRCA(ML_tree, nodes[[i]])
  ML_tree <- bind.tip(ML_tree,  nodes_n[i], edge.length=0, where = n)
  sim_dat <- rbind(sim_dat, c(.25,0,.25,0,.25,.25,0,0))
}

# bind the final node, which has a different fossilization
n <- getMRCA(ML_tree, nodes[[6]])
ML_tree <- bind.tip(ML_tree,  nodes_n[6], edge.length=0, where = n)
sim_dat <- rbind(sim_dat, c(.25, .25, 0, 0, .25, 0, .25, 0))
rownames(sim_dat) <- c(rownames(sim_dat[1:(dim(sim_dat)[1]-length(nodes_n)),]), nodes_n)


#make simmaps
sim <- make.simmap(ML_tree, sim_dat, nsim = 100,
                 Q = Q_full, state_labels = eight_state_state, 
                 pi= c(.25, 0, .25, 0, .25, .25, 0, 0))

sim_best <- make.simmap(ML_tree, sim_dat, nsim = 100,
                   Q = Q, state_labels = eight_state_state, 
                   pi= c(.25, 0, .25, 0, .25, .25, 0, 0))

#drop the additional tips for clearer visualizations
sim <- lapply(sim,drop.tip.simmap,tip=nodes_n)
sim_best <- lapply(sim,drop.tip.simmap,tip=nodes_n)
class(sim)<-"multiPhylo"
class(sim_best)<-"multiPhylo"


#fit ML ancestral reconstructions, using the Q matrix estimated from the
# BayesTraits posterior summary
ancr_sim <- ancr(fitMk(ML_tree, sim_dat, fixedQ = Q_full))
ancr_sim_b <- ancr(fitMk(ML_tree, sim_dat, fixedQ = Q))


# this method closely follows a suggestion by Liam Revell in a phytools blog (http://blog.phytools.org/2020/06/mapping-multi-state-discrete-character.html)
#plot black and red tree for figure 2 

{png(file = "figures/figure_2_phylo.png", bg = "transparent", units = "in", res = 1000, width = 4.2, height = 4.2)
plotTree(sim[[1]],ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=2.3, type = "fan")

par(fg="transparent",lend=1)
plotTree(sim[[1]],ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
         lwd=1.3,color="white",add=TRUE, type = "fan")
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency

for(i in 1:100) {
  plot(sim[[i]], colors=sapply(cols1, make.transparent,alpha=0.03),
      add=TRUE, lwd=1.3, ftype="off", fsize=0.5, offset=0.5, type = "fan")
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}

dev.off()}


#trim down tree for figure 3
families <- unique(data$Family)
indicator <- rep(NA, length(families))
names(indicator) <- families

for(i in 1:length(families)){
  if(length(which(data$Family == families[i] & data$lung == 0)) > 0){
    indicator[i] <- 1
  }else{indicator[i] <- 0}
}

fams <- families[which(indicator == 0)]
rm(families)

lungless_data <- data[!(data$Family %in% fams),]


genera <- unique(lungless_data$Genus)
indicator <- rep(NA, length(genera))
names(indicator) <- genera

for(i in 1:length(genera)){
  if(length(which(lungless_data$Genus == genera[i] & lungless_data$lung == 0)) > 0){
    indicator[i] <- 1
  }else{indicator[i] <- 0}
}


lungless_genera <- genera[which(indicator == 1)]
lunged_genera <- genera[which(indicator == 0)]

tips <- rep(NA, (length(fams) + length(lungless_genera) + length(lunged_genera)))
names(tips) <- c(fams, lungless_genera, lunged_genera)


for(i in 1:length(tips)){
  if(names(tips[i]) %in% fams){
    taxa <- data$Taxa[which(data$Family == names(tips[i]))]
    ecology <- data$eight_state[which(data$Family == names(tips[i]))]
    eco_mode <- as.numeric(names(sort(table(ecology), decreasing = T))[1])
    tips[i] <- taxa[which(ecology == eco_mode)][1]}
  if(names(tips[i]) %in% lunged_genera){
    taxa <- data$Taxa[which(data$Genus == names(tips[i]))]
    ecology <- data$eight_state[which(data$Genus == names(tips[i]))]
    eco_mode <- as.numeric(names(sort(table(ecology), decreasing = T))[1])
    tips[i] <- taxa[which(ecology == eco_mode)][1]}
  if(names(tips[i]) %in% lungless_genera){
    taxa <- data$Taxa[which(data$Genus == names(tips[i]) & data$lung == 0)]
    ecology <- data$eight_state[which(data$Genus == names(tips[i]) & data$lung == 0)]
    eco_mode <- as.numeric(names(sort(table(ecology), decreasing = T))[1])
    tips[i] <- taxa[which(ecology == eco_mode)][1]}
}

tips <- tips[!(tips %in% data$Taxa[which(data$Family == "Microhylidae")])]

tips[which(names(tips) == "Crinia")] <- "Crinia_riparia"

tips <- c(tips, "Litoria_subglandulosa", "Rana_sauteri", "Thoropa_miliaris", "Kalophrynus_sinensis",
          "Otophryne_robusta", "Phrynomantis_bifasciatus")


sim_trim <- lapply(sim,keep.tip.simmap,tip=tips)
class(sim_trim)<-"multiPhylo"


# this method closely follows a suggestion by Liam Revell in a phytools blog (http://blog.phytools.org/2020/06/mapping-multi-state-discrete-character.html)

{png(file = "figures/figure_4_phylo.png", units = "in", width = 7.5, height = 2.5, res = 1000, bg = "transparent")
  plotTree(sim_trim[[1]], pts=F, ftype = "off", fsize = .3, direction = "upwards",
           offset=1, lwd=2.75, mar = c(.1,.1,.5,.1))
  par(fg="transparent",lend=1)
  plotTree(sim_trim[[1]], pts=F, ftype = "off", direction = "upwards", offset=1,
           lwd=1.5, mar = c(.1,.1,.5,.1), color="white", add=TRUE)
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) {
    plot(sim_trim[[i]], pts=F, colors=sapply(cols2, make.transparent,alpha=0.03), offset=1,
         direction = "upwards", add=TRUE, lwd=1.5, ftype="off",  mar = c(.1,.1,.5,.1))
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  segments(160,20,160,70)
  points(rep(160,6), seq(20,70,10), pch = 16, cex = .5)
dev.off()}





 ############# Suppl. Figs for eight-state models ######################################
cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#B9FFB4",
                          "#8B5E3C", "#CC6600"), 
                          c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                            "X_terr", "Lu_terr"))

{jpeg(file = "figures/supp_fig_full_8state.jpg", bg = "transparent", units = "in", res = 1000, width = 10, height = 10)
  plotTree(sim[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(sim[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")

  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) { 
    plot(sim[[i]], colors=sapply(cols2, make.transparent,alpha=0.03),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  par(fg="black",lend=1)
  nodelabels("", pie = ancr_sim[1]$ace, piecol = cols2, cex = .3)
  dev.off()}

{jpeg(file = "figures/supp_fig_best_8state.jpg", bg = "transparent", units = "in", res = 1000, width = 10, height = 10)
  plotTree(sim[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(sim[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")
  
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) { 
    plot(sim[[i]], colors=sapply(cols2, make.transparent,alpha=0.03),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  par(fg="black",lend=1)
  nodelabels("", pie = ancr_sim_b[1]$ace, piecol = cols2, cex = .3)  
  dev.off()}
rm(list = ls())









######################################################################################
# supplemental figures for the four-state models
rm(list = ls())

dep_allzero_2 <- read.csv("bayestraits/output/processed_logs/dep_allzero_2.txt", sep = "\t")
indep_1 <- read.csv("bayestraits/output/processed_logs/indep_1.txt", sep = "\t")
dep_1 <- read.csv("bayestraits/output/processed_logs/dep_1.txt", sep = "\t")


ML_tree_4state <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_aqu.nex") 

data_dep <- read.csv(file = "lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")
data_dep <- data_dep[data_dep$Taxa %in% ML_tree_4state$tip.label,]

data_indep <- read.csv(file = "lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")
data_indep <- data_dep[data_dep$Taxa %in% ML_tree_4state$tip.label,]

source("lung_loss_git/scripts/functions/summary_posterior_function.R")

rate_sum_dep <- summary_posterior(dep_1, grep("q", colnames(dep_1)))
rate_sum_best <- summary_posterior(dep_allzero_2, grep("q", colnames(dep_allzero_2)))
rate_sum_indep <- summary_posterior(indep_1, match(c("alpha1", "beta1", "alpha2", "beta2"), 
                                                       colnames(indep_1)))

state_labels <- four_state_state <- c("X_S", "Lu_S", "X_P", "Lu_P")
trait_col_names <- four_state_col <- c("ecology", "lung")
mode = "dependent"
model = "dependent"


source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")
cols <- setNames(c("purple", "blue", "red", "#B9FFB4"), four_state_state)

Q_dep <- custom_Q(rate_sum_dep[,2], model = "dependent", scale = .001)
Q_best <- custom_Q(rate_sum_best[,2], model = "dependent", scale = .001)
Q_indep <- custom_Q(rate_sum_indep[,2], model = "independent", scale = .001)

sim_dat_dep <- make_simmap_data(data_dep, ML_tree_4state, mode = "double binary", 
                                four_state_state, four_state_col)

A <- make.simmap(ML_tree_4state, sim_dat_dep, nsim = 100,
                 Q = Q_dep, state_labels = four_state_state, pi=rep(1/4,4))
B <- make.simmap(ML_tree_4state, sim_dat_dep, nsim = 100,
                 Q = Q_best, state_labels = four_state_state, pi=rep(1/4,4))
C <- make.simmap(ML_tree_4state, sim_dat_dep, nsim = 100,
                 Q = Q_indep, state_labels = four_state_state, pi=rep(1/4,4))


ancr_dep <- ancr(fitMk(ML_tree_4state, sim_dat_dep, fixedQ = Q_dep))
ancr_best <- ancr(fitMk(ML_tree_4state, sim_dat_dep, fixedQ = Q_best))
ancr_indep <- ancr(fitMk(ML_tree_4state, sim_dat_dep, fixedQ = Q_indep))


plotSimmap(A[[2]], cols, pts=F, fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
nodelabels(pie = ancr_dep[1]$ace, piecol = cols, cex = .2)  




plotSimmap(B[[2]], cols, pts=F, fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
plotSimmap(C[[2]], cols, pts=F, fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)


{jpeg(file = "figures/supp_fig_4dep.jpg", bg = "transparent", units = "in", res = 1200, width = 10, height = 10)
  
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
  par(fg="black",lend=1)
  nodelabels(pie = ancr_dep[1]$ace, piecol = cols, cex = .2)  
  dev.off()}


{jpeg(file = "figures/supp_fig_4dep_best.jpg", bg = "transparent", units = "in", res = 1200, width = 10, height = 10)
  
  plotTree(B[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(B[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) {
    plot(B[[i]], colors=sapply(cols, make.transparent,alpha=0.02),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  par(fg="black",lend=1)
  nodelabels(pie = ancr_best[1]$ace, piecol = cols, cex = .2)  
  dev.off()}

{jpeg(file = "figures/supp_fig_4indep.jpg", bg = "transparent", units = "in", res = 1200, width = 10, height = 10)
  
  plotTree(C[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(C[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) {
    plot(C[[i]], colors=sapply(cols, make.transparent,alpha=0.02),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  par(fg="black",lend=1)
  nodelabels(pie = ancr_indep[1]$ace, piecol = cols, cex = .2)  
  
  dev.off()}


cols <- setNames(c("purple", "blue", "red", "#B9FFB4"), four_state_state)
cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#B9FFB4",
                          "brown", "#CC6600"), 
                          c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                            "X_terr", "Lu_terr"))

{png(file = "figures/supp_fig_labels.png", bg = "transparent", units = "in", res = 1200, width = 5, height = 5)
  par(mar = c(0,0,0,0))
  plot(1:25, 1:25, ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  points(rep(2,12), seq(2,24,length.out = 12), bg = c(cols2,cols), col = "black", pch = 21, cex = 2)
  text(rep(3,12), seq(2,24,length.out = 12), 
       labels = c("lotic, lungless (not suctorial or fossorial)", 
                  "lotic, lungless (suctorial/fossorial)",
                  "lotic, lunged (not suctorial or fossorial)", 
                  "lotic, lunged (suctorial/fossorial)",
                  "lentic, lungless", "lentic, lunged",
                  "terrestrial, lungless",
                  "terrestrial, lunged", "lotic, lungless", "lotic, lunged", 
                  "lentic, lungless", "lentic, lunged"), 
       col = "black", adj = c(0,.5), cex = 1)
  
  dev.off()}




