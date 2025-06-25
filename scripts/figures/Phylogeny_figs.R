setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2024")



rm(list = ls())
library(ape)
library(phytools)
help(package="phytools")
help(package="ape")

eight_bestmod_1 <- read.csv("bayestraits/output/processed_logs/eight_state_bestmod_1.txt", sep = "\t")
eight_state_res_3 <- read.csv("bayestraits/output/processed_logs/eight_state_res_3.txt", sep = "\t")
eight_state_full_2 <- read.csv("bayestraits/output/processed_logs/eight_state_full_2.txt", sep = "\t")

source("lung_loss_git/scripts/functions/summary_posterior_function.R")


rate_sum <- summary_posterior(eight_bestmod_1, grep("q", colnames(eight_bestmod_1)))
rate_sum_full <- summary_posterior(eight_state_full_2, grep("q", colnames(eight_bestmod_1)))

ML_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/MaxLh_tree_full.nex") 
state_labels <- eight_state_state <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P", "X_terr", "Lu_terr")
trait_col_names <- eight_state_col <- c("ecology", "lung", "Spec_lotic", "terrestrial")
mode = "eight_state"
model = "eight_state"

Nruns <- 10
data <- read.csv(file = "lung_loss_git/processed_data/lung_data/full_data.csv")
data <- data[data$Taxa %in% ML_tree$tip.label,]

source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")

cols1 <- setNames(c("gold", "purple", "white", "white", "red", "white",
                         "#8B5E3C", "white"), 
                         c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                           "X_terr", "Lu_terr"))
cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#B9FFB4",
                          "#8B5E3C", "#CC6600"), 
                          c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                            "X_terr", "Lu_terr"))

Q <- custom_Q(rate_sum[,2], model = "eight_state", scale = .01)
Q_full  <- custom_Q(rate_sum_full[,2], model = "eight_state", scale = .01)

sim_dat <- make_simmap_data(data, ML_tree, mode = "eight_state", 
                            eight_state_state, eight_state_col)

{
tip_1 <- getMRCA(ML_tree, c("Ascaphus_truei", "Rana_sauteri"))
tip_2 <- getMRCA(ML_tree, c("Heleophryne_regis", "Rana_sauteri"))
tip_3 <- getMRCA(ML_tree, c("Oophaga_pumilio", "Allobates_nidicola"))
tip_4 <- getMRCA(ML_tree, c("Leptodactylus_fuscus", "Adenomera_marmorata"))
tip_5 <- getMRCA(ML_tree, c("Thoropa_miliaris", "Alsodes_monticola"))
}

tips_n <- c("tip_1", "tip_2", "tip_3", "tip_4", "tip_5")
tips <- c(tip_1, tip_2, tip_3, tip_4, tip_5)


for(i in 1:length(tips)){
ML_tree <- bind.tip(ML_tree,  tips_n[i], edge.length=0, where = tips[i])
sim_dat <- rbind(sim_dat, c(.25,0,.25,0,.25,.25,0,0))
}

rownames(sim_dat) <- c(rownames(sim_dat[1:(dim(sim_dat)[1]-length(tips)),]), tips_n)
sim_dat[500:520,]

A <- make.simmap(ML_tree, sim_dat, nsim = 100,
                 Q = Q, state_labels = eight_state_state, 
                 pi= c(.25, 0, .25, 0, .25, .25, 0, 0))
B <- make.simmap(ML_tree, sim_dat, nsim = 100,
                 Q = Q_full, state_labels = eight_state_state, 
                 pi= c(.25, 0, .25, 0, .25, .25, 0, 0))

A <- lapply(A,drop.tip.simmap,tip=tips_n)
B <- lapply(B,drop.tip.simmap,tip=tips_n)

class(A)<-"multiPhylo"
class(B)<-"multiPhylo"

ancr_A <- ancr(fitMk(ML_tree, sim_dat, fixedQ = Q))
ancr_B <- ancr(fitMk(ML_tree, sim_dat, fixedQ = Q_full))
?ancr
node_1 <- getMRCA(A[[7]], c("Bufo_bufo", "Melanophryniscus_pachyrhynus"))
node_2 <- getMRCA(A[[7]], c("Crinia_riparia", "Taudactylus_diurnus"))

nodes <- match(c(node_1, node_2), rownames(ancr_A[1]$ace))


#test out whether it worked and makes sense

plot(A[[1]], colors=cols2, ftype = "off", )


plot(A[[7]], colors=cols2, lwd=1.3, ftype="off", fsize=0.5, offset=0.5, type = "fan")



nodelabels("", pie = ancr_A[1]$ace, piecol = cols2, cex = .3)




plot(B[[7]], colors=cols1, lwd=1.3, ftype="off", fsize=0.5, offset=0.5, type = "fan")

# this method closely follows a suggestion by Liam Revell in a phytools blog (http://blog.phytools.org/2020/06/mapping-multi-state-discrete-character.html)


#plot black and red tree for figure 2 

{png(file = "figures/figure_2_phylo.png", bg = "transparent", units = "in", res = 1000, width = 4.2, height = 4.2)
plotTree(A[[1]],ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=2.3, type = "fan")

par(fg="transparent",lend=1)
plotTree(A[[1]],ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
         lwd=1.3,color="white",add=TRUE, type = "fan")
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency

for(i in 1:100) {
  plot(A[[i]], colors=sapply(cols1, make.transparent,alpha=0.03),
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



A_trim <- lapply(A,keep.tip.simmap,tip=tips)
class(A_trim)<-"multiPhylo"

{png(file = "figures/figure_4_phylo.png", width = 7.5, height = 3, 
      units = "in", res = 1000, bg = "transparent")
par(xpd = TRUE)

plotSimmap(A_trim[[2]], cols2, pts=F, ftype = "off", fsize = .3, direction = "upwards",
           outline = TRUE, add=FALSE, lwd=1.75, offset=1, mar = c(.1,.1,.5,.1))
segments(160,20,160,70)
points(rep(160,6), seq(20,70,10), pch = 16, cex = .5)
dev.off()}

max(ML_tree$edge.length)

?plotSimmap

# this method closely follows a suggestion by Liam Revell in a phytools blog (http://blog.phytools.org/2020/06/mapping-multi-state-discrete-character.html)

{png(file = "figures/figure_4_phylo.png", units = "in", width = 7.5, height = 2.5, res = 1000, bg = "transparent")
  plotTree(A_trim[[1]], pts=F, ftype = "off", fsize = .3, direction = "upwards",
           offset=1, lwd=2.75, mar = c(.1,.1,.5,.1))
  par(fg="transparent",lend=1)
  plotTree(A_trim[[1]], pts=F, ftype = "off", direction = "upwards", offset=1,
           lwd=1.5, mar = c(.1,.1,.5,.1), color="white", add=TRUE)
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) {
    plot(A_trim[[i]], pts=F, colors=sapply(cols2, make.transparent,alpha=0.03), offset=1,
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

{jpeg(file = "figures/supp_fig_best_8state.jpg", bg = "transparent", units = "in", res = 1000, width = 10, height = 10)
  plotTree(A[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(A[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")

  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) { 
    plot(A[[i]], colors=sapply(cols2, make.transparent,alpha=0.03),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  par(fg="black",lend=1)
  nodelabels("", pie = ancr_A[1]$ace, piecol = cols2, cex = .3)
  dev.off()}

{jpeg(file = "figures/supp_fig_full_8state.jpg", bg = "transparent", units = "in", res = 1000, width = 10, height = 10)
  plotTree(B[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(B[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")
  
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) { 
    plot(B[[i]], colors=sapply(cols2, make.transparent,alpha=0.03),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  par(fg="black",lend=1)
  nodelabels("", pie = ancr_B[1]$ace, piecol = cols2, cex = .3)  
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

############# Presentation Figures ##################

eco_cols <- setNames(c("blue", "blue", "green", "green"), four_state_state)
lung_cols <- setNames(c("red", "white", "red", "white"), four_state_state)



{png(file = "figures/eco_simm_presentation.png", bg = "transparent", units = "in", res = 2000, width = 5, height = 12)
par(fg="transparent",lend=1)
plotSimmap(B[[1]], color="white", pts=F, fsize = .2, ftype="off",
           outline = F, type = "phylogram", add=FALSE, lwd=1.4, offset=0.5)

for(i in 1:100) {
  plot(B[[i]], colors=sapply(eco_cols, make.transparent,alpha=0.02),
       add=TRUE, lwd=1.4, ftype="off", fsize=0.35, offset=.5, type = "phylogram")
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}
dev.off()
}


{png(file = "figures/lung_simm_presentation.png", bg = "transparent", units = "in", res = 2000, width = 5, height = 12)
  par(fg="transparent",lend=1)
  plotSimmap(B[[1]], color="white", pts=F, fsize = .2, ftype="off",
             outline = F, type = "phylogram", add=FALSE, lwd=1.4, offset=0.5)
  
  for(i in 1:100) {
    plot(B[[i]], colors=sapply(lung_cols, make.transparent,alpha=0.02),
         add=TRUE, lwd=1.4, ftype="off", fsize=0.35, offset=.5, type = "phylogram")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  dev.off()
}


 plotSimmap(B[[1]], eco_cols, pts=F, fsize = .2, ftype="off",
           outline = F, type = "phylogram", add=TRUE, lwd=2, offset=0.5)




plotSimmap(B[[1]], eco_cols, pts=F, fsize = .2, ftype="off",
           outline = F, type = "phylogram", add=FALSE, lwd=2, offset=0.5)


?plotSimmap

plot(B[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
         lwd=5, type = "fan")





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






######################################################
# phylo code for Moey Rojas by Jack Phillips 5-12-2024

library(ape)
library(phytools)

rm(list = ls())


# load in your tree (any file type is fine - .tre, .nex, .phy)
ML_tree_4state <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_aqu.nex") 

# load in your dataset, which needs a column with taxon names that match the tree 
# and a column with 0s and 1s for a trait of interest
data_dep <- read.csv(file = "lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")
data_dep <- data_dep[data_dep$Taxa %in% ML_tree_4state$tip.label,]
# (second line will filter out any taxa in the dataset not in your tree)


eco <-  cbind(rep(0, dim(data_dep)[1]), rep(0, dim(data_dep)[1]))
# this will be your simulation dataset, so I'm just matching dimensions

rownames(eco) <- data_dep$Taxa
# this is key, it won't work if the taxon names don't match the tree
#to check, you can do some of this kinda thing:
table(data_dep$Taxa %in% ML_tree_4state$tip.label) 
table(ML_tree_4state$tip.label %in% data_dep$Taxa) 
#should be only trues, no falses.


colnames(eco) <- c("lotic", "lentic") #change to match what your states are

#populate your simulation dataset
eco[which(data_dep$ecology == 0), 1] <- 1
eco[which(data_dep$ecology == 1), 2] <- 1
eco[which(is.na(data_dep$ecology)), c(1,2)] <- c(.5, .5)
# last line allows for unknown tip states, if they are NAs in your dataset. 
# we are essentially just giving them a flat prior and letting the data decide.



#lungs <-  cbind(rep(0, dim(data_dep)[1]), rep(0, dim(data_dep)[1]))
#rownames(lungs) <- data_dep$Taxa
#colnames(lungs) <- c("lungless", "lunged")
#lungs[which(data_dep$lung == 0), 1] <- 1
#lungs[which(data_dep$lung == 1), 2] <- 1
#lungs[which(is.na(data_dep$lung)), c(1,2)] <- c(.5, .5)

head(eco) 
#head(lungs) 

# here you are actually doing the simulation, with an "all rates different (ARD)" model
# you can also try the "SYM" model, and in a paper you would want to compare them, 
# but looking at your data, the SYM model probably won't even run it is so bad.
# I usually do 100 simulations, but you can pick. Remember to adjust your state names again here.
eco_simmap <- make.simmap(ML_tree_4state, eco, model = "ARD", nsim = 100, state_labels = c("lotic", "lentic"))

# set up some colors for plotting, change the state names to match what you want to see
cols <- setNames(c("blue", "green"), c("lotic", "lentic"))

# do a tester on a single simulation run:
plot(eco_simmap[[1]], cols, ftype = "off")

#you can turn on tip labels like this: (on big trees I turn them off, but yours might be small enough)
plot(eco_simmap[[1]], cols, ftype = "i")


# to build an actual figure, you can run this block, which will produce a graded color figure that 
# smooths out a lot of the randomness. I adjusted it to only plot 10 layers, but feel free to increase
# that resolution by adding more iterations to the for-loop, just make sure you also lower the 
# transparency (alpha) value when you do. 100 iterations works well with an alpha of .02


{#first go through starting on the line starting with "plotTree.." to make sure it looks good,
  # then run the whole thing from the bracket to make sure it works. If it is being weird, 
  # run dev.off() until it says "null device 1", then run the whole block from the bracket
  # and it will pop up in your directory.
jpeg(file = "phylo_plot.jpg", bg = "transparent", units = "in", res = 1200, width = 4, height = 10)
  
plotTree(eco_simmap[[1]],ftype="off",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
         lwd=3.5, type = "phylogram")

par(fg="transparent",lend=1)
plotTree(eco_simmap[[1]],ftype="off",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
         lwd=1.5,color="white",add=TRUE, type = "phylogram")

## now plot our 100 stochastic map trees pulled from simmap object (same lwd as white tree)
## with 99% transparency
for(i in 1:10) { #here's where you increase iterations
  plot(eco_simmap[[i]], colors=sapply(cols, make.transparent,alpha=0.2), #alpha on this line
       add=TRUE, lwd=1.5, ftype="off", fsize=0.35, offset=3, type = "phylogram")
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}
par(fg="black",lend=1)
nodelabels(pie = ancr_best[1]$ace, piecol = cols, cex = .2)  
dev.off() 

}


# sorry if anything is over-explained, I just always default to assuming people don't 
# know code stuff :)







plot(eco_simmap[[2]], cols)








lung_simmap <- make.simmap(ML_tree_4state, lungs, nsim = 10, state_labels = c("lungless", "lunged"))







sim_dat_dep <- make_simmap_data(data_dep, ML_tree_4state, mode = "double binary", 
                                four_state_state, four_state_col)


eco

data_dep$ecology







names(eco) <- data_dep$Taxa

A <- make.simmap(ML_tree_4state, eco, nsim = 100, pi=rep(1/4,4))


#        , state_labels = c("lotic", "lentic")



state_labels <- four_state_state <- c("X_S", "Lu_S", "X_P", "Lu_P")
trait_col_names <- four_state_col <- c("ecology", "lung")



which(data_dep$ecology == 0)
      
      

source("lung_loss_git/scripts/functions/make_simmap_data_function.R")
sim_dat_dep <- make_simmap_data(data_dep, ML_tree_4state, mode = "double binary", 
                                four_state_state, four_state_col)



source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")
cols <- setNames(c("purple", "blue", "red", "#B9FFB4"), four_state_state)

Q <- custom_Q(rate_sum_dep[,2], model = "dependent", scale = .001)

source("lung_loss_git/scripts/functions/make_simmap_data_function.R")
sim_dat_dep <- make_simmap_data(data_dep, ML_tree_4state, mode = "double binary", 
                                four_state_state, four_state_col)

A <- make.simmap(ML_tree_4state, sim_dat_dep, nsim = 100,
                 Q = Q, state_labels = four_state_state, pi=rep(1/4,4))




A <- make.simmap()

data <- read.csv(file = "lung_loss_git/processed_data/lung_data/full_data.csv")
data <- data[data$Taxa %in% ML_tree$tip.label,]





























 
A <- BayesTraits_simmap(eight_bestmod_1, 100, data, model = "eight_state",
                        trees = trees, trait_col_names = eight_state_col, 
                        state_labels = eight_state_state, pis = c(.25, 0, .25,
                                                                  0, .25, .25, 
                                                                  0, 0))



Q <- custom_Q(rate_sum[,2], model = "eight_state", scale = .01)
#Q <- custom_Q(rates, model = "eight_state", scale = .01)

sim_dat <- make_simmap_data(data, ML_tree, mode = "eight_state", 
                            eight_state_state, eight_state_col)


tip <- getMRCA(ML_tree, c("Ascaphus_truei", "Leiopelma_hochstetteri"))
ML_tree <- bind.tip(ML_tree, "Leio_anc" ,edge.length=0,where=tip)
sim_dat2 <- rbind(sim_dat, c(0,0,.5,0,0,.5,0,0))
rownames(sim_dat2) <- c(rownames(sim_dat), "Leio_anc")
sim_dat <- sim_dat2; rm(sim_dat2)

A <- make.simmap(ML_tree, sim_dat, nsim = 100,
                 Q = Q, state_labels = eight_state_state, 
                 pi= c(.25, 0, .25, 0, .25, .25, 0, 0))

A <- lapply(A,drop.tip.simmap,tip="Leio_anc")
class(A)<-"multiPhylo"

