#######################################################################################################
###############################    Ancestral Reconstructions   ########################################

#######################################################################################################
# load in data first
# load in bayestraits run with ancestral state estimations (ASE)

dep_ase_4062 <- read.csv("Bayestraits/results/dep_ASE.txt",sep = "\t")
dep_ase_4062 <- dep_ase_4062[round(dim(dep_ase_4062)[1]*.25) : dim(dep_ase_4062)[1],]    # and burnin

# and master run (3 converged runs without ancestral estimations)

load("Bayestraits/results/maxLH_dep_master_export.Rdata")

# check convergence (that the ASE run is representative of other runs)

{require(vioplot)
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                 
labels <- c("q12 (gain lungs in S)", "q13 (S to P no lungs)", "q21 (lose lungs in S)",
              "q24 (S to P w/ lungs)", "q31 (P to S no lungs)", "q34 (gain lungs in P)",
              "q42 (P to S w/ lungs)", "q43 (lose lungs in P)")
par(mfrow=c(1,2))
par(xpd = TRUE)
par(mar = c(2.1,4.1,2.1,0.1))
vioplot(dep_ase_4062[,8:15], lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
        colMed = "transparent", col = colorBlindGrey8, main = "ASE Rate Values",
        ylab = "Rate value")
vioplot(master_export_maxLH_dep[,8:15], lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
        colMed = "transparent", col = colorBlindGrey8, main = "Master Rate Values",
        ylab = "")}

rm(list = c("colorBlindGrey8", "labels", "i", "iteration"))

#######################################################################################################

# formatting ASE from bayestraits into a usable format

library(ape)

BT_ML_tree <- read.tree(file = "Trees/edited_trees/maxLH_BT_vis_tree.tre")
load(file = "bayestraits/data/nodes_to_keep.Rdata")      # save this for later plotting, since it takes a while to calculate


nodes <- as.data.frame(matrix(nrow = length(nodes_to_keep), ncol = 8))     # matrix of reconstructions. Each row corresponds to an internal node
                                                             # named columns give the estimated probability that node is a given state

colnames(nodes) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)", "lungless", "lunged", "lotic", "lentic")
rownames(nodes) <- nodes_to_keep

averages <- sapply(dep_ase_4062[,20:dim(dep_ase_4062)[2]], "mean")  #average the probability each node is in each given state
                                                          #across the entire posterior

for(i in 1:N_nodes){                                        #fill in the nodes matrix with each state and combined states,
  iteration <- 1 + (i-1)*4                                  #such as "lunged" or "lotic"
  nodes[i,1:4] <- averages[iteration : (iteration + 3)]
  nodes[i,5] <- sum(nodes[i,1], nodes[i,3], na.rm = TRUE)
  nodes[i,6] <- 1 - nodes[i,5]
  nodes[i,7] <- sum(nodes[i,1], nodes[i,2], na.rm = TRUE)
  nodes[i,8] <- 1 - nodes[i,7]
}

save(nodes, file = "data/ASE_R_formatted.R")

#############################################################################################
# tip data for visualization
library(geiger)

data_no_endo <- read.csv(file = "data/no_endo_lung_data.csv")

data_ASE <- data_no_endo[data_no_endo$Taxa %in% BT_ML_tree$tip.label,]
rownames(data_ASE) <- data_4062$Taxa
data_4062$state <- rep(NA)

for(i in 1:length(data_4062$state)){
  if(is.na(data_4062$ecology[i])){
    data_4062$state[i] <- NA}
  if(!is.na(data_4062$ecology[i])){
    if(data_4062$ecology[i] == 0 & data_4062$lung[i] == 0){
      data_4062$state[i] <- 1}
    if(data_4062$ecology[i] == 0 & data_4062$lung[i] == 1){
      data_4062$state[i] <- 2}
    if(data_4062$ecology[i] == 1 & data_4062$lung[i] == 0){
      data_4062$state[i] <- 3}
    if(data_4062$ecology[i] == 1 & data_4062$lung[i] == 1){
      data_4062$state[i] <- 4}}}

treeWData <- treedata(BT_4062_tree, data_4062, sort = T)
state <- as.factor(treeWData$data[,dim(treeWData$data)[2]])        #This automatically created a factor w/ four levels

######################################################################################################################
##############                  plotting BT reconstructions

par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plot(treeWData$phy, show.tip.label = FALSE, type = "f")
#nodelabels(node = 512, text = " ", bg = "#FFFF0050", cex = 2, frame = "circle")
#nodelabels(node = 385, text = " ", bg = "#FFFF0050", cex = 1, frame = "circle")
reconCol <- c("magenta", "blue", "red", "green")
nodelabels(node = ((tips+1):(2*tips-1)), pie = as.matrix(nodes[,1:4]), cex = .45, piecol = c("magenta", "blue", "red", "green"))
tiplabels(text = rep("", length(treeWData$phy$tip.label)), tip = c(1:tips), cex = .2, 
          frame = "circle", bg = reconCol[state]) #Plot tips


#######################################################################################################
######      Stochastic Character Mapping in Phytools using Bayestraits Rate Matrices      #############

# bayestraits does not have a native stochastic mapping tool, so we will use phytools to do the mapping,
# using the BT posterior to generate Q_matrices while incorporating uncertainty

#              Draw rate matrices from 100 different runs of posterior
load("Bayestraits/results/maxLH_dep_master_export.Rdata")
data_no_endo <- read.csv(file = "data/no_endo_lung_data.csv")

set.seed(123)

random_runs <- sample(1:(dim(master_export_maxLH_dep)[1]), 100, replace = FALSE)  # get a sample of 100 random runs from the posterior

rate_values_max_lh <- vector(mode = "list", length = 100)             # empty list for the rate values for each random run

for(i in 0:99){
  start <- which(colnames(master_export_maxLH_dep) == "q12")
  run <- random_runs[i+1]
  rate_values_max_lh[i+1] <- list(master_export_maxLH_dep[(run),start:(start+7)])
}

# phytools does not have a native dependent model, but this can easily be dealt with by treating the 
# dependent model as a constrained 4 state model, with each combination of the two binary states treated
# as one of four states (lungless stream; lunged, stream; lungless, pond; lunged, pond)
rm(i, run, start, N_nodes, tips, iteration, averages)

source("scripts/extra_scripts/custom_Q.R")

custom_Q(rate_values_max_lh[1])  #example Q matrix from one run of the posterior

############### make simmap objects ####################
library(phytools)

##### set up simmap for pseudo-dependent analysis

# data_4062 loaded above.
# simmap requires a specific way data are input, which we create below with a quick function

source("scripts/extra_scripts/make_simmap_data_function.R")

lung_eco_matrix_max_LH <- make_simmap_data(data_no_endo, BT_4062_tree)
lung_eco_matrix_portik <- make_simmap_data(data_no_endo, BT_portik_tree)

##########################################################################
############ states for plotting later

{states <- matrix(nrow = dim(lung_eco_matrix_max_LH)[1], ncol = 2)
states[,1]  <- BT_4062_tree$tip.label
for(i in 1:dim(lung_eco_matrix_max_LH)[1]){
  name <- states[i,1]
  if(length(which(lung_eco_matrix_max_LH[which(names(lung_eco_matrix_max_LH[,1]) == name),] == 1))>0){
    states[i,2] <- which(lung_eco_matrix_max_LH[which(names(lung_eco_matrix_max_LH[,1]) == name),] == 1)
  }else{states[i,2] <- 5}
}
states[,2] <- as.numeric(states[,2])}
rm(i, name)

###############################
# make a dummy simmap with 100 sims to fill later (doesn't matter what the values are) (may take a couple minutes)
master_simmap <- make.simmap(BT_4062_tree,lung_eco_matrix_max_LH, nsim=100, Q = custom_Q(rate_values_max_lh[1]), pi=rep(.25,4))

#now fill each index of the simmap object with different Q matrices
for(i in 1:100){
  master_simmap[[i]] <- make.simmap(BT_4062_tree,lung_eco_matrix_max_LH, nsim=1, Q = custom_Q(rate_values_max_lh[i]), pi=rep(.25,4))
  print(100-i) # count down for sanity
}

save(master_simmap, file = "data/4062_simmap_100.Rdata")

####################### plotting ###############################

cols2 <- setNames(c("#990099", "#3399FF", "#FF3333","#EBFFE1", "grey"), c("X_S", "Lu_S", "X_P", "Lu_P"))


#tester of one run to see if colors look good
plotSimmap(master_simmap[[2]], cols2, pts=F, ftype="off", fsize = .2, type = "fan", add=TRUE, lwd=2, offset=0.5)

# this method closely follows a suggestion by Liam Revell in a phytools blog (http://blog.phytools.org/2020/06/mapping-multi-state-discrete-character.html)

plotTree(BT_4062_tree,ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=5, type = "fan")

par(fg="transparent",lend=1)
plotTree(BT_4062_tree,ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
         lwd=2.5,color="white",add=TRUE, type = "fan")
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency

for(i in 1:100) {plot(master_simmap[[i]],
                      colors=sapply(cols2, make.transparent,alpha=0.01),
                      add=TRUE, lwd=2, ftype="off", fsize=0.5, offset=0.5, type = "fan")
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}
par(fg="black")

legend(x="bottomleft",c("lungless in streams", "lunged in streams", "lungless in ponds", "lunged in ponds"),pch=22,
       pt.bg=cols2,pt.cex=1.5,bty="n",cex=0.7)
nodelabels(node = ((tips+1):(2*tips-1)), pie = as.matrix(nodes[,1:4]), cex = .45, piecol = c("magenta", "blue", "red", "green"))
tiplabels(text = rep("", length(treeWData$phy$tip.label)), tip = c(1:tips), cex = .2, 
          frame = "circle", bg = reconCol[state]) #Plot tips

#######################################################################################################
########### try with Portik tree

BT_portik_tree <- read.nexus(file = "Bayestraits/tree/portik_tree.nex")

master_simmap <- make.simmap(BT_portik_tree,lung_eco_matrix_portik, nsim=100, Q = custom_Q(rate_values[1]), pi=c(.25,.25,.25,25))

#now fill each index of the simmap object with different Q matrices
for(i in 1:100){
  master_simmap[[i]] <- make.simmap(BT_portik_tree,lung_eco_matrix_portik, nsim=1, Q = custom_Q(rate_values[i]), pi=c(.25,.25,.25,25))
  print(100-i) # count down for sanity
}

plotTree(BT_portik_tree,ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=5, type = "fan")

par(fg="transparent",lend=1)
plotTree(BT_portik_tree,ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
         lwd=2.5,color="white",add=TRUE, type = "fan")
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency

for(i in 1:100) {plot(master_simmap[[i]],
                      colors=sapply(cols2, make.transparent,alpha=0.01),
                      add=TRUE, lwd=2, ftype="off", fsize=0.5, offset=0.5, type = "fan")
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}
par(fg="black")

legend(x="bottomleft",c("lungless in streams", "lunged in streams", "lungless in ponds", "lunged in ponds"),pch=22,
       pt.bg=cols2,pt.cex=1.5,bty="n",cex=0.7)
nodelabels(node = ((tips+1):(2*tips-1)), pie = as.matrix(nodes[,1:4]), cex = .45, piecol = c("magenta", "blue", "red", "green"))
tiplabels(text = rep("", length(treeWData$phy$tip.label)), tip = c(1:tips), cex = .2, 
          frame = "circle", bg = reconCol[state]) #Plot tips












