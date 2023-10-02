#######################################################################################################
###############################    Ancestral Reconstructions   ########################################

#######################################################################################################
require(ape)
require(phytools)
data <- read.csv("git/data/no_endo_lung_data.csv")                              # load in tip data
BT_ML_tree <- read.tree(file = "git/trees/edited_trees/maxLH_BT_vis_tree.tre")  # load in ML_tree
treedata <- data[data$Taxa %in% BT_ML_tree$tip.label,]; rm(data)                # trim tip data to tree and remove extra file

cols <- c("magenta", "blue", "red", "green")

load(file = "git/data/ASE_R_formatted.R")                                       # load in ASE for reconstructed nodes (from setup/BT_data_set)

######################################################################################################################
##############                  plotting BT reconstructions

par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plot(BT_ML_tree, show.tip.label = TRUE, cex = .2, tip.color = cols[treedata$state[match(BT_ML_tree$tip.label, treedata$Taxa)]], type = "f")
nodelabels(node = as.numeric(rownames(nodes)), pie = as.matrix(nodes[,1:4]), cex = .2, piecol = c("magenta", "blue", "red", "green"))
tiplabels(cex = .5, pch = 21, bg = cols[treedata$state[match(BT_ML_tree$tip.label, treedata$Taxa)]])

#######################################################################################################
######      Stochastic Character Mapping in Phytools using Bayestraits Rate Matrices      #############

# bayestraits does not have a native stochastic mapping tool, so we will use phytools to do the mapping,
# using the BT posterior to generate Q_matrices while incorporating uncertainty

#              Draw rate matrices from 100 different runs of posterior
data <- read.csv(file = "dep_mcmc_test.txt", sep = "\t")
data_no_endo <- read.csv(file = "git/data/no_endo_lung_data.csv")

set.seed(123)

random_runs <- sample(1:(dim(data)[1]), 100, replace = FALSE)  # get a sample of 100 random runs from the posterior

rate_values <- vector(mode = "list", length = 100)             # empty list for the rate values for each random run

for(i in 0:99){
  start <- which(colnames(data) == "q12")
  run <- random_runs[i+1]
  rate_values[i+1] <- list(data[(run),start:(start+7)])
}

# phytools does not have a native dependent model, but this can easily be dealt with by treating the 
# dependent model as a constrained 4 state model, with each combination of the two binary states treated
# as one of four states (lungless stream; lunged, stream; lungless, pond; lunged, pond)
rm(i, run, start)

source("git/scripts/extra_scripts/custom_Q.R")

custom_Q(rate_values[1])  #example Q matrix from one run of the posterior

############### make simmap objects ####################

##### set up simmap for pseudo-dependent analysis

# data_4062 loaded above.
# simmap requires a specific way data are input, which we create below with a quick function

source("git/scripts/extra_scripts/make_simmap_data_function.R")

lung_eco_matrix <- make_simmap_data(data_no_endo, BT_ML_tree)

###############################
# make a dummy simmap with 100 sims to fill later (doesn't matter what the values are) (may take a couple minutes)
master_simmap <- make.simmap(BT_ML_tree,lung_eco_matrix, nsim=100, Q = custom_Q(rate_values[1]), pi=rep(.25,4))

#now fill each index of the simmap object with different Q matrices
for(i in 1:100){
  master_simmap[[i]] <- make.simmap(BT_ML_tree,lung_eco_matrix, nsim=1, Q = custom_Q(rate_values[i]), pi=rep(.25,4))
  print(100-i) # count down for sanity
}

save(master_simmap, file = "git/data/simmap_100.Rdata")

####################### plotting ###############################

cols2 <- setNames(c("#990099", "#3399FF", "#FF3333","#EBFFE1", "grey"), c("X_S", "Lu_S", "X_P", "Lu_P"))


#tester of one run to see if colors look good
plotSimmap(master_simmap[[1]], cols2, pts=F, ftype="off", fsize = .2, 
           outlin = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)


# this method closely follows a suggestion by Liam Revell in a phytools blog (http://blog.phytools.org/2020/06/mapping-multi-state-discrete-character.html)

plotTree(BT_ML_tree,ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=5, type = "fan")

par(fg="transparent",lend=1)
plotTree(BT_ML_tree,ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
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
nodelabels(node = as.numeric(rownames(nodes)), pie = as.matrix(nodes[,1:4]), cex = .2, piecol = c("#990099", "#3399FF", "#FF3333","#EBFFE1"))
#cols <- c("#990099", "#3399FF", "#FF3333","#EBFFE1")
#tiplabels(cex = .5, pch = 21, bg = cols[treedata$state[match(BT_ML_tree$tip.label, treedata$Taxa)]])




#
trees <- read.nexus(file = "git/bayestraits/trees/tree_set.nex")

master_simmap_ind <- make.simmap(BT_ML_tree,lung_eco_matrix, nsim=100, Q = custom_Q(rate_values[1]), pi=rep(.25,4))


#now fill each index of the simmap object with different Q matrices
for(i in 1:100){
  master_simmap_ind[[i]] <- make.simmap(trees[[data$Tree.No[random_runs[i]]]],lung_eco_matrix, nsim=1, Q = custom_Q(rate_values[i]), pi=rep(.25,4))
  print(100-i) # count down for sanity
}

genera <- unique(treedata$Genus)
genus_tips <- rep(NA, length(genera))

for(i in 1:length(genera)){
  A <- treedata[which(treedata$Genus == genera[i]),]
  genus_tips[i] <- A$Taxa[sample(dim(A)[1], 1)]
}

simmap_trim <- lapply(master_simmap_ind, keep.tip.simmap, tip = genus_tips)
class(simmap_trim)<-"multiPhylo"

simmap_trim
cols2 <- setNames(c("#990099", "#3399FF", "#FF3333","#EBFFE1", "grey"), c("X_S", "Lu_S", "X_P", "Lu_P"))
plotSimmap(simmap_trim[[1]], cols2, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)


par(mfrow = c(3,3))
plotSimmap(sample(simmap_trim,9), cols2, pts=F, ftype="off", fsize = .2, 
           outlin = TRUE, type = "fan", add=FALSE, lwd=1, offset=0.5)

par(mfrow = c(1,1))



# plotSimmap(master_simmap_ind[[1]])





















