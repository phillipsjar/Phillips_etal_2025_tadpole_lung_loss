data <- read.csv(file = "dep_mcmc_test.txt", sep = "\t")

require(ape)

BT_ML_tree <- read.tree(file = "git/trees/edited_trees/maxLH_BT_vis_tree.tre")
load(file = "git/bayestraits/data/nodes_to_keep.Rdata")      # nodes actually estimated in BayesTraits


nodes <- as.data.frame(matrix(nrow = length(nodes_to_keep), ncol = 8))     # matrix of reconstructions. Each row corresponds to an internal node
# named columns give the estimated probability that node is a given state

colnames(nodes) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)", "lungless", "lunged", "lotic", "lentic")
rownames(nodes) <- nodes_to_keep

node_cols <- rep(NA, length = length(nodes_to_keep)*4)
for(i in 1:length(nodes_to_keep)){
  iter <- 1 + (i-1)*4
node_cols[iter] <-paste("X",nodes_to_keep[i],"node...P.0.0.", sep = "")
node_cols[iter+1] <- paste("X",nodes_to_keep[i],"node...P.0.1.", sep = "")
node_cols[iter+2] <- paste("X",nodes_to_keep[i],"node...P.1.0.", sep = "")
node_cols[iter+3] <- paste("X",nodes_to_keep[i],"node...P.1.1.", sep = "")
}


averages <- sapply(data[,match(node_cols, colnames(data))], "mean")  #average the probability each sampled node is in each given state
#across the entire posterior

for(i in 1:length(nodes_to_keep)){                                        #fill in the nodes matrix with each state and combined states,
  iteration <- 1 + (i-1)*4                                  #such as "lunged" or "lotic"
  nodes[i,1:4] <- averages[iteration : (iteration + 3)]
  nodes[i,5] <- sum(nodes[i,1], nodes[i,3], na.rm = TRUE)
  nodes[i,6] <- 1 - nodes[i,5]
  nodes[i,7] <- sum(nodes[i,1], nodes[i,2], na.rm = TRUE)
  nodes[i,8] <- 1 - nodes[i,7]
}

save(nodes, file = "git/data/ASE_R_formatted.R")

rm(list = c("data", "averages", "iter", "i", "iteration", "node_cols"))

#check and make sure it makes sense
data <- read.csv("git/data/no_endo_lung_data.csv")
treedata <- data[data$Taxa %in% BT_ML_tree$tip.label,]; rm(data)
cols <- c("magenta", "blue", "red", "green")

par(mar = c(1,1,1,1))
plot(BT_ML_tree, show.tip.label = TRUE, cex = .25, 
     tip.color = cols[treedata$state[match(BT_ML_tree$tip.label, treedata$Taxa)]], type = "f")
nodelabels(node = nodes_to_keep, pie = as.matrix(nodes[,1:4]), cex = .3, piecol = c("magenta", "blue", "red", "green"))
par(mar = c(5.1, 4.1, 4.1, 2.1))

rm(list = c("BT_ML_tree", "nodes", "treedata", "cols", "nodes_to_keep"))






