# script currently non-functional


require(stringr)
require(ape)
require(phytools)

data <- read.csv(file = "dep_mcmc_test.txt", sep = "\t")
trees <- read.nexus(file = "git/bayestraits/trees/tree_set.nex")

run <- sample(dim(data)[1], 1)
run <- data[run,]
tree <- trees[[run$Tree.No]]

rm(list = c("data", "trees"))
load(file = "git/bayestraits/data/nodes_to_keep.Rdata")      # nodes actually estimated in BayesTraits
data <- read.csv("git/data/no_endo_lung_data.csv")
treedata <- data[data$Taxa %in% tree$tip.label,]; rm(data)
cols <- c("magenta", "blue", "red", "green")

plot(tree, show.tip.label = TRUE, cex = .25, 
     tip.color = cols[treedata$state[match(tree$tip.label, treedata$Taxa)]], type = "f")
nodelabels(node = nodes_to_keep, frame = "circle", cex = .2)

rm(data)


treedata$state[match(tree$tip.label, treedata$Taxa)]

par(mar = c(1,1,1,1))
plot(tree, show.tip.label = TRUE, cex = .25, 
     tip.color = cols[treedata$state[match(tree$tip.label, treedata$Taxa)]], type = "f")


run_nodes <- run[grep("node...P.....$", names(run))]


recon_nodes <- gsub("X","", names(run_nodes))
recon_nodes <- gsub("node...P.....$","", recon_nodes)
nodes <- unique(recon_nodes) #nodes reconstructed in BayesTraits


node_cols <- rep(NA, length = length(nodes)*4)
for(i in 1:length(nodes)){
  iter <- 1 + (i-1)*4
  node_cols[iter]   <- paste("X",nodes[i],"node...P.0.0.", sep = "")
  node_cols[iter+1] <- paste("X",nodes[i],"node...P.0.1.", sep = "")
  node_cols[iter+2] <- paste("X",nodes[i],"node...P.1.0.", sep = "")
  node_cols[iter+3] <- paste("X",nodes[i],"node...P.1.1.", sep = "")
}

match(node_cols, names(run_nodes))

node_values <- as.data.frame(matrix(nrow = length(nodes), ncol = 8))

colnames(node_values) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)", "lungless", "lunged", "lotic", "lentic")
rownames(node_values) <- nodes


run_nodes[iteration : (iteration + 3)]


i = 1
iteration <- 1 + (i-1)*4     

for(i in 1:length(nodes)){  
  node_values[i,1:4] <- run_nodes[grep(rownames(node_values)[i], names(run_nodes))]
}

str_extract(str1, "\\d{10}")




for(i in 1:length(nodes)){                                        #fill in the nodes matrix with each state and combined states,
  iteration <- 1 + (i-1)*4                                  #such as "lunged" or "lotic"
  node_values[i,1:4] <- run_nodes[iteration : (iteration + 3)]
  node_values[i,5] <- sum(node_values[i,1], node_values[i,3], na.rm = TRUE)
  node_values[i,6] <- 1 - node_values[i,5]
  node_values[i,7] <- sum(node_values[i,1], node_values[i,2], na.rm = TRUE)
  node_values[i,8] <- 1 - node_values[i,7]
}

data <- read.csv("git/data/no_endo_lung_data.csv")
treedata <- data[data$Taxa %in% tree$tip.label,]; rm(data)
cols <- c("magenta", "blue", "red", "green")



treedata$state[match(tree$tip.label, treedata$Taxa)]

par(mar = c(1,1,1,1))
plot(tree, show.tip.label = TRUE, cex = .25, 
     tip.color = cols[treedata$state[match(tree$tip.label, treedata$Taxa)]], type = "f")
nodelabels(node = as.numeric(rownames(node_values)), pie = as.matrix(node_values[,1:4]), cex = .45, piecol = cols)
par(mar = c(5.1, 4.1, 4.1, 2.1))











