#Load in combined Bayestraits posterior

load(file = "bayestraits_exports/master_dep_export.Rdata")    #dependent
load(file = "bayestraits_exports/master_indep_export.Rdata")  #independent

data_no_endo <- read.csv(file = "lung_loss_git/processed_data/lung_data/no_endo_lung_data.csv")


require(ape)

BT_ML_tree <- read.tree(file = "lung_loss_git/trees/edited_trees/maxLH_BT_vis_tree.tre")

node_cols <- colnames(master_indep_export)[grep("node...P.", colnames(master_indep_export))] #nodes estimated in BT

dep_nodes <- as.data.frame(matrix(nrow = length(node_cols)/4, ncol = 8))     # matrix of reconstructions. Each row corresponds to an internal node
# named columns give the estimated probability that node is a given state
colnames(dep_nodes) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)", "lungless", "lunged", "lotic", "lentic")
rownames(dep_nodes) <- unique(gsub("X", "", gsub("node.........", "", node_cols)))


dep_averages <- sapply(master_dep_export[,match(node_cols, colnames(master_dep_export))], "mean")  #average the probability each sampled node is in each given state
#across the entire posterior


indep_nodes <- as.data.frame(matrix(nrow = length(node_cols)/4, ncol = 8))     # matrix of reconstructions. Each row corresponds to an internal node
# named columns give the estimated probability that node is a given state
colnames(indep_nodes) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)", "lungless", "lunged", "lotic", "lentic")
rownames(indep_nodes) <- unique(gsub("X", "", gsub("node.........", "", node_cols)))

indep_averages <- sapply(master_indep_export[,match(node_cols, colnames(master_indep_export))], "mean")  #average the probability each sampled node is in each given state



for(i in 1:dim(dep_nodes)[1]){                                        #fill in the dep_nodes matrix with each state and combined states,
  iteration <- 1 + (i-1)*4                                  #such as "lunged" or "lotic"
  dep_nodes[i,1:4] <- dep_averages[iteration : (iteration + 3)]
  dep_nodes[i,5] <- sum(dep_nodes[i,1], dep_nodes[i,3], na.rm = TRUE)
  dep_nodes[i,6] <- 1 - dep_nodes[i,5]
  dep_nodes[i,7] <- sum(dep_nodes[i,1], dep_nodes[i,2], na.rm = TRUE)
  dep_nodes[i,8] <- 1 - dep_nodes[i,7]
}
for(i in 1:dim(indep_nodes)[1]){                                        #fill in the indep_nodes matrix with each state and combined states,
  iteration <- 1 + (i-1)*4                                  #such as "lunged" or "lotic"
  indep_nodes[i,1:4] <- indep_averages[iteration : (iteration + 3)]
  indep_nodes[i,5] <- sum(indep_nodes[i,1], indep_nodes[i,3], na.rm = TRUE)
  indep_nodes[i,6] <- 1 - indep_nodes[i,5]
  indep_nodes[i,7] <- sum(indep_nodes[i,1], indep_nodes[i,2], na.rm = TRUE)
  indep_nodes[i,8] <- 1 - indep_nodes[i,7]
}

save(dep_nodes, file = "lung_loss_git/processed_data/BT_output/ASE_R_formatted.R")
save(indep_nodes, file = "lung_loss_git/processed_data/BT_output/ASE_R_formatted.R")

rm(list = c("data", "dep_averages", "indep_averages", 
            "iter", "i", "iteration", "node_cols"))

#check and make sure it makes sense
data <- read.csv("git/data/no_endo_lung_data.csv")
treedata <- data[data$Taxa %in% BT_ML_tree$tip.label,]; rm(data)
cols <- c("magenta", "blue", "red", "green")

par(mar = c(1,1,1,1))
plot(BT_ML_tree, show.tip.label = TRUE, cex = .25, 
     tip.color = cols[treedata$state[match(BT_ML_tree$tip.label, treedata$Taxa)]], type = "f")
nodelabels(node = nodes_to_keep, pie = as.matrix(nodes[,1:4]), cex = .3, 
           piecol = c("magenta", "blue", "red", "green"))
par(mar = c(5.1, 4.1, 4.1, 2.1))

rm(list = c("BT_ML_tree", "dep_nodes", "indep_nodes", "treedata", "cols", "nodes_to_keep"))


#############################################################################
#                 Now for the six-state models

load(file = "bayestraits_exports/master_six_state_dep_export.Rdata")    #dependent
load(file = "bayestraits_exports/master_six_state_indep_export.Rdata")    #dependent

data_no_endo <- read.csv(file = "lung_loss_git/processed_data/lung_data/no_endo_lung_data.csv")


BT_ML_tree <- read.tree(file = "lung_loss_git/trees/edited_trees/maxLH_BT_vis_tree.tre")


node_cols <- colnames(master_six_state_dep_export)[grep("node.P.", 
                      colnames(master_six_state_dep_export))] #nodes estimated in BT

node_cols <- node_cols[7:length(node_cols)]
  
  
six_dep_nodes <- as.data.frame(matrix(nrow = length(node_cols)/6, ncol = 8))     # matrix of reconstructions. Each row corresponds to an internal node
# named columns give the estimated probability that node is a given state
colnames(six_dep_nodes) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)", "lungless", "lunged", "lotic", "lentic")
rownames(six_dep_nodes) <- unique(gsub("X", "", gsub("node.....", "", node_cols)))



six_dep_averages <- sapply(master_six_state_dep_export[,match(node_cols, 
                                        colnames(master_six_state_dep_export))], "mean")  #average the probability each sampled node is in each given state


#across the entire posterior






















