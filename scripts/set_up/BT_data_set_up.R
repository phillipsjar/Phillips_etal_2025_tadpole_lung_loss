# This script creates data matrices for Bayestraits and automates listing out taxa for
# Bayestraits ancestral reconstruction.
rm(list = ls())
setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2024")

#Requires having run through the set_up script to generate the data_no_endo files and correct trees
source("lung_loss_git/scripts/set_up/Initial_set_up.R")

require(ape)

data_aqu <- read.csv("lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")
aqu_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_aqu.nex")

aqu_tree

data_full <- read.csv("lung_loss_git/processed_data/lung_data/full_data.csv")

table(data_full$Genus)
length(unique(data_full$Genus))


full_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_full.nex")

full_tree

tree_set_aqu <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/tree_set_aqu.nex")
tree_set_full <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/tree_set_full.nex")
tree_set_six <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/tree_set_six.nex")

###########################################################################
################## HOW TO Make Bayestraits data matrices ##################
###########################################################################


# little function for brevity
BT_data_mat = function(data, tree, traits){
  temp_data <- data[data$Taxa %in% tree$tip.label,] #only include taxa in tree
  mat <- matrix(nrow = dim(temp_data)[1], ncol = (length(traits)+1))
  for(i in 1:length(traits)){
    mat[,1] <- temp_data$Taxa
    mat[,(1+i)] <- temp_data[,match(traits[i], colnames(temp_data))]}
  mat[which(is.na(mat[,2])),2] <- rep("-", length(which(is.na(mat[,2])))) #allows ecology to be unknown
  return(mat)}

data_six <- data_full[!(is.na(data_full$six_state)),]
data_eight <- data_full[!(is.na(data_full$eight_state)),]

BT_aqu_mat <- BT_data_mat(data_aqu, aqu_tree, c("ecology", "lung"))
six_state_mat <- BT_data_mat(data_six, full_tree, "six_state")  #for six state data (including specialized and non-specialized lotic as two different states)
eight_state_mat <- BT_data_mat(data_eight, full_tree, "eight_state")  #for six state data (including specialized and non-specialized lotic as two different states)

table(six_state_mat[,2])
table(eight_state_mat[,2])

dim(BT_aqu_mat) # doesn't include terrestrial species 
dim(eight_state_mat) # doesn't include species missing any ecology or lung info

#export tables in BT-friendly format
write.table(BT_aqu_mat,row.names=F, col.names=F, file = "lung_loss_git/bayestraits_trees_data/data/data_matrix.txt", sep = "\t",quote = FALSE)
write.table(six_state_mat,row.names=F, col.names=F, file = "lung_loss_git/bayestraits_trees_data/data/6state_data_matrix.txt", sep = "\t",quote = FALSE)
write.table(eight_state_mat,row.names=F, col.names=F, file = "lung_loss_git/bayestraits_trees_data/data/8state_data_matrix.txt", sep = "\t",quote = FALSE)


rm(list = c("eight_state_mat", "six_state_mat", "BT_aqu_mat", "data_aqu", "data_eight",
            "data_six", "BT_data_mat", "data_full"))

#########################################################################
## HOW TO FINISH BT COMMAND FILES with tags and nodes ##
#########################################################################

# Bayestraits performs ancestral reconstructions on a node by node basis
# you must list taxa in nodes, which this code automatically creates for
# all taxa in tree. Ancestral reconstructions will only work if nodes are present and monophyletic across trees
require(ape)
require(phytools)

ML_nodes_aqu <- vector(mode = "list", length = aqu_tree$Nnode)             # empty list to be filled with the descendents of each node of ML tree
ML_nodes_full <- vector(mode = "list", length = full_tree$Nnode)             # empty list to be filled with the descendents of each node of ML tree


for(i in 1:length(ML_nodes_aqu)){
  A <- getDescendants(aqu_tree, aqu_tree$Nnode + 1 + i)              # descendents of ML_tree node
  B <-  aqu_tree$tip.label[A]                                          # convert to names (alternatively, sort by number to only include tips)
  ML_nodes_aqu[[i]] <- B[!(is.na(B))]                                        # remove NAs, leaving only tips
}

for(i in 1:length(ML_nodes_full)){
  A <- getDescendants(full_tree, full_tree$Nnode + 1 + i)              # descendents of ML_tree node
  B <-  full_tree$tip.label[A]                                          # convert to names (alternatively, sort by number to only include tips)
  ML_nodes_full[[i]] <- B[!(is.na(B))]                                        # remove NAs, leaving only tips
}

ML_node_check_aqu  <- rep(NA, aqu_tree$Nnode)                               # vector to record whether node is present in tree set as well
ML_node_check_full <- rep(NA, full_tree$Nnode)                               # vector to record whether node is present in tree set as well

#################################################
#   Next step takes a couple minutes
#################################################
for(i in 1:100){                                                         # iterate across 100 trees
    for(j in 1:aqu_tree$Nnode){                                        # iterate across all nodes of ML tree
      if(!(is.monophyletic(tree_set_aqu[[i]], ML_nodes_aqu[[j]]))){           # if a clade in ML_tree is not monophyletic in at least 1 tree in the treeset,
        ML_node_check_aqu[j] <- 1                                            # then mark it as such (can take a while, a lot of comparisons, especially for big trees)
      }
  }
}

#################################################
#   Next step takes a couple minutes
#################################################
for(i in 1:100){                                                         # iterate across 100 trees
  for(j in 1:full_tree$Nnode){                                        # iterate across all nodes of ML tree
    if(!(is.monophyletic(tree_set_full[[i]], ML_nodes_full[[j]]))){           # if a clade in ML_tree is not monophyletic in at least 1 tree in the treeset,
      ML_node_check_full[j] <- 1                                            # then mark it as such (can take a while, a lot of comparisons, especially for big trees)
    }
  }
}

nodes_to_keep_aqu <- aqu_tree$Nnode + 1 + which(is.na(ML_node_check_aqu))     # change to reflect the way R keeps track of nodes vs. tips (assumes tree is bifurcating)
nodes_to_keep_full <- full_tree$Nnode + 1 + which(is.na(ML_node_check_full))     # change to reflect the way R keeps track of nodes vs. tips (assumes tree is bifurcating)

save(nodes_to_keep_aqu, file = "lung_loss_git/bayestraits_trees_data/data/nodes_to_keep_aqu.Rdata")      # save this for later plotting, since it takes a while to calculate
save(nodes_to_keep_full, file = "lung_loss_git/bayestraits_trees_data/data/nodes_to_keep_full.Rdata")      # save this for later plotting, since it takes a while to calculate
#load(file = "lung_loss_git/bayestraits_trees_data/data/nodes_to_keep_aqu.Rdata")
#load(file = "lung_loss_git/bayestraits_trees_data/data/nodes_to_keep_full.Rdata")

aqu_data <- read.csv("lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")
aqu_data <- aqu_data[aqu_data$Taxa %in% aqu_tree$tip.label,]
cols <- c("magenta", "blue", "red", "green")
par(xpd = TRUE)
plot(aqu_tree, show.tip.label = TRUE, cex = .25, 
     tip.color = cols[aqu_data$state[match(aqu_tree$tip.label, aqu_data$Taxa)]], type = "f")


nodelabels(node = nodes_to_keep_aqu, frame = "circle", cex = .2)


full_data <- read.csv("lung_loss_git/processed_data/lung_data/full_data.csv")
full_data <- full_data[full_data$Taxa %in% full_tree$tip.label,]

cols2 <- c("gold", "purple", "#99CCFF", "blue", "red", "lightgreen", "#844200", "#D28231")

plot(full_tree, show.tip.label = TRUE, cex = .25, 
     tip.color = cols2[full_data$eight_state[match(full_tree$tip.label, full_data$Taxa)]], type = "f")

nodelabels(node = nodes_to_keep_full, frame = "circle", cex = .2)


rm(list = c("A", "B", "i", "j", "ML_node_check_aqu", "ML_nodes_aqu", "aqu_data", 
            "ML_node_check_full", "ML_nodes_full", "cols", "cols2", "full_data"))





# these functions actually write out the taxa names as needed by bayestraits for tags and nodes

namelist <- function(x,nodes) {
  for (i in 1:length(x)){
    strlist <- c("AddNode ",nodes[i],"node ",nodes[i],"tag")
    cat(paste(strlist, collapse=''),"\n")
  }
}

fnlist <- function(x) {
  z <- deparse(substitute(x))
  # cat(z, "\n")
  nams=names(x) 
  for (i in 1:length(x)){
    names <- cat(nams[i],  x[[i]], "\n", "\n")
    cat(paste(names, collapse=''))
  } 
}

# Get list of nodes and total number of tips 
# Create empty list
desc_list_aqu <- vector(mode = "list", length = length(nodes_to_keep_aqu))  #output list of descendents for each node
desc_list_full <- vector(mode = "list", length = length(nodes_to_keep_full))  #output list of descendents for each node

for(i in 1:length(nodes_to_keep_aqu)){
  a <- aqu_tree$tip.label[getDescendants(aqu_tree, nodes_to_keep_aqu[i])]   #get descendents for each node 
  b <- a[which(a != "NA")]                                       #remove nodes
  desc_list_aqu[[i]] <- b                                 #add to output as list
  names(desc_list_aqu)[i] <- paste("AddTag ", nodes_to_keep_aqu[i], "tag", sep = "")
}

for(i in 1:length(nodes_to_keep_full)){
  a <- full_tree$tip.label[getDescendants(full_tree, nodes_to_keep_full[i])]   #get descendents for each node 
  b <- a[which(a != "NA")]                                       #remove nodes
  desc_list_full[[i]] <- b                                 #add to output as list
  names(desc_list_full)[i] <- paste("AddTag ", nodes_to_keep_full[i], "tag", sep = "")
}

#preview what they look like
fnlist(desc_list_aqu)
namelist(desc_list_aqu,nodes_to_keep_aqu)


### creates and writes portion of command file into a new file named 4062_node_list.txt

sink("lung_loss_git/bayestraits_trees_data/aqu_node_list.txt")
cat("\n","\n")
fnlist(desc_list_aqu)
cat("\n","\n")
namelist(desc_list_aqu,nodes_to_keep_aqu)
cat("\n","run", sep = "")
sink()

sink("lung_loss_git/bayestraits_trees_data/full_node_list.txt")
cat("\n","\n")
fnlist(desc_list_full)
cat("\n","\n")
namelist(desc_list_full,nodes_to_keep_full)
cat("\n","run", sep = "")
sink()

#clean up
rm(list = c("a", "b", "i", "desc_list_aqu", "fnlist", "namelist", "nodes_to_keep_aqu",
            "aqu_tree", "full_tree", "tree_set_aqu", "tree_set_full", "desc_list_full",
            "nodes_to_keep_full"))







