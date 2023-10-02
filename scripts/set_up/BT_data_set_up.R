# This script creates data matrices for Bayestraits and automates listing out taxa for
# Bayestraits ancestral reconstruction.

#Requires having run through the set_up script to generate the data_no_endo files and correct trees
require(ape)

data_no_endo <- read.csv("git/data/no_endo_lung_data.csv")

BT_ML_tree <- read.nexus(file = "git/bayestraits/trees/maxLH_tree.nex")
BT_tree_set <- read.nexus(file = "git/bayestraits/trees/tree_set.nex")
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


BT_mat <- BT_data_mat(data_no_endo, BT_ML_tree, c("ecology", "lung"))

six_state_mat <- BT_data_mat(data_no_endo, BT_ML_tree, "six_state")  #for six state data (including specialized and non-specialized lotic as two different states)

dim(BT_mat) # check that the numbers match set up 
dim(six_state_mat)

#export tables in BT-friendly format
write.table(BT_mat,row.names=F, col.names=F, file = "git/bayestraits/data/data_matrix.txt", sep = "\t",quote = FALSE)
write.table(six_state_mat,row.names=F, col.names=F, file = "git/bayestraits/data/6state_data_matrix.txt", sep = "\t",quote = FALSE)


rm(list = c("BT_mat", "six_state_mat", "BT_data_mat", "data_no_endo"))

#########################################################################
## HOW TO FINISH BT COMMAND FILES with tags and nodes ##
#########################################################################

# Bayestraits performs ancestral reconstructions on a node by node basis
# you must list taxa in nodes, which this code automatically creates for
# all taxa in tree. Ancestral reconstructions will only work if nodes are present and monophyletic across trees
require(ape)
require(phytools)

ML_nodes <- vector(mode = "list", length = BT_ML_tree$Nnode)             # empty list to be filled with the descendents of each node of ML tree


for(i in 1:length(ML_nodes)){
  A <- getDescendants(BT_ML_tree, BT_ML_tree$Nnode + 1 + i)              # descendents of ML_tree node
  B <-  BT_ML_tree$tip.label[A]                                          # convert to names (alternatively, sort by number to only include tips)
  ML_nodes[[i]] <- B[!(is.na(B))]                                        # remove NAs, leaving only tips
}

ML_node_check <- rep(NA, BT_ML_tree$Nnode)                               # vector to record whether node is present in tree set as well

#################################################
#   Next step takes a couple minutes
#################################################
for(i in 1:100){                                                         # iterate across 100 trees
    for(j in 1:BT_ML_tree$Nnode){                                        # iterate across all nodes of ML tree
      if(!(is.monophyletic(BT_tree_set[[i]], ML_nodes[[j]]))){           # if a group in the ML_tree is not monophyletic in at least 1 tree in the treeset,
        ML_node_check[j] <- 1                                            # then mark it as such (can take a while, a lot of comparisons, especially for big trees)
      }
  }
}

nodes_to_keep <- BT_ML_tree$Nnode + 1 + which(is.na(ML_node_check))     # change to reflect the way R keeps track of nodes vs. tips (assumes tree is bifurcating)
save(nodes_to_keep, file = "git/bayestraits/data/nodes_to_keep.Rdata")      # save this for later plotting, since it takes a while to calculate
#load(file = "git/bayestraits/data/nodes_to_keep.Rdata")

data <- read.csv("git/data/no_endo_lung_data.csv")
treedata <- data[data$Taxa %in% tree$tip.label,]; rm(data)
cols <- c("magenta", "blue", "red", "green")

plot(BT_ML_tree, show.tip.label = TRUE, cex = .25, 
     tip.color = cols[treedata$state[match(BT_ML_tree$tip.label, treedata$Taxa)]], type = "f")

nodelabels(node = nodes_to_keep, frame = "circle", cex = .2)

rm(list = c("A", "B", "i", "j", "ML_node_check", "ML_nodes", "data", "treedata"))

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
desc_list <- vector(mode = "list", length = length(nodes_to_keep))  #output list of descendents for each node

for(i in 1:length(nodes_to_keep)){
  a <- BT_ML_tree$tip.label[getDescendants(BT_ML_tree, nodes_to_keep[i])]   #get descendents for each node 
  b <- a[which(a != "NA")]                                       #remove nodes
  desc_list[[i]] <- b                                 #add to output as list
  names(desc_list)[i] <- paste("AddTag ", nodes_to_keep[i], "tag", sep = "")
}

#preview what they look like
fnlist(desc_list)
namelist(desc_list,nodes_to_keep)


### creates and writes portion of command file into a new file named 4062_node_list.txt

sink("git/Bayestraits/node_list.txt")
cat("\n","\n")
fnlist(desc_list)
cat("\n","\n")
namelist(desc_list,nodes_to_keep)
cat("\n","run", sep = "")
sink()

#clean up
rm(list = c("a", "b", "i", "desc_list", "fnlist", "namelist", "nodes_to_keep",
            "BT_ML_tree", "BT_tree_set"))







