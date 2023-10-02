
namelist <- function(x,tips) {
  for (i in 1:length(x)){
    strlist <- c("AddNode ",i+tips,"node ",i+tips,"tag")
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
desc_list <- vector(mode = "list", length = BT_ML_tree$Nnode)  #output list of descendents for each node
nodes <- BT_ML_tree$Nnode
tips <- length(BT_ML_tree$tip.label)

for(i in (tips+1):(tips+nodes)){
  a <- BT_ML_tree$tip.label[getDescendants(BT_ML_tree, i)]   #get descendents for each node 
  b <- a[which(a != "NA")]                                       #remove nodes
  desc_list[i-tips] <- list(b)                                   #add to output as list
  names(desc_list)[i-tips] <- paste("AddTag ", i, "tag", sep = "")
}

#preview what they look like
fnlist(desc_list)
namelist(desc_list,tips)


### creates and writes portion of command file into a new file named 4062_node_list.txt

sink("Bayestraits/node_list.txt")
cat("\n","\n")
fnlist(desc_list)
cat("\n","\n")
namelist(desc_list,tips)
cat("\n","run", sep = "")
sink()

#clean up
rm(list = c("a", "b", "i", "nodes", "tips", "desc_list", "fnlist", "namelist", 
            "BT_ML_tree"))

