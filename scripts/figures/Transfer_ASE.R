convert_ancestral_state = function(tree1, tree2, ASE){
  nodes <- as.data.frame(ASE)
  nodes$tree2_node_number <- rep(NA)
  if(dim(nodes)[1] != tree1$Nnode){
    print("ASE dataframe must have the same number of rows as nodes present in tree 1");
    stop()}
  
  tree1_nodes <- vector(mode = "list", length = tree1$Nnode) #blank list of lists with a list for each node of tree 1
  
  names(tree1_nodes) <- (length(tree1$tip.label)+1): (length(tree1$tip.label)+tree1$Nnode) #name each list with the node number for tree 1
  for(i in 1:length(tree1_nodes)){
    tree1_nodes[i] <- list(tree1$tip.label[getDescendants(tree1, names(tree1_nodes)[i])]) #fill each list the descendents in tree 1
  }
  tree1_nodes <- lapply(tree1_nodes, function(x) x[!is.na(x)])                            # remove NAs (internal nodes)
  missing <- tree1$tip.label[!(tree1$tip.label %in% tree2$tip.label)]                      # vector of the taxa missing in tree 2
  tree_nodes_fixed <- lapply(tree1_nodes, function(x) x[!(x %in% missing)])               # remove the missing taxa from the lists
  
  for(i in 1:dim(nodes)[1]){
    if(length(getMRCA(tree2, unlist(tree_nodes_fixed[i])) > 0)){                         # use the new lists of taxa (all present in tree 2)
      nodes$tree2_node_number[i] <- getMRCA(tree2, unlist(tree_nodes_fixed[i]))          # to get the associated nodes in tree 2
    }
  }
  return(nodes)
}
