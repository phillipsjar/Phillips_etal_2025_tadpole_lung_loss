

# convert a given tree to a time tree, based on an already calibrated time tree.
# the topology can be different, and is extremely fast. 
time.calibrate = function(tree1, cal_tree){
  require(dispRity)
  require(ape)
  require(phytools)
  node_list_tree1 <- vector(mode = "list", length = tree1$Nnode)     #output list of descendents for each node of our tree
  df <- as.data.frame(matrix(nrow = tree1$Nnode, ncol = 3))
  colnames(df) <- c("node1", "node2", "node_age")
  df$node1 <- seq(1:tree1$Nnode)
  
  cal_tree_ages <- tree.age(cal_tree) #requires dispRity, outputs a nice easy table of ages for each node in cal_tree
  
  #This next bit finds all the descendents for each node in our target tree (tree1), removes internal nodes, and outputs that as a list of taxa
  
  for(i in 1:tree1$Nnode){
    A <- getDescendants(tree1,(length(tree1$tip.label)+i))
    B <- A[which(A<=(length(tree1$tip.label)))]
    node_list_tree1[i] <- list(tree1$tip.label[B])}
  for(i in 1:tree1$Nnode){
    A <- unlist(node_list_tree1[i])
    B <- getMRCA(cal_tree, A)
    df$node2[i] <- B
    df$node_age[i] <- cal_tree_ages$ages[B]}
    rm(A)
    rm(B)
    rm(i)
  
    df2 <- df[(which(!(duplicated(df$node2)))),]  # remove all duplicates, leaving only the lowest number
    rm(df)
    
    cal <- as.data.frame(matrix(nrow = dim(df2)[1], ncol = 4))  # creates a calibration table in the style chronos uses
    colnames(cal) <- colnames(makeChronosCalib(tree1))
    cal$node <- (df2$node1+length(tree1$tip.label))
    cal$soft.bounds <- rep(TRUE)
    for(i in 1:dim(cal)[1]){
      cal$age.min[i] <- round(df2$node_age[i]*.99)              # .99 means the lower bound must be within 1% of the age in cal_tree
      cal$age.max[i] <- round(df2$node_age[i]*1.01)             # this can be widened, at the expense of calculation time.
    }
    rm(i)

  time_tree <- chronos(tree1, cal = cal)
  return(time_tree)
}