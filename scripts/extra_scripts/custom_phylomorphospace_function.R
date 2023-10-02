
# Custom phylomorphospace script function specifically for two binary traits (easily modifiable for true 4 state data)

#tree <- read.nexus(file = "Bayestraits/tree/maxLH_tree.nex") #reload full tree (with some missing data) for plotting ASE
#
#data <- read.csv(file = "data/no_endo_lung_data.csv")
#
#tips <- cbind(data$ecology, data$lung)
#row.names(tips) <- data$Taxa
#
#load(file = "data/4062_ASE_R_formatted.R")
#nodes <- as.data.frame(nodes)
#nodes <- cbind(nodes$lentic, nodes$lunged)
#rownames(nodes) <- 1:dim(nodes)[1]
##
#cols <- c("magenta", "skyblue", "red", "lightgreen")
#jitter = .025
#rm(data)
#tree is phylogenetic tree
#tips is an N_tips x 2 matrix of tip data (NAs acceptable and will be assumed equal probability of each state)
#nodes is an N_nodes x 2 matrix of node data (must be provided)
#cols is a vector of four colors for the four quadrants of the plot. I recommend avoiding greyscale.

custom_phylomorphospace = function(tree, tips, nodes, cols = c("magenta", "skyblue", "red", "lightgreen"), jitter = 0){
  require(ape)
  if(tree$Nnode != dim(nodes)[1]){
    print("Node data dimensions must match number of nodes in supplied tree")
    stop()
  }
  
  if(length(row.names(tips) %in% tree$tip.label) > 0){
    true_tips <- tips[row.names(tips) %in% tree$tip.label,] #trim taxa not present in the tree if needed
    print("Warning: number of tips in tree differs from tip data. Data trimmed to accommodate tree")
  } else{true_tips <- tips}
  if(missing(jitter)){
    jitter <- 0         #default value
  }
  cols <- c(cols, "grey")
  true_tips <- true_tips[match(tree$tip.label, row.names(true_tips)),] #sort tips as found in tree, for indexing
  
  states <- rep(NA, dim(true_tips)[1])
  names(states) <- tree$tip.label
  for(i in 1:length(states)){
    if(!(is.na(true_tips[i,1])) & !(is.na(true_tips[i,2]))){
    if(true_tips[i,1] == 0 & true_tips[i,2] == 0){states[i] <- 1}
    if(true_tips[i,1] == 0 & true_tips[i,2] == 1){states[i] <- 2}
    if(true_tips[i,1] == 1 & true_tips[i,2] == 0){states[i] <- 3}
    if(true_tips[i,1] == 1 & true_tips[i,2] == 1){states[i] <- 4}}
    if(true_tips[i,1] == 1 & is.na(true_tips[i,2])){states[i] <- 5} #for the four possible unknown combinations
    if(true_tips[i,1] == 0 & is.na(true_tips[i,2])){states[i] <- 5} #for the four possible unknown combinations
    if(is.na(true_tips[i,1]) & true_tips[i,2] == 1){states[i] <- 5} #for the four possible unknown combinations
    if(is.na(true_tips[i,1]) & true_tips[i,2] == 0){states[i] <- 5} #for the four possible unknown combinations
  }
  
  
  for(i in 1:dim(true_tips)[1]){       #add some randomly distributed noise for aesthetically pleasing jitter
    true_tips[i,1] <-   true_tips[i,1] + rnorm(1,0,jitter)
    true_tips[i,2] <-   true_tips[i,2] + rnorm(1,0,jitter)}
  
  true_tips[which(is.na(true_tips[,1])),1] <- .5
  true_tips[which(is.na(true_tips[,1])),2] <- round(true_tips[which(is.na(true_tips[,1])),2])
  true_tips[which(is.na(true_tips[,2])),1] <- round(true_tips[which(is.na(true_tips[,1])),1])
  true_tips[which(is.na(true_tips[,2])),2] <- .5
   
  pal1 <- colorRampPalette(c(cols[1], cols[2]))(101)
  pal2 <- colorRampPalette(c(cols[1], cols[3]))(101)
  pal3 <- colorRampPalette(c(cols[2], cols[4]))(101)
  pal4 <- colorRampPalette(c(cols[3], cols[4]))(101)  

  
  ################################## build a vector to connect tips to nodes to deeper nodes
  require(phytools)
  
  N_tips <- length(tree$tip.label)
  N_nodes <- tree$Nnode
  node_desc <- vector(mode = "list", length = N_nodes)
  node_dir_anc <- vector(mode = "list", length = N_nodes)
  tip_dir_anc <- vector(mode = "list", length = N_tips)
  

  for(i in 1:N_nodes){
    a <- getDescendants(tree, (N_tips + i))
    node_desc[i] <- list(a) }

  for(i in (N_tips+2):(N_nodes+N_tips)){  #tree tips first, then nodes, then skip the root, which has no ancestor
    a <- which(sapply(node_desc, function(x) i %in% x))                      # Index positions of matches
    b <- sort(a, decreasing = TRUE)[1]
    node_dir_anc[(i-N_tips)] <- list(b)
  }
  
  for(i in 1:N_tips){
    a <- which(sapply(node_desc, function(x) i %in% x))                      # Index positions of matches
    b <- sort(a, decreasing = TRUE)[1]
    tip_dir_anc[(i)] <- list(b)
  }
  
  for(i in 2:(N_nodes)){
    segments(nodes[i,1], nodes[i,2], nodes[unlist(node_dir_anc[i]),1], nodes[unlist(node_dir_anc[i]),2])
  }
  for(i in 1:(N_tips)){
    segments(true_tips[i,1], true_tips[i,2], nodes[unlist(tip_dir_anc[i]),1], nodes[unlist(tip_dir_anc[i]),2])
  } 
  
points(nodes[,1], nodes[,2], bg = "grey", pch = 21)

points(nodes[which(nodes[,1] > .67),1] , nodes[which(nodes[,1] > .67),2],
       bg = pal4[round(nodes[which(nodes[,1] > .67),2], 1)*100+1], pch = 21)

points(nodes[which(nodes[,2] > .67),1] , nodes[which(nodes[,2] > .67),2],
       bg = pal3[round(nodes[which(nodes[,2] > .67),1], 1)*100+1], pch = 21)
points(nodes[which(nodes[,1] < .33),1] , nodes[which(nodes[,1] < .33),2],
       bg = pal1[round(nodes[which(nodes[,1] < .33),2], 1)*100+1], pch = 21)
points(nodes[which(nodes[,2] < .33),1] , nodes[which(nodes[,2]  < .33),2],
       bg = pal2[round(nodes[which(nodes[,2]  < .33),1], 1)*100+1], pch = 21)
points(true_tips[,1], true_tips[,2], bg = cols[states], pch = 21, cex = .5)
}

#plot(tips[1], tips[2], ylim = c(-0.1,1.1), xlim = c(-0.1,1.1),
#     xlab = "Probability of lentic habitat", pch = 21,
#     ylab = "Probability of lungs", lwd = 1.5, cex = 1.6, type = "n")
#custom_phylomorphospace(tree,tips,nodes, jitter = .02)
