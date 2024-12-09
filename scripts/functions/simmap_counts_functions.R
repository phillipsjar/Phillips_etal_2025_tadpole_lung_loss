avg_simmap_counts = function(simmap){
  require(phytools)
  N <- length(simmap)
  list <- vector(mode = "list", length = N)             
  for(i in 1:N){
    list[i] <- list(countSimmap(simmap[[i]])$Tr)
  }
  Y <- do.call(cbind, list)
  Y <- array(Y, dim=c(dim(list[[1]]), length(list)))
  counts <- round(apply(Y, c(1,2), mean, na.rm = TRUE),1)
  rownames(counts) <- rownames(countSimmap(simmap[[1]])$Tr)
  colnames(counts) <- colnames(countSimmap(simmap[[1]])$Tr)
  return(counts)
}

sd_simmap_counts = function(simmap){
  require(phytools)
  N <- length(simmap)
  list <- vector(mode = "list", length = N)             
  for(i in 1:N){
    list[i] <- list(countSimmap(simmap[[i]])$Tr)
  }
  Y <- do.call(cbind, list)
  Y <- array(Y, dim=c(dim(list[[1]]), length(list)))
  counts <- round(apply(Y, c(1,2), sd, na.rm = TRUE),1)
  rownames(counts) <- rownames(countSimmap(simmap[[1]])$Tr)
  colnames(counts) <- colnames(countSimmap(simmap[[1]])$Tr)
  return(counts)
}