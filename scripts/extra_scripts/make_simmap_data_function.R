## function to create a dataset phytools is comfortable with two binary traits allowing missing data

make_simmap_data = function(data, tree){
  data <- data[data$Taxa %in% tree$tip.label,] #trim data to tree
  lung <- (data$lung)                                          #lung data
  eco <- (data$ecology)                                        #ecology data binary    
  lung_eco <- cbind(lung, eco)                                              #stick 'em together    
  lung_eco_matrix <- matrix(0, nrow = dim(lung_eco)[1], ncol = 4)           #empty matrix (for simmap with all zeroes)
  rownames(lung_eco_matrix) <- data$Taxa                               #rownames for later
  colnames(lung_eco_matrix) <- c("X_S", "Lu_S", "X_P", "Lu_P")              #col names will reflect states
  
  # this loop populates the matrix to indicate which state each taxon is in, allowing for unknowns
  
  for(i in 1:dim(lung_eco)[1]){
    if(!anyNA(lung_eco[i,])){ 
      if(lung_eco[i,1] == 0 & lung_eco[i,2] == 0){
        lung_eco_matrix[i,1] <- 1}    
      if(lung_eco[i,1] == 1 & lung_eco[i,2] == 0){
        lung_eco_matrix[i,2] <- 1}
      if(lung_eco[i,1] == 0 & lung_eco[i,2] == 1){
        lung_eco_matrix[i,3] <- 1}
      if(lung_eco[i,1] == 1 & lung_eco[i,2] == 1){
        lung_eco_matrix[i,4] <- 1}}
    if(is.na(lung_eco[i,2])){
      dummy_lung_status <- lung_eco[i,1]
      if(dummy_lung_status == 0){
        lung_eco_matrix[i,1] <- .5
        lung_eco_matrix[i,3] <- .5}
      if(dummy_lung_status == 1){
        lung_eco_matrix[i,2] <- .5
        lung_eco_matrix[i,4] <- .5}}
    if(is.na(lung_eco[i,1])){
      dummy_eco_status <- lung_eco[i,2]
      if(dummy_eco_status == 0){
        lung_eco_matrix[i,1] <- .5
        lung_eco_matrix[i,2] <- .5}
      if(dummy_eco_status == 1){
        lung_eco_matrix[i,3] <- .5
        lung_eco_matrix[i,4] <- .5}}
    
    
  } # loop to manually make a broken up matrix so simmap is happy
  return(lung_eco_matrix)
}
