## function to create a dataset phytools is comfortable with more than one binary traits, allowing missing data

make_simmap_data_six_state = function(data, tree, state_labels, trait_col_names){
    
  data <- data[rownames(data) %in% tree$tip.label,] #trim data to tree
  if(dim(data)[1] < 1){
    stop("Rownames do not match tip labels")
  }
  data <- data[,match(trait_col_names, colnames(data))] # keep only relevant columns
  
  N_states <- length(state_labels)
  raw_states  <- data  
  
  state_matrix <- matrix(0, nrow = dim(raw_states)[1], ncol = N_states)           #empty matrix (for simmap with all zeroes)
  rownames(state_matrix) <- rownames(data)                               #rownames for later
  colnames(state_matrix) <- state_labels                                 #col names will reflect states
  
  # this loop populates the matrix to indicate which state each taxon is in, allowing for unknowns
  
  for(i in 1:dim(raw_states)[1]){
    if(!anyNA(raw_states[i,])){      
      if(raw_states[i,1] == 0 & raw_states[i,2] == 0 & raw_states[i,3] == 0){
        state_matrix[i,1] <- 1}    
      if(raw_states[i,1] == 0 & raw_states[i,2] == 0 & raw_states[i,3] == 1){
        state_matrix[i,2] <- 1}
      if(raw_states[i,1] == 0 & raw_states[i,2] == 1 & raw_states[i,3] == 0){
        state_matrix[i,3] <- 1}
      if(raw_states[i,1] == 0 & raw_states[i,2] == 1 & raw_states[i,3] == 1){
        state_matrix[i,4] <- 1}
      if(raw_states[i,1] == 1 & raw_states[i,2] == 0 & raw_states[i,3] == 0){
        state_matrix[i,5] <- 1}      
      if(raw_states[i,1] == 1 & raw_states[i,2] == 1 & raw_states[i,3] == 0){
        state_matrix[i,6] <- 1}}
    
    if(is.na(raw_states[i,1])){
      dummy_lung_status <- raw_states[i,2]
      if(dummy_lung_status == 0){
        state_matrix[i,1] <- .5
        state_matrix[i,5] <- .5}
      if(dummy_lung_status == 1){
        state_matrix[i,3] <- .5
        state_matrix[i,6] <- .5}}
 
    if(is.na(raw_states[i,2])){
      dummy_eco_status <- raw_states[i,c(1,3)]
      if(dummy_eco_status[1] == 0 & dummy_eco_status[2] == 0){
        state_matrix[i,1] <- .5
        state_matrix[i,3] <- .5}
      if(dummy_eco_status[1] == 0 & dummy_eco_status[2] == 1){
        state_matrix[i,2] <- .5
        state_matrix[i,4] <- .5}
      if(dummy_eco_status[1] == 1 & dummy_eco_status[2] == 0){
        state_matrix[i,5] <- .5
        state_matrix[i,6] <- .5}}
  } # loop to manually make a broken up matrix so simmap is happy
  return(state_matrix)
}





  
  
  
  
  
  
  
  
  
  
  
  