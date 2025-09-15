## function to create a dataset phytools is comfortable with more than one binary traits, allowing missing data

make_simmap_data = function(data, tree, mode = "double binary",
                            state_labels, trait_col_names){
  
  rownames(data) <- data$Taxa

  data <- data[rownames(data) %in% tree$tip.label, , drop = F] #trim data to tree
  names <- data$Taxa
  if(dim(data)[1] < 1){
    stop("Rownames do not match tip labels")
  }
  
  N_states <- length(state_labels)
  
  
  if(mode == "single"){
    if(length(trait_col_names) != 1){
      stop("single mode used, but number of traits does not
           equal one")}
    
    data <- data[,match(trait_col_names, colnames(data)), drop = FALSE] # keep only relevant columns
    
    state_matrix <- matrix(0, nrow = dim(data)[1], ncol = N_states)           #empty matrix (for simmap with all zeroes)
    rownames(state_matrix) <- names                             #rownames for later
    colnames(state_matrix) <- state_labels                                 #col names will reflect states
    
    
    for(i in 1:dim(state_matrix)[1]){
    for(j in 1:N_states){
      if(!anyNA(data[i])){
        if(data[i] == j){
          state_matrix[i,j] <- 1}}
      if(anyNA(data[i])){
        state_matrix[i,j] <- 1/N_states}
    }}
    
    
    
  }
  
if(mode != "single"){
  data <- data[,match(trait_col_names, colnames(data))] # keep only relevant columns
  
  if(mode == "double binary"){
    if(N_states != 4){
      stop("double binary mode used (default), but number of states does not
           equal four")}
    if(length(trait_col_names) != 2){
      stop("double binary mode used (default), but number of traits does not
           equal two")}
  
    state_matrix <- matrix(0, nrow = dim(data)[1], ncol = N_states)           #empty matrix (for simmap with all zeroes)
    rownames(state_matrix) <- rownames(data)                               #rownames for later
    colnames(state_matrix) <- state_labels                                 #col names will reflect states
    
    for(i in 1:dim(state_matrix)[1]){
      if(!anyNA(data[i,])){
        if(data[i,1] == 0 & data[i,2] == 0){
          state_matrix[i,1] <- 1}
        if(data[i,1] == 1 & data[i,2] == 0){
          state_matrix[i,3] <- 1}
        if(data[i,1] == 0 & data[i,2] == 1){
          state_matrix[i,2] <- 1}
        if(data[i,1] == 1 & data[i,2] == 1){
          state_matrix[i,4] <- 1}
      }
      if(is.na(data[i,1]) & data[i,2] == 0){
        state_matrix[i,1] <- .5
        state_matrix[i,3] <- .5
        }
      if(is.na(data[i,1]) & data[i,2] == 1){
        state_matrix[i,2] <- .5
        state_matrix[i,4] <- .5
      }
      if(data[i,2] == "-" & data[i,1] == 0 | is.na(data[i,2]) & data[i,1] == 0){
        state_matrix[i,1] <- .5
        state_matrix[i,2] <- .5
      }
      if(data[i,2] == "-" & data[i,1] == 1 | is.na(data[i,2]) & data[i,1] == 1){
        state_matrix[i,3] <- .5
        state_matrix[i,4] <- .5
      }
      if(is.na(data[i,1]) & data[i,2] == "-" | is.na(data[i,2])){
        state_matrix[i,1:4] <- (1/N_states)
      }}}
  
  
  if(mode == "six_state"){
    if(N_states != 6){
      stop("six state mode used, which was developed for a very specific use-case")}
    if(length(trait_col_names) != 3){
      stop("three binary traits needed")}
    
    state_matrix <- matrix(0, nrow = dim(data)[1], ncol = N_states)           #empty matrix (for simmap with all zeroes)
    rownames(state_matrix) <- rownames(data)                               #rownames for later
    colnames(state_matrix) <- state_labels                                 #col names will reflect states
  
    for(i in 1:dim(data)[1]){
      if(!anyNA(data[i,])){      
        if(data[i,1] == 0 & data[i,2] == 0 & data[i,3] == 0){
          state_matrix[i,1] <- 1}    
        if(data[i,1] == 0 & data[i,2] == 0 & data[i,3] == 1){
          state_matrix[i,2] <- 1}
        if(data[i,1] == 0 & data[i,2] == 1 & data[i,3] == 0){
          state_matrix[i,3] <- 1}
        if(data[i,1] == 0 & data[i,2] == 1 & data[i,3] == 1){
          state_matrix[i,4] <- 1}
        if(data[i,1] == 1 & data[i,2] == 0 & data[i,3] == 0){
          state_matrix[i,5] <- 1}      
        if(data[i,1] == 1 & data[i,2] == 1 & data[i,3] == 0){
          state_matrix[i,6] <- 1}}
      
      if(is.na(data[i,1]) & data[i,2] == 0 & data[i,3] == 0){
        state_matrix[i,1] <- .5
        state_matrix[i,5] <- .5}
      if(is.na(data[i,1]) & data[i,2] == 1 & data[i,3] == 0){
        state_matrix[i,3] <- .5
        state_matrix[i,6] <- .5}
      if(data[i,2] == "-" & data[i,1] == 0 & data[i,3] == 0){
        state_matrix[i,1] <- .5
        state_matrix[i,3] <- .5}
      if(data[i,2] == "-" & data[i,1] == 0 & data[i,3] == 1){
        state_matrix[i,2] <- .5
        state_matrix[i,4] <- .5}
      if(data[i,2] == "-" & data[i,1] == 1 & data[i,3] == 0){
        state_matrix[i,5] <- .5
        state_matrix[i,6] <- .5}
    }
    }
  
  if(mode == "eight_state"){
    if(N_states != 8){
      stop("eight state mode used, which was developed for a very specific use-case")}
    if(length(trait_col_names) != 4){
      stop("four binary traits needed")}
    
    state_matrix <- matrix(0, nrow = dim(data)[1], ncol = N_states)           #empty matrix (for simmap with all zeroes)
    rownames(state_matrix) <- rownames(data)                               #rownames for later
    colnames(state_matrix) <- state_labels                                 #col names will reflect states
    
    for(i in 1:dim(data)[1]){
      if(!anyNA(data[i,])){      
        if(data[i,1] == 0 & data[i,2] == 0 & data[i,3] == 0 & data[i,4] == 0){
          state_matrix[i,1] <- 1}    
        if(data[i,1] == 0 & data[i,2] == 0 & data[i,3] == 1 & data[i,4] == 0){
          state_matrix[i,2] <- 1}
        if(data[i,1] == 0 & data[i,2] == 1 & data[i,3] == 0 & data[i,4] == 0){
          state_matrix[i,3] <- 1}
        if(data[i,1] == 0 & data[i,2] == 1 & data[i,3] == 1 & data[i,4] == 0){
          state_matrix[i,4] <- 1}
        if(data[i,1] == 1 & data[i,2] == 0 & data[i,3] == 0 & data[i,4] == 0){
          state_matrix[i,5] <- 1}      
        if(data[i,1] == 1 & data[i,2] == 1 & data[i,3] == 0 & data[i,4] == 0){
          state_matrix[i,6] <- 1}
        if(data[i,2] == 0 & data[i,4] == 1){
          state_matrix[i,7] <- 1}       
        if(data[i,2] == 1 & data[i,4] == 1){
          state_matrix[i,8] <- 1}}}}}
  
  return(state_matrix)
}

















  
  
  
  
  
  
  
  
  
  
  
  