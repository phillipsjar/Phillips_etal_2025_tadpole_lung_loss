custom_Q = function(rates, state_labels = c("X_S", "Lu_S", "X_P", "Lu_P"), model = "dependent"){

  if(model == "dependent" | model == "independent"){Q_dim <- 4}
  if(model == "six_state"){Q_dim <- 6}
  
  rate_matrix <- matrix(0, nrow = Q_dim, ncol = Q_dim)         # matrix is populated with zeroes to begin
  colnames(rate_matrix) <- state_labels
  rownames(rate_matrix) <- state_labels
  
  if(model == "dependent"){
  rate_matrix[1,2] <- as.numeric(rates[1]*.001)                   # rates are multiplied by .001 to counteract scaling in BT
  rate_matrix[1,3] <- as.numeric(rates[2]*.001)
  rate_matrix[2,1] <- as.numeric(rates[3]*.001)
  rate_matrix[2,4] <- as.numeric(rates[4]*.001)
  rate_matrix[3,1] <- as.numeric(rates[5]*.001)
  rate_matrix[3,4] <- as.numeric(rates[6]*.001)
  rate_matrix[4,2] <- as.numeric(rates[7]*.001)
  rate_matrix[4,3] <- as.numeric(rates[8]*.001)
  }
  
  if(model == "independent"){
    rate_matrix[1,2] <- as.numeric(rates[3]*.001)                   # rates are multiplied by .001 to counteract scaling in BT
    rate_matrix[1,3] <- as.numeric(rates[1]*.001)
    rate_matrix[2,1] <- as.numeric(rates[4]*.001)
    rate_matrix[2,4] <- as.numeric(rates[1]*.001)
    rate_matrix[3,1] <- as.numeric(rates[2]*.001)
    rate_matrix[3,4] <- as.numeric(rates[3]*.001)
    rate_matrix[4,2] <- as.numeric(rates[2]*.001)
    rate_matrix[4,3] <- as.numeric(rates[4]*.001)
  }
  
  if(model == "six_state"){
    rate_matrix[1,2] <- as.numeric(rates[1]*.01)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[1,3] <- as.numeric(rates[2]*.01)
    rate_matrix[1,5] <- as.numeric(rates[3]*.01)
    rate_matrix[2,1] <- as.numeric(rates[4]*.01)
    rate_matrix[2,4] <- as.numeric(rates[5]*.01)
    rate_matrix[3,1] <- as.numeric(rates[6]*.01)
    rate_matrix[3,4] <- as.numeric(rates[7]*.01)
    rate_matrix[3,6] <- as.numeric(rates[8]*.01)
    rate_matrix[4,2] <- as.numeric(rates[9]*.01)
    rate_matrix[4,3] <- as.numeric(rates[10]*.01)
    rate_matrix[5,1] <- as.numeric(rates[11]*.01)
    rate_matrix[5,6] <- as.numeric(rates[12]*.01)
    rate_matrix[6,3] <- as.numeric(rates[13]*.01)
    rate_matrix[6,5] <- as.numeric(rates[14]*.01)
  }
  
  for(i in 1:Q_dim){
    rate_matrix[i,i] <- -(sum(rate_matrix[i,]))}       # the rates along the diagonal should sum each row to zero
  return(rate_matrix)
}