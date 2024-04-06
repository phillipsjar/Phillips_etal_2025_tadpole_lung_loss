custom_Q = function(rates, state_labels = c("X_S", "Lu_S", "X_P", "Lu_P"), model = "dependent"){

  if(model == "dependent" | model == "independent"){Q_dim <- 4}
  if(model == "six_state"){
    Q_dim <- 6
    state_labels = c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")}
  if(model == "eight_state"){
    Q_dim <- 8
    state_labels = c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P", "X_terr", "Lu_terr")
    }

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
    rate_matrix[1,5] <- as.numeric(rates[4]*.01)
    rate_matrix[2,1] <- as.numeric(rates[6]*.01)
    rate_matrix[2,4] <- as.numeric(rates[8]*.01)
    rate_matrix[3,1] <- as.numeric(rates[11]*.01)
    rate_matrix[3,4] <- as.numeric(rates[13]*.01)
    rate_matrix[3,6] <- as.numeric(rates[15]*.01)
    rate_matrix[4,2] <- as.numeric(rates[17]*.01)
    rate_matrix[4,3] <- as.numeric(rates[18]*.01)
    rate_matrix[5,1] <- as.numeric(rates[21]*.01)
    rate_matrix[5,6] <- as.numeric(rates[26]*.01)
    rate_matrix[6,3] <- as.numeric(rates[28]*.01)
    rate_matrix[6,5] <- as.numeric(rates[30]*.01)
  }
  
  
  if(model == "eight_state"){
    rate_matrix[1,2] <- as.numeric(rates[1]*.001)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[1,3] <- as.numeric(rates[2]*.001)
    rate_matrix[1,4] <- as.numeric(rates[3]*.001)
    rate_matrix[1,5] <- as.numeric(rates[4]*.001)
    rate_matrix[1,6] <- as.numeric(rates[5]*.001)
    rate_matrix[1,7] <- as.numeric(rates[6]*.001)
    rate_matrix[1,8] <- as.numeric(rates[7]*.001)
  
    rate_matrix[2,1] <-  as.numeric(rates[8]*.001)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[2,3] <-  as.numeric(rates[9]*.001)
    rate_matrix[2,4] <- as.numeric(rates[10]*.001)
    rate_matrix[2,5] <- as.numeric(rates[11]*.001)
    rate_matrix[2,6] <- as.numeric(rates[12]*.001)
    rate_matrix[2,7] <- as.numeric(rates[13]*.001)
    rate_matrix[2,8] <- as.numeric(rates[14]*.001)
    
    rate_matrix[3,1] <- as.numeric(rates[15]*.001)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[3,2] <- as.numeric(rates[16]*.001)
    rate_matrix[3,4] <- as.numeric(rates[17]*.001)
    rate_matrix[3,5] <- as.numeric(rates[18]*.001)
    rate_matrix[3,6] <- as.numeric(rates[19]*.001)
    rate_matrix[3,7] <- as.numeric(rates[20]*.001)
    rate_matrix[3,8] <- as.numeric(rates[21]*.001)
    
    rate_matrix[4,1] <- as.numeric(rates[22]*.001)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[4,2] <- as.numeric(rates[23]*.001)
    rate_matrix[4,3] <- as.numeric(rates[24]*.001)
    rate_matrix[4,5] <- as.numeric(rates[25]*.001)
    rate_matrix[4,6] <- as.numeric(rates[26]*.001)
    rate_matrix[4,7] <- as.numeric(rates[27]*.001)
    rate_matrix[4,8] <- as.numeric(rates[28]*.001)
    
    rate_matrix[5,1] <- as.numeric(rates[29]*.001)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[5,2] <- as.numeric(rates[30]*.001)
    rate_matrix[5,3] <- as.numeric(rates[31]*.001)
    rate_matrix[5,4] <- as.numeric(rates[32]*.001)
    rate_matrix[5,6] <- as.numeric(rates[33]*.001)
    rate_matrix[5,7] <- as.numeric(rates[34]*.001)
    rate_matrix[5,8] <- as.numeric(rates[35]*.001)
    
    rate_matrix[6,1] <- as.numeric(rates[36]*.001)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[6,2] <- as.numeric(rates[37]*.001)
    rate_matrix[6,3] <- as.numeric(rates[38]*.001)
    rate_matrix[6,4] <- as.numeric(rates[39]*.001)
    rate_matrix[6,5] <- as.numeric(rates[40]*.001)
    rate_matrix[6,7] <- as.numeric(rates[41]*.001)
    rate_matrix[6,8] <- as.numeric(rates[42]*.001)
    
    rate_matrix[7,1] <- as.numeric(rates[43]*.001)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[7,2] <- as.numeric(rates[44]*.001)
    rate_matrix[7,3] <- as.numeric(rates[45]*.001)
    rate_matrix[7,4] <- as.numeric(rates[46]*.001)
    rate_matrix[7,5] <- as.numeric(rates[47]*.001)
    rate_matrix[7,6] <- as.numeric(rates[48]*.001)
    rate_matrix[7,8] <- as.numeric(rates[49]*.001)
    
    rate_matrix[8,1] <- as.numeric(rates[50]*.001)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[8,2] <- as.numeric(rates[51]*.001)
    rate_matrix[8,3] <- as.numeric(rates[52]*.001)
    rate_matrix[8,4] <- as.numeric(rates[53]*.001)
    rate_matrix[8,5] <- as.numeric(rates[54]*.001)
    rate_matrix[8,6] <- as.numeric(rates[55]*.001)
    rate_matrix[8,7] <- as.numeric(rates[56]*.001)}
  
  for(i in 1:Q_dim){
    rate_matrix[i,i] <- -(sum(rate_matrix[i,]))}       # the rates along the diagonal should sum each row to zero
  return(rate_matrix)
}