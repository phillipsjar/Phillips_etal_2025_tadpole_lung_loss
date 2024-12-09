custom_Q = function(rates, state_labels = c("X_S", "Lu_S", "X_P", "Lu_P"), model = "dependent", scale = 1){

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
  rate_matrix[1,2] <- as.numeric(rates[1]*scale)                   # rates are multiplied by .001 to counteract scaling in BT
  rate_matrix[1,3] <- as.numeric(rates[2]*scale)
  rate_matrix[2,1] <- as.numeric(rates[3]*scale)
  rate_matrix[2,4] <- as.numeric(rates[4]*scale)
  rate_matrix[3,1] <- as.numeric(rates[5]*scale)
  rate_matrix[3,4] <- as.numeric(rates[6]*scale)
  rate_matrix[4,2] <- as.numeric(rates[7]*scale)
  rate_matrix[4,3] <- as.numeric(rates[8]*scale)
  }
  
  if(model == "independent"){
    rate_matrix[1,2] <- as.numeric(rates[3]*scale)                   # rates are multiplied by .001 to counteract scaling in BT
    rate_matrix[1,3] <- as.numeric(rates[1]*scale)
    rate_matrix[2,1] <- as.numeric(rates[4]*scale)
    rate_matrix[2,4] <- as.numeric(rates[1]*scale)
    rate_matrix[3,1] <- as.numeric(rates[2]*scale)
    rate_matrix[3,4] <- as.numeric(rates[3]*scale)
    rate_matrix[4,2] <- as.numeric(rates[2]*scale)
    rate_matrix[4,3] <- as.numeric(rates[4]*scale)
  }
  
  
  if(model == "six_state"){
    rate_matrix[1,2] <- as.numeric(rates[1]*scale)                   # rates are multiplied by .01 to counteract scaling in BT
    rate_matrix[1,3] <- as.numeric(rates[2]*scale)
    rate_matrix[1,5] <- as.numeric(rates[4]*scale)
    rate_matrix[2,1] <- as.numeric(rates[6]*scale)
    rate_matrix[2,4] <- as.numeric(rates[8]*scale)
    rate_matrix[3,1] <- as.numeric(rates[11]*scale)
    rate_matrix[3,4] <- as.numeric(rates[13]*scale)
    rate_matrix[3,6] <- as.numeric(rates[15]*scale)
    rate_matrix[4,2] <- as.numeric(rates[17]*scale)
    rate_matrix[4,3] <- as.numeric(rates[18]*scale)
    rate_matrix[5,1] <- as.numeric(rates[21]*scale)
    rate_matrix[5,6] <- as.numeric(rates[26]*scale)
    rate_matrix[6,3] <- as.numeric(rates[28]*scale)
    rate_matrix[6,5] <- as.numeric(rates[30]*scale)
  }
  
  
  if(model == "eight_state"){
    rate_matrix[1,2] <- as.numeric(rates[1]*scale)                   
    rate_matrix[1,3] <- as.numeric(rates[2]*scale)
    rate_matrix[1,4] <- as.numeric(rates[3]*scale)
    rate_matrix[1,5] <- as.numeric(rates[4]*scale)
    rate_matrix[1,6] <- as.numeric(rates[5]*scale)
    rate_matrix[1,7] <- as.numeric(rates[6]*scale)
    rate_matrix[1,8] <- as.numeric(rates[7]*scale)
  
    rate_matrix[2,1] <-  as.numeric(rates[8]*scale)                   
    rate_matrix[2,3] <-  as.numeric(rates[9]*scale)
    rate_matrix[2,4] <- as.numeric(rates[10]*scale)
    rate_matrix[2,5] <- as.numeric(rates[11]*scale)
    rate_matrix[2,6] <- as.numeric(rates[12]*scale)
    rate_matrix[2,7] <- as.numeric(rates[13]*scale)
    rate_matrix[2,8] <- as.numeric(rates[14]*scale)
    
    rate_matrix[3,1] <- as.numeric(rates[15]*scale)                   
    rate_matrix[3,2] <- as.numeric(rates[16]*scale)
    rate_matrix[3,4] <- as.numeric(rates[17]*scale)
    rate_matrix[3,5] <- as.numeric(rates[18]*scale)
    rate_matrix[3,6] <- as.numeric(rates[19]*scale)
    rate_matrix[3,7] <- as.numeric(rates[20]*scale)
    rate_matrix[3,8] <- as.numeric(rates[21]*scale)
    rate_matrix[4,1] <- as.numeric(rates[22]*scale)                   
    rate_matrix[4,2] <- as.numeric(rates[23]*scale)
    rate_matrix[4,3] <- as.numeric(rates[24]*scale)
    rate_matrix[4,5] <- as.numeric(rates[25]*scale)
    rate_matrix[4,6] <- as.numeric(rates[26]*scale)
    rate_matrix[4,7] <- as.numeric(rates[27]*scale)
    rate_matrix[4,8] <- as.numeric(rates[28]*scale)
    rate_matrix[5,1] <- as.numeric(rates[29]*scale)                   
    rate_matrix[5,2] <- as.numeric(rates[30]*scale)
    rate_matrix[5,3] <- as.numeric(rates[31]*scale)
    rate_matrix[5,4] <- as.numeric(rates[32]*scale)
    rate_matrix[5,6] <- as.numeric(rates[33]*scale)
    rate_matrix[5,7] <- as.numeric(rates[34]*scale)
    rate_matrix[5,8] <- as.numeric(rates[35]*scale)
    rate_matrix[6,1] <- as.numeric(rates[36]*scale)                   
    rate_matrix[6,2] <- as.numeric(rates[37]*scale)
    rate_matrix[6,3] <- as.numeric(rates[38]*scale)
    rate_matrix[6,4] <- as.numeric(rates[39]*scale)
    rate_matrix[6,5] <- as.numeric(rates[40]*scale)
    rate_matrix[6,7] <- as.numeric(rates[41]*scale)
    rate_matrix[6,8] <- as.numeric(rates[42]*scale)
    rate_matrix[7,1] <- as.numeric(rates[43]*scale)                  
    rate_matrix[7,2] <- as.numeric(rates[44]*scale)
    rate_matrix[7,3] <- as.numeric(rates[45]*scale)
    rate_matrix[7,4] <- as.numeric(rates[46]*scale)
    rate_matrix[7,5] <- as.numeric(rates[47]*scale)
    rate_matrix[7,6] <- as.numeric(rates[48]*scale)
    rate_matrix[7,8] <- as.numeric(rates[49]*scale)
    rate_matrix[8,1] <- as.numeric(rates[50]*scale)                   
    rate_matrix[8,2] <- as.numeric(rates[51]*scale)
    rate_matrix[8,3] <- as.numeric(rates[52]*scale)
    rate_matrix[8,4] <- as.numeric(rates[53]*scale)
    rate_matrix[8,5] <- as.numeric(rates[54]*scale)
    rate_matrix[8,6] <- as.numeric(rates[55]*scale)
    rate_matrix[8,7] <- as.numeric(rates[56]*scale)}
  
  for(i in 1:Q_dim){
    rate_matrix[i,i] <- -(sum(rate_matrix[i,]))}       # the rates along the diagonal should sum each row to zero
  return(rate_matrix)
}