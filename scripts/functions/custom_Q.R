custom_Q = function(rates){
  rates <- unlist(rates)
  rate_matrix <- matrix(0, nrow = 4, ncol = 4)         # matrix is populated with zeroes to begin
  colnames(rate_matrix) <- c("X_S", "Lu_S", "X_P", "Lu_P")
  rownames(rate_matrix) <- c("X_S", "Lu_S", "X_P", "Lu_P")
  rate_matrix[1,2] <- rates[1]*.001                   # rates are multiplied by .0001 to counteract scaling in BT
  rate_matrix[1,3] <- rates[2]*.001
  rate_matrix[2,1] <- rates[3]*.001
  rate_matrix[2,4] <- rates[4]*.001
  rate_matrix[3,1] <- rates[5]*.001
  rate_matrix[3,4] <- rates[6]*.001
  rate_matrix[4,2] <- rates[7]*.001
  rate_matrix[4,3] <- rates[8]*.001
  for(i in 1:4){
    rate_matrix[i,i] <- -(sum(rate_matrix[i,]))}       # the rates along the diagonal should sum each row to zero
  return(rate_matrix)
}