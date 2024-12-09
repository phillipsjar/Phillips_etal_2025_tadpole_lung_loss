combine_posteriors = function(data1, data2, data3){
  master <- as.data.frame(matrix(ncol = dim(data1)[2], nrow = (dim(data1)[1] + dim(data2)[1] + dim(data3)[1])))
  colnames(master) <- colnames(data1)
  master$run_number <- c(rep(1,dim(data1)[1]), rep(2,dim(data2)[1]), rep(3,dim(data3)[1]))
  
  master[which(master$run_number == 1),1:dim(data1)[2]] <- data1[1:dim(data1)[2]]
  master[which(master$run_number == 2),1:dim(data2)[2]] <- data2[1:dim(data2)[2]]
  master[which(master$run_number == 3),1:dim(data3)[2]] <- data3[1:dim(data3)[2]]
  master <- as.data.frame(master)
  return(master)
}