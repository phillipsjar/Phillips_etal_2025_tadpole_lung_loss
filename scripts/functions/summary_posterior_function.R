
summary_posterior = function(data, cols){
require(HDInterval)

summary <- data[, cols]

rate_sum <- matrix(NA, ncol = 6, nrow = length(cols))
rownames(rate_sum) <- colnames(summary)
colnames(rate_sum) <- c("mean", "median", "SD", "lower_interval", "upper_interval", "count")

i = 1
for(i in 1:dim(rate_sum)[1]){
  rate_sum[i,1] <- mean(summary[,match(rownames(rate_sum)[i], colnames(summary))])
  rate_sum[i,2] <- median(summary[,match(rownames(rate_sum)[i], colnames(summary))])
  rate_sum[i,3] <- sd(summary[,match(rownames(rate_sum)[i], colnames(summary))])
  rate_sum[i,4] <- hdi(summary[,match(rownames(rate_sum)[i], colnames(summary))])[1]
  rate_sum[i,5] <- hdi(summary[,match(rownames(rate_sum)[i], colnames(summary))])[2]
  rate_sum[i,6] <- i}
return(rate_sum)
}
