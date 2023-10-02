#check scales of bayestraits dependent and independent models


dep_scale <- read.csv("git/bayestraits/output/processed_logs/scale_test_dep.txt",sep = "\t")
indep_scale <- read.csv("git/bayestraits/output/processed_logs/scale_test_indep.txt",sep = "\t")

dep_scale <- dep_scale[round(dim(dep_scale)[1]*.25):dim(dep_scale)[1],]
indep_scale <- indep_scale[round(dim(indep_scale)[1]*.25):dim(indep_scale)[1],]

summary = function(data,column.names){
  if (sum(dim(data)[2]) != 0){
    cols <- match(column.names, colnames(data))
    par(mfrow = c(dim(data[cols])[2]/4, 4))
    for(i in 1:dim(data[,cols])[2]){
      smoothScatter(data$Iteration, data[,cols[i]], ylim = c(0,max(data[,cols])),
                    cex=4,nr=500, xlab = "iter", 
                    ylab = colnames(data)[cols[i]])}} else{
                      par(mfrow = c(1,1))
                      smoothScatter(data$Iteration, data[,cols],cex=4,nr=500, xlab = "iter", ylab = colnames(rows))
                    }
}

summary(dep_scale, column.names = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"))
summary(indep_scale, column.names = c("alpha1", "alpha2", "beta1", "beta2"))

plot(dep_scale$Iteration, dep_scale$RJRates...Mean)
plot(indep_scale$Iteration, indep_scale$RJRates...Mean)

