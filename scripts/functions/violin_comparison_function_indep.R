violin_comparison_dep = function(data1, data2, data3, name){
  cols   <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
                       
  columns <- c("alpha1", "beta1", "alpha2", "beta2")
  labels <- c("alpha1 (S to P)", "beta1 (P to S)", "alpha2 (gain lungs)", "beta2 (lose lungs)")
  require(vioplot)
  par(mfrow=c(3,1))
  par(xpd = TRUE)
  par(mar = (c(2,4.1,1.1,2.1)))
  vioplot(data1[,match(columns, colnames(data1))], lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
          colMed = "transparent", col = cols, main = paste(name,"1"),
          ylab = "Rate value")
  vioplot(data2[,match(columns, colnames(data1))], lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
          colMed = "transparent", col = cols, main = paste(name,"2"),
          ylab = "Rate value")
  par(mar = (c(4.1,4.1,1.1,2.1)))
  vioplot(data3[,match(columns, colnames(data1))], lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
          colMed = "transparent", col = cols, main = paste(name,"3"),
          ylab = "Rate value")
  legend(0,((par("usr")[3]) - (par("usr")[3])*.025), legend = labels, ncol = 4, fill = cols, bty = "n", cex = 1)
}
