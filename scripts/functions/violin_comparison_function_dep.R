violin_comparison_dep = function(data1, data2, data3, name){
  cols   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                       
  labels <- c("q12 (gain lungs in S)", "q13 (S to P no lungs)", "q21 (lose lungs in S)",
              "q24 (S to P w/ lungs)", "q31 (P to S no lungs)", "q34 (gain lungs in P)",
              "q42 (P to S w/ lungs)", "q43 (lose lungs in P)")
  require(vioplot)
  par(mfrow=c(3,1))
  par(xpd = TRUE)
  par(mar = (c(2,4.1,1.1,2.1)))
  vioplot(data1[,8:15], lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
          colMed = "transparent", col = cols, main = paste(name,"1"),
          ylab = "Rate value")
  vioplot(data2[,8:15], lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
          colMed = "transparent", col = cols, main = paste(name,"2"),
          ylab = "Rate value")
  par(mar = (c(4.1,4.1,1.1,2.1)))
  vioplot(data3[,8:15], lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
          colMed = "transparent", col = cols, main = paste(name,"3"),
          ylab = "Rate value")
  legend(0,((par("usr")[3]) - (par("usr")[3])*.025), legend = labels, ncol = 4, fill = cols, bty = "n", cex = 1)
}
