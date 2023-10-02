####                visualize posterior for figure 4A
library(HDInterval)
library(vioplot)
# calculate highest density intervals of the master posterior
#BT_data <- master_maxLH_dep
#states_OI <- c("q12", "q13", "q15", "q21", "q24", "q31", "q34", "q36", "q42", "q43", "q51", "q56", "q63", "q65")
states_OI <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
BT_data <- master_dep_export[,match(states_OI, colnames(master_dep_export))]*.001

posterior_intervals <- hdi(BT_data, credMass = .95)

#lungless

{pdf(file = "figures/rate_violin_plots.pdf", bg = "transparent", width = 5, height = 3.5)
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                  
  labels <- c("gain lungs (lotic)", "lotic to lentic (no lungs)", "lose lungs (lotic)",
              "lotic to lentic (lunged)", "lentic to lotic (no lungs)", "gain lungs (lentic)",
              "lentic to lotic (lunged)", "lose lungs (lentic)")
  par(mar = c(5.1,4.1,4.1,2.1)) #reset default margins
  par(mfrow=c(1,1))
  par(xpd = TRUE)
  vioplot(BT_data, lineCol = "transparent", rectCol = "transparent",  xaxt = "n", yaxt = "n",
          colMed = "transparent", col = "darkgrey", main = "",
          ylab = "")
  #segments(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3], col = "white", lwd = 2)  
  text(x = (1:length(labels))+((par("usr")[2]-par("usr")[1])*.05), y = par("usr")[3]-(par("usr")[4]-par("usr")[3])*.05,
       labels = labels, xpd = NA, adj = 1, srt = 25, cex = .65)
  
  axis(2,cex.axis=.75)
  
  #legend(0,((par("usr")[3]) - (par("usr")[3])*.025), legend = labels, ncol = 3, fill = "darkgray", bty = "n", cex = .75)
  for(i in 1:dim(posterior_intervals)[2]){
    #segments(i, 0, i, par("usr")[3]-(par("usr")[4]-par("usr")[3])*.05, lty = "dashed")
    segments(i, posterior_intervals[2,i], i,posterior_intervals[1,i], lwd = 1.75, col = "white")
    segments(i, posterior_intervals[2,i], i,posterior_intervals[1,i], lwd = 1.5, col = "black")
    segments(i - ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[1,i], 
             i + ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[1,i], lwd = 1.75, col = "white")
    segments(i - ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[2,i], 
             i + ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[2,i], lwd = 1.75, col = "white")
    segments(i - ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[1,i], 
             i + ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[1,i], lwd = 1.5, col = "black")
    segments(i - ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[2,i], 
             i + ((par("usr")[2] - par("usr")[1])/(length(labels)*4)), posterior_intervals[2,i], lwd = 1.5, col = "black")
    points(i, median(BT_data[,i]), pch = 21, bg = "white", cex = 1.25)
  }
  dev.off()}


avgs_dep <- sapply(master_dep_export[,grep("q12", colnames(master_dep_export)):(grep("q12", colnames(master_dep_export))+7)], "mean");
avgs_indep <- sapply(master_indep_export[,grep("alpha1", colnames(master_indep_export)):(grep("alpha1", colnames(master_indep_export))+3)], "mean");




require(scales)
avgs_scaled <- round(rescale(avgs_maxLH_dep,c(.5,7.5),c(0,5)),2)

round(avgs_scaled,2)