rm(list = ls())

{
BT_output <- read.csv("bayestraits/output/processed_logs/dep_allzero_1.txt", sep = "\t")
summary <- BT_output[, 1:max(grep("q", colnames(BT_output)))]; rm(BT_output)
losses <- c("q43", "q21")
gains <- c("q34", "q12")
cols <- match(c(losses, gains), colnames(summary))
}

{BT_output2 <- read.csv("bayestraits/output/processed_logs/eight_state_res_2.txt", sep = "\t")
summary2 <- BT_output2[, 1:max(grep("q", colnames(BT_output2)))]; rm(BT_output2)
losses <- c("q65", "q31", "q42", "q87")
gains <- c("q56", "q13", "q24", "q78")
cols2 <- match(c(losses, gains), colnames(summary2))
}





par(mfrow = c(1,2))

library(vioplot)
{pdf(file = "figures/Fig_4_A.pdf", width = 2.195, height = 2.6)
  par(mar = c(1,1,1,1))
  m <- max(max(summary[,cols]), max(summary2[,cols2]))
  vioplot(log(summary[,cols]*.1+1), xaxt = "n", yaxt = "n", colMed = "black",
          ylim = c(0,log(1+m+.2)))
  segments(.5,0,4.5,0, col = "red", lty = 2)
  points(1:4, lapply(log(summary[,cols]*.1+1), "median"), pch = 21, cex = 1.4, bg = "black")
  points(1:4, lapply(log(summary[,cols]*.1+1), "median"), pch = 21, cex = 1.2, bg = "white")
  axis(side = 4, at = log(c(1,2,6,11,21,51)), labels = c(0,1,5,10,20,50),
       tck = .025, lwd.ticks = 2, lwd = 2)
  axis(side = 4, at = log(c(1,2,6,11,21,51)), labels = c(0,1,5,10,20,50),
       tck = -.01, lwd.ticks = 2, lwd = 1)
  dev.off()}

{pdf(file = "figures/Fig_4_B.pdf", width = 4, height = 2.6)
par(mar = c(1,1,1,1))
m <- max(max(summary[,cols]), max(summary2[,cols2]))
vioplot(log(summary2[,cols2]+1), xaxt = "n", yaxt = "n", colMed = "black",
        ylim = c(0,log(1+m+.2)))
segments(.5,0,8.5,0, col = "red", lty = 2)
points(1:8, lapply(log(summary2[,cols2]+1), "median"), pch = 21, cex = 1.2, bg = "white")
axis(side = 2, at = log(c(1,2,6,11,21,51)), labels = c(0,1,5,10,20,50),
     tck = .01, lwd.ticks = 2, lwd = 2)
axis(side = 2, at = log(c(1,2,6,11,21,51)), labels = c(0,.01,.05,.10,20,50),
     tck = -.025, lwd.ticks = 2, lwd = 1)

dev.off()}


library(scales)
source("lung_loss_git/scripts/functions/summary_posterior_function.R")
#a <- summary_posterior(summary, grep("q", colnames(summary)))
#a[,2]
#1.5/(rescale(a[,1], to = c(.25, 12), from = range(a[,1])))*100
#
#rescale(a[,1], to = c(.25, 12), from = range(a[,1]))



a <- summary_posterior(summary, grep("q", colnames(summary)))
b <- summary_posterior(summary2, grep("q", colnames(summary2)))

a[,1] <- a[,1]*.1

l <- rescale(c(a[,1],b[,1]), to = c(.25, 12), from = range(c(a[,1],b[,1])))

1.5/l*100
















