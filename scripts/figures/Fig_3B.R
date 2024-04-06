rm(list = ls())

{BT_output <- read.csv("bayestraits/output/processed_logs/dep_mcmc_1.txt", sep = "\t")

summary <- BT_output[, 1:max(grep("q", colnames(BT_output)))]; rm(BT_output)
summary <- summary[50:dim(summary)[1],]

losses <- c("q43", "q21")
gains <- c("q34", "q12")

cols <- match(c(losses, gains), colnames(summary))

library(vioplot)

{pdf(file = "figures/Fig_3_A.pdf", width = 2.195, height = 2.6)
  par(mar = c(1,1,1,1))
  vioplot(summary[,cols], xaxt = "n", yaxt = "n", colMed = "black",
          ylim = c(0,12.5))
  segments(.5,0,4.5,0, col = "red", lty = 2)
  points(1:4, lapply(summary[,cols], "median"), pch = 21, cex = 1.2, bg = "white")
  axis(side = 2, at = seq(0,12,2), labels = rep("",7))
  dev.off()}}

{BT_output <- read.csv("bayestraits/output/processed_logs/eight_state_1.txt", sep = "\t")

summary <- BT_output[, 1:max(grep("q", colnames(BT_output)))]; rm(BT_output)


losses <- c("q65", "q31", "q42", "q87")
gains <- c("q56", "q13", "q24", "q78")

cols <- match(c(losses, gains), colnames(summary))

library(vioplot)

{pdf(file = "figures/Fig_3_B.pdf", width = 4, height = 2.6)
par(mar = c(1,1,1,1))
vioplot(summary[,cols], xaxt = "n", yaxt = "n", colMed = "black", 
        ylim = c(0,10))
segments(.5,0,8.5,0, col = "red", lty = 2)
points(1:8, lapply(summary[,cols], "median"), pch = 21, cex = 1.2, bg = "white")
axis(side = 4, at = seq(0,8,2), labels = rep("",5))
dev.off()}}

rm(list = ls())

