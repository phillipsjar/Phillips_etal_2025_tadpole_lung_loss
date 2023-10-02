### proportions of lunged and lungless genera in different habitats for figure 5

{lunged_lentic <- unique(data$Genus[which(data$lung == 1 & data$ecology == 1)])
lungless_lentic_bufo <- unique(data$Genus[which(data$lung == 0 & data$ecology == 1 & data$Family == "Bufonidae")])
lungless_lentic_one_bufo <- c(unique(data$Genus[which(data$lung == 0 & data$ecology == 1 & data$Family != "Bufonidae")]), "Bufonidae")

lunged_lotic_fos <- unique(data$Genus[which(data$lung == 1 & data$guild == "lotic, fossorial")])
lungless_lotic_fos <- unique(data$Genus[which(data$lung == 0 & data$guild == "lotic, fossorial")])

lunged_lotic_suc <- unique(data$Genus[which(data$lung == 1 & data$guild == "lotic, suctorial" | data$lung == 1 & data$guild == "lotic, gastromyzophorous")])
lungless_lotic_suc_bufo <- unique(data$Genus[which(data$lung == 0 & data$guild == "lotic, suctorial" & data$Family == "Bufonidae" | data$lung == 0 & data$guild == "lotic, gastromyzophorous" & data$Family == "Bufonidae")])
lungless_lotic_suc_no_bufo <- unique(data$Genus[which(data$lung == 0 & data$guild == "lotic, suctorial" & data$Family != "Bufonidae" | data$lung == 0 & data$guild == "lotic, gastromyzophorous" & data$Family != "Bufonidae" )])


lunged_lotic_semi <- unique(data$Genus[which(data$lung == 1 & data$guild == "lotic, semiterrestrial")])
lungless_lotic_semi <- unique(data$Genus[which(data$lung == 0 & data$guild == "lotic, semiterrestrial")])

lunged_lotic_gen <- unique(data$Genus[which(data$lung == 1 & data$ecology == 0 & data$guild != "lotic, fossorial" & 
                                              data$guild != "lotic, suctorial" & data$guild != "lotic, gastromyzophorous" & 
                                              data$guild != "lotic, semiterrestrial")])
lungless_lotic_gen_bufo <- unique(data$Genus[which(data$lung == 0 & data$ecology == 0 & data$guild != "lotic, fossorial" & 
                                                     data$guild != "lotic, suctorial" & data$guild != "lotic, gastromyzophorous" & 
                                                     data$guild != "lotic, semiterrestrial" & data$Family == "Bufonidae")])
lungless_lotic_gen_no_bufo <- unique(data$Genus[which(data$lung == 0 & data$ecology == 0 & data$guild != "lotic, fossorial" & 
                                                        data$guild != "lotic, suctorial" & data$guild != "lotic, gastromyzophorous" & 
                                                        data$guild != "lotic, semiterrestrial" & data$Family != "Bufonidae")])

lunged_endo <- unique(data$Genus[which(data$lung == 1 & data$ecology == 2)])
lungless_endo_bufo <- unique(data$Genus[which(data$lung == 0 & data$ecology == 2 & data$Family == "Bufonidae")])
lungless_endo_no_bufo <- unique(data$Genus[which(data$lung == 0 & data$ecology == 2 & data$Family != "Bufonidae")])}


{bar_matrix <- matrix(nrow = 6, ncol = 6)
  rownames(bar_matrix) <- c("lunged count", "lungless count_bufo", "no_bufo_count", "lunged percentage", "lungless percentage bufo", "no_bufo_percentage")
  colnames(bar_matrix) <- c("Lentic", "Lotic, non-specialized", "Lotic, suctorial", "Lotic, fossorial", "Lotic, semiterrestrial", "Endotroph")
  
  bar_matrix[1,1] <- as.numeric(length(lunged_lentic))
  bar_matrix[2,1] <- as.numeric(length(lungless_lentic_bufo))
  bar_matrix[3,1] <- as.numeric(length(lungless_lentic_one_bufo))
  bar_matrix[1,2] <- as.numeric(length(lunged_lotic_gen))
  bar_matrix[2,2] <- as.numeric(length(lungless_lotic_gen_bufo))
  bar_matrix[3,2] <- as.numeric(length(lungless_lotic_gen_no_bufo))
  bar_matrix[1,3] <- as.numeric(length(lunged_lotic_suc))
  bar_matrix[2,3] <- as.numeric(length(lungless_lotic_suc_bufo))
  bar_matrix[3,3] <- as.numeric(length(lungless_lotic_suc_no_bufo))
  bar_matrix[1,4] <- as.numeric(length(lunged_lotic_fos))
  bar_matrix[2,4] <- as.numeric(length(lungless_lotic_fos))
  bar_matrix[3,4] <- as.numeric(0)
  bar_matrix[1,5] <- as.numeric(length(lunged_lotic_semi))
  bar_matrix[2,5] <- as.numeric(length(lungless_lotic_semi))
  bar_matrix[3,5] <- as.numeric(0)
  bar_matrix[1,6] <- as.numeric(length(lunged_endo))
  bar_matrix[2,6] <- as.numeric(length(lungless_endo_bufo))
  bar_matrix[3,6] <- as.numeric(length(lungless_endo_no_bufo))
}

for(i in 1:6){
  n <- sum(as.numeric(bar_matrix[1:3,i]))
  bar_matrix[4,i] <- as.numeric(bar_matrix[1,i])/n
  bar_matrix[5,i] <- as.numeric(bar_matrix[2,i])/n
  bar_matrix[6,i] <- as.numeric(bar_matrix[3,i])/n
}


#{pdf(file = "figures/stacked_plots_fig_5.pdf", width = 1.5, height = 4)
#par(mar = c(.1,.1,.1,.1))
barplot(bar_matrix[4:6,], col = c("grey", "red", "red"), xaxt = "n", yaxt = "n", 
        horiz = T, space = .5)
#dev.off()}





######################################################
############# time calibrate maxLH tree ##############
######################################################

# custom function built for this purpose:      (see attached script in extra_scripts for details)
#   time.calibrate = function(tree1, cal_tree){
#     require(dispRity)
#     require(ape)
#     require(phytools)

source("scripts/extra_scripts/time_calibrate_function.R")

#time_tree_4062 <- time.calibrate(tree_4062, consensus_tree_7239)
#write.tree(time_tree_4062, file = "trees/edited_trees/maxLH_time_tree.tre")

rm(list = c("time.calibrate"))


time_tree_4062 <- read.tree("trees/edited_trees/maxLH_time_tree.tre")


######## Figure 3D

{pdf(file = "figures/RJ_proportion_plot.pdf", bg = "transparent", width = 6.5, height = 3.5)
  par(mar = c(5.1,4.1,4.1,1.1))
  par(xpd = TRUE)
  cols <- c("#E8E8E8", "grey", "#393939")
  BF <- c(2, 5, 10)
  p = BF/(1+BF)
  bp <- barplot(master_RJ, ylim = c(0,1), col = cols, beside = TRUE)
  segments(par("usr")[1], p[1], par("usr")[2], p[1], lty = 2, col = "red")
  segments(par("usr")[1], p[2], par("usr")[2], p[2], lty = 3, col = "red")
  segments(par("usr")[1], p[3], par("usr")[2], p[3], lty = 4, col = "red")
  barplot(master_RJ, ylim = c(0,1), col = cols, beside = TRUE, add = TRUE)
  text(bp,(master_RJ+.05), round(master_RJ, 3), cex = .75)
  #text(par("usr")[1] + .1, p[1]+.015, "weak support", adj = c(0,0), cex = .75)
  #text(par("usr")[1] + .1, p[2]+.015, "support", adj = c(0,0), cex = .75)
  #text(par("usr")[1] + .1, p[3]+.015, "strong support", adj = c(0,0), cex = .75)
  
  dev.off()}

