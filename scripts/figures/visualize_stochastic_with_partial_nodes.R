{BT_4062_tree <- read.tree(file = "Trees/edited_trees/maxLH_BT_vis_tree.tre")
data <- read.csv(file = "data/no_endo_lung_data.csv")
genera <- unique(data$Genus)
data2 <- data[1:(length(genera)*3),]
data <- data[data$Jetz_name %in% BT_4062_tree$tip.label,]
for(i in 1:length(genera)){
  iter <- (i*3)-2
  data2[c(iter, (iter+1), (iter+2)), ] <- data[(which(data$Genus == genera[i])[1:3]),]}
rm(list = c("data", "genera", "i", "iter"))
data2 <- data2[!(is.na(data2$Taxa)),]}

sub_tree <- read.tree(file = "Trees/edited_trees/maxLH_sub_vis_tree.tre")
load(file = "data/4062_ASE_R_formatted.R")


source("scripts//figures/Transfer_ASE.R")

nodes2 <- convert_ancestral_state(BT_4062_tree, sub_tree, nodes)

par(mar = c(1,1,1,1))
par(mfrow = c(1,1))
plot(BT_4062_tree, type = "fan", cex = .4)
tips <- length(BT_4062_tree$tip.label)
nodelabels(node = ((tips+1):(2*tips-1)), pie = as.matrix(nodes[,1:4]), cex = .45, piecol = c("magenta", "blue", "red", "green"))
rm(tips)

plot(sub_tree, type = "fan", cex = .4)
nodelabels(node = nodes2$tree2_node_number[!(is.na(nodes2$tree2_node_number))], pie = as.matrix(nodes2[!(is.na(nodes2$tree2_node_number)),1:4]), cex = .45, piecol = c("magenta", "blue", "red", "green"))

uncertain <- which(nodes2$lungless > .25 & nodes2$lungless < .75 | nodes2$lotic > .25 & nodes2$lotic < .75)

plot(sub_tree, type = "fan", cex = .4)
nodelabels(node = nodes2$tree2_node_number[uncertain], pie = as.matrix(nodes2[uncertain,1:4]), cex = .45, piecol = c("magenta", "blue", "red", "green"))




load(file = "data/4062_simmap_100.Rdata")

uncertain1 <- which(nodes$lungless > .25 & nodes$lungless < .75 | nodes$lotic > .25 & nodes$lotic < .75)

cols2 <- setNames(c("purple", "skyblue", "red", "white", "grey"), c("X_S", "Lu_S", "X_P", "Lu_P"))
plot(master_simmap[[1]], cols2, pts=F, ftype="off", fsize = .2, type = "fan", lwd=3, outline = TRUE)
tips <- length(BT_4062_tree$tip.label)
nodelabels(node = as.numeric(rownames(nodes[uncertain1,])) + tips, pie = as.matrix(nodes[uncertain1,1:4]), cex = .45, piecol = c("purple", "skyblue", "red", "white"))

























