# build a subset of MaxLH tree
BT_4062_tree <- read.tree(file = "Trees/edited_trees/maxLH_BT_vis_tree.tre")

data <- read.csv(file = "data/no_endo_lung_data.csv")
data <- data[data$Taxa %in% BT_4062_tree$tip.label,]

genera <- unique(data$Genus)
data2 <- data[1:(length(genera)*3),]

for(i in 1:length(genera)){
  iter <- (i*3)-2
  data2[c(iter, (iter+1), (iter+2)), ] <- data[(which(data$Genus == genera[i])[1:3]),]
}

rm(list = c("data", "genera", "i", "iter"))

data2 <- data2[!(is.na(data2$Taxa)),]

tree <- BT_4062_tree
data <- data2

new_tree = function(tree, data){
  require(geiger)
  require(ape)
  rownames(data)<-data$Taxa
  name.check(tree,data)->overlap
  drop.tip(tree,c(overlap$tree_not_data, overlap$data_not_tree)) -> trimmed_tree
  trimmed_tree$tip.label[match(data$Taxa, trimmed_tree$tip.label)] <- data$Taxa
  return(trimmed_tree)
}

sub_tree <- new_tree(BT_4062_tree, data2)




write.tree(sub_tree, file = "Trees/edited_trees/maxLH_sub_vis_tree.tre")



