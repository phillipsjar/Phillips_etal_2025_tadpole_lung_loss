#Code for supplementary figure 1

rm(list = ls()) #clean up

data_original <- read.csv("raw_data/Full_data_matrix.csv")
keep_cols <- c("Family", "Genus", "Taxa", "Portik", "guild", "ecology", "terrestrial", "lung", "Spec_lotic")

data_original <- data_original[,match(keep_cols, colnames(data_original))] # keep only relevant columns

data_lung <- read.csv("lung_loss_git/processed_data/lung_data/full_data.csv")

data_full <- data_original[which(data_original$guild != "unknown" & data_original$guild != "viviparous" & 
        data_original$guild != "direct development" & data_original$guild != 
        "gastric brooding" & data_original$guild != "mouth brooding"),]
ML_tree <- read.tree("lung_loss_git/trees/original_trees/TreePL/Rooted_Anura_bestTree.tre")


new_tree = function(tree, data){
  require(geiger)
  require(ape)
  rownames(data)<-data$tree_names
  name.check(tree,data)->overlap
  drop.tip(tree,c(overlap$tree_not_data, overlap$data_not_tree)) -> trimmed_tree
  trimmed_tree$tip.label[match(data$tree_names, trimmed_tree$tip.label)] <- data$Taxa
  return(trimmed_tree)
}

tree_table <- matrix(nrow = length(ML_tree$tip.label),
                     ncol = 4)
rownames(tree_table) <- ML_tree$tip.label
colnames(tree_table) <- c("frogs", "DDs", "lung_data", "plot")

tree_table[,1] <- ML_tree$tip.label

for(i in 1:dim(tree_table)[1]){
if(gsub("_.*", "", rownames(tree_table)[i]) %in% data_original$Genus){
  tree_table[i,1] = 1}else{tree_table[i,1] = 0}
  if(gsub("_.*", "", rownames(tree_table)[i]) %in% data_full$Genus){
    tree_table[i,2] = 1}else{tree_table[i,2] = 0}
  if(rownames(tree_table)[i] %in% data_lung$tree_names){
    tree_table[i,3] = 1}else{tree_table[i,3] = 0}
  if(tree_table[i,1] == 1 &
     tree_table[i,2] == 1){tree_table[i,4] = 1}
  if(tree_table[i,1] == 1 &
     tree_table[i,2] == 0){tree_table[i,4] = 2}
  if(tree_table[i,3] == 1){tree_table[i,4] = 3}
  }


names(which(tree_table[,1] == 1))

unique(data_original$Genus)

library(phytools)
#frog_data <- data_original[which(tree_table[,1] == 1),]

frog_tree <- drop.tip(ML_tree, which(tree_table[,1] == 0))
nodes <- rep(NA, frog_tree$Nnode)

for(i in 1:length(nodes)){
desc <- getDescendants(frog_tree, (length(frog_tree$tip.label)+i))
desc <- frog_tree$tip.label[desc[which(desc <= length(frog_tree$tip.label))]]
if("1" %in% tree_table[match(desc, rownames(tree_table)),2]){
  nodes[i] <- 0
}else{nodes[i] <- 1}}

#data_lung$Taxa[which(!(data_lung$Taxa %in% frog_data$Taxa[which(!(is.na(frog_data$lung)))]))]


edges <- rep("black", Nedge(frog_tree))

edges[which(frog_tree$tip.label[frog_tree$edge[,2]] %in% 
              names(which(tree_table[,2] == 0)))] <- "orange"

edges[which(frog_tree$tip.label[frog_tree$edge[,2]] %in% 
                names(which(tree_table[,3] == 1)))] <- "red"




{jpeg(file = "figures/suppl_1.jpeg", width = 7, height = 7, units = "in", res = 800)
par(mar = rep(0,4))
plot(frog_tree, type = "fan", show.tip.label = F, lwd = 5,
     edge.color = edges)
tips1 <- match(names(which(tree_table[,4] == 1)), frog_tree$tip.label)
tiplabels(tip = tips1, pch = 16, cex = .4, col = "black")
tips2 <- match(names(which(tree_table[,4] == 2)), frog_tree$tip.label)
tiplabels(tip = tips2, pch = 16, cex = .4, col = "orange")
tips3 <- match(names(which(tree_table[,4] == 3)), frog_tree$tip.label)
tiplabels(tip = tips3, pch = 16, cex = .4, col = "red")
dev.off()}



length(tips1)
length(tips2)
length(tips3)

length(tips1) +length(tips2) +length(tips3)
















                          
                          
                          