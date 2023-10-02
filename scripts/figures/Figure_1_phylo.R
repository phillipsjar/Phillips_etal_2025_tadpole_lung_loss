library(phytools)
library(ape)

{data_original <- read.csv("data/Full_matrix.csv")                         #full dataset with all columns
keep_cols <- c("Family", "Genus", "Taxa", "Jetz_name", "guild", "ecology", "lung")
data <- data_original[!(is.na(data_original$lung)),]  #full dataset with lung data
data <- data[,match(keep_cols, colnames(data))] # keep only relevant columns
rm(list = c("data_original", "keep_cols"))}

consensus_tree_7239 <- read.tree("trees/Jetz_pyron_Consensus_7239.tre")
time_tree_4062 <- read.tree(file = "Trees/maxLH_time_tree.tre")

for(i in 1:dim(data)[1]){ 
  if(data$Taxa[i] %in% time_tree_4062$tip.label & data$Jetz_name[i] == ""){      
    data$Jetz_name[i] <- data$Taxa[i]}   # if the real name is in the smaller tree, use it
  if(!(data$Taxa[i] %in% time_tree_4062$tip.label) & data$Taxa[i] %in% consensus_tree_7239$tip.label & data$Jetz_name[i] == ""){
    data$Jetz_name[i] <- data$Taxa[i]    # if the real name is not in the smaller tree, but is in the bigger tree, use that
  }}

data_4062_full <- data[data$Jetz_name %in% time_tree_4062$tip.label,] # now only include taxa found in the maxLH tree

new_tree = function(tree, data){
  require(geiger)
  require(ape)
  rownames(data)<-data$Jetz_name
  name.check(tree,data)->overlap
  drop.tip(tree,c(overlap$tree_not_data, overlap$data_not_tree)) -> trimmed_tree
  trimmed_tree$tip.label[match(data$Jetz_name, trimmed_tree$tip.label)] <- data$Taxa
  return(trimmed_tree)
}

plot_4062_tree <- new_tree(time_tree_4062, data_4062_full )

simmap_mat <- matrix(0, ncol = 2, nrow = (dim(data_4062_full)[1]+1))
colnames(simmap_mat) <- c("lungless", "lunged")
rownames(simmap_mat) <- c(data_4062_full$Taxa, "leio_anc")
simmap_mat[(dim(data_4062_full)[1]+1),] <- c(0,1)
for(i in 1:dim(data_4062_full)[1]){
  if(data_4062_full$lung[i] == 0){simmap_mat[i,1] <- 1}
  if(data_4062_full$lung[i] == 1){simmap_mat[i,2] <- 1}}

plot_4062_tree <- bind.tip(plot_4062_tree, "leio_anc", edge.length=0, where=getMRCA(plot_4062_tree, c("Ascaphus_truei", "Leiopelma_hochstetteri")), position=0)



test <- make.simmap(plot_4062_tree, simmap_mat, nsim=100, model = "ARD", pi=c(0,1))

test2 <- densityMap(test, plot=FALSE)
plot(test2, type = "fan", ftype="off", lwd=2.5, 
     legend = TRUE, outline = TRUE)

n <- length(test2$cols)
pal1 <- colorRampPalette(c("grey", "red"))(n)
test2$cols[1:n]<-pal1

pdf(file = "figures/lung_tree_full.pdf", height = 7, width = 7, bg = "transparent")
plot(test2, type = "fan", ftype="off", lwd=2.5, 
     legend = FALSE, outline = TRUE)
dev.off()





