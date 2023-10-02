library(ape)
library(geiger)
library(phytools)

setwd("/Users/jack/desktop/Research/lunglessness/2021_dataspace/final_dataspace_for_paper")

# general workflow and idea
#
#load Portik tree and other trees for data wrangling
{portik_tree <- read.tree("trees/original_trees/Portik_genus_tree.tre")
tree_4062 <- read.tree("trees/original_trees/Jetz_pyron_maxLH_4062.tre")
consensus_tree_7239 <- read.tree("trees/original_trees/Jetz_pyron_Consensus_7239.tre")}

#load data (same code as in main file)
{data_original <- read.csv("data/Full_matrix.csv")                         #full dataset with all columns
keep_cols <- c("Family", "Genus", "Taxa", "Jetz_name", "guild", "ecology", "lung")
data <- data_original[!(is.na(data_original$lung)),]  #full dataset with lung data
data <- data[,match(keep_cols, colnames(data))] # keep only relevant columns
rm(list = c("data_original", "keep_cols"))

  for(i in 1:dim(data)[1]){ 
    if(data$Taxa[i] %in% tree_4062$tip.label & data$Jetz_name[i] == ""){      
      data$Jetz_name[i] <- data$Taxa[i]}   # if the real name is in the smaller tree, use it
    if(!(data$Taxa[i] %in% tree_4062$tip.label) & data$Taxa[i] %in% consensus_tree_7239$tip.label & data$Jetz_name[i] == ""){
      data$Jetz_name[i] <- data$Taxa[i]    # if the real name is not in the smaller tree, but is in the bigger tree, use that
    }}
  data_no_endo <- data[c(which(data$ecology != 2), which(is.na(data$ecology))),]       #remove endotrophs but keep taxa with unknown ecology
  data_in_tree <- data_no_endo[data_no_endo$Jetz_name %in% consensus_tree_7239$tip.label,] # double check we only include taxa in the least restricted tree
  rm(list = c("data_no_endo", "data"))
  }

dim(data_in_tree)

{
data_genera <- as.data.frame(matrix(nrow = length(portik_tree$tip.label), ncol = 3))
colnames(data_genera) <- c("Genus", "Jetz_name", "Taxa")
rownames(data_genera) <- portik_tree$tip.label
data_genera$Genus <- portik_tree$tip.label
}

for(i in 1:dim(data_genera)[1]){
# some outgroup taxa are shared between Portik and consensus tree, so just pull those over as is
if(data_genera$Genus[i] %in% consensus_tree_7239$tip.label){
  data_genera$Jetz_name[i] <- data_genera$Taxa[i] <- data_genera$Genus[i]
}
in_set <- grep(data_genera$Genus[i], data_in_tree$Genus)
in_tree <- grep(data_genera$Genus[i], consensus_tree_7239$tip.label)

if(length(in_set) > 0){
  sample <- sample(1:length(in_set), 1, replace = FALSE)  # get a random number between 1 and N taxa in my dataset   to index with
  data_genera$Jetz_name[i] <- data_in_tree$Jetz_name[in_set][sample]
  data_genera$Taxa[i] <- data_in_tree$Taxa[in_set][sample]}

if(length(in_set) == 0 && length(in_tree) > 0){
  sample <- sample(1:length(in_tree), 1, replace = FALSE)  # get a random number between 1 and N taxa in my dataset   to index with
  data_genera$Jetz_name[i] <- consensus_tree_7239$tip.label[in_tree][sample]
  data_genera$Taxa[i] <- NA}}




keep <- data_genera[!(is.na(data_genera$Jetz_name)),]

trimmed_portik_tree <- keep.tip(portik_tree, keep$Genus)

trimmed_portik_tree$tip.label[match(keep$Genus, trimmed_portik_tree$tip.label)] <- keep$Jetz_name

rm(list = c("tree_4062", "i", "in_tree", "sample", "portik_tree", "in_set", "keep", "data_genera"))

plot(trimmed_portik_tree, cex = .2)

#################################################################################

# time calibrate
# tree is trimmed as little possible to preserve tips for time calibration

#custom function for time_calibration - see script for details
source("scripts/extra_scripts/time_calibrate_function.R")

# this may take a couple minutes, but should be relatively fast
portik_time_tree <- time.calibrate(trimmed_portik_tree, consensus_tree_7239)
write.tree(portik_time_tree, file = "trees/edited_trees/portik_time_tree_full.tre")

rm(list = c("consensus_tree_7239", "trimmed_portik_tree", "time.calibrate"))

#################################################################################

portik_time_tree <- read.tree(file = "trees/edited_trees/portik_time_tree_full.tre")
# now we trim the tree to only taxa we have data for
# not elegant at all, but other solutions were buggy

length(portik_time_tree$tip.label)
#start with 400 tips


list <- rep(NA, length(portik_time_tree$tip.label))
names(list) <- portik_time_tree$tip.label
for(i in 1:length(list)){
  if(portik_time_tree$tip.label[i] %in% data_in_tree$Taxa){ #if taxa are in both: keep em
    list[i] <- 1}
  if(!(portik_time_tree$tip.label[i] %in% data_in_tree$Taxa) & portik_time_tree$tip.label[i] %in% data_in_tree$Jetz_name){
    list[i] <- 2}
}



par(mfrow = c(1,2))
trimmed_portik_time_tree <- keep.tip(portik_time_tree, names(c(list[which(list == 1)], list[which(list == 2)]))) 
plot(trimmed_portik_time_tree, cex = .2)

trimmed_portik_time_tree$tip.label <- data_in_tree$Taxa[match(trimmed_portik_time_tree$tip.label, data_in_tree$Jetz_name)] #change tip labels
plot(trimmed_portik_time_tree, cex = .2, direction = "leftwards")


rm(list = c("portik_time_tree", "data_in_tree", "i", "list"))
#################################################################################

write.tree(trimmed_portik_time_tree, file = "trees/edited_trees/portik_tree_edit.tre")

#rm(trimmed_portik_time_tree)
par(mfrow = c(1,1))
plot(portik_time_tree, type = "f", cex = .4, font = 3)









