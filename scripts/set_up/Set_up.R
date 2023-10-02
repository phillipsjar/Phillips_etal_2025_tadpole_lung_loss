############### full script for all analyses and visualization in Phillips et al. 2023 - Phylogenetics of Lung loss in anuran tadpoles
### R version 

# set working directory
setwd("/Users/jack/desktop/Research/lunglessness/2021_dataspace/final_dataspace_for_paper")
setwd("/Users/WomackLab/desktop/Jack/final_lungloss_dataspace")


# uploading data
{
  data_original <- read.csv("git/data/Full_data_matrix.csv")                         #full dataset with all columns
  keep_cols <- c("Family", "Genus", "Taxa", "Portik", "guild", "ecology", "lung", "Spec_lotic")
  data_original <- data_original[,match(keep_cols, colnames(data_original))] # keep only relevant columns
  
####################################################

data_original$six_state <- rep(NA)
data_original$six_state[which(data_original$ecology == 0 & data_original$lung == 0 & data_original$Spec_lotic == 0)] <- 1
data_original$six_state[which(data_original$ecology == 0 & data_original$lung == 0 & data_original$Spec_lotic == 1)] <- 2
data_original$six_state[which(data_original$ecology == 0 & data_original$lung == 1 & data_original$Spec_lotic == 0)] <- 3
data_original$six_state[which(data_original$ecology == 0 & data_original$lung == 1 & data_original$Spec_lotic == 1)] <- 4
data_original$six_state[which(data_original$ecology == 1 & data_original$lung == 0)] <- 5
data_original$six_state[which(data_original$ecology == 1 & data_original$lung == 1)] <- 6

data_original$state <- rep(NA)
data_original$state[which(data_original$ecology == 0 & data_original$lung == 0)] <- 1
data_original$state[which(data_original$ecology == 0 & data_original$lung == 1)] <- 2
data_original$state[which(data_original$ecology == 1 & data_original$lung == 0)] <- 3
data_original$state[which(data_original$ecology == 1 & data_original$lung == 1)] <- 4

data <- data_original[which(data_original$lung %in% 0:1),]  #full dataset with lung data
#data_vis <- data_original[which(data_original$lung %in% c(0,1,"?")),]  
data_vis$lung[which(data_vis$lung == "?")] <- NA
  
rm(list = c("data_original", "keep_cols"))}
  
#data <- data_vis

######################################################
################# summarize data #####################
######################################################

# number of taxa with lung data:
#dim(data)[1]
# 513 species

# 47 families represented
#length(unique(data$Family))
# taxa distribution within families
#table(data$Family)

# 249 genera represented
#length(unique(data$Genus))
# taxa distribution within genera
#table(data$Genus)


######################################################
################## upload trees ######################
######################################################
require(ape)
# MaxLH tree from Portik et al (2022) 
ML_tree <- read.tree("git/trees/original_trees/TreePL/Rooted_Anura_bestTree.tre")

# set of 100 trees from the posterior distribution of tree space in Portik et al (2022)
tree_set <- read.tree("git/trees/original_trees/TreePL/TreePL-Rooted_Anura_bootstraps.tre")

######################################################
############## match datasets to trees ###############
######################################################

data$tree_names <- rep(NA)

for(i in 1:dim(data)[1]){ 
  if(data$Taxa[i] %in% ML_tree$tip.label){      
    data$tree_names[i] <- data$Taxa[i]}   # if taxa are present in the tree, include them
  if(!(data$Taxa[i] %in% ML_tree$tip.label) & data$Portik[i] %in% ML_tree$tip.label){
    data$tree_names[i] <- data$Portik[i]}}
    
#remove 10 endotroph taxa
data_no_endo <- data[c(which(data$ecology != 2), which(is.na(data$ecology))),]       #remove endotrophs but keep taxa with unknown ecology

#data_no_endo_vis <- data_vis[c(which(data_vis$ecology != 2), which(is.na(data_vis$ecology))),]       #remove endotrophs but keep taxa with unknown ecology


##############################################################
#             Trim trees to those datasets
##############################################################

## A quick function for brevity - trims trees to the dataset using tree_names from above and then changes the tip labels 
## back to true name of sampled taxa.
  
new_tree = function(tree, data){
  require(geiger)
  require(ape)
  rownames(data)<-data$tree_names
  name.check(tree,data)->overlap
  drop.tip(tree,c(overlap$tree_not_data, overlap$data_not_tree)) -> trimmed_tree
  trimmed_tree$tip.label[match(data$tree_names, trimmed_tree$tip.label)] <- data$Taxa
  return(trimmed_tree)
}

### Trim trees for different purposes
# Visualization tree of all sampled taxa that can be placed in tree
Full_sampled_tree <- new_tree(ML_tree, data[data$tree_names %in% ML_tree$tip.label,] )

### Bayestraits trees must only include taxa with binary data (no endotrophs)

# BT tree used for ancestral reconstructions (ML tree with no endotrophs)
BT_ML_tree <- new_tree(ML_tree, data_no_endo[data_no_endo$tree_names %in% ML_tree$tip.label,] )

#tree_set_posterior_100[[1]]

# BT tree set of 100 trees (no endotrophs)
{BT_tree_set_trimmed <- tree_set #dummy tree set to be written over
for(i in 1:100){
  BT_tree_set_trimmed[[i]] <- new_tree(tree_set[[i]], data_no_endo[data_no_endo$tree_names %in% tree_set[[i]]$tip.label,])
  BT_tree_set_trimmed[[i]]$node.label <- NULL
  if((i/10) == round((i/10))){print((100 - i)*.1)}} # counter down to zero as loop finishes
}

#write nexus files to bayestraits folder 

write.nexus(BT_ML_tree, file = "git/bayestraits/trees/maxLH_tree.nex", translate = TRUE)
write.nexus(BT_tree_set_trimmed, file = "git/bayestraits/trees/tree_set.nex", translate = TRUE)

#write tree files (with node labels) to another folder for later visualization
write.tree(Full_sampled_tree, file = "git/trees/edited_trees/All_taxa_vis_tree.tre") #includes endotrophs
write.tree(BT_ML_tree, file = "git/trees/edited_trees/maxLH_BT_vis_tree.tre")        #removes endotrophs



rm(list = c("tree_set", "ML_tree", "i", "new_tree", "BT_tree_set_trimmed"))

write.csv(data, file = "git/data/full_lung_data.csv")
write.csv(data_no_endo, file = "git/data/no_endo_lung_data.csv")

rm(list = c("BT_ML_tree", "data", "data_no_endo", "Full_sampled_tree"))


####################################################################################################
# to automatically generate Bayestraits data files and ancestral reconstruction code for 
# command files, see BT_data_set.R in extra_scripts

####################################################################################################

# this script completes set up, creating trees for bayestraits analyses and organizing data





