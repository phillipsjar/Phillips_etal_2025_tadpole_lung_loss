############### full script for all analyses and visualization in Phillips et al. 2023 - Phylogenetics of Lung loss in anuran tadpoles
### R version 

# set working directory
setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2024")
#setwd("/Users/WomackLab/desktop/Jack/final_lungloss_dataspace")

# uploading data
{
  data_original <- read.csv("raw_data/Full_data_matrix.csv")                         #full dataset with all columns
  keep_cols <- c("Family", "Genus", "Taxa", "Portik", "guild", "ecology", "terrestrial", "lung", "Spec_lotic")
  data_original <- data_original[,match(keep_cols, colnames(data_original))] # keep only relevant columns
  
####################################################

data_original$eight_state <- rep(NA)
data_original$eight_state[which(data_original$ecology == 0 & data_original$lung == 0 & data_original$Spec_lotic == 0)] <- 1
data_original$eight_state[which(data_original$ecology == 0 & data_original$lung == 0 & data_original$Spec_lotic == 1)] <- 2
data_original$eight_state[which(data_original$ecology == 0 & data_original$lung == 1 & data_original$Spec_lotic == 0)] <- 3
data_original$eight_state[which(data_original$ecology == 0 & data_original$lung == 1 & data_original$Spec_lotic == 1)] <- 4
data_original$eight_state[which(data_original$ecology == 1 & data_original$lung == 0)] <- 5
data_original$eight_state[which(data_original$ecology == 1 & data_original$lung == 1)] <- 6
data_original$eight_state[which(data_original$terrestrial == 1 & data_original$lung == 0)] <- 7
data_original$eight_state[which(data_original$terrestrial == 1 & data_original$lung == 1)] <- 8


data_original$state <- rep(NA)
data_original$state[which(data_original$ecology == 0 & data_original$lung == 0)] <- 1
data_original$state[which(data_original$ecology == 0 & data_original$lung == 1)] <- 2
data_original$state[which(data_original$ecology == 1 & data_original$lung == 0)] <- 3
data_original$state[which(data_original$ecology == 1 & data_original$lung == 1)] <- 4

data <- data_original[which(data_original$lung %in% c(0,1,"-")),]  

rm(list = c("data_original", "keep_cols"))}

table(data$eight_state)  

data[which(data$eight_state == 8),]

#data <- data_vis

######################################################
################# summarize data #####################
######################################################

# number of taxa with lung data:
dim(data)[1]
# 513 species

# 51 families represented
length(unique(data$Family))
# taxa distribution within families
table(data$Family)

# 249 genera represented
length(unique(data$Genus))
# taxa distribution within genera
table(data$Genus)


######################################################
################## upload trees ######################
######################################################
require(ape)
# MaxLH tree from Portik et al (2022) 
ML_tree <- read.tree("lung_loss_git/trees/original_trees/TreePL/Rooted_Anura_bestTree.tre")

# set of 100 trees from the posterior distribution of tree space in Portik et al (2022)
tree_set <- read.tree("lung_loss_git/trees/original_trees/TreePL/TreePL-Rooted_Anura_bootstraps.tre")

######################################################
############## match datasets to trees ###############
######################################################

data$tree_names <- rep(NA)

for(i in 1:dim(data)[1]){ 
  if(data$Taxa[i] %in% ML_tree$tip.label){      
    data$tree_names[i] <- data$Taxa[i]}   # if taxa are present in the tree, include them
  if(!(data$Taxa[i] %in% ML_tree$tip.label) & data$Portik[i] %in% ML_tree$tip.label){
    data$tree_names[i] <- data$Portik[i]}}
    
#data[is.na(data$tree_names),c(1,2,3,4)]

dim(data)
#dim(data[which(data$terrestrial == 1),c(3,6,7,8)])
#remove 18 terrestrial taxa
data_aqu <- data[which(data$terrestrial != 1),]  #remove terrestrial taxa but keep taxa with unknown ecology
dim(data_aqu)


data_terr <- data[!(is.na(data$eight_state)),]

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
Full_tree <- new_tree(ML_tree, data_terr[data_terr$tree_names %in% ML_tree$tip.label,] )
#lose taxa that can't fit in tree

### Bayestraits analysis I - only no terrestrial species for a binary comparison

aqu_tree <- new_tree(ML_tree, data_aqu[data_aqu$tree_names %in% ML_tree$tip.label,] )
#no terrestrial taxa and lose taxa that can't be placed in tree


# BT analysis I - set of 100 trees (no terrestrial)
{tree_set_aqu <- tree_set #dummy tree set to be written over
  for(i in 1:100){
    tree_set_aqu[[i]] <- new_tree(tree_set[[i]], data_aqu[data_aqu$tree_names %in% tree_set[[i]]$tip.label,])
    tree_set_aqu[[i]]$node.label <- NULL
    if((i/10) == round((i/10))){print((100 - i)*.1)}} # counter down to zero as loop finishes
}

# full tree set of 100 trees (all possible taxa)
{tree_set_full <- tree_set #dummy tree set to be written over
for(i in 1:100){
  tree_set_full[[i]] <- new_tree(tree_set[[i]], data_terr[data_terr$tree_names %in% tree_set[[i]]$tip.label,])
  tree_set_full[[i]]$node.label <- NULL
  if((i/10) == round((i/10))){print((100 - i)*.1)}} # counter down to zero as loop finishes
}

#write nexus files to bayestraits folder 

write.nexus(Full_tree, file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_full.nex", translate = TRUE)
write.nexus(aqu_tree, file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_aqu.nex", translate = TRUE)
write.nexus(tree_set_full, file = "lung_loss_git/bayestraits_trees_data/trees/tree_set_full.nex", translate = TRUE)
write.nexus(tree_set_aqu, file = "lung_loss_git/bayestraits_trees_data/trees/tree_set_aqu.nex", translate = TRUE)

#write tree files (with node labels) to another folder for later visualization
write.tree(Full_tree, file = "lung_loss_git/trees/edited_trees/All_taxa_vis_tree.tre") #includes endotrophs
write.tree(aqu_tree, file = "lung_loss_git/trees/edited_trees/maxLH_BT_vis_tree.tre")        #removes endotrophs



rm(list = c("tree_set", "ML_tree", "i", "new_tree", 
            "tree_set_aqu", "tree_set_full", "Full_tree", "aqu_tree"))

write.csv(data, file = "lung_loss_git/processed_data/lung_data/full_data.csv")
write.csv(data_aqu, file = "lung_loss_git/processed_data/lung_data/aqu_lung_data.csv")

rm(list = c("data", "data_aqu"))


####################################################################################################
# to automatically generate Bayestraits data files and ancestral reconstruction code for 
# command files, see BT_data_set.R in extra_scripts

####################################################################################################

# this script completes set up, creating trees for bayestraits analyses and organizing data
