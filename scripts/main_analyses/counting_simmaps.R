rm(list = ls())

eight_bestmod_1 <- read.csv("bayestraits/output/processed_logs/eight_state_bestmod_1.txt", sep = "\t")
eight_res_3 <- read.csv("bayestraits/output/processed_logs/eight_state_res_3.txt", sep = "\t")

ML_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/MaxLh_tree_full.nex") 
trees <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/tree_set_full.nex")

state_labels <- eight_state_state <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P", "X_terr", "Lu_terr")
trait_col_names <- eight_state_col <- c("ecology", "lung", "Spec_lotic", "terrestrial")
mode = "eight_state"
model = "eight_state"

Nruns <- 10
data <- read.csv(file = "lung_loss_git/processed_data/lung_data/full_data.csv")
data <- data[data$Taxa %in% ML_tree$tip.label,]

source("lung_loss_git/scripts/functions/Bayestraits_simmap.R")
source("lung_loss_git/scripts/functions/simmap_counts_functions.R")

A <- BayesTraits_simmap(eight_bestmod_1, 1000, data, model = "eight_state",
                        trees = trees, trait_col_names = eight_state_col, 
                        state_labels = eight_state_state, pis = c(.25, 0, .25,
                                                                  0, .25, .25, 
                                                                  0, 0))


a <- avg_simmap_counts(A)

B <- BayesTraits_simmap(eight_res_3, 1000, data, model = "eight_state",
                        trees = trees, trait_col_names = eight_state_col, 
                        state_labels = eight_state_state, pis = c(.25, 0, .25,
                                                                  0, .25, .25, 
                                                                  0, 0))


b <- avg_simmap_counts(B)

a
b





