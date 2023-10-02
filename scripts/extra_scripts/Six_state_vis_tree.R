## function to create a dataset phytools is comfortable with two binary traits allowing missing data

vis_tree <- read.tree("Trees/edited_trees/MaxLH_vis_tree.tre")
dataset <- read.csv(file = "data/no_endo_lung_data_vis.csv")
labels <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")

load(file = "bayestraits/results/tree_set_dep_sixstate_master_export.Rdata")
states_OI <- c("q12", "q13", "q15", "q21", "q24", "q31", "q34", "q36", "q42", "q43", "q51", "q56", "q63", "q65")
six_state_tree_set <- master_export_6state_tree_set_dep[,match(states_OI, colnames(master_export_6state_tree_set_dep))]

six_state <- six_state_tree_set

random_runs <- sample(1:(dim(six_state)[1]), 100, replace = FALSE)  # get a sample of 100 random runs from the posterior

rate_values <- vector(mode = "list", length = 100)             # empty list for the rate values for each random run

for(i in 1:100){
  rate_values[i] <- list(six_state[random_runs,][i,])
}


custom_Q_six = function(rates, state_labels = c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")){
  rates <- unlist(rates)
  rate_matrix <- matrix(0, nrow = 6, ncol = 6)         # matrix is populated with zeroes to begin
  colnames(rate_matrix) <- state_labels
  rownames(rate_matrix) <- state_labels
  rate_matrix[1,2] <- rates[1]*.0001                   # rates are multiplied by .0001 to counteract scaling in BT
  rate_matrix[1,3] <- rates[2]*.0001
  rate_matrix[1,5] <- rates[3]*.0001
  rate_matrix[2,1] <- rates[4]*.0001
  rate_matrix[2,4] <- rates[5]*.0001
  rate_matrix[3,1] <- rates[6]*.0001
  rate_matrix[3,4] <- rates[7]*.0001
  rate_matrix[3,6] <- rates[8]*.0001
  rate_matrix[4,2] <- rates[9]*.0001
  rate_matrix[4,3] <- rates[10]*.0001
  rate_matrix[5,1] <- rates[11]*.0001
  rate_matrix[5,6] <- rates[12]*.0001
  rate_matrix[6,3] <- rates[13]*.0001
  rate_matrix[6,5] <- rates[14]*.0001
  for(i in 1:6){
    rate_matrix[i,i] <- -(sum(rate_matrix[i,]))}       # the rates along the diagonal should sum each row to zero
  return(rate_matrix)
}
rm(list = c("master_export_6state_tree_set_dep", "states_OI", "i", "labels", "six_state", "six_state_tree_set"))


source("scripts/extra_scripts/make_simmap_data_six_state_function.R")

trait_col_names <- c("ecology", "lung", "Spec_lotic")
state_labels <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")
rownames(dataset) <- dataset$Taxa

simmap_matrix_six_state_ML <- make_simmap_data_six_state(dataset, vis_tree, state_labels, trait_col_names)

Q1 <- custom_Q_six(rate_values[1], state_labels)

Q1[1,3] <- 0
Q1[2,4] <- 0
Q1[1,3] <- 0
Q1[6,3] <- 0.001624238

master_simmap <- make.simmap(vis_tree,simmap_matrix_six_state_ML, nsim=10, Q = Q1, pi=c(0,0,0,0,0,1))
cols <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "white"), 
                 c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P"))

plotTree(vis_tree,ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=4, type = "fan", color = "black", part = .995)

par(fg="transparent",lend=1)
plotTree(vis_tree,ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
         lwd=2.5,color="white",add=TRUE, type = "fan", part = .995)
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency

for(i in 1:10) {plot(master_simmap[[i]], part = .995,
                      colors=sapply(cols, make.transparent,alpha=0.1),
                      add=TRUE, lwd=2.5, ftype="off", fsize=0.5, offset=0.5, type = "fan")
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}

tiplabels(text = rownames(simmap_matrix_six_state_ML), 
          tip = match(vis_tree$tip.label, rownames(simmap_matrix_six_state_ML)))





rownames(simmap_matrix_six_state_ML)

rownames(simmap_matrix_six_state_ML)[match(vis_tree$tip.label, rownames(simmap_matrix_six_state_ML))]


cols <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#DEFADC"), 
                 c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P"))


plotSimmap(master_simmap[[3]], colors = cols, outline = TRUE, type = "fan", 
           lwd = 1.5, fsize = .2, part = .995)



labels <- c("q12 (specialize in S no lungs)", "q13 (gain lungs in S unspec)", "q15 (S to P no lungs)",
            "q21 (de-specialize in S no lungs)", "q24 (gain lungs in S spec)", 
            "q31 (lose lungs in S unspec)", "q34 (specialize in S lunged)", "q36 (S to P with lungs)",
            "q42 (lose lungs in S spec)", "q43 (de-specialize in S with lungs)",
            "q51 (P to S no lungs)", "q56 (gain lungs in P)", "q63 (P to S with lungs)", "q65 (lose lungs in P)")


?plotSimmap


sort(tree_4062$tip.label[grep("Rhacophorus", tree_4062$tip.label)])



