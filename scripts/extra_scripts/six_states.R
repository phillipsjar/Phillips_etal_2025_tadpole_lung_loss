library(ape)
library(vioplot)
library(HDInterval)

setwd("/Users/jack/desktop/Research/lunglessness/2021_dataspace/final_dataspace_for_paper")

mod <- matrix(c(0, 1, 2, 0, 3, 0,
         4, 0, 0, 5, 0, 0, 
         6, 0, 0, 7, 0, 8, 
         0, 9, 10, 0, 0, 0, 
         11, 0, 0, 0, 0, 12, 
         0, 0, 13, 0, 14, 0), 6)




{
BT_4062_tree <- read.tree(file = "Trees/edited_trees/maxLH_BT_vis_tree.tre")
BT_portik_tree <- read.tree(file = "Trees/edited_trees/portik_BT_vis_tree.tre")
data_no_endo <- read.csv(file = "data/no_endo_lung_data.csv")
data_4062 <- data_no_endo[data_no_endo$Taxa %in% BT_4062_tree$tip.label,]
data_portik <- data_no_endo[data_no_endo$Taxa %in% BT_portik_tree$tip.label,]
}


#six_state <- read.csv("Bayestraits/results/six_state_test_dep.txt",sep = "\t")
#six_state <- six_state[round(dim(six_state)[1]*.25) : dim(six_state)[1],]

load(file = "bayestraits/results/maxLH_dep_sixstate_master_export.Rdata")
load(file = "bayestraits/results/tree_set_dep_sixstate_master_export.Rdata")
load(file = "bayestraits/results/portik_dep_sixstate_master_export.Rdata")





states_OI <- c("q12", "q13", "q15", "q21", "q24", "q31", "q34", "q36", "q42", "q43", "q51", "q56", "q63", "q65")
six_state_ML <- master_export_6state_maxLH_dep[,match(states_OI, colnames(master_export_6state_maxLH_dep))]
six_state_tree_set <- master_export_6state_tree_set_dep[,match(states_OI, colnames(master_export_6state_tree_set_dep))]
six_state_portik <- master_export_6state_portik_dep[,match(states_OI, colnames(master_export_6state_portik_dep))]




rm(list = c("master_export_6state_maxLH_dep", "master_export_6state_tree_set_dep", "master_export_6state_portik_dep"))
  
  
labels <- c("q12 (specialize in S no lungs)", "q13 (gain lungs in S unspec)", "q15 (S to P no lungs)",
            "q21 (de-specialize in S no lungs)", "q24 (gain lungs in S spec)", 
            "q31 (lose lungs in S unspec)", "q34 (specialize in S lunged)", "q36 (S to P with lungs)",
            "q42 (lose lungs in S spec)", "q43 (de-specialize in S with lungs)",
            "q51 (P to S no lungs)", "q56 (gain lungs in P)", "q63 (P to S with lungs)", "q65 (lose lungs in P)")



{six_state <- six_state_tree_set
posterior_intervals <- hdi(six_state, credMass = .95)
par(xpd = TRUE)
vioplot(six_state, lineCol = "transparent", rectCol = "transparent",  xaxt = "n", 
        colMed = "transparent", main = "Estimated Rate Values",
        ylab = "Rate value")
text(x = (1:length(labels))+((par("usr")[2]-par("usr")[1])*.05), y = par("usr")[3]-(par("usr")[4]-par("usr")[3])*.05,
     labels = labels, xpd = NA, adj = 1, srt = 25, cex = .65)
for(i in 1:dim(posterior_intervals)[2]){
  if(posterior_intervals[2,i] == 0){
    arrows(i, 0, i,-((par("usr")[4] -  par("usr")[2])/750), angle = 90, lwd = 1.75, col = "black")
    arrows(i, 0, i,-((par("usr")[4] -  par("usr")[2])/750), angle = 90, lwd = 1.5, col = "white", lty = "dotted")}
  if(median(six_state[,i]) != posterior_intervals[1,i]){
    arrows(i, median(six_state[,i]), i,posterior_intervals[1,i], angle = 90, lwd = 1.75, col = "black")}
  if(median(six_state[,i]) != posterior_intervals[2,i]){
    arrows(i, median(six_state[,i]), i,posterior_intervals[2,i], angle = 90, lwd = 1.75, col = "black")}  
  if(median(six_state[,i]) != posterior_intervals[1,i]){
    arrows(i, median(six_state[,i]), i,posterior_intervals[1,i], angle = 90, lwd = 1.5, col = "white", lty = "dotted")}
  if(median(six_state[,i]) != posterior_intervals[2,i]){
    arrows(i, median(six_state[,i]), i,posterior_intervals[2,i], angle = 90, lwd = 1.5, col = "white", lty = "dotted")}  
  points(i, median(six_state[,i]), pch = 21, bg = "white", cex = 1.5)
}}

rm(list = c("posterior_intervals", "six_state", "data_no_endo", "i", "labels", "states_OI"))
############


set.seed(123)


six_state <- six_state_tree_set
  
random_runs <- sample(1:(dim(six_state)[1]), 100, replace = FALSE)  # get a sample of 100 random runs from the posterior

rate_values <- vector(mode = "list", length = 100)             # empty list for the rate values for each random run

for(i in 1:100){
  rate_values[i] <- list(six_state[random_runs,][i,])
}



# phytools does not have a native dependent model, but this can easily be dealt with by treating the 
# dependent model as a constrained 4 state model, with each combination of the two binary states treated
# as one of four states (lungless stream; lunged, stream; lungless, pond; lunged, pond)

#rate_values

states_OI <- c("q12", "q13", "q15", "q21", "q24", "q31", "q34", "q36", "q42", "q43", "q51", "q56", "q63", "q65")

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

custom_Q_six(rate_values[1])  #example Q matrix from one run of the posterior

############### make simmap objects ####################
library(phytools)

##### set up simmap for pseudo-dependent analysis

# data_4062 loaded above.
# simmap requires a specific way data are input, which we create below with a quick function

make_simmap_data_onetrait = function(data, tree, state_labels){
  if(is.matrix(data) | is.data.frame(data)){
    print("only supply one vector of state data, not full dataset")
    stop()
  }
  n_states <- length(unique(data[!(is.na(data))]))
  trait_matrix <- matrix(0, nrow = length(data), ncol = n_states)           #empty matrix (for simmap with all zeroes)
  if(is.null(names(data))){
    names(data) <- tree$tip.label
    print("No taxon names supplied, so order is assumed to be the same as phy$tip.labels")
  }
  data <- data[names(data) %in% tree$tip.label] #trim data to tree
  
  names(trait_matrix) <- names(data)
  colnames(trait_matrix) <- state_labels
  # this loop populates the matrix to indicate which state each taxon is in, allowing for unknown habitat but not lung status
  for(i in 1:dim(trait_matrix)[1]){
    for(j in 1:n_states){
      if(!(is.na(data[i]))){
      if(data[i] == j){trait_matrix[i,j] <- 1}}
      if(is.na(data[i])){trait_matrix[i,1:n_states] <- 1/n_states}
    }}
  return(trait_matrix)}

labels <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")



tree <- read.tree("Trees/edited_trees/MaxLH_vis_tree.tre")
data <- read.csv(file = "data/no_endo_lung_data_vis.csv")
rownames(data) <- data$Taxa
trait_col_names <- c("ecology", "lung", "Spec_lotic")
state_labels <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")


source("scripts/extra_scripts/make_simmap_data_six_state_function.R")



simmap_matrix_six_state_ML <- make_simmap_data_six_state(data, tree, state_labels, trait_col_names)




trait_data_ML <- data_4062$six_state
trait_data_portik <- data_portik$six_state
names(trait_data_ML) <- data_4062$Taxa
names(trait_data_portik) <- data_portik$Taxa
simmap_matrix_six_state_ML <- make_simmap_data_onetrait(trait_data_ML, BT_4062_tree, labels)
simmap_matrix_six_state_portik <- make_simmap_data_onetrait(trait_data_portik, BT_portik_tree, labels)

rm(list = c("data_4062", "data_portik", "six_state", "six_state_ML", 
            "six_state_portik", "six_state_tree_set", "i", "random_runs", "states_OI", 
            "trait_data_ML", "trait_data_portik", "make_simmap_data_onetrait"))

#test_simmap <- make.simmap(BT_4062_tree,simmap_matrix_six_state, nsim=100, pi=rep((1/6),6))
#master_simmap <- test_simmap
#rm(test_simmap)



master_simmap_ML <- make.simmap(BT_4062_tree,simmap_matrix_six_state_ML, nsim=1, Q = custom_Q_six(rate_values[1], labels), pi=c(0,0,0,0,0,1))



master_simmap_portik <- make.simmap(BT_portik_tree,simmap_matrix_six_state_portik, nsim=100, Q = custom_Q_six(rate_values[1], labels), pi=c(0,0,0,0,0,1))

#now fill each index of the simmap object with different Q matrices
for(i in 1:100){
  master_simmap_ML[[i]] <- make.simmap(BT_4062_tree,simmap_matrix_six_state_ML, nsim=1, Q = custom_Q_six(rate_values[i], labels), pi=rep((1/6),6))
  print(100-i) # count down for sanity
}
#now for the portik tree
for(i in 1:100){
  master_simmap_portik[[i]] <- make.simmap(BT_portik_tree,simmap_matrix_six_state_portik, nsim=1, Q = custom_Q_six(rate_values[i], labels), pi=rep((1/6),6))
  print(100-i) # count down for sanity
}

cols <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "white"), 
                 c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P"))


plot(master_simmap_ML, colors = cols, outline = TRUE, type = "fan", fsize = .2)
plot(master_simmap_portik[[1]], colors = cols, outline = TRUE, type = "fan", fsize = .2)

?plot.simmap

plotTree(BT_4062_tree,ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=5, type = "fan", color = "black")

par(fg="transparent",lend=1)
plotTree(BT_4062_tree,ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
         lwd=2.5,color="white",add=TRUE, type = "fan")
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency

for(i in 1:100) {plot(master_simmap[[i]],
                      colors=sapply(cols, make.transparent,alpha=0.01),
                      add=TRUE, lwd=2.5, ftype="off", fsize=0.5, offset=0.5, type = "fan")
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}






par(fg="black")

legend(x="bottomleft",c("lungless generalist in streams", "lunged generalist in streams", "lungless specialist in streams", "lunged specialist in streams", "lungless in ponds", "lunged in ponds"),pch=22,
       pt.bg=cols,pt.cex=1.5,bty="n",cex=0.7)

par(mfrow=c(3,3))

for(i in 1:9){
  plot(master_simmap[[i]], colors = cols, ftype = "off")
}

plot(master_simmap[[1:9]], colors = cols, ftype = "off")



final_transitions <- matrix(NA, nrow = 100, ncol = 14)
colnames(final_transitions) <- c("q12 (specialize in S no lungs)", "q13 (gain lungs in S unspec)", "q15 (S to P no lungs)",
                                 "q21 (de-specialize in S no lungs)", "q24 (gain lungs in S spec)", 
                                 "q31 (lose lungs in S unspec)", "q34 (specialize in S lunged)", "q36 (S to P with lungs)",
                                 "q42 (lose lungs in S spec)", "q43 (de-specialize in S with lungs)",
                                 "q51 (P to S no lungs)", "q56 (gain lungs in P)", "q63 (P to S with lungs)", "q65 (lose lungs in P)")
for(i in 1:100){
    sim <- make.simmap(BT_4062_tree,simmap_matrix_six_state, nsim=100, Q = custom_Q_six(rate_values[i], labels), pi=rep((1/6),6))
    a <- as.data.frame(countSimmap(sim)$Tr)
    a <- a[c(3, 4, 6, 8, 11, 14, 17, 19, 21, 22, 26, 31, 34, 36)]
    final_transitions[i,] <- sapply(a[1:dim(a)[2]], "mean")
print(paste(100-i, "remaining"))}  
  
a
  
  
head(final_transitions)








