#######################################################################################################
######      Stochastic Character Mapping in Phytools using Bayestraits Rate Matrices      #############

# bayestraits does not have a native stochastic mapping tool, so we will use phytools to do the mapping,
# using the BT posterior to generate Q_matrices while incorporating uncertainty

#              Draw rate matrices from 100 different runs of posterior
BayesTraits_simmap = function(BT_output, Nruns, Nsim, data, model, trees, state_labels, trait_col_names){
  require(phytools)
if(class(trees) == "phylo"){
  stop("currently implemented for Bayestraits runs that use multiple trees")
}
  
random_runs <- sample(1:(dim(BT_output)[1]), Nruns, replace = FALSE)  # get a sample of 100 random runs from the posterior
tree_numbers <- BT_output$Tree.No[random_runs]

rate_values <- vector(mode = "list", length = Nruns)             # empty list for the rate values for each random run

if(model == "dependent" & length(which(colnames(BT_output) == "q12")) == 0){
  stop("dependent (default) model chosen, but q12 rate not available")}

if(model == "dependent"){
for(i in 0:(Nruns-1)){
  start <- which(colnames(BT_output) == "q12")
  run <- random_runs[i+1]
  rate_values[i+1] <- list(BT_output[(run),start:(start+7)])
}
}

if(model == "independent"){
  for(i in 0:(Nruns-1)){
    start <- which(colnames(BT_output) == "alpha1")
    run <- random_runs[i+1]
    rate_values[i+1] <- list(BT_output[(run),start:(start+3)])
  }
}

if(model == "six_state"){
  for(i in 0:(Nruns-1)){
    start <- which(colnames(BT_output) == "q12")
    run <- random_runs[i+1]
    rate_values[i+1] <- list(BT_output[(run),start:(start+29)])
  }
}

source("lung_loss_git/scripts/functions/custom_Q.R")
#rate_values[[1]]
#custom_Q(rate_values[[1]], state_labels, model)

source("lung_loss_git/scripts/functions/make_simmap_data_function.R")


if(model == "dependent" | model == "independent"){mode = "double binary"}
if(model == "six_state"){mode = model}

sim_data <- make_simmap_data(data = data, tree = trees[[1]], mode = mode, state_labels = state_labels, trait_col_names = trait_col_names)


master_simmap <- vector(mode = "list", length = Nruns*Nsim)             # empty list to be filled with the descendents of each node of ML tree
class(master_simmap) <- "Multiphylo"

#now fill each index of the simmap object with different Q matrices
if(model == "dependent" | model == "independent"){
  for(j in 1:(Nruns*Nsim)){
    for(i in 1:Nruns){
      master_simmap[[j]] <- make.simmap(trees[[tree_numbers[i]]], sim_data, nsim = 1,
                                   Q = custom_Q(rate_values[[i]], model = model, 
                                   state_labels = state_labels), pi=rep(.25,4))
    }}
  print(Nruns-i) # count down for sanity
}



if(model == "six_state"){
  for(j in 1:(Nruns*Nsim)){
    for(i in 1:Nruns){
      master_simmap[[j]] <- make.simmap(trees[[tree_numbers[i]]], sim_data, nsim = 1,
                            Q = custom_Q(rate_values[[i]], model = model, 
                              state_labels = state_labels), pi=rep((1/6),6))}
    print(Nruns-i) # count down for sanity
}}


return(master_simmap)
}
