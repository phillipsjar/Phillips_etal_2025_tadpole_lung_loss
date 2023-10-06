
source("lung_loss_git/scripts/functions/Bayestraits_simmap.R")

load(file = "bayestraits_exports/master_dep_export.Rdata")
load(file = "bayestraits_exports/master_indep_export.Rdata")
load(file = "bayestraits_exports/master_six_state_dep_export.Rdata")
load(file = "bayestraits_exports/master_six_state_indep_export.Rdata")
Nsim <- 10
Nruns <- 100
data <- read.csv(file = "lung_loss_git/processed_data/lung_data/no_endo_lung_data.csv")
trees <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/tree_set.nex")

four_state_state = c("X_S", "Lu_S", "X_P", "Lu_P")
four_state_col = c("ecology", "lung")

six_state_state = c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")
six_state_col = c("ecology", "lung", "Spec_lotic")

set.seed(123)


{t_indep <- BayesTraits_simmap(master_indep, Nruns, Nsim, data, model = "independent", 
                         trees, state_labels = four_state_state, 
                         trait_col_names = four_state_col)
t_dep <- BayesTraits_simmap(master_dep, Nruns, Nsim, data, model = "dependent", 
                         trees, state_labels = four_state_state, 
                         trait_col_names = four_state_col)
t_six_indep <- BayesTraits_simmap(master_six_state_indep_export, Nruns, Nsim, data, 
                         model = "six_state", trees, state_labels = six_state_state, 
                         trait_col_names = six_state_col)
t_six_dep <- BayesTraits_simmap(master_six_state_dep_export, Nruns, Nsim, data, 
                         model = "six_state", trees, state_labels = six_state_state, 
                         trait_col_names = six_state_col)}



rm(list = c("four_state_col", "four_state_state", "six_state_col", "six_state_state", 
"Nsim", "BayesTraits_simmap", "custom_Q", "make_simmap_data", "master_dep",
"master_indep", "master_six_state_dep_export", "master_six_state_indep_export"))






cols <- setNames(c("#990099", "#3399FF", "#FF3333","#EBFFE1", "grey"), c("X_S", "Lu_S", "X_P", "Lu_P"))
cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#EBFFE1"), 
                  c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P"))
require(phytools)
#tester of one run to see if colors look good
{plotSimmap(t_indep[[1]], cols, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
title(main="independent model",font.main=3,
      line=-1)}

{plotSimmap(t_dep[[1]], cols, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
title(main="dependent model",font.main=3,
      line=-1)}

{plotSimmap(t_six_indep[[1]], cols2, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
title(main="six_independent model",font.main=3,
      line=-1)}

{plotSimmap(t_six_dep[[2]], cols2, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
title(main="six_dependent model",font.main=3,
      line=-1)}



t_indep_count <- countSimmap(t_indep[[3]])
t_dep_count <- countSimmap(t_dep[[1]])
t_six_indep_count <- countSimmap(t_six_indep[[1]])
t_six_dep_count <- countSimmap(t_six_dep[[1]])



avg_simmap_counts = function(simmap){
  require(phytools)
  N <- length(simmap)
  list <- vector(mode = "list", length = N)             
  for(i in 1:N){
    list[i] <- list(countSimmap(simmap[[i]])$Tr)
  }
  Y <- do.call(cbind, list)
  Y <- array(Y, dim=c(dim(list[[1]]), length(list)))
  counts <- round(apply(Y, c(1,2), mean, na.rm = TRUE),1)
  rownames(counts) <- rownames(countSimmap(simmap[[1]])$Tr)
  colnames(counts) <- colnames(countSimmap(simmap[[1]])$Tr)
  return(counts)
}

sd_simmap_counts = function(simmap){
  require(phytools)
  N <- length(simmap)
  list <- vector(mode = "list", length = N)             
  for(i in 1:N){
    list[i] <- list(countSimmap(simmap[[i]])$Tr)
  }
  Y <- do.call(cbind, list)
  Y <- array(Y, dim=c(dim(list[[1]]), length(list)))
  counts <- round(apply(Y, c(1,2), sd, na.rm = TRUE),1)
  rownames(counts) <- rownames(countSimmap(simmap[[1]])$Tr)
  colnames(counts) <- colnames(countSimmap(simmap[[1]])$Tr)
  return(counts)
}

avg_simmap_counts(t_indep)
sd_simmap_counts(t_indep)

avg_simmap_counts(t_dep)
avg_simmap_counts(t_six_indep)
avg_simmap_counts(t_six_dep)












#     [,1] [,2]
#[1,]    4    5
#[2,]    6    7

















