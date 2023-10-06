
source("lung_loss_git/scripts/functions/Bayestraits_simmap.R")

load(file = "bayestraits_exports/master_dep_export.Rdata")
load(file = "bayestraits_exports/master_indep_export.Rdata")
load(file = "bayestraits_exports/master_six_state_dep_export.Rdata")
load(file = "bayestraits_exports/master_six_state_indep_export.Rdata")
Nsim <- 10
data <- read.csv(file = "lung_loss_git/processed_data/lung_data/no_endo_lung_data.csv")
trees <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/tree_set.nex")

four_state_state = c("X_S", "Lu_S", "X_P", "Lu_P")
four_state_col = c("ecology", "lung")

six_state_state = c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")
six_state_col = c("ecology", "lung", "Spec_lotic")

set.seed(123)


t1 <- BayesTraits_simmap(master_dep, Nsim, data, model = "dependent", 
                         trees, state_labels = four_state_state, 
                         trait_col_names = four_state_col)
t2 <- BayesTraits_simmap(master_indep, Nsim, data, model = "independent", 
                         trees, state_labels = four_state_state, 
                         trait_col_names = four_state_col)
t3 <- BayesTraits_simmap(master_six_state_dep_export, Nsim, data, 
                         model = "six_state", trees, state_labels = six_state_state, 
                         trait_col_names = six_state_col)
t4 <- BayesTraits_simmap(master_six_state_indep_export, Nsim, data, 
                         model = "six_state", trees, state_labels = six_state_state, 
                         trait_col_names = six_state_col)


rm(list = c("four_state_col", "four_state_state", "six_state_col", "six_state_state", 
"Nsim", "BayesTraits_simmap", "custom_Q", "make_simmap_data", "master_dep",
"master_indep", "master_six_state_dep_export", "master_six_state_indep_export"))






cols <- setNames(c("#990099", "#3399FF", "#FF3333","#EBFFE1", "grey"), c("X_S", "Lu_S", "X_P", "Lu_P"))
cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#EBFFE1"), 
                  c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P"))

#tester of one run to see if colors look good
plotSimmap(t1[[1]], cols, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
plotSimmap(t2[[1]], cols, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
plotSimmap(t3[[1]], cols2, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
plotSimmap(t4[[2]], cols2, pts=F, ftype="off", fsize = .2, 
           outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
