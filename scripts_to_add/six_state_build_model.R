data <- read.csv("git/data/no_endo_lung_data.csv")                              # load in tip data
BT_ML_tree <- read.tree(file = "git/trees/edited_trees/maxLH_BT_vis_tree.tre")  # load in ML_tree
treedata <- data[data$Taxa %in% BT_ML_tree$tip.label,]; rm(data)                # trim tip data to tree and remove extra file
rownames(treedata) <- treedata$Taxa
trait_col_names <- c("ecology", "lung", "Spec_lotic")
state_labels <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")

cols <- c("gold", "purple", "#99CCFF", "blue", "red", "#EBFFE1")

par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plot(BT_ML_tree, show.tip.label = TRUE, cex = .2, tip.color = cols[treedata$six_state[match(BT_ML_tree$tip.label, treedata$Taxa)]], type = "f")


source("git/scripts/extra_scripts/make_simmap_data_six_state_function.R")

mod <- matrix(c(0, 1, 2, 0, 3, 0,
                4, 0, 0, 5, 0, 0, 
                6, 0, 0, 7, 0, 8, 
                0, 9, 10, 0, 0, 0, 
                11, 0, 0, 0, 0, 12, 
                0, 0, 13, 0, 14, 0), 6)

simmap_matrix_six_state <- make_simmap_data_six_state(treedata, BT_ML_tree, state_labels, trait_col_names)


t <- make.simmap(BT_ML_tree, simmap_matrix_six_state, model = mod, nsim = 10)

cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "white"), 
                 c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P"))


plotSimmap(t[[3]], cols2, pts=F, ftype="off", fsize = .2, 
           outlin = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)





