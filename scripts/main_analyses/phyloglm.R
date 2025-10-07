setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2024")

# upload Bayestraits trees (used in 4-state analysis, so only binary lotic/lentic)
require(ape)

# upload lung and ecology data
data_full <- read.csv("lung_loss_git/processed_data/lung_data/full_data.csv")
full_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_full.nex")



data_full <- data_full[data_full$Taxa %in% full_tree$tip.label,]
rownames(data_full) <- data_full$Taxa

data_aqua <- data_full[which(data_full$ecology != 2),]
data_full$ecology[which(data_full$ecology == 2)] <- NA




#simple phyloglm:
require(phylolm)
packageVersion("phylolm")

mod <- phyloglm(lung ~ ecology, data_aqua, full_tree, method = c("logistic_MPLE"),
                 start.beta=NULL, start.alpha=NULL,
                 boot = 1000, full.matrix = TRUE)

summary(mod)
# ecology recovered as statistically significant (p = 0.03325) and estimated effect size = 2.6677e


mod2 <- phyloglm(lung ~ Spec_lotic, data_aqua, full_tree, method = c("logistic_MPLE"),
                start.beta=NULL, start.alpha=NULL,
                boot = 1000, full.matrix = TRUE)


summary(mod2)



mod3 <- phyloglm(lung ~ ecology+Spec_lotic+terrestrial, data_full, full_tree, method = c("logistic_MPLE"),
                 start.beta=NULL, start.alpha=NULL,
                 boot = 1000, full.matrix = TRUE)

summary(mod3)

mod4 <- phyloglm(lung ~ terrestrial, data_full, full_tree, method = c("logistic_MPLE"),
                 start.beta=NULL, start.alpha=NULL,
                 boot = 1000, full.matrix = TRUE)

summary(mod4)


