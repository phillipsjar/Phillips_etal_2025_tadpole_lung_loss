setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2025")

library(geiger)

rm(list = ls())

data_full <- read.csv("lung_loss_git/processed_data/lung_data/full_data.csv")
full_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_full.nex")
data_corHMM <- data_full[data_full$Taxa %in% full_tree$tip.label,]
dataset_corHMM <- as.data.frame(cbind(data_corHMM$Taxa, data_corHMM$lung, 
                                      data_corHMM$ecology, data_corHMM$Spec_lotic,
                                      data_corHMM$terrestrial)); rm(data_corHMM)

colnames(dataset_corHMM) <- c("Genus_sp", "lung", "ecology", "Spec", "Terr")
dataset_corHMM$ecology[which(dataset_corHMM$Terr == 1)] <- rep(0)
dataset_corHMM$Spec[which(dataset_corHMM$Terr == 1)] <- rep(0)

{
  tip_1 <- getMRCA(full_tree, c("Ascaphus_truei", "Rana_sauteri"))
  tip_2 <- getMRCA(full_tree, c("Heleophryne_regis", "Rana_sauteri"))
  tip_3 <- getMRCA(full_tree, c("Oophaga_pumilio", "Allobates_nidicola"))
  tip_4 <- getMRCA(full_tree, c("Leptodactylus_fuscus", "Adenomera_marmorata"))
  tip_5 <- getMRCA(full_tree, c("Thoropa_miliaris", "Alsodes_monticola"))
}
tips_n <- c("tip_1", "tip_2", "tip_3", "tip_4", "tip_5")
tips <- c(tip_1, tip_2, tip_3, tip_4, tip_5)

for(i in 1:length(tips)){
  full_tree <- bind.tip(full_tree,  tips_n[i], edge.length=0, where = tips[i])
  row <- nrow(dataset_corHMM)+1
  dataset_corHMM[row,] <- c(tips_n[i], "?", "?",0,0)
}

rownames(dataset_corHMM) <- dataset_corHMM$Genus_sp
name.check(full_tree, dataset_corHMM)



write.csv(dataset_corHMM, "test_data.csv", row.names = FALSE)
write.tree(full_tree, "test_tree.tre")









tail(dataset_corHMM)

# Copy your dataset
dataset_corHMM_Q <- dataset_corHMM[,1:3]

# Replace "?" with proper NA
dataset_corHMM_Q[dataset_corHMM_Q == "?"] <- NA

# Convert trait columns to numeric (everything except first column)
for (col in names(dataset_corHMM_Q)[-1]) {
  dataset_corHMM_Q[[col]] <- as.numeric(dataset_corHMM_Q[[col]])
}

# Drop rows with any missing values
dataset_corHMM_Q <- na.omit(dataset_corHMM_Q)

any(is.na(dataset_corHMM_Q))

state_strings <- apply(dataset_corHMM_Q[, -1], 1, paste, collapse = "_")
table(state_strings)
length(unique(state_strings))

str(dataset_corHMM_Q)
head(dataset_corHMM_Q)




# Now this should succeed
auto_matrix <- getStateMat4Dat(dataset_corHMM_Q, collapse = FALSE)





dataset_corHMM_Q <- dataset_corHMM[,1:4]

auto_matrix <- getStateMat4Dat(dataset_corHMM_Q, collapse = T)
auto_matrix <- getStateMat4Dat(dataset_corHMM_Q, collapse = F)




scheme <- c("0_0_0_0", "0_0_1_0", "1_0_0_0", "1_0_1_0", 
            "0_1_0_0", "1_1_0_0", "0_0_0_1", "1_0_0_1")
Q <- auto_matrix$rate.mat

# allow transitions to and from terrestriality

{lungless_terr <- match("0_0_0_1", auto_matrix$legend)
lunged_terr <-    match("1_0_0_1", auto_matrix$legend)
lungless_lent <-  match("0_1_0_0", auto_matrix$legend)
lungless_lot <-   match("0_0_0_0", auto_matrix$legend)
lungless_spec <-  match("0_0_1_0", auto_matrix$legend)
lunged_lent <-    match("1_1_0_0", auto_matrix$legend)
lunged_lot <-     match("1_0_0_0", auto_matrix$legend)
lunged_spec <-    match("1_0_1_0", auto_matrix$legend)}

m <- max(Q)

Q[lungless_terr,lungless_lent] <- m + 1
Q[lungless_lent,lungless_terr] <- m + 2

Q[lungless_terr,lungless_spec] <- m + 3
Q[lungless_spec,lungless_terr] <- m + 4

Q[lunged_terr,lunged_lent] <- m + 5
Q[lunged_lent,lunged_terr] <- m + 6

Q[lunged_terr,lunged_spec] <- m + 7
Q[lunged_spec,lunged_terr] <- m + 8


#set up a basic corhmm model that estimates each parameter independently
full_model <- corHMM(full_tree, dataset_corHMM, rate.cat = 1, 
                          rate.mat = Q, root.p = "maddfitz",
                          collapse = F, lower.bound = 0, node.states = "marginal",
                          nstarts =5, n.cores = 6)




#rename stuff for easier visualization
dimnames(full_model$solution) <- list(auto_matrix$legend, auto_matrix$legend)

#reorder Q matrix to match the scheme for consistency
full_model$solution <- full_model$solution[match(scheme, rownames(full_model$solution)), 
                                    match(scheme, colnames(full_model$solution))]

#plot the model to take a look
plotMKmodel(full_model$solution*100, rate.cat = 1, arrow.scale = 3)


scheme

combined <- paste(dataset_corHMM$lung, dataset_corHMM$ecology, 
                  dataset_corHMM$Spec, dataset_corHMM$Terr, sep = "_")
names(combined) <- dataset_corHMM$Genus_sp

state_matrix <- matrix(0, nrow = length(combined), ncol = 8)
rownames(state_matrix) <- names(combined)                   
colnames(state_matrix) <- scheme                                 

split_state <- function(s) unlist(strsplit(s, "_"))

for(i in seq_along(combined)){
  input <- combined[i]
  
  if(input %in% scheme){
    # Exact match — assign 1 to that column
    state_matrix[i, match(input, scheme)] <- 1
  } else if(grepl("\\?", input)) {
    # Partial / uncertain state
    input_bits <- split_state(input)
    
    # Compare to each scheme entry
    compatible <- sapply(scheme, function(sch) {
      scheme_bits <- split_state(sch)
      all(mapply(function(x, y) x == "?" || x == y, input_bits, scheme_bits))
    })
    
    n_compatible <- sum(compatible)
    if(n_compatible > 0){
      state_matrix[i, compatible] <- 1 / n_compatible
    } else {
      warning(paste("No matches for:", input))
    }
  } else {
    warning(paste("Mismatch or unknown format:", input))
  }
}


Q_simmap <- full_model$solution
for(i in 1:dim(Q_simmap)[1]){
  Q_simmap[i,i] <- -(sum(Q_simmap[i,], na.rm = T))
}
Q_simmap[which(is.na(Q_simmap))] <- rep(0)

Q_simmap

simmap <- make.simmap(full_model$phy, state_matrix, nsim = 10, Q = Q_simmap,
                 state_labels = scheme, 
                 pi= c(.25, 0, .25, 0, .25, .25, 0, 0))

simmap <- lapply(simmap,drop.tip.simmap,tip=tips_n)
class(simmap)<-"multiPhylo"

cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#B9FFB4",
                    "#8B5E3C", "#CC6600"), scheme)

plotSimmap(simmap[[3]], col = cols2, ftype = "off", fsize = .1, 
           type = "fan")

  
col_anc <- c("gold", "#8B5E3C", "purple", "red", "#99CCFF", "#CC6600", "blue", "#B9FFB4")

nodelabels("", pie = full_model$states, piecol = col_anc, cex = .2, frame = "n")




eight_bestmod_1 <- read.csv("bayestraits/output/processed_logs/eight_state_bestmod_1.txt", sep = "\t")
source("lung_loss_git/scripts/functions/summary_posterior_function.R")
rate_sum <- summary_posterior(eight_bestmod_1, grep("q", colnames(eight_bestmod_1)))


scheme <- c("0_0_0_0", "0_0_1_0", "1_0_0_0", "1_0_1_0", 
            "0_1_0_0", "1_1_0_0", "0_0_0_1", "1_0_0_1")
auto <- auto_matrix$legend

Q

# Use which() to get non-zero indices
idx <- which(Q != 0, arr.ind = TRUE)

# Construct labels
Q_labels <- sort(paste0("q", idx[,1], idx[,2]))
n <- gsub("q", "", Q_labels)

result <- do.call(rbind, strsplit(n, ""))
rownames(result) <- n

#these are the output positions

match(auto[as.numeric(result[i,1])], scheme)


for(i in 1:nrow(result)){
  result[i,1] <- match(auto[as.numeric(result[i,1])], scheme)
  result[i,2] <- match(auto[as.numeric(result[i,2])], scheme)
}

result <- apply(result, 1, function(x) paste(x, collapse = ""))

result

BT_solution <- rate_sum[,2]

na <- gsub("q", "", names(BT_solution))

match(result, na)


BT_solution


BT_solution <- BT_solution[match(result, na)]

full_model$solution[which(is.na(full_model$solution))] <- rep(0)












#set up a basic corhmm model that estimates each parameter independently
jump_started <- corHMM(full_tree, dataset_corHMM, rate.cat = 1, 
                     rate.mat = Q, ip = full_model$solution,
                     collapse = TRUE, lower.bound = 0, node.states = "marginal",
                     nstarts =5, n.cores = 6)





?corHMM

#rename stuff for easier visualization
dimnames(full_model$solution) <- list(auto_matrix$legend, auto_matrix$legend)

#reorder Q matrix to match the scheme for consistency
full_model$solution <- full_model$solution[match(scheme, rownames(full_model$solution)), 
                                           match(scheme, colnames(full_model$solution))]

#plot the model to take a look
plotMKmodel(full_model$solution*100, rate.cat = 1, arrow.scale = 3)












?makeSimmap







Q_dollo <- Q

#    lungless_terr
#    lunged_terr
#    lungless_lent 
#    lungless_lot
#    lungless_spec
#    lunged_lent
#    lunged_lot
#    lunged_spec
    
    
r <- Q_dollo[lungless_lot,  lunged_lot]
Q_dollo[lungless_lent, lunged_lent] <- r
Q_dollo[lungless_spec, lunged_spec] <- r
Q_dollo[lungless_terr, lunged_terr] <- r

dollo_model <- corHMM(full_tree, dataset_corHMM, rate.cat = 1, 
                     rate.mat = Q_dollo, 
                     collapse = TRUE, lower.bound = 0, node.states = "none",
                     nstarts =5, n.cores = 6)
#rename stuff for easier visualization
dimnames(dollo_model$solution) <- list(auto_matrix$legend, auto_matrix$legend)

dollo_model$solution <- dollo_model$solution[match(scheme, rownames(dollo_model$solution)), 
                                           match(scheme, colnames(dollo_model$solution))]

#plot the model to take a look
plotMKmodel(dollo_model$solution*100, rate.cat = 1, arrow.scale = 3)









try <- 





