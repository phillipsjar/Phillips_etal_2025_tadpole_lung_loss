
rm(list = ls())
library(ape)
library(phytools)

BT_output <- read.csv("bayestraits/output/processed_logs/eight_state_1.txt", sep = "\t")

{summary <- BT_output[, 1:max(grep("q", colnames(BT_output)))]
  
  rate_sum <- matrix(NA, ncol = 4, nrow = length(grep("q", colnames(summary))))
  rownames(rate_sum) <- colnames(summary)[grep("q", colnames(summary))]
  colnames(rate_sum) <- c("mean", "median", "SD", "count")
  
  
  for(i in 1:dim(rate_sum)[1]){
    rate_sum[i,1] <- mean(summary[,match(rownames(rate_sum)[i], colnames(summary))])
    rate_sum[i,2] <- median(summary[,match(rownames(rate_sum)[i], colnames(summary))])
    rate_sum[i,3] <- sd(summary[,match(rownames(rate_sum)[i], colnames(summary))])
    rate_sum[i,4] <- i}}

rm(BT_output)

ML_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/MaxLh_tree_full.nex") 
state_labels <- eight_state_state <- c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P", "X_terr", "Lu_terr")
trait_col_names <- eight_state_col <- c("ecology", "lung", "Spec_lotic", "terrestrial")
mode = "eight_state"
model = "eight_state"

Nruns <- 10
data <- read.csv(file = "lung_loss_git/processed_data/lung_data/full_data.csv")
data <- data[data$Taxa %in% ML_tree$tip.label,]

source("lung_loss_git/scripts/functions/custom_Q.R")
source("lung_loss_git/scripts/functions/make_simmap_data_function.R")

cols1 <- setNames(c("red", "purple", "white", "white", "red", "white",
                         "brown", "white"), 
                         c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                           "X_terr", "Lu_terr"))
cols2 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#EBFFE1",
                          "brown", "#CC6600"), 
                          c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                            "X_terr", "Lu_terr"))

Q <- custom_Q(rate_sum[,2], model = "eight_state")
sim_dat <- make_simmap_data(data, ML_tree, mode = "eight_state", 
                            eight_state_state, eight_state_col)


tip <- getMRCA(ML_tree, c("Ascaphus_truei", "Leiopelma_hochstetteri"))
ML_tree <- bind.tip(ML_tree, "Leio_anc" ,edge.length=0,where=tip)
sim_dat2 <- rbind(sim_dat, c(0,0,.5,0,0,.5,0,0))
rownames(sim_dat2) <- c(rownames(sim_dat), "Leio_anc")
sim_dat <- sim_dat2; rm(sim_dat2)

A <- make.simmap(ML_tree, sim_dat, nsim = 100,
                 Q = Q, state_labels = eight_state_state, pi=rep(1/8,8))

A <- lapply(A,drop.tip.simmap,tip="Leio_anc")
class(A)<-"multiPhylo"


plot(A[[1]], colors=cols2, lwd=1.3, ftype="off", fsize=0.5, offset=0.5, type = "fan")

# this method closely follows a suggestion by Liam Revell in a phytools blog (http://blog.phytools.org/2020/06/mapping-multi-state-discrete-character.html)

{png(file = "figures/figure_2_phylo.png", bg = "transparent", units = "in", res = 1000, width = 4.2, height = 4.2)
plotTree(A[[1]],ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=2.3, type = "fan")

par(fg="transparent",lend=1)
plotTree(A[[1]],ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
         lwd=1.3,color="white",add=TRUE, type = "fan")
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency

for(i in 1:100) {
  plot(A[[i]], colors=sapply(cols1, make.transparent,alpha=0.03),
      add=TRUE, lwd=1.3, ftype="off", fsize=0.5, offset=0.5, type = "fan")
  print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
}

dev.off()}

master_count <- matrix(0, nrow = 8, ncol = 8)
colnames(master_count) <- colnames(countSimmap(A[[1]])$Tr)
rownames(master_count) <- rownames(countSimmap(A[[1]])$Tr)

for(i in 1:100){
  master_count[1,2] <- countSimmap(A[[i]])$Tr[1,2] + master_count[1,2]
  master_count[1,3] <- countSimmap(A[[i]])$Tr[1,3] + master_count[1,3]
  master_count[1,5] <- countSimmap(A[[i]])$Tr[1,5] + master_count[1,5]
  master_count[1,7] <- countSimmap(A[[i]])$Tr[1,7] + master_count[1,7]
  master_count[2,1] <- countSimmap(A[[i]])$Tr[2,1] + master_count[2,1]
  master_count[2,4] <- countSimmap(A[[i]])$Tr[2,4] + master_count[2,4]
  master_count[2,7] <- countSimmap(A[[i]])$Tr[2,7] + master_count[2,7]
  master_count[3,1] <- countSimmap(A[[i]])$Tr[3,1] + master_count[3,1]
  master_count[3,4] <- countSimmap(A[[i]])$Tr[3,4] + master_count[3,4]
  master_count[3,6] <- countSimmap(A[[i]])$Tr[3,6] + master_count[3,6]
  master_count[3,8] <- countSimmap(A[[i]])$Tr[3,8] + master_count[3,8]
  master_count[4,2] <- countSimmap(A[[i]])$Tr[4,2] + master_count[4,2]
  master_count[4,3] <- countSimmap(A[[i]])$Tr[4,3] + master_count[4,3]
  master_count[4,8] <- countSimmap(A[[i]])$Tr[4,8] + master_count[4,8]
  master_count[5,1] <- countSimmap(A[[i]])$Tr[5,1] + master_count[5,1]
  master_count[5,6] <- countSimmap(A[[i]])$Tr[5,6] + master_count[5,6]
  master_count[5,7] <- countSimmap(A[[i]])$Tr[5,7] + master_count[5,7]
  master_count[6,3] <- countSimmap(A[[i]])$Tr[6,3] + master_count[6,3]
  master_count[6,5] <- countSimmap(A[[i]])$Tr[6,5] + master_count[6,5]
  master_count[6,8] <- countSimmap(A[[i]])$Tr[6,8] + master_count[6,8]
  master_count[7,1] <- countSimmap(A[[i]])$Tr[7,1] + master_count[7,1]
  master_count[7,2] <- countSimmap(A[[i]])$Tr[7,2] + master_count[7,2]
  master_count[7,5] <- countSimmap(A[[i]])$Tr[7,5] + master_count[7,5]
  master_count[7,8] <- countSimmap(A[[i]])$Tr[7,8] + master_count[7,8]
  master_count[8,3] <- countSimmap(A[[i]])$Tr[8,3] + master_count[8,3]
  master_count[8,4] <- countSimmap(A[[i]])$Tr[8,4] + master_count[8,4]
  master_count[8,6] <- countSimmap(A[[i]])$Tr[8,6] + master_count[8,6]
  master_count[8,7] <- countSimmap(A[[i]])$Tr[8,7] + master_count[8,7]
}

master_count/100





families <- unique(data$Family)
indicator <- rep(NA, length(families))
names(indicator) <- families

for(i in 1:length(families)){
  if(length(which(data$Family == families[i] & data$lung == 0)) > 0){
    indicator[i] <- 1
  }else{indicator[i] <- 0}
}

fams <- families[which(indicator == 0)]
rm(families)

lungless_data <- data[!(data$Family %in% fams),]


genera <- unique(lungless_data$Genus)
indicator <- rep(NA, length(genera))
names(indicator) <- genera

for(i in 1:length(genera)){
  if(length(which(lungless_data$Genus == genera[i] & lungless_data$lung == 0)) > 0){
    indicator[i] <- 1
  }else{indicator[i] <- 0}
}


lungless_genera <- genera[which(indicator == 1)]
lunged_genera <- genera[which(indicator == 0)]

tips <- rep(NA, (length(fams) + length(lungless_genera) + length(lunged_genera)))
names(tips) <- c(fams, lungless_genera, lunged_genera)


for(i in 1:length(tips)){
  if(names(tips[i]) %in% fams){
    taxa <- data$Taxa[which(data$Family == names(tips[i]))]
    ecology <- data$eight_state[which(data$Family == names(tips[i]))]
    eco_mode <- as.numeric(names(sort(table(ecology), decreasing = T))[1])
    tips[i] <- taxa[which(ecology == eco_mode)][1]}
  if(names(tips[i]) %in% lunged_genera){
    taxa <- data$Taxa[which(data$Genus == names(tips[i]))]
    ecology <- data$eight_state[which(data$Genus == names(tips[i]))]
    eco_mode <- as.numeric(names(sort(table(ecology), decreasing = T))[1])
    tips[i] <- taxa[which(ecology == eco_mode)][1]}
  if(names(tips[i]) %in% lungless_genera){
    taxa <- data$Taxa[which(data$Genus == names(tips[i]) & data$lung == 0)]
    ecology <- data$eight_state[which(data$Genus == names(tips[i]) & data$lung == 0)]
    eco_mode <- as.numeric(names(sort(table(ecology), decreasing = T))[1])
    tips[i] <- taxa[which(ecology == eco_mode)][1]}
}

tips <- tips[!(tips %in% data$Taxa[which(data$Family == "Microhylidae")])]

tips[which(names(tips) == "Crinia")] <- "Crinia_riparia"

tips <- c(tips, "Litoria_subglandulosa", "Rana_sauteri", "Thoropa_miliaris", "Kalophrynus_sinensis",
          "Otophryne_robusta", "Phrynomantis_bifasciatus")



A_trim <- lapply(A,keep.tip.simmap,tip=tips)
class(A_trim)<-"multiPhylo"

{pdf(file = "figures/figure_4_phylo.pdf", width = 7.5, height = 3, bg = "transparent")
par(xpd = TRUE)
plotSimmap(A_trim[[3]], cols2, pts=F, ftype = "off", fsize = .3, direction = "upwards",
           outline = TRUE, add=FALSE, lwd=1.75, offset=1, mar = c(.1,.1,.5,.1))
dev.off()}




############# Suppl. Fig. 3
cols3 <- setNames(c("gold", "purple", "#99CCFF", "blue", "red", "#E2FFE2",
                          "brown", "#CC6600"), 
                          c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P",
                            "X_terr", "Lu_terr"))

{jpeg(file = "figures/supp_fig_3.jpg", bg = "transparent", units = "in", res = 1000, width = 10, height = 10)
  plotTree(A[[1]],ftype="i",fsize=0.35,offset=3, #black tree with wide edges to be outlines later
           lwd=5, type = "fan")
  
  par(fg="transparent",lend=1)
  plotTree(A[[1]],ftype="i",fsize=0.35,offset=3,     #white tree to create the outline (smaller lwd)
           lwd=3,color="white",add=TRUE, type = "fan")

  ?plotSimmap
  
  ## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
  ## with 99% transparency
  
  for(i in 1:100) { 
    plot(A[[i]], colors=sapply(cols3, make.transparent,alpha=0.03),
         add=TRUE, lwd=3, ftype="i", fsize=0.35, offset=3, type = "fan")
    print(100-i) #just added so you can see your progress, counts down to zero (can take a while)
  }
  
  dev.off()}



