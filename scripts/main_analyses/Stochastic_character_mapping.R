setwd("/Users/jack/desktop/Research/lunglessness/lung_loss_Phillips_etal_2024")

library(phytools)

source("lung_loss_git/scripts/set_up/Initial_set_up.R")
source("lung_loss_git/scripts/functions/Bayestraits_simmap.R")

load(file = "bayestraits_exports/master_dep_export.Rdata")
load(file = "bayestraits_exports/master_indep_export.Rdata")
load(file = "bayestraits_exports/master_six_state_dep_export.Rdata")
load(file = "bayestraits_exports/master_six_state_indep_export.Rdata")


Nruns <- 10
data <- read.csv(file = "lung_loss_git/processed_data/lung_data/no_endo_lung_data.csv")
trees <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/tree_set.nex")
data <- data[data$Taxa %in% trees[[1]]$tip.label,]

four_state_state = c("X_S", "Lu_S", "X_P", "Lu_P")
four_state_col = c("ecology", "lung")

six_state_state = c("X_S_unsp", "X_S_sp", "Lu_S_unsp", "Lu_S_sp", "X_P", "Lu_P")
six_state_col = c("ecology", "lung", "Spec_lotic")

set.seed(123)


{t_indep <- BayesTraits_simmap(master_indep, Nruns, data, model = "independent", 
                         trees, state_labels = four_state_state, 
                         trait_col_names = four_state_col)
t_dep <- BayesTraits_simmap(master_dep, Nruns, data, model = "dependent", 
                         trees, state_labels = four_state_state, 
                         trait_col_names = four_state_col)
t_six_indep <- BayesTraits_simmap(master_six_state_indep_export, Nruns, data, 
                         model = "six_state", trees, state_labels = six_state_state, 
                         trait_col_names = six_state_col)
t_six_dep <- BayesTraits_simmap(master_six_state_dep_export, Nruns, data, 
                         model = "six_state", trees, state_labels = six_state_state, 
                         trait_col_names = six_state_col)}


t_six_dep <- BayesTraits_simmap(master_six_state_dep_export, Nruns, data, 
                                model = "six_state", trees, state_labels = six_state_state, 
                                trait_col_names = six_state_col)


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


{plotSimmap(t_six_dep[[1]], cols2, pts=F, ftype="off", fsize = .2, 
            outline = TRUE, type = "fan", add=FALSE, lwd=2, offset=0.5)
  title(main="six_dependent model",font.main=3,
        line=-1)}



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
    ecology <- data$ecology[which(data$Family == names(tips[i]))]
    if(mean(ecology, na.rm = TRUE) > .5){
      tips[i] <- taxa[which(ecology == 1)][1]}
    if(mean(ecology, na.rm = TRUE) < .5){
      tips[i] <- taxa[which(ecology == 0)][1]}
    if(mean(ecology, na.rm = TRUE) == .5){
      tips[i] <- taxa[sample(length(taxa), 1)]}
  }
  if(names(tips[i]) %in% lunged_genera){
    taxa <- data$Taxa[which(data$Genus == names(tips[i]))]
    ecology <- data$ecology[which(data$Genus == names(tips[i]))]
    if(mean(ecology, na.rm = TRUE) > .5){
      tips[i] <- taxa[which(ecology == 1)][1]}
    if(mean(ecology, na.rm = TRUE) < .5){
      tips[i] <- taxa[which(ecology == 0)][1]}
    if(mean(ecology, na.rm = TRUE) == .5){
      tips[i] <- taxa[sample(length(taxa), 1)]}
  }
  if(names(tips[i]) %in% lungless_genera){
    taxa <- data$Taxa[which(data$Genus == names(tips[i]) & data$lung == 0)]
    ecology <- data$ecology[which(data$Genus == names(tips[i]) & data$lung == 0)]
    if(mean(ecology, na.rm = TRUE) > .5){
      tips[i] <- taxa[which(ecology == 1)][1]}
    if(mean(ecology, na.rm = TRUE) < .5){
      tips[i] <- taxa[which(ecology == 0)][1]}
    if(mean(ecology, na.rm = TRUE) == .5){
      tips[i] <- taxa[sample(length(taxa), 1)]}
  }
}

i = which(names(tips) == "Clinotarsus")

tips[which(names(tips) == "Crinia")] <- "Crinia_riparia"

tips <- c(tips, "Litoria_subglandulosa", "Rana_sauteri", "Thoropa_miliaris")


plotSimmap(t_six_dep[[2]], cols2, pts=F, fsize = .3, direction = "upwards",
           outline = TRUE, add=FALSE, lwd=1.75, offset=2)

t_six_dep_trim <- lapply(t_six_dep,keep.tip.simmap,tip=tips)
class(t_six_dep_trim)<-"multiPhylo"

plotSimmap(t_six_dep_trim[[2]], cols2, pts=F, fsize = .3, direction = "upwards",
           outline = TRUE, add=FALSE, lwd=1.75, offset=2)

plot_sim <- t_six_dep_trim[[2]]
plot_sim$tip.label <- names(tips)[match(plot_sim$tip.label, tips)]

labels <- rep(NA, length(plot_sim$tip.label))

for(i in 1:length(labels)){
  if(plot_sim$tip.label[i] %in% fams){
    labels[i] <- plot_sim$tip.label[i]
  }else{ labels[i] <- ""}
}

labels <- c("Rhacophoridae", "Dicroglossidae", "Phrynobatrachidae", "Alytidae",
            "Bombinatoridae", "Rhinophrynidae", "Pipidae", "Scaphiopodidae", 
            "Pelodytidae", "Pelobatidae", "Batrachylidae", "Alsodidae", "Hylodidae", 
            "Ceratophryidae", "Hemiphractidae", "Dendrobatidae", "Odontophrynidae", 
            "Allophrynidae", "Leptodactylidae", "Calyptocephalellidae", "Limnodynastinae", 
            "Hemisotidae", "Hyperoliidae", "Pyxicephalidae", "Conrauidae", "Petropedetidae", 
            "Ptychadenidae", "Rhacophoridae")


ntips <- length(plot_sim$tip.label)

pdf(file = "figures/figure_4_phylo.pdf", width = 7.5, height = 3, bg = "transparent")
par(xpd = TRUE)
plotSimmap(plot_sim, cols2, pts=F, ftype = "off", fsize = .3, direction = "upwards",
            outline = TRUE, add=FALSE, lwd=1.75, offset=1, mar = c(.1,.1,.5,.1))
#tiplabels(labels, tip = match(labels, plot_sim$tip.label), frame = "none", adj = c(0,.5), offset = 4, cex = .4)
dev.off()




?plotSimmap

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
sd_simmap_counts(t_dep)

avg_simmap_counts(t_six_indep)
sd_simmap_counts(t_six_indep)

avg_simmap_counts(t_six_dep)
sd_simmap_counts(t_six_dep)












#     [,1] [,2]
#[1,]    4    5
#[2,]    6    7

















