data_full <- read.csv("lung_loss_git/processed_data/lung_data/full_data.csv")
full_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/maxLH_tree_full.nex")

tree_data <- data_full[data_full$Taxa %in% full_tree$tip.label,]
row.names(tree_data) <- tree_data$Taxa

x <- tree_data$lung
names(x) <- tree_data$Taxa

library(phytools)

t <- make.simmap(full_tree, x, model = "ARD", nsim = 10, pi = c(0,1))
cols <- setNames(c("red", "black"), c(0,1))

c <- countSimmap(t)
range(c$Tr[,4])

par(mar = c(4.2, 2.2, 2.2, 2.2))
plot(density(c$Tr[,4], bw = .4))


t2 <- densityMap(t, plot=FALSE)
plot(t2, type = "fan", ftype="off", lwd=2.5, offset = .5,
     legend = FALSE, outline = FALSE)

n <- length(t2$cols)
pal1 <- colorRampPalette(c("red", "black"))(n)
t2$cols[1:n]<-pal1

plot(t2, type = "fan", ftype="off", lwd=2.5, offset = .5,
     legend = FALSE, outline = FALSE)

########### to fix ancestor of leiopelmatoidea to lunged

tip <- getMRCA(full_tree, c("Ascaphus_truei", "Leiopelma_hochstetteri"))
tree2 <- bind.tip(full_tree, "Leio_anc" ,edge.length=0,where=tip)
x2 <- c(x,1)
names(x2) <- c(names(x), "Leio_anc")

t <- make.simmap(tree2, x2, model = "ARD", nsim = 10)
c <- countSimmap(t)
range(c$Tr[,4])

par(mar = c(4.2, 2.2, 2.2, 2.2))
hist(c$Tr[,4])
plot(density(c$Tr[,4], bw = .4))

trees <- lapply(t,drop.tip.simmap,tip="Leio_anc")
class(trees)<-"multiPhylo"

cols <- setNames(c("red", "black"), c(0,1))
plot(trees[[3]], type = "fan", ftype = "off", colors = cols)
plot(trees[[3]], type = "fan", colors = cols, fsize = .2)


?plot.simmap

#########################################################

countSimmap(trees[[3]])

tree_data$state











