BT_4062_tree <- read.tree(file = "Trees/edited_trees/maxLH_BT_vis_tree.tre")
data_no_endo <- read.csv(file = "data/no_endo_lung_data.csv")
load("Bayestraits/results/tree_set_dep_master_export.Rdata")

rates <- master_export_tree_set_dep
rm(master_export_tree_set_dep)

set.seed(123)

random_runs <- sample(1:(dim(rates)[1]), 100, replace = FALSE)  # get a sample of 100 random runs from the posterior

rate_values_max_lh <- vector(mode = "list", length = 100)             # empty list for the rate values for each random run

for(i in 0:99){
  start <- which(colnames(rates) == "q12")
  run <- random_runs[i+1]
  rate_values_max_lh[i+1] <- list(rates[(run),start:(start+7)])
}

# phytools does not have a native dependent model, but this can easily be dealt with by treating the 
# dependent model as a constrained 4 state model, with each combination of the two binary states treated
# as one of four states (lungless stream; lunged, stream; lungless, pond; lunged, pond)

custom_Q = function(rates){
  rates <- unlist(rates)
  rate_matrix <- matrix(0, nrow = 4, ncol = 4)         # matrix is populated with zeroes to begin
  colnames(rate_matrix) <- c("X_S", "Lu_S", "X_P", "Lu_P")
  rownames(rate_matrix) <- c("X_S", "Lu_S", "X_P", "Lu_P")
  rate_matrix[1,2] <- rates[1]*.0001                   # rates are multiplied by .0001 to counteract scaling in BT
  rate_matrix[1,3] <- rates[2]*.0001
  rate_matrix[2,1] <- rates[3]*.0001
  rate_matrix[2,4] <- rates[4]*.0001
  rate_matrix[3,1] <- rates[5]*.0001
  rate_matrix[3,4] <- rates[6]*.0001
  rate_matrix[4,2] <- rates[7]*.0001
  rate_matrix[4,3] <- rates[8]*.0001
  for(i in 1:4){
    rate_matrix[i,i] <- -(sum(rate_matrix[i,]))}       # the rates along the diagonal should sum each row to zero
  return(rate_matrix)
}
rm(i, run, start, N_nodes, tips, iteration, averages)
custom_Q(rate_values_max_lh[2])



library(phytools)

##### set up simmap for pseudo-dependent analysis

# data_4062 loaded above.
# simmap requires a specific way data are input, which we create below with a quick function

source("scripts/extra_scripts/make_simmap_data_function.R")

lung_eco_matrix_max_LH <- make_simmap_data(data_no_endo, BT_4062_tree)

master_simmap <- make.simmap(BT_4062_tree,lung_eco_matrix_max_LH, nsim=100, Q = custom_Q(rate_values_max_lh[2]), pi=rep(.25,4))

for(i in 1:100){
  master_simmap[[i]] <- make.simmap(BT_4062_tree,lung_eco_matrix_max_LH, nsim=1, Q = custom_Q(rate_values_max_lh[i]), pi=rep(.25,4))
  print(100-i) # count down for sanity
}

cols2 <- setNames(c("#990099", "#3399FF", "#FF3333","lightgreen", "grey"), c("X_S", "Lu_S", "X_P", "Lu_P"))
cols2 <- setNames(c("magenta", "blue", "#FF3333","#99FF99", "grey"), c("X_S", "Lu_S", "X_P", "Lu_P"))




png(file = "figures/test_lend2.png", bg = "black", width = 2400 , height =  2400)
par(lend = 0)
plotSimmap(master_simmap[[2]], colors=sapply(cols2, make.transparent,alpha=0.5), pts=F, 
           ftype="off", fsize = .2, type = "fan", lwd=12, offset=0.5)
dev.off()


#tester of one run to see if colors look good
png(file = "figures/conclusion_slide.png", bg = "transparent", width = 300 , height =  1000)
plotSimmap(master_simmap[[2]], cols2, pts=F, ftype="off", fsize = .2, type = "phylogram", lwd=2, offset=0.5)
dev.off()
# this method closely follows a suggestion by Liam Revell in a phytools blog (http://blog.phytools.org/2020/06/mapping-multi-state-discrete-character.html)

plotTree(BT_4062_tree,ftype="off",fsize=0.5,offset=0.5, #black tree with wide edges to be outlines later
         lwd=5, type = "fan")

par(fg="transparent",lend=1)
plotTree(BT_4062_tree,ftype="off",fsize=0.5,offset=0.5,     #white tree to create the outline (smaller lwd)
         lwd=2.75,color="white",add=TRUE, type = "fan")
## now plot our 100 stochastic map trees pulled from master (same lwd as white tree)
## with 99% transparency






for(i in 1:10) {plot(master_simmap[[i*10]],
                      colors=sapply(cols2, make.transparent,alpha=0.1),
                      add=TRUE, lwd=2.75, ftype="off", fsize=0.5, offset=0.5, type = "fan")
  print(10-i) #just added so you can see your progress, counts down to zero (can take a while)
}
par(fg="black")

legend(x="bottomleft",c("lungless in streams", "lunged in streams", "lungless in ponds", "lunged in ponds"),pch=22,
       pt.bg=cols2,pt.cex=1.5,bty="n",cex=0.7)



for(i in 1:100){
iter <- c(1:50,seq(55,100,5))
if(i %in% iter){
  png(file = paste("figures/title_gif",i,".png", sep = "_"), bg = "transparent", width = 2400 , height =  2400)
  plotTree(BT_4062_tree,ftype="off",fsize=0.5,offset=0.5,
           lwd=13.5, color="black", type = "fan")
  
  for(j in 1:i){ plot(master_simmap[[j]], type = "fan",
                      colors=sapply(cols2,make.transparent,alpha=0.0225),
                      add=TRUE,lwd=13.5,ftype="off",fsize=0.5,offset=0.5)}
dev.off()}}



{i = 100
png(file = paste("figures/title_gif",i+1,".png", sep = "_"), bg = "transparent", width = 2400 , height =  2400)
plotTree(BT_4062_tree,ftype="off",fsize=0.5,offset=0.5,
         lwd=12, color="black", type = "fan")

for(j in 1:i){ plot(master_simmap[[j]], type = "fan",
                    colors=sapply(cols2,make.transparent,alpha=0.025),
                    add=TRUE,lwd=12,ftype="off",fsize=0.5,offset=0.5)}
dev.off()}





800*6

dev.off()

?make.transparent

#"#FF00FF0D" "#00FFFF0D" "#FF33330D" "#99FF990D" "#BEBEBE0D" 


sapply(cols2, make.transparent,alpha=0.05)

"#FF00FF19"







