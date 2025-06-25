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



install.packages("phylolm")
#simple phyloglm:
require(phylolm)
packageVersion("phylolm")

mod <- phyloglm(lung ~ ecology, data_aqua, full_tree, method = c("logistic_MPLE"),
                 start.beta=NULL, start.alpha=NULL,
                 boot = 1000, full.matrix = TRUE)

summary(mod)
# ecology recovered as statistically significant (p = 0.04128) and estimated effect size = 0.269086.


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

mod$aic
mod2$aic
mod3$aic





anova(mod)



car::Anova(mod, type = "3")

table(data_full$ecology)

table(data_full$Spec_lotic)
table(data_aqua$ecology, data_aqua$lung)
table(data_aqua$ecology)
table(data_aqua$lung)

data_aqua


table(data_aqua$ecology[which(data_aqua$Family == "Bufonidae")])


table(data_full$lung, data_full$ecology, data_full$Spec_lotic, data_full$terrestrial)


data_full$Spec_lotic















# visualization of the binary data being input into the analyses

# custom function for a phylomorphospace with flexibility for certain aesthetic properties
source("scripts/final_scripts/extra_scripts/custom_phylomorphospace_function.R")

{ML_tree <- read.nexus(file = "lung_loss_git/bayestraits_trees_data/trees/MaxLh_tree_aqu.nex") 

tips <- as.data.frame(cbind(data_4062$ecology, data_4062$lung)) #add ecology and lung to an nx2 matrix
rownames(tips) <- data_4062$Taxa}

#Bayestraits Ancestral State Estimates (see Ancestral_reconstructions.R) 
{load(file = "data/4062_ASE_R_formatted.R")
nodes <- as.data.frame(nodes)
nodes <- cbind(nodes$lentic, nodes$lunged)}

########################################################################################################
# figure 2 plotting

{pdf(file = "figures/phylomorphospace.pdf", bg = "transparent", width = 4.25 , height =  3)
par(xpd = TRUE)
par(mar=c(2.1, 2.1, 1.1, 1.1))
plot(1, 1, ylim = c(-0.1,1.1), xlim = c(-0.1,1.1), xlab = "", pch = 21, ylab = "", lwd = 1.5, cex = 1.6, type = "n",
     xaxt = "n", yaxt = "n" )
axis(2,cex.axis=.75)
axis(1,cex.axis=.75)

cc <- summary(mod_4062)$coefficients
x_sim <- seq(0, 1, length = 100)
polygon(c(x_sim, sort(x_sim, decreasing = TRUE)), c((exp(cc[1,4] + cc[2,4]*x_sim))/(1 + (exp(cc[1,4] + cc[2,4]*x_sim))), 
       (exp(cc[1,5] + cc[2,5]*sort(x_sim, decreasing = TRUE)))/(1 + (exp(cc[1,5] + cc[2,5]*sort(x_sim, decreasing = TRUE))))),
        col = "#EEEEEE")
curve(plogis(cc[1,1]+cc[2,1]*x),col="black", lwd = 2, add=TRUE, xlim = c(0,1))

custom_phylomorphospace(BT_4062_tree,tips,nodes, jitter = .02) #using default colors

load("Bayestraits/data/master_export.Rdata")

avgs <- sapply(master_export[,grep("q12", colnames(master_export)):(grep("q12", colnames(master_export))+7)], "mean");
require(scales)
avgs_scaled <- rescale(avgs,c(0,20),c(0,10))

#plot(1, 1, ylim = c(-0.1,1.1), xlim = c(-0.1,1.1), xlab = "", pch = 21, ylab = "", lwd = 1.5, cex = 1.6, type = "n")
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                
arrows(x0 = -.05, x1 = -.05, col = colorBlindGrey8[1],
       y0 = .25, y1 = .75, lwd = 2.5,
       length = 0.1,
       angle = 20) 
arrows(x0 = -.11, x1 = -.11, col = colorBlindGrey8[3],
       y0 = .75, y1 = .25, lwd = 2.5,
       length = 0.1,
       angle = 20) 
arrows(x0 = .25, x1 = .75,  col = colorBlindGrey8[2],
       y0 = -.05, y1 = -.05, lwd = 2.5,
       length = 0.1,
       angle = 20) 
arrows(x0 = .75, x1 = .25, col = colorBlindGrey8[4],
       y0 = -.11, y1 = -.11, lwd = 2.5,
       length = 0.1,
       angle = 20) 
arrows(x0 = .75, x1 = .25, col = colorBlindGrey8[5],
       y0 = 1.11, y1 = 1.11, lwd = 2.5,
       length = 0.1,
       angle = 20) 
arrows(x0 = .25, x1 = .75, col = colorBlindGrey8[7],
       y0 = 1.05, y1 = 1.05, lwd = 2.5,
       length = 0.1,
       angle = 20) 
arrows(x0 = 1.05, x1 = 1.05, col = colorBlindGrey8[6],
       y0 = .25, y1 = .75, lwd = 2.5,
       length = 0.1,
       angle = 20) 
arrows(x0 = 1.11, x1 = 1.11, col = colorBlindGrey8[8],
       y0 = .75, y1 = .25, lwd = 2.5,
       length = 0.1,
       angle = 20) 

  text(x = -.05, y = .775, labels = round(avgs[1], 2), cex = .55, srt = 90, adj = c(0,.5))
  arrows(x0 = -.05, x1 = -.05, 
        y0 = .25, y1 = .75, lwd = .5,
        length = 0.1,
        angle = 20) 
  text(x = -.11, y = .225, labels = round(avgs[3], 2), cex = .55, srt = 90, adj = c(1,.5))
  arrows(x0 = -.11, x1 = -.11,
         y0 = .75, y1 = .25, lwd = .5 ,
         length = 0.1,
         angle = 20) 
  
  text(x = .775, y = -.05, labels = round(avgs[2], 2), cex = .55, adj = c(0,.5))
  arrows(x0 = .25, x1 = .75, 
         y0 = -.05, y1 = -.05, lwd = .5,
         length = 0.1,
         angle = 20) 
  text(x = .225, y = -.1, labels = round(avgs[5], 2), cex = .55, adj = c(1,.5))
  arrows(x0 = .75, x1 = .25,
         y0 = -.11, y1 = -.11, lwd = .5,
         length = 0.1,
         angle = 20) 

  
  text(x = .225, y = 1.11, labels = round(avgs[7], 2), cex = .55, adj = c(1,.5))
  arrows(x0 = .75, x1 = .25,
         y0 = 1.11, y1 = 1.11, lwd = .5,
         length = 0.1,
         angle = 20) 
  text(x = .775, y = 1.05, labels = round(avgs[4], 2), cex = .55, adj = c(0,.5))
  arrows(x0 = .25, x1 = .75,
         y0 = 1.05, y1 = 1.05, lwd = .5,
         length = 0.1,
         angle = 20) 
  text(x = 1.05, y = .775, labels = round(avgs[6], 3), srt = 90, cex = .55, adj = c(0,.5))
  arrows(x0 = 1.05, x1 = 1.05,
         y0 = .25, y1 = .75, lwd = .5,
         length = 0.1,
         angle = 20) 
  text(x = 1.11, y = .225, labels = round(avgs[8], 5), srt  = 90, cex = .55, adj = c(1,.5))
  arrows(x0 = 1.11, x1 = 1.11,
         y0 = .75, y1 = .25, lwd = .5,
         length = 0.1,
         angle = 20)
    dev.off()}


                                


















 par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plot(BT_4062_tree, show.tip.label = FALSE, type = "f")
reconCol <- c("magenta", "blue", "red", "green")
tips <- length(BT_4062_tree$tip.label)          # number of tips in tree
nodelabels(node = ((tips+1):(2*tips-1)), pie = as.matrix(nodes[,1:4]), cex = .45, piecol = c("magenta", "blue", "red", "green"))

par(mar=c(4.1,4.1,2,2))
phylomorphospace(BT_4062_tree, X = PMS_tips, label = "off")

rm(nodes, data_4062)



A = PMS_nodes

?phylomorphospace



length(BT_4062_tree$tip.label)
dim(PMS_tips)

plot(jitter(data_4062$lung, .5) ~ jitter(data_4062$ecology, .5))
plot(jitter(data_portik$lung, .5) ~ jitter(data_portik$ecology, .5))




tree <- tree_4062_partial
simple_lung_full <- data_4062_partial
rownames(simple_lung_full)<-simple_lung_full$Taxa;


simple_lung_nobufo <- simple_lung_full[which(simple_lung_full$Family != "Bufonidae"),]
rownames(simple_lung_nobufo)<-simple_lung_nobufo$Taxa;
name.check(tree,simple_lung_nobufo) -> overlap
drop.tip(tree,overlap$tree_not_data) -> tree_no_bufo


simple_lung_full <- na.omit(simple_lung_full)

?na.omit

simple_lung_full <- read.csv("data/412_matrix.csv")
simple_lung_full <- simple_lung_full[simple_lung_full$Taxa %in% tree_396$tip.label,]
#simple_lung_full$lung <-
simple_lung_full <- na.omit(simple_lung_full)

simple_lung_full <- subset(simple_lung_full, ecology == "1" | ecology == "0")
rownames(simple_lung_full)<-simple_lung_full$Taxa;
name.check(tree_396,simple_lung_full) -> overlap
drop.tip(tree_396,overlap$tree_not_data) -> tree_375

complex_nobufo[which(complex_nobufo$lung == 4),]

complex_nobufo$lung[which(complex_nobufo$ecology == 1)]

simple_lung_nobufo <- simple_lung_full[which(simple_lung_full$family != "bufonidae"),]
rownames(simple_lung_nobufo)<-simple_lung_nobufo$Taxa;
name.check(tree_396,simple_lung_nobufo) -> overlap
drop.tip(tree_396,overlap$tree_not_data) -> tree_302


complex_full <- read.csv("phylolm/138_matrix_detailedlung.csv")
complex_nobufo <- complex_full[which(complex_full$family != "bufonidae"),]
rownames(complex_full)<-complex_full$Taxa;
rownames(complex_nobufo)<-complex_nobufo$Taxa;

name.check(tree_396,complex_full) -> overlap
drop.tip(tree_396,overlap$tree_not_data) -> tree_138

name.check(tree_396,complex_nobufo) -> overlap
drop.tip(tree_396,overlap$tree_not_data) -> tree_106






#simplest phylolm:

mod1 <- phyloglm(lung ~ ecology, simple_lung_full, tree, method = c("logistic_MPLE"),
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)

summary(mod1)


#remove bufonids

mod2 <- phyloglm(lung ~ ecology, simple_lung_nobufo, tree_no_bufo, method = c("logistic_MPLE"),
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod2)


mod3 <- glm(lung ~ ecology, data = simple_lung_full, family = binomial(link = "logit"))
summary(mod3)


#complex lung data (logistic phyloglm)

mod3 <- phyloglm(ecology ~ lung, complex_full, tree_138, method = c("logistic_IG10"),
                      start.beta=NULL, start.alpha=NULL,
                      boot = 0, full.matrix = TRUE)
summary(mod3)


#complex lung data and remove bufonids

mod4 <- phyloglm(ecology ~ lung, complex_nobufo, tree_106, method = c("logistic_IG10"),
                 start.beta=NULL, start.alpha=NULL,
                 boot = 0, full.matrix = TRUE)
summary(mod4)


#add in complex ecology data

mod <- phylolm(lung ~ ecology_2, complex_full, tree_138, 
               start.beta=NULL, start.alpha=NULL,
               boot = 0, full.matrix = TRUE)

mod <- phyloglm(ecology_2 ~ lung, complex_full, tree_138, method = c("poisson_GEE"),
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)

mod <- phyloglm(ecology_2 ~ lung, complex_nobufo, tree_106, method = c("poisson_GEE"),
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)
summary(mod)




dev.off()
plot(as.factor(complex_full$ecology_2) ~ complex_full$lung)
plot(as.factor(complex_nobufo$ecology_2) ~ complex_nobufo$lung)


plotdata <- matrix(ncol = 5, nrow = 2)
rownames(plotdata) <- c("lentic", "lotic")
colnames(plotdata) <- c("buds", "reduced", "average", "large", "greatly enlarged")

for(i in 1:5){
  plotdata[1,i] <- length(which(complex_full$ecology == 1 & complex_full$lung == i-1))    #lunged, lentic
  plotdata[2,i] <- length(which(complex_full$ecology == 0 & complex_full$lung == i-1))    #lunged, lentic
}
data_percentage <- apply(plotdata, 2, function(x){x*100/sum(x,na.rm=T)})

pdf(file="figures/barplot_binary.pdf",width=7,height=4)
barplot(data_percentage, 
        col= c("forest green", "deep sky blue") , 
        border= "black" , 
        space=0.04, 
        font.axis=2, 
        ylab = "percentage",
        xlab = "lung state")
dev.off()

plotdata2 <- matrix(ncol = 5, nrow = 2)
rownames(plotdata2) <- c("lentic", "lotic")
colnames(plotdata2) <- c("buds", "reduced", "average", "large", "greatly enlarged")
for(i in 1:5){
  plotdata2[1,i] <- length(which(complex_nobufo$ecology == 1 & complex_nobufo$lung == i-1))    #lunged, lentic
  plotdata2[2,i] <- length(which(complex_nobufo$ecology == 0 & complex_nobufo$lung == i-1))    #lunged, lentic
}
data_percentage2 <- apply(plotdata2, 2, function(x){x*100/sum(x,na.rm=T)})

pdf(file="figures/barplot_binary_nobufo.pdf",width=7,height=4)
barplot(data_percentage2, 
        col= c("forest green", "deep sky blue") , 
        border= "black" , 
        space=0.04, 
        font.axis=2, 
        ylab = "number of species")
dev.off()

plotdata3 <- matrix(ncol = 5, nrow = 4)
rownames(plotdata3) <- c("specialized lotic", "lotic", "lentic", "specialized lentic")
colnames(plotdata3) <- c("buds", "reduced", "average", "large", "greatly enlarged")
for(i in 1:5){
  plotdata3[1,i] <- length(which(complex_full$ecology_2 == 0 & complex_full$lung == i-1))    #
  plotdata3[2,i] <- length(which(complex_full$ecology_2 == 1 & complex_full$lung == i-1))    #
  plotdata3[3,i] <- length(which(complex_full$ecology_2 == 2 & complex_full$lung == i-1))    #
  plotdata3[4,i] <- length(which(complex_full$ecology_2 == 3 & complex_full$lung == i-1))    #
}
data_percentage3 <- apply(plotdata3, 2, function(x){x*100/sum(x,na.rm=T)})


pdf(file="figures/barplot_full.pdf",width=7,height=4)
barplot(data_percentage3, 
        col = c("dark blue", "deep sky blue", "green", "forest green"),
        border= "black" , 
        space=0.04, 
        font.axis=2, 
        ylab = "number of species")
dev.off()

plotdata4 <- matrix(ncol = 5, nrow = 4)
rownames(plotdata4) <- c("specialized lotic", "lotic", "lentic", "specialized lentic")
colnames(plotdata4) <- c("buds", "reduced", "average", "large", "greatly enlarged")
for(i in 1:5){
  plotdata4[1,i] <- length(which(complex_nobufo$ecology_2 == 0 & complex_nobufo$lung == i-1))    #
  plotdata4[2,i] <- length(which(complex_nobufo$ecology_2 == 1 & complex_nobufo$lung == i-1))    #
  plotdata4[3,i] <- length(which(complex_nobufo$ecology_2 == 2 & complex_nobufo$lung == i-1))    #
  plotdata4[4,i] <- length(which(complex_nobufo$ecology_2 == 3 & complex_nobufo$lung == i-1))    #
}
data_percentage4 <- apply(plotdata4, 2, function(x){x*100/sum(x,na.rm=T)})

pdf(file="figures/barplot_full_nobufo.pdf",width=7,height=4)
barplot(data_percentage4, 
        col = c("dark blue", "deep sky blue", "green", "forest green"),
        border= "black" , 
        space=0.04, 
        font.axis=2, 
        ylab = "number of species")
dev.off()


pdf(file="figures/barplot_full_nobufo_labelled.pdf",width=7,height=4)
barplot(data_percentage4, 
        col = c("dark blue", "deep sky blue", "green", "forest green"),
        border= "black" , 
        space=0.04, 
        font.axis=2, 
        ylab = "number of species")

int <- 1.04
for(i in 1:5){
  if(plotdata4[1,i] != 0){
    text(.54+(i-1)*int, (data_percentage4[1,i])/2 , labels=round(plotdata4[1,i]), col="white", lwd = 2, cex = 1)
  }
  if (plotdata4[2,i] != 0){
    text(.54+(i-1)*int, data_percentage4[1,i]+data_percentage4[2,i]/2 , labels=round(plotdata4[2,i]), col="black", lwd = 2, cex = 1)
  }
  if(plotdata4[3,i] != 0){
    text(.54+(i-1)*int, data_percentage4[1,i]+data_percentage4[2,i]+data_percentage4[3,i]/2 , labels=round(plotdata4[3,i]), col="black", lwd = 2, cex = 1)
  }
  if(plotdata4[4,i] != 0){
    text(.54+(i-1)*int, data_percentage4[1,i]+data_percentage4[2,i]+data_percentage4[3,i]+data_percentage4[4,i]/2 , labels=round(plotdata4[4,i]), col="white", lwd = 2, cex = 1)
  }
  else{}
}
dev.off()







lung_bin_comp <- replace(complex_full$lung, complex_full$lung > 0, 1)
lung_bin_comp <- to.matrix(x = lung_bin_comp,c(0,1))
rownames(lung_bin_comp) <- complex_full$Taxa
colnames(lung_bin_comp) <- c("lungless", "lunged")

eco_comp <- to.matrix(x = complex_full$ecology,c(0,1))
rownames(eco_comp) <- complex_full$Taxa
colnames(eco_comp) <- c("lotic", "lentic")

lung_comp <- to.matrix(x = complex_full$lung,c(0,1,2,3,4))
lung_comp <- as.matrix(complex_full$lung)
rownames(lung_comp) <- complex_full$Taxa
colnames(lung_comp) <- c("buds", "reduced", "average", "large", "greatly enlarged")


lung_tree1 <-make.simmap(tree_138,lung_bin_comp, nsim=100, model = "ER", pi="estimated")
cols <- setNames(c("black","red"), c("lunged", "lungless"))
plotSimmap(lung_tree1[[1]],cols,pts=F,ftype="off")

lung_tree2 <-make.simmap(tree_138,eco_comp, nsim=100, model = "ER", pi="estimated")
cols2 <- setNames(c("blue","green"), c("lotic", "lentic"))
plotSimmap(lung_tree2[[1]],cols2,pts=F,ftype="off")

obj<-contMap(tree_138,lung_comp)
plot(obj,lwd=7,xlim=c(-0.2,3.6))

lung_tree3 <-make.simmap(tree_138,lung_comp, nsim=100, model = "ER", pi="estimated")
cols3 <- setNames(c("red","light grey","grey","dark grey","black"), c("buds", "reduced", "average", "large", "greatly enlarged"))


plotSimmap(lung_tree3[[1]],cols3,pts=F,ftype="off")


pdf(file="figures/138treewithnames.pdf",width=6,height=8)
par(mar=c(1,1,1,1))
plot(tree_138, cex = .4)
table(data2$lung)
table(data$lung)
dev.off()




data2 <- data[which(data1$family != "bufonidae"),]
rownames(data2)<-data2$Taxa;
name.check(tree_396,data2) -> overlap

drop.tip(tree_396,overlap$tree_not_data) -> tree_nobufo

#mod <- phylolm( ~ ecology, data, tree_315, 
#                start.beta=NULL, start.alpha=NULL,
#                boot = 0, full.matrix = TRUE)

mod <- phyloglm(ecology ~ lung, data1, tree_138, method = c("logistic_IG10"),
                start.beta=NULL, start.alpha=NULL,
                boot = 0, full.matrix = TRUE)

plot( as.factor(data1$ecology) ~ data1$lung)
plot( as.factor(data1$ecology) ~ data1$lung)

?phyloglm


summary(mod)

modnobufo <- phyloglm(ecology ~ lung, data2, tree_nobufo, method = c("logistic_IG10"),
                      start.beta=NULL, start.alpha=NULL,
                      boot = 0, full.matrix = TRUE)

summary(modnobufo)

plot(as.factor(data2$ecology) ~ (data2$lung))



