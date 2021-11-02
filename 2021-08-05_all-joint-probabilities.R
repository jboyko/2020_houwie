setwd("~/2020_houwie/")
require(geiger)
require(corHMM)
require(OUwie)
require(partitions)
require(expm)
require(MASS)
require(phytools)
require(ggtree)
require(ggplot2)
require(reshape2)
require(gridExtra)
require(ggplotify)
source("hOUwieNode.R")
source("Utils.R")


# balenced_tree <- read.tree(text = "((S1:0.5,S2:0.5):0.5,(S3:0.5,S4:0.5):0.5);")
pectinate_tree <- read.tree(text = "(((S1:0.33, S2:0.33):0.33,S3:0.66):0.34,S4:1);")

plot(pectinate_tree)

# generate all possible discrete datasets
nTip <- 4
possible_discrete <- sapply(1:16, function(x) as.numeric(intToBits(x)[1:4]))
possible_continuous <- (possible_discrete * 5) + 5
possible_combinations <- expand.grid(Disc=1:16, Cont=1:16)

# parameters for estimation (the MLE is an alpha that goes to infinity and sigma goes to 0)
time_slice <- 1.1
rate.cat <- 1
rate <- 1
alpha = c(4,4,4,4)
sigma.sq = c(1,1,1,1)
theta0 = 5
theta = c(5,5,10,10)
root <- c(1,0)
Q1 <- matrix(c(-rate[1], rate[1],rate[1],-rate[1]),2,2)
Q2 <- equateStateMatPars(getFullMat(list(getRateCatMat(2), getRateCatMat(2)), getRateCatMat(2)), 1:6)
Q2[is.na(Q2)] <- 0
Q2[Q2 > 0] <- rate
diag(Q2) <- -rowSums(Q2)

# a function for examing joint probs for all possible datasets
getJoint <- function(poss_disc){
  data <- data.frame(sp = c("S1", "S2", "S3", "S4"), reg = poss_disc, x = c(5,10)[poss_disc + 1])
  out_CD <- getAllJointProbs(pectinate_tree, data, 1, time_slice, Q1, alpha[1:2], sigma.sq[1:2], theta[2:3])
  return(out_CD)
}

full_joint_res <- apply(possible_discrete[,1:14], 2, getJoint)
all_joint_probs <- as.data.frame(do.call(rbind, lapply(full_joint_res, function(x) x$joint_probability_table)))
all_joint_probs$dataset <- as.factor(c("A", "B")[c(rep(1, 3*8), rep(2, 1*8), rep(1, 10*8))])

### results figure (C)
head(all_joint_probs)
C <- ggplot(all_joint_probs, aes(x = disc, y = cont, color = total)) +
  scale_color_gradient(low="#542788", high="#f1a340")+
  labs(x= "Discrete ln(Likelihood)", y="Continuous ln(Likelihood)") +
  geom_smooth(method=lm,  linetype="dashed", color="black", fill="grey") +
  geom_point(shape = 19, size = 3) +
  geom_point(data = all_joint_probs[all_joint_probs$dataset == "B",], shape = 1,size = 5, colour = "black") +
  theme_classic()


### concept figure (A)
out_CD <- full_joint_res[[3]]
poss_disc <- possible_discrete[,3]
data <- data.frame(sp = c("S1", "S2", "S3", "S4"), reg = poss_disc, x = c(5,10)[poss_disc + 1])
# out_CD <- getAllJointProbs(pectinate_tree, data, 1, time_slice, Q1, alpha[1:2], sigma.sq[1:2], theta[2:3])
# out_CID <- getAllJointProbs(pectinate_tree, data, 1, time_slice, Q1, alpha[1:2], sigma.sq[1:2], c(7.5,7.5))
# out_CIDx <- getAllJointProbs(pectinate_tree, data, 2, time_slice, Q2, alpha, sigma.sq, theta)
# log(sum(exp(out_CD$joint_probability_table[,3])))
# log(sum(exp(out_CID$joint_probability_table[,3])))
# log(sum(exp(out_CIDx$joint_probability_table[,3])))


out <- vector("list", 8)

i <- 1
for(i in 1:8){
  stochastic_map <- out_CD$simmap_list[[i]]
  
  cols<-setNames(RColorBrewer::brewer.pal(6, "BrBG")[c(1,3)], c("1","2"))
  g <- as.grob(~plotSimmap(stochastic_map, color=cols, fsize = 1e-10, lwd = 5, ylim=c(-0.25,4.5), outline = TRUE))
  data$reg <- as.factor(data$reg)
  data$x <- as.numeric(data$x)
  if(i < 5){
    xlab = ""
  }else{
    xlab = "Cont. Value"
  }
  
  p <- ggplot(data, aes(x=sp, y=x, fill = reg)) + 
    geom_bar(stat="identity", colour="black") + 
    coord_flip() + 
    theme_classic() +
    theme(axis.line.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),legend.position="none") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(6, "BrBG")[c(6,4)]) +
    scale_y_continuous(n.breaks = 3) +
    labs(y=xlab)
  
  g_phy <- grid.arrange(g, p, ncol=2, widths=c(2,1))
  
  if(i < 5){
    ylab = ""
  }else{
    ylab = "Total -Ln(Likelihood)"
  }
  
  stacked_dat <- melt(-out_CD$joint_probability_table[,c(1,2)])
  stacked_dat$Var1 <- as.factor(stacked_dat$Var1)
  stacked_dat$Var2 <- factor(stacked_dat$Var2, levels = c("cont", "disc"))
  new_dat <- stacked_dat[stacked_dat$Var1 == i,]
  g_bar <- ggplot(new_dat, aes(fill=Var2, y=value, x = Var1)) +
    geom_bar(position="stack", stat="identity", color="black") +
    coord_flip() +
    scale_fill_manual(values = RColorBrewer::brewer.pal(6, "BrBG")[c(5,2)]) +
    ylab(ylab) + 
    xlab("") +
    ylim(c(0,50)) +
    theme_classic() +
    theme(axis.line.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),legend.position="none")
  
  out[[i]] <- grid.arrange(g_phy, g_bar, ncol=1, heights=c(3,1))
}

A <- grid.arrange(out[[1]], out[[2]], out[[3]], out[[4]], 
             out[[5]], out[[6]], out[[7]], out[[8]],
             ncol=4)

#### figure b

out_CD <- full_joint_res[[5]]
poss_disc <- possible_discrete[,5]
data <- data.frame(sp = c("S1", "S2", "S3", "S4"), reg = poss_disc, x = c(5,10)[poss_disc + 1])
out <- vector("list", 8)

for(i in 1:8){
  stochastic_map <- out_CD$simmap_list[[i]]
  
  cols<-setNames(RColorBrewer::brewer.pal(6, "BrBG")[c(1,3)], c("1","2"))
  g <- as.grob(~plotSimmap(stochastic_map, color=cols, fsize = 1e-10, lwd = 5, ylim=c(-0.25,4.5), outline = TRUE))
  data$reg <- as.factor(data$reg)
  data$x <- as.numeric(data$x)
  if(i < 5){
    xlab = ""
  }else{
    xlab = "Cont. Value"
  }
  
  p <- ggplot(data, aes(x=sp, y=x, fill = reg)) + 
    geom_bar(stat="identity", colour="black") + 
    coord_flip() + 
    theme_classic() +
    theme(axis.line.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),legend.position="none") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(6, "BrBG")[c(6,4)]) +
    scale_y_continuous(n.breaks = 3) +
    labs(y=xlab)
  
  g_phy <- grid.arrange(g, p, ncol=2, widths=c(2,1))
  
  if(i < 5){
    ylab = ""
  }else{
    ylab = "Total -Ln(Likelihood)"
  }
  
  stacked_dat <- melt(-out_CD$joint_probability_table[,c(1,2)])
  stacked_dat$Var1 <- as.factor(stacked_dat$Var1)
  stacked_dat$Var2 <- factor(stacked_dat$Var2, levels = c("cont", "disc"))
  new_dat <- stacked_dat[stacked_dat$Var1 == i,]
  g_bar <- ggplot(new_dat, aes(fill=Var2, y=value, x = Var1)) +
    geom_bar(position="stack", stat="identity", color="black") +
    coord_flip() +
    scale_fill_manual(values = RColorBrewer::brewer.pal(6, "BrBG")[c(5,2)]) +
    ylab(ylab) + 
    xlab("") +
    ylim(c(0,50)) +
    theme_classic() +
    theme(axis.line.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),legend.position="none")
  
  out[[i]] <- grid.arrange(g_phy, g_bar, ncol=1, heights=c(3,1))
}

B <- grid.arrange(out[[1]], out[[2]], out[[3]], out[[4]], 
                  out[[5]], out[[6]], out[[7]], out[[8]],
                  ncol=4)



A <- grid.arrange(A, B, ncol=2)
grid.arrange(A, C, ncol=2)
grid.arrange(B, C, ncol=2)



