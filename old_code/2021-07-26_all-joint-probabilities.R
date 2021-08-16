# fit the OU likelihood to the stochastic mapping
setwd("~/2020_houwie/")
require(geiger)
require(corHMM)
require(OUwie)
require(partitions)
require(expm)
require(MASS)
require(phytools)
source("hOUwieNode.R")


balenced_tree <- read.tree(text = "((S1:0.5,S2:0.5):0.5,(S3:0.5,S4:0.5):0.5);")
pectinate_tree <- read.tree(text = "(((S1:0.33, S2:0.33):0.33,S3:0.66):0.34,S4:1);")

plot(balenced_tree)
plot(pectinate_tree)

nTip <- 4
possible_datasets <- sapply(1:16, function(x) as.numeric(intToBits(x)[1:4]))

time_slice <- 1.1
rate.cat <- 1
rate <- 1
alpha = c(4,4)
sigma.sq = c(1,1)
theta0 = 5
theta = c(5,10)
root <- c(1,0)
Q <- matrix(c(-rate[1], rate[1],rate[1],-rate[1]),2,2)
# Q <- equateStateMatPars(getFullMat(list(getRateCatMat(2), getRateCatMat(2)), getRateCatMat(2)), 1:6)
# Q[Q > 0] <- rate
# diag(Q) <- -rowSums(Q)
data <- data.frame(sp = c("S1", "S2", "S3", "S4"), reg = possible_datasets[,1]+1, x = c(5, 5, 10, 10))
out <- getAllJointProbs(balenced_tree, data, rate.cat, time_slice, Q, alpha, sigma.sq, theta)
# model <- equateStateMatPars(getFullMat(list(getRateCatMat(2), getRateCatMat(2)), getRateCatMat(2)), 1:6)
corHMM(phy, data[,c(1,2)], 1, model = "ARD", p = rate[c(2,1)], root.p = c(0.5,0.5))
log(sum(exp(out[,1])))

plot(data$disc, data$cont, cex=25+data[,3], pch=20, col = rgb(0,0,0,0.5))
#text(label=rownames(data), data[,1], data[,2], col="red")

# library(plyr)
# library(ggplot2)
# data <- as.data.frame(out)
# data$disc_prob <- exp(data$disc)
# data$disc_weight <- data$disc_prob/sum(data$disc_prob)
# results <- data.frame()
# for (nmodels in sequence(nrow(data)-1)) {
#   for (rep in sequence(20)) {
#     for (type in sequence(2)) {
#       if(type==1) {
#         focal_rows <- sample(sequence(nrow(data)), size=nmodels, prob=data$disc_weight, replace=FALSE)
#       } else {
#         focal_rows <- sample(sequence(nrow(data)), size=nmodels, prob=NULL, replace=FALSE)
#       }
#       results <- plyr::rbind.fill(results, data.frame(nmodel=nmodels, type=ifelse(type==1,"discrete-weighted", "flat-weighted"), total_lnl=log(sum(exp(data$total[focal_rows])))))
#     }
#   }
# }
# ggplot(results, aes(x=nmodel, y=total_lnl)) + geom_point(aes(colour=factor(type), alpha=0.3)) + facet_wrap(~type)

out_CD <- getAllJointProbs(phy, data, rate.cat, time_slice, Q, alpha, sigma.sq, theta)
out_CID_1 <- getAllJointProbs(phy, data, rate.cat, time_slice, Q, alpha, sigma.sq, c(7.5,7.5,7.5,7.5))
out_CID_2 <- getAllJointProbs(phy, data, rate.cat, time_slice, Q, alpha, sigma.sq, theta[c(1,1,2,2)])

log(sum(exp(out_CD[,3])))
log(sum(exp(out_CID_1[,3])))
log(sum(exp(out_CID_2[,3])))

true_index <- sort(out_true[,3], index.return = TRUE)$ix
out_true <- data.frame(out_true[true_index,])
out_slow <- data.frame(out_slow[true_index,])
out_fast <- data.frame(out_fast[true_index,])
out_true$disc_weight <- exp(out_true$disc)/(exp(out_true$disc) + exp(out_true$cont))
out_slow$disc_weight <- exp(out_slow$disc)/(exp(out_slow$disc) + exp(out_slow$cont))
out_fast$disc_weight <- exp(out_fast$disc)/(exp(out_fast$disc) + exp(out_fast$cont))
out_true$type <- "true"
out_slow$type <- "slow"
out_fast$type <- "fast"
n_maps <- dim(out_true)[1]
out_true$map_no <- 1:n_maps
out_slow$map_no <- 1:n_maps
out_fast$map_no <- 1:n_maps

log(sum(exp(out_true$total)))
log(sum(exp(out_fast$total)))
log(sum(exp(out_slow$total)))

out <- rbind(out_true, out_slow, out_fast)
head(out)

ggplot(out, aes(x = map_no, y = total)) + geom_point(aes(colour = disc_weight, alpha=0.2)) + facet_wrap(~type)





