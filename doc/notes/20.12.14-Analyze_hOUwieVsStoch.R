# this script will analyze whether hOUwie can preform better than independent corHMM > stoch Map > OUwie

require(viridis)
require(ggplot2)

setwd("~/2020_hOUwie/LiamTestResults/")

load("20_12_09-Test_q0.5_hOUwieVsLiam.R")
df.Sig2_1 <- df.Q1_Sig2_1
df.Sig2_2 <- df.Q1_Sig2_2
load("20_12_10-Test_q1_hOUwieVsLiam.R")
df.Sig2_1 <- rbind(df.Sig2_1, df.Q1_Sig2_1)
df.Sig2_2 <- rbind(df.Sig2_2, df.Q1_Sig2_2)
load("20_12_10-Test_q2_hOUwieVsLiam.R")
df.Sig2_1 <- rbind(df.Sig2_1, df.Q1_Sig2_1)
df.Sig2_2 <- rbind(df.Sig2_2, df.Q1_Sig2_2)
load("20_12_11-Test_q4_hOUwieVsLiam.R")
df.Sig2_1 <- rbind(df.Sig2_1, df.Q1_Sig2_1)
df.Sig2_2 <- rbind(df.Sig2_2, df.Q1_Sig2_2)
load("20_12_12-Test_q8_hOUwieVsLiam.R")
df.Sig2_1 <- rbind(df.Sig2_1, df.Q1_Sig2_1)
df.Sig2_2 <- rbind(df.Sig2_2, df.Q1_Sig2_2)
df.Sig2_ratio <- df.Sig2_2
df.Sig2_ratio[,3] <- df.Sig2_2[,3]/df.Sig2_1[,3]
df.Sig2_ratio <- df.Sig2_ratio[df.Sig2_ratio[,3] < 50,]

cols <- viridis(4)

ggplot(df.Sig2_1, aes(x=par, y=val, fill=fit)) + 
  labs(x = "Transition Rate", y = expression('estimated '*sigma[1]^2)) +
  scale_fill_manual(values=cols) + 
  theme(text = element_text(size = 20)) + 
  geom_boxplot()

ggplot(df.Sig2_2, aes(x=par, y=val, fill=fit)) + 
  labs(x = "Transition Rate", y = expression('estimated '*sigma[2]^2)) +
  scale_fill_manual(values=cols) + 
  theme(text = element_text(size = 20)) + 
  geom_boxplot()

jpeg(filename = "~/2020_hOUwie/doc/hOUwieVsStoch.jpg")
ggplot(df.Sig2_ratio, aes(x=par, y=val, fill=fit)) + 
  labs(x = "Transition Rate", y = expression('estimated '*sigma[2]^2/sigma[1]^2)) +
  scale_fill_manual(values=cols) + 
  theme(text = element_text(size = 20)) + 
  geom_boxplot()
dev.off()

Qlist <- list()
Qrates <- unique(df.Sig2_ratio[,1])
Models <- unique(df.Sig2_ratio[,2])
SDRatio <- MeanRatio <- matrix(0, length(Qrates), length(Models), dimnames = list(Qrates, Models))
for(i in 1:length(unique(Qrates))){
  Qlist_i <- df.Sig2_ratio[df.Sig2_ratio[,1] == Qrates[i],]
  for(j in 1:length(Models)){
    Ratio_j <- Qlist_i[Qlist_i[,2] == Models[j], 3]
    MeanRatio[i, j] <- median(Ratio_j)
    SDRatio[i, j] <- max(Ratio_j)
  }
}

MedianTablePercent <- as.data.frame(round(apply(MeanRatio, 1, function(x) (x-10)/10*100), 2))
write.table(MedianTablePercent, "~/2020_hOUwie/doc/MedianTablePercent.csv")

MedianTablePercent <- as.data.frame(round(apply(MeanRatio, 1, function(x) (x-10)/10*100), 2))









### for the grant
head(df.Sig2_ratio)
grant.df <- df.Sig2_ratio
grant.df <- grant.df[!grant.df$par == "q_8", ]
grant.df <- grant.df[!grant.df$par == "q_1", ]
grant.df <- grant.df[!grant.df$fit == "hOUwTru",]
grant.df$fit[grant.df$fit == "TrueMap"] <- "True Model"
grant.df$fit[grant.df$fit == "StocMap"] <- "Existing Methodology"
grant.df$fit[grant.df$fit == "hOuwWt"] <- "hOUwie Estimate"
grant.df$par[grant.df$par == "q_0.5"] <- "Low"
grant.df$par[grant.df$par == "q_2"] <- "Medium"
grant.df$par[grant.df$par == "q_4"] <- "High"
grant.df$par <- factor(grant.df$par, levels = c("Low", "Medium", "High"), ordered = TRUE)
grant.df$fit <- factor(grant.df$fit, levels = c("True Model", "Existing Methodology", "hOUwie Estimate"), ordered = TRUE)
cols <- viridis(4)[2:4]
colnames(grant.df)[2] <- "Model"

pdf(file = "~/GetMoney/grants/2021-GPSC-hOUwie/hOUwieFig1.pdf", width = 7.5, height = 7.5)
ggplot(grant.df, aes(x=par, y=val, fill=Model)) + 
  labs(x = "Simulated Rate of \nDiscrete Variable Change", y = "Estimated Rate of \nContinuous Variable Change") +
  scale_fill_manual(values=cols, labels=c("\nTrue Model\n", "\nExisting\nMethodology\n", "\nhOUwie\nEstimiate\n")) + 
  theme(text = element_text(size = 18)) + 
  geom_boxplot()
dev.off()

