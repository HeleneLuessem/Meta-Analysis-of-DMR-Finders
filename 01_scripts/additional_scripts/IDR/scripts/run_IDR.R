

if (!require(idr)) install.packages('idr')
library(idr)

#data <- read.table("temp/Merged_Regions_Scores.tsv")
#data$V1 <- NULL
#data$V2 <- NULL
#data$V3 <- NULL
#ranks <- rev(seq(-nrow(data), -1, by=1))
#data$V4 <- (ranks * as.numeric(data$V4 == 0)) + data$V4
#data$V5 <- (ranks * as.numeric(data$V5 == 0)) + data$V5
#data$V6 <- (ranks * as.numeric(data$V6 == 0)) + data$V6
#data$V7 <- (ranks * as.numeric(data$V7 == 0)) + data$V7
#data$V8 <- (ranks * as.numeric(data$V8 == 0)) + data$V8

#write.table(data, "IDR_Input_2.tsv", sep="\t", col.names=F, row.names=F, quote=F)

#mu <- 1 #2.6
#sigma <- 0.1 #1.3
#rho <- 0.9#0.8
#p <- 0.7#0.7
#idr.out <- est.IDR(data, mu, sigma, rho, p, eps=0.001, max.ite=200)



p <- 0.7
rho <- 0.8
sigma <- 1.3
mu <- 2.6

data <- read.table("temp/Merged_Regions_Scores.tsv")
ranks1 <- -sample(nrow(data))
ranks2 <- -sample(nrow(data))
ranks3 <- -sample(nrow(data))
ranks4 <- -sample(nrow(data))
ranks5 <- -sample(nrow(data))
data$V4 <- (ranks1 * as.numeric(data$V4 == 0)) + data$V4
data$V5 <- (ranks2 * as.numeric(data$V5 == 0)) + data$V5
data$V6 <- (ranks3 * as.numeric(data$V6 == 0)) + data$V6
data$V7 <- (ranks4 * as.numeric(data$V7 == 0)) + data$V7
data$V8 <- (ranks5 * as.numeric(data$V8 == 0)) + data$V8
data$V1 <- NULL
data$V2 <- NULL
data$V3 <- NULL
idr <- est.IDR(data, mu, sigma, rho, p, eps=0.001, max.ite=600)
data$IDR <- idr$IDR
data$idr <- idr$idr
data$NUM <- c(1:nrow(data))

idr$para


current <- Sys.time()
#fileName <- paste("IDR_out-",current,".tsv",sep="")
fileName <- "IDR_out.tsv"
write.table(data, fileName, sep='\t', row.names=F, quote=F)

