
path_to_file_500<-"F:/Masterarbeit/Results Sets/Extension/WGBS/Parameter Optimization/Window Size/Plots/Small_in_Large/500_CF/MethylKit_DMRs_std_500_CF.tsv"
path_to_file_1000<-"F:/Masterarbeit/Results Sets/Extension/WGBS/Parameter Optimization/Window Size/Plots/Small_in_Large/1000_CF/MethylKit_DMRs_std_1000_CF.tsv"

data_500 <- read.table(path_to_file_500)
data_1000 <- read.table(path_to_file_1000)

data_500$V4 <- NULL
data_500$V5 <- NULL
data_500$V6 <- NULL
data_500$V7 <- NULL
data_500$V8 <- NULL
#------------------
data_1000$V4 <- NULL
data_1000$V5 <- NULL
data_1000$V6 <- NULL
data_1000$V7 <- NULL
data_1000$V8 <- NULL

colnames(data_500) <- c("Chr", "Start_500", "End_500")
colnames(data_1000) <- c("Chr", "Start_1000", "End_1000")

data_500$Start  <- (data_500$Start - data_500$Start %% 1000 + 1)
data_500$End    <- (data_500$End + data_500$End %% 1000)
data_1000$Start <- data_1000$Start_1000
data_1000$End   <- data_1000$End_1000

# Join Data
merged_500_1000_inner <- merge(data_500,data_1000,by=c("Chr", "Start","End"))
merged_500_1000_left  <- merge(data_500,data_1000,by=c("Chr", "Start","End"), all.x=TRUE) 
merged_500_1000_right  <- merge(data_500,data_1000,by=c("Chr", "Start","End"), all.y=TRUE) 
merged_500_1000_full  <- merge(data_500,data_1000,by=c("Chr", "Start","End"), all=TRUE)

# A
tpm_1000 <- merged_500_1000_full
tpm_1000$Start <- NULL
tpm_1000$End <- NULL
tpm_1000$Start_500 <- NULL
tpm_1000$End_500 <- NULL
tpm_1000$NA.x <- NULL
tpm_1000$NA.y <- NULL
# Remove NA rows
tpm_1000 <- na.omit(tpm_1000, cols=c("Start_1000", "End_1000"))
# Mark Duplicates
tpm_1000$Dup <- duplicated(tpm_1000)
A <- nrow(tpm_1000[tpm_1000[, "Dup"] == TRUE,])
# 1000 as well as 2 500 found (A)
print("1000 coverd by 2 500")
nrow(tpm_1000[tpm_1000[, "Dup"] == TRUE,])

# 1000 and 1 500 found (B)
B <- nrow(merged_500_1000_right)-2*nrow(tpm_1000[tpm_1000[, "Dup"] == TRUE,])-nrow(subset(merged_500_1000_right, is.na(Start_500)))
print("1000 coverd by 1 500")
nrow(merged_500_1000_right)-2*nrow(tpm_1000[tpm_1000[, "Dup"] == TRUE,])-nrow(subset(merged_500_1000_right, is.na(Start_500)))

# 1000 found, but no 500 (C)
C <- nrow(subset(merged_500_1000_right, is.na(Start_500)))
print("1000 coverd by 0 500")
nrow(subset(merged_500_1000_right, is.na(Start_500)))

# 500 found, but no 1000 (D)
D <- nrow(subset(merged_500_1000_left, is.na(Start_1000)))
print("500 coverd by 0 1000")
nrow(subset(merged_500_1000_left, is.na(Start_1000)))

print("1000:")
print(A+B+C)

print("500:")
print(2*A+B+D)
  

