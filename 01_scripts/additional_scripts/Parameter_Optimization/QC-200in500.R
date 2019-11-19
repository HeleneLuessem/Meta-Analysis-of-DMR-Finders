path_to_file_small<-"F:/Masterarbeit/Results Sets/Extension/WGBS/Parameter Optimization/Window Size/Plots/Small_in_Large/200_CF/RnBeads_DMRs_std_200_CF.tsv"
path_to_file_big<-"F:/Masterarbeit/Results Sets/Extension/WGBS/Parameter Optimization/Window Size/Plots/Small_in_Large/500_CF/RnBeads_DMRs_std_500_CF.tsv"



small<-200
big<-500


data_small <- read.table(path_to_file_small)
data_big <- read.table(path_to_file_big)

data_small$V4 <- NULL
data_small$V5 <- NULL
data_small$V6 <- NULL
data_small$V7 <- NULL
data_small$V8 <- NULL
#------------------
data_big$V4 <- NULL
data_big$V5 <- NULL
data_big$V6 <- NULL
data_big$V7 <- NULL
data_big$V8 <- NULL

colnames(data_small) <- c("Chr", "Start_small", "End_small")
colnames(data_big) <- c("Chr", "Start_big", "End_big")

# Set small windows to big ones
data_small$Start  <- (data_small$Start_small - data_small$Start_small %% big + 1)

small_end_200 <- data_small[data_small$End_small %% 500 == 200,]
small_end_200$End <- (small_end_200$End_small + 300)

small_end_400 <- data_small[data_small$End_small %% 500 == 400,]
small_end_400$End <- (small_end_400$End_small + 100)

small_end_100 <- data_small[data_small$End_small %% 500 == 100,]
small_end_100$End <- (small_end_100$End_small - 100)

small_end_300 <- data_small[data_small$End_small %% 500 == 300,]
small_end_300$End <- (small_end_300$End_small + 200)

small_end_0 <- data_small[data_small$End_small %% 500 == 0,]
small_end_0$End <- (small_end_0$End_small)

data_small <- do.call("rbind", list(small_end_200, small_end_400, small_end_100, small_end_300, small_end_0))


data_big$Start <- data_big$Start_big
data_big$End   <- data_big$End_big

# Join Data
merged_500_1000_inner <- merge(data_small,data_big,by=c("Chr", "Start","End"))
merged_500_1000_left  <- merge(data_small,data_big,by=c("Chr", "Start","End"), all.x=TRUE) 
merged_500_1000_right  <- merge(data_small,data_big,by=c("Chr", "Start","End"), all.y=TRUE) 
merged_500_1000_full  <- merge(data_small,data_big,by=c("Chr", "Start","End"), all=TRUE)

# A
tpm_1000 <- merged_500_1000_full
tpm_1000$Start <- NULL
tpm_1000$End <- NULL
tpm_1000$Start_small <- NULL
tpm_1000$End_small <- NULL
tpm_1000$NA.x <- NULL
tpm_1000$NA.y <- NULL
# Remove NA rows
tpm_1000 <- na.omit(tpm_1000, cols=c("Start_big", "End_big"))
# Mark Duplicates
tpm_1000$Dup <- duplicated(tpm_1000)

dup = tpm_1000[tpm_1000$Dup == T,]
triples <- nrow(tpm_1000[tpm_1000$Dup == T,]) - nrow(unique(dup))
doubles <- nrow(unique(dup))
# 500 as well as 2 200 found (A.1)
A_1 <- doubles
print("500 coverd by 2 200")
print(doubles)

# 500 as well as 2 200 found (A.2)
A_2 <- triples
print("500 coverd by 3 200")
print(triples)

# 500 and 1 200 found (B)
B <- nrow(merged_500_1000_right)-2*doubles-3*triples-nrow(subset(merged_500_1000_right, is.na(Start_small)))
print("500 coverd by 1 200")
nrow(merged_500_1000_right)-2*doubles-3*triples-nrow(subset(merged_500_1000_right, is.na(Start_small)))

# 500 found, but no 200 (C)
C <- nrow(subset(merged_500_1000_right, is.na(Start_small)))+A_2
print("500 coverd by 0 200")
nrow(subset(merged_500_1000_right, is.na(Start_small)))+A_2

# 200 found, but no 500 (D)
D <- nrow(subset(merged_500_1000_left, is.na(Start_big)))
print("200 coverd by 0 500")
nrow(subset(merged_500_1000_left, is.na(Start_big)))

print("500:")
print(A_1+A_2+B+C)

print("200:")
print(2*A_1+3*A_2+B+D)
