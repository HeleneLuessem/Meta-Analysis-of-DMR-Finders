
args <- commandArgs(TRUE)

path <- as.character(args[1])

data <- read.table(path)

# Creates rank vector of the negative scores (-> low score gets low rank)

# Sort data
data1 <- data[order(data$V8),]

# Write to temp file
write.table(data1, "temp_ordered", sep="\t", col.names=F, row.names=F, quote=F)

# Read from file
data2 <- read.table("temp_ordered")

# Get row num and add as column
id <- rownames(data2)
data2$V9 <- NULL
data2$V8 <- id

write.table(data2, args[2], sep="\t", col.names=F, row.names=F, quote=F)

