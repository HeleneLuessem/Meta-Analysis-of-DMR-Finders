suppressPackageStartupMessages(library(bsseq))

path_to_config_file <- "config.tsv"

# Agruments are:
# [1] File with list of files from group A
# [2] File with list of files from group B
# [3] Outputname (File will contain DMRs)

args <- commandArgs(TRUE)

# Check if number of arguments is correct
if (length(args)!=3) {
        stop("Please select exactly 2 files (one for each group) containing the list of input files and provide an output name!", call.=FALSE)
} else if (length(args)==3) {
        groupA <- as.character(args[1])
        groupB <- as.character(args[2])
        outputPath <- as.character(args[3])
}

# Check if input files exist
if(!file.exists(groupA)){
        stop(paste("File with list of input files from group A", groupA, "does not exist! Please pass a valid file."), call.=FALSE)
}
if(!file.exists(groupB)){
        stop(paste("File with list of input files from group B", groupB, "does not exist! Please pass a valid file."), call.=FALSE)
}

# Set Parameters to Default Values
# BSmooth
ns	<- 70 		# Minimum number of methylation loci in a smoothing window
h 	<- 1000		# Minimum smoothing window, in bases
maxGapL	<- 100000000	# maximum gap between two methylation loci, before the smoothing is broken across the gap
mc.cores <- 1

# tstat
local.correct <- TRUE

# dmrFinder
maxGapD <- 300

# Overwrite Default values
parameters <- read.delim(path_to_config_file, header= FALSE, sep="\t", comment.char="#")
colnames(parameters) <- c("parameter", "value")

options(digits=5)
ns		<- strtoi(parameters[parameters$parameter == "Min Number of Loci in Window", ]$value)
h		<- strtoi(parameters[parameters$parameter == "Min Smoothing Window", ]$value)
maxGapL		<- strtoi(parameters[parameters$parameter == "Maximum Gap Loci", ]$value)
local.correct	<- as.logical(parameters[parameters$parameter == "Local Correction", ]$value)
maxGapD		<- strtoi(parameters[parameters$parameter == "Maximum Distance between two CpGs in one DMR", ]$value)
mc.cores	<- strtoi(parameters[parameters$parameter == "Number of Threads", ]$value)

cat(sprintf("\nMin Number of Loci in Window: %s\n", ns))
cat(sprintf("\nMin Smoothing Window %s\n", h))
cat(sprintf("\nMaximum Gap Loci %s\n", maxGapL))
cat(sprintf("\nLocal Correction %s\n", local.correct))
cat(sprintf("\nMaximum Gap DMRs %s\n", maxGapD))
cat(sprintf("\nCores %s\n", mc.cores))

# Iterate over file to read out input files for a given group
readListOfFiles = function(filepath, group, sample.names) {
	i = 0
	f = file(filepath, "r")
	while( TRUE ) {
		print(line)
		line = readLines(f, n = 1)
		if ( length(line) > 0 ) {
			if (i >= 1){
				i = i + 1
				dat <- read.table(line, skip = 0, row.names = NULL, col.names = c("chr", "from", "to", "methPerc", "coverage", "strand", "NA", "NA","NA","NA","NA"), colClasses = c("character", "integer", "integer", "character", "character", "character", "character", "character","character","character","character"))[,c('chr','from','strand','methPerc','coverage')]
				dat <- transform(dat,methPerc= as.integer(0.5 + as.double(methPerc) * as.double(coverage)))
				dat$strand <- "+"
				dat <- transform(dat, coverage=as.integer(coverage))
				BS.new <- BSseq(pos = dat$from, chr = dat$chr, M = as.matrix(dat$methPerc, ncol = 1), Cov = as.matrix(dat$coverage, ncol = 1), sampleNames = "all")
				sampleNames(BS.new) <- paste(group, as.character(i), sep="")
				sample.names[i] <- paste(group, as.character(i), sep="")
				BS.samples <- combine(BS.samples, BS.new)
			} else {
				dat <- read.table(line, skip = 0, row.names = NULL, col.names = c("chr", "from", "to", "methPerc", "coverage", "strand", "NA", "NA","NA","NA","NA"), colClasses = c("character", "integer", "integer", "character", "character", "character", "character", "character","character","character","character"))[,c('chr','from','strand','methPerc','coverage')]
				dat <- transform(dat,methPerc= as.integer(0.5 + as.double(methPerc) * as.double(coverage)))
				dat$strand <- "+"
				dat <- transform(dat, coverage=as.integer(coverage))	
				BS.samples <- BSseq(pos = dat$from, chr = dat$chr, M = as.matrix(dat$methPerc, ncol = 1), Cov = as.matrix(dat$coverage, ncol = 1), sampleNames = "all")
				i = i + 1
				sampleNames(BS.samples) <- paste(group, as.character(i), sep="")
				sample.names[i] <- paste(group, as.character(i), sep="")
			}
		} else {
                        break
                }
	}
	close(f)
	validObject(BS.samples)
	return(list(BS.samples,sample.names))
}

# Sample Names
samples_A <- vector()
samples_B <- vector()

# Read in samples in BS objects
read_A <- readListOfFiles(groupA, "A", samples_A)
read_B <- readListOfFiles(groupB, "B", samples_B)
print(3)
BS.groupA <- read_A[[1]]
BS.groupB <- read_B[[1]]
samples_A <- read_A[[2]]
samples_B <- read_B[[2]]
print(4)
# Merge into one object
print(BS.groupA)
print(BS.groupB)
BS.samples <- combine(BS.groupA, BS.groupB)
print(BS.samples)

#Smooth Data
BS.smooth <- BSmooth(BSseq=BS.samples, ns=ns, h=h, maxGap=maxGapL, mc.cores=mc.cores) #BPPARAM=MulticoreParam(workers=mc.cores))

# Calculate t statistics
BS.tstat <- BSmooth.tstat(	BSseq=BS.smooth, 
				group1=samples_A, #c("A1", "A2", "A2"), #samples_A, 
				group2=samples_B, #c("B1", "B2", "B3")) #samples_B)
				local.correct=local.correct,
				mc.cores=mc.cores)
#print("Stats")
#print(getStats(BS.tstat))

# Finding DMRs
dmrs <- dmrFinder(bstat=BS.tstat, maxGap=maxGapD)
#print(dmrs)

# Store in file
write.table(dmrs, file = outputPath, sep = "\t", col.names = NA, quote = FALSE)

#head(granges(BS.samples), n = 4)
#save(BS.samples, file = "BS.samples.rda")


# --------------------------------------------------------------------------------------
# FIRST TRY - READ IN DATA MANUALLY
## Group A
#BS.r1 <- read.bed("s1")
#sampleNames(BS.r1) <- "C1"
#BS.r2 <- read.bed("s2")
#sampleNames(BS.r2) <- "C2"
#BS.r3 <- read.bed("s3")
#sampleNames(BS.r3) <- "C3"
#BS.samplesA <- combine(combine(BS.r1, BS.r2), BS.r3)
#print(BS.samplesA)
#pData(BS.samplesA)$Rep <- c("sample_A_1", "sample_A_2", "sample_A_3")
#validObject(BS.samplesA)

## Group B
#BS.r4 <- read.bed("s4")
#sampleNames(BS.r4) <- "E1"
#BS.r5 <- read.bed("s5")
#sampleNames(BS.r5) <- "E2"
#BS.r6 <- read.bed("s6")
#sampleNames(BS.r6) <- "E3"
#BS.samplesB <- combine(combine(BS.r4, BS.r5), BS.r6)
#print(BS.samplesB)
#pData(BS.samplesB)$Rep <- c("sample_B_1", "sample_B_2", "sample_B_3")
#validObject(BS.samplesB)

#pData(BS.samplesA)
#pData(BS.samplesB)
#BS.samples <- combine(BS.samplesA, BS.samplesB)
#print(BS.samples)

