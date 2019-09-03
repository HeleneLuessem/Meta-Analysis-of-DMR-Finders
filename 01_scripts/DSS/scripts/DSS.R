.libPaths("/home/users/s9helues/miniconda2/envs/DSS/lib/R/library")
suppressPackageStartupMessages(library(DSS))

path_to_config_file<-"config.tsv"

# Set to Default
# dmlTest
equal.disp <- FALSE
smoothing <- TRUE
smoothing.span <- 500

#CallDML
#delta1 <- 0 # Default
p.threshold1 <- 0.001

# CallDMR
#delta2 <- 0.00001 # Default 0.1
p.threshold2 <- 0.00001 #0.01
minlen <- 50
minCG <- 3
dis.merge <- 100#50
pct.sig <- 0.5

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

# Read in parameter choices from the config file
parameters <- read.delim(path_to_config_file, header= FALSE, sep="\t", comment.char="#")
colnames(parameters) <- c("parameter", "value")

# Overwrite Default Values of Parameters
options(digits=5)
equal.disp	<- as.logical(parameters[parameters$parameter == "Equal Dispersion", ]$value)
smoothing	<- as.logical(parameters[parameters$parameter == "Smoothing", ]$value)
smoothing.span	<- strtoi(parameters[parameters$parameter == "Size of Smoothing Window [bp]", ]$value)
#delta1		<- as.double(toString(parameters[parameters$parameter == "DML Delta", ]$value))
p.threshold1	<- as.double(toString(parameters[parameters$parameter == "DML p threshold", ]$value))
#delta2		<- as.double(toString(parameters[parameters$parameter == "DMR Delta", ]$value))
p.threshol2 	<- as.double(toString(parameters[parameters$parameter == "DMR p threshold", ]$value))
minlen		<- strtoi(parameters[parameters$parameter == "Minimum DMR Length", ]$value)
minCG 		<- strtoi(parameters[parameters$parameter == "Minimum Number of CpGs", ]$value)
dis.merge	<- strtoi(parameters[parameters$parameter == "Maximum Distance between two CpGs in one DMR", ]$value)
pct.sig		<- as.double(toString(parameters[parameters$parameter == "Percentage of CG sites with significant p-values per DMR", ]$value))

cat(sprintf("\nEqual Dispersion: %s\n", equal.disp))
cat(sprintf("Smoothing: %s\n", smoothing))
cat(sprintf("Size of Smoothing Window [bp]: %s\n", smoothing.span))
#cat(sprintf("DML Delta: %s\n", delta1))
cat(sprintf("DML p threshold: %s\n", p.threshold1))
#cat(sprintf("DMR Delta: %s\n", delta2))
cat(sprintf("DMR p threshold: %s\n", p.threshol2))
cat(sprintf("Minimum Length of a DMR [bp]: %s\n", minlen))
cat(sprintf("Minimum Number of CpGs: %s\n", minCG))
cat(sprintf("Maximum Distance between two CpGs in one DMR: %s\n", dis.merge))
cat(sprintf("Percentage of CG sites with significant p-values per DMR: %s\n", pct.sig))


# Iterate over file to read out input files for a given group
readListOfFiles = function(filepath, group) {
	files <- list()
	sampleNames <- c()
	i = 1
	f = file(filepath, "r")
	while ( TRUE ) {
		line = readLines(f, n = 1)
		if ( length(line) > 0 ) {
	     		files[[i]] <- read.table(line, header=TRUE)
			sampleNames <- c(sampleNames, paste(group, i, sep = ""))
                        i = i + 1
	   	} else { 
			break
		}
	}
	close(f)
	return(list(files, sampleNames))
}

dataA <- readListOfFiles(groupA, "A")
dataB <- readListOfFiles(groupB, "B")

#cat("\n\n Input Files:\n")
inputFilesAB <- c(dataA[[1]], dataB[[1]])
#print(inputFilesAB)
#cat("\n\n Sample Names:\n")
sampleNamesAB <- c(dataA[[2]], dataB[[2]])
#print(sampleNamesAB)

BSobj <- makeBSseqData(inputFilesAB, sampleNamesAB )#[1:1000,]
#BSobj
#sampleNames(BSobj)
# Perform statistical test for DML. 3 Steps are performed
# (1) Estimate mean methylation levels for all CpG site (smoothing is possible)
# (2) Estimate dispersions at each CpG sites
# (3) Conduct Wald test
#cat("\nPerform statistical Test for DML\n")
#dmlTest <- DMLtest(BSobj, group1=dataA[[2]], group2=dataB[[2]], equal.disp=equal.disp, smoothing=smoothing, smoothing.span=smoothing.span)

zz <- file("all.Rout", open = "wt")
sink(zz)
sink(zz, type = "message")

dmlTest <- DMLtest(BSobj, group1=dataA[[2]], group2=dataB[[2]], equal.disp=equal.disp, smoothing=smoothing, smoothing.span=smoothing.span)

sink(type = "message")
sink()

#head(dmlTest)
#cat("Done\n\n")
#cat("Call DMLs\n")
dmls <- callDML(dmlTest, p.threshold=p.threshold1)
#head(dmls)
#cat("Done\n\n")

#cat("Call DMRs\n")
#dmrs <- callDMR(dmlTest, p.threshold=p.threshold2, minlen=minlen, minCG=minCG, dis.merge=dis.merge, pct.sig=pct.sig)
dmrs <- callDMR(dmlTest, p.threshold=p.threshold2, minlen=minlen, minCG=minCG, dis.merge=dis.merge, pct.sig=pct.sig)
#cat("called")
#head(dmrs)
write.table(dmrs, file = outputPath, sep = "\t", col.names = NA, quote = FALSE)
#cat(paste("Done. Results written to ", outputPath, "\n\n", sep = ""))

#Standards
#makeBSseqData( list(dat1.1, dat1.2, dat2.1, dat2.2), c("C1","C2", "N1", "N2") )[1:1000,]
#DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"))
#DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"), smoothing=TRUE)
#callDML(dmlTest, p.threshold=0.001)
#callDMR(dmlTest, p.threshold=0.01)
