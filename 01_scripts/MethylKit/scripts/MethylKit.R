suppressPackageStartupMessages(library(methylKit))

path_to_config_file<-"config.tsv"

# Agruments are:
# [1] Number of samples in Group A
# [2] Number of samples in Group B
# [3] Outputname (File will contain DMRs)
# [4]..[n+1] Input Files 

#############################
# Read in Script Parameters #
#############################
#zz <- file("all.Rout", open = "wt")
#sink(zz)
#sink(zz, type = "message")

args <- commandArgs(TRUE)

# Check if number of arguments is correct
if (length(args)<5) {
        stop("Please enter exactly 1 output path to write the results and input files", call.=FALSE)
} else if (length(args)>=5) {
	num_A <- as.integer(args[1])
        num_B <- as.integer(args[2])	
       	outputPath <- as.character(args[3])
	args_files <- args
	args <- as.list(args)
	inputFiles <- args[4:length(args)]
}

##################################
# Set Tool Parameters to Default #
##################################

# ---- methRead ----
assembly<-"hg38"
dbtype<-NA
header=TRUE
skip=0
sep="\t"
context="CpG"
resolution="region" # Default "base"
treatment=c(1,1,1,0,0,0)
mincov=10

# ---- unite ----
destrand=FALSE
min.per.group=3L #Depends on numer of samples
chunk.size=1000000
mc.cores=1

# ---- tileMethylCounts ----
win.size=1000
step.size=1000
nCpG=0
###############################################
# Overwrite Default values if Tool Parameters #
###############################################
parameters <- read.delim(path_to_config_file, header= FALSE, sep="\t", comment.char="#")
colnames(parameters) <- c("parameter", "value")

options(digits=5)

# Create Treatment Vector and Sample ID Vector
treat <- vector("integer")
contr <- vector("integer")
sample.A <- list("character")
sample.B <- list("character")
for (i in 1:num_A){
        treat[i] = 1
	sample.A[i] = paste("test", i, sep="_")
}
for (i in 1:num_B){
        contr[i] = 0
        sample.B[i] = paste("ctrl", i, sep="_")
}
treatment <- c(treat, contr)
sample.id <- c(sample.A, sample.B)
# Read in other parameters from config file
mincov		<- strtoi(parameters[parameters$parameter == "Minimum Read Coverage", ]$value)
min.per.group	<- strtoi(parameters[parameters$parameter == "Minimum Sample Number per Group", ]$value)
mc.cores	<- strtoi(parameters[parameters$parameter == "Number of Threads", ]$value)
win.size	<- strtoi(parameters[parameters$parameter == "Tiling Window Size", ]$value)
step.size	<- strtoi(parameters[parameters$parameter == "Tiling Window Step Size", ]$value)
nCpG 		<- strtoi(parameters[parameters$parameter == "Minimum Number of CpGs", ]$value)
message("Minimum Sample Number per Group: ", 	min.per.group)
message("Cores: ", 				mc.cores)
message("Tiling Window Size: ", 		win.size)
message("Tiling Window Step Size: ", 		step.size)
message("Minimum Number of CpGs: ",		nCpG)

print(inputFiles)

myobj=methRead(	location=inputFiles,
		sample.id=sample.id,
           	assembly=assembly,
		dbtype=dbtype,
		header=header,
		skip=skip,
		sep=sep,
		mincov=mincov,
		context=context,
		resolution=resolution,
           	treatment=treatment)


filtered.myobj=filterByCoverage(myobj,lo.count=mincov,lo.perc=NULL,hi.count=NULL,hi.perc=100)


# Generate Plots and Statistics
#getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
#getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
#getMethylationStats(myobj[[3]],plot=TRUE,both.strands=FALSE)
#getMethylationStats(myobj[[4]],plot=TRUE,both.strands=FALSE)
#getMethylationStats(myobj[[5]],plot=TRUE,both.strands=FALSE)
#getMethylationStats(myobj[[6]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[3]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[4]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[5]],plot=TRUE,both.strands=FALSE)
#getCoverageStats(myobj[[6]],plot=TRUE,both.strands=FALSE)

# Merge samples
meth=unite(	object=filtered.myobj, 
		destand=destrand, 
		min.per.group=min.per.group, 
		chunk.size=chunk.size,
		mc.cores=mc.cores)

# Generate Plots and Statistics
getCorrelation( meth, plot=TRUE)
clusterSamples( meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)

# Create Tiles of length 1000
tiles=tileMethylCounts(	object=meth, 
			win.size=win.size, 
			step.size=step.size,
			cov.bases=nCpG,
			mc.cores=mc.cores)
#te<-getData(tiles)
#write.table(tiles, file = "tiles", sep = "\t", col.names = NA, quote = FALSE)
#cat("TILES\n")

# Try unite here
#tilesu=unite(   object=tiles,
#                destand=destrand,
#                min.per.group=min.per.group,
#                chunk.size=chunk.size,
#                mc.cores=mc.cores)

# Finding differentially methylated bases
diffMeth=calculateDiffMeth(tiles, mc.cores=mc.cores)

# Only needed to get subset of data
#dmrs=getMethylDiff(diffMeth)

results<-getData(diffMeth)
write.table(results, file = outputPath, sep = "\t", col.names = NA, quote = FALSE)


#sink(type = "message")
#sink()

