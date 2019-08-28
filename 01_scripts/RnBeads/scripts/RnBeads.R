cat("Load library RnBeads\n")
suppressPackageStartupMessages(library(RnBeads))
suppressPackageStartupMessages(library(RnBeads.hg38))
cat("Done\n\n")

# Path to file containing the tool's parameters
path_to_config_file<-"config.tsv"

# Set RnBeads parameters to defauls
coverage <- 10
threshold <- 1

# Directory where the bed files and the sample annotation file is located
bed.dir <- "02_data/RnBeads/converted"

# Path to sample annotation file
sample.annotation <- "02_data/RnBeads/sample_annotation.csv"

# Define Data Source
data.source <- c(bed.dir, sample.annotation)

# Directory where the output should be written to
results.dir <- "03_results/RnBeads"

# Directory where the analysis and report files should be written to
analysis.dir <- file.path(results.dir, "analysis")
report.dir <- file.path(results.dir, "reports")

# Set analysis options
rnb.options(import.default.data.type="bed.dir")

# Define input format (column order,...)
rnb.options(import.bed.style="BisSNP")
rnb.options(differential.comparison.columns="Sample_Group")

# Choose genome. Default is hg19
rnb.options(assembly="hg38")

# Overwrite Default Values of Parameters
parameters <- read.delim(path_to_config_file, header= FALSE, sep="\t", comment.char="#")
colnames(parameters) <- c("parameter", "value")

options(digits=5)
coverage	<- strtoi(parameters[parameters$parameter == "Minimum Read Coverage", ]$value)
threshold	<- as.double(toString(parameters[parameters$parameter == "Max Quantile of NAs per Site", ]$value))

cat(sprintf("\nCoverage: %s\n", coverage))
cat(sprintf("\nThreshold: %s\n", threshold))

# Default tiles have a size of 5000 kb. Here 1000 / 500 bps are used
REGION_SET <- c("tiling1kb", "tiling500bp", "tiling200bp")
ASSEMBLY <- "hg38"
rnb.load.annotation.from.db(REGION_SET, assembly=ASSEMBLY)
# Here you can choose which output files you want (region.types are the 4 default types: tiling=5000, genes, promotors, cpgislands)
#rnb.options(region.types=union(rnb.getOption("region.types"), REGION_SET))
rnb.getOption("region.types")

rnb.options(region.types=c("tiling", "genes", "promotors", "cpgislands", "tiling1kb", "tiling500bp", "tiling200bp"))
print("PEGION TYPES")
rnb.getOption("region.types")

# Import data (less runtime as no reports)
imp <- rnb.execute.import(data.source=data.source, data.type="bs.bed.dir",  dry.run=FALSE, verbose=TRUE)

print("# Lines before Coverage Filter")
nrow(meth(imp))

# Try more cores
#setModuleNumCores(imp, 10L)

# Replace low coverage sites with NA
imp.filtered <- rnb.execute.low.coverage.masking(imp, coverage)$dataset

# More cores, geht nicht mehr seit neuer conda env
#setModuleNumCores(imp.filtered, 10L)

# Remove NA entries
filtered.set.noNA <- rnb.execute.na.removal(imp.filtered, threshold)$dataset

print("# Sites after Coverage Filter")
nrow(meth(filtered.set.noNA))

rnb.set <- filtered.set.noNA

# Try more cores
#setModuleNumCores(rnb.set, 10L)

# Run differential methylation analysis (dont use rnb.executeDiffMeth as no results are written to folder)
rnb.run.differential(rnb.set=rnb.set, dir.reports=report.dir)

# Catch warnings
cat(sapply(warnings(), toString) , file="rnbeads_warnings", sep="\n")



# RUN 1, not working (groups not found)
# Import Data
#imp <- rnb.run.import(data.source=data.source, dir.reports=report.dir, data.type="bs.bed.dir")
#rnb.set <- imp$rnb.set
# Run differential methylation analysis alternative (less runtime as no reports)
#rnb.execute.computeDiffMeth(rnb.set, "Sample_Group", c("A", "B"))

# RUN 2 Working, runtime to high
#rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation, data.dir=bed.dir, data.type="bs.bed.dir")

# No effect
#rnb.options(filtering.coverage.threshold=0)
