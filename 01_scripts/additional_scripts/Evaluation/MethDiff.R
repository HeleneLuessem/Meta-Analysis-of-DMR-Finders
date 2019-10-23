args <- commandArgs(TRUE)

path_to_sig_file 	<- as.character(args[1])
path_to_file_bsmooth 	<- as.character(args[2])
path_to_file_dss 	<- as.character(args[3])
path_to_file_methylkit 	<- as.character(args[4])
path_to_file_metilene 	<- as.character(args[5])
path_to_file_rnbeads 	<- as.character(args[6])

sig 		<- read.table(path_to_sig_file)
bsmooth 	<- read.table(path_to_file_bsmooth)
dss 		<- read.table(path_to_file_dss)
methylkit 	<- read.table(path_to_file_methylkit)
metilene 	<- read.table(path_to_file_metilene)
rnbeads 	<- read.table(path_to_file_rnbeads)


sig$V10 <- NULL
colnames(sig) <- c("MR_chr", "MR_start", "MR_end", "MR_BS", "MR_DS", "MR_MH", "MR_MT", "MR_RB", "MR_IDR")
colnames(bsmooth)   <- c("bsmooth_chr", "bsmooth_start", "bsmooth_end", "bsmooth_cpgs", "bsmooth_methA", "bsmooth_mathB", "bsmooth_methDiff", "bsmooth_QM")
colnames(dss)       <- c("dss_chr", "dss_start", "dss_end", "dss_cpgs", "dss_methA", "dss_mathB", "dss_methDiff", "dss_QM")
colnames(methylkit) <- c("methylkit_chr", "methylkit_start", "methylkit_end", "methylkit_cpgs", "methylkit_methA", "methylkit_mathB", "methylkit_methDiff", "methylkit_QM")
colnames(metilene)  <- c("metilene_chr", "metilene_start", "metilene_end", "metilene_cpgs", "metilene_methA", "metilene_mathB", "metilene_methDiff", "metilene_QM")
colnames(rnbeads)   <- c("rnbeads_chr", "rnbeads_start", "rnbeads_end", "rnbeads_cpgs", "rnbeads_methA", "rnbeads_mathB", "rnbeads_methDiff", "rnbeads_QM")

sig$bsmooth   <- NA
sig$dss       <- NA
sig$methylkit <- NA
sig$metilene  <- NA
sig$rnbeads   <- NA


for(i in 1:nrow(sig)){
	if (i %% 1000 == 0){
		print(i)
	}
	chr   <- sig[i, ]$MR_chr
	start <- sig[i, ]$MR_start
	end   <- sig[i, ]$MR_end
 
       	if (sig[i, ]$MR_BS > 0){
		match_bsmooth <- bsmooth[bsmooth$bsmooth_chr == chr & bsmooth$bsmooth_start >= start & bsmooth$bsmooth_end <= end, ]
		methDiff <- NA 
		# Only 1 match
		if (nrow(match_bsmooth) == 1){
			methDiff <- match_bsmooth$bsmooth_methDiff
		# Multiple Matches --> aggregate 
		} else if (nrow(match_bsmooth) > 1){
			methDiff <- mean(match_bsmooth$bsmooth_methDiff)
		} 
		sig[i, ]$bsmooth <- methDiff
	}
	if (sig[i, ]$MR_DS > 0){
		match_dss <- dss[dss$dss_chr == chr & dss$dss_start >= start & dss$dss_end <= end, ]
		methDiff <- NA 
		# Only 1 match
		if (nrow(match_dss) == 1){
			methDiff <- match_dss$dss_methDiff
		# Multiple Matches --> aggregate 
		} else if (nrow(match_dss) > 1){
			methDiff <- mean(match_dss$dss_methDiff)
		}
		sig[i, ]$dss <- methDiff
	}
	if (sig[i, ]$MR_MH > 0){
		match_methylkit <- methylkit[methylkit$methylkit_chr == chr & methylkit$methylkit_start >= start & methylkit$methylkit_end <= end, ]
		methDiff <- NA 
		# Only 1 match
		if (nrow(match_methylkit) == 1){
			methDiff <- match_methylkit$methylkit_methDiff
		# Multiple Matches --> aggregate 
		} else if (nrow(match_methylkit) > 1){
			methDiff <- mean(match_methylkit$methylkit_methDiff)
		}
		sig[i, ]$methylkit <- methDiff
 	}
	if (sig[i, ]$MR_MT > 0){
		match_metilene <- metilene[metilene$metilene_chr == chr & metilene$metilene_start >= start & metilene$metilene_end <= end, ]
		methDiff <- NA 
		# Only 1 match
		if (nrow(match_metilene) == 1){
			methDiff <- match_metilene$metilene_methDiff
		# Multiple Matches --> aggregate 
		} else if (nrow(match_metilene) > 1){
			methDiff <- mean(match_metilene$metilene_methDiff)
		}
		sig[i, ]$metilene <- methDiff
	}
	if (sig[i, ]$MR_RB > 0){	      	
		match_rnbeads <- rnbeads[rnbeads$rnbeads_chr == chr & rnbeads$rnbeads_start >= start & rnbeads$rnbeads_end <= end, ]
		methDiff <- NA 
		# Only 1 match
		if (nrow(match_rnbeads) == 1){
			methDiff <- match_rnbeads$rnbeads_methDiff
		# Multiple Matches --> aggregate 
		} else if (nrow(match_rnbeads) > 1){
			methDiff <- mean(match_rnbeads$rnbeads_methDiff)
		}
		sig[i, ]$rnbeads <- methDiff
	}
}

sig$methylkit <- sig$methylkit/100

# Calculate mean from all available MethDiffs
sig$MeanDifference <- rowMeans(subset(sig, select = c(bsmooth, dss, methylkit, metilene, rnbeads)), na.rm = T)

# Write to file
write.table(sig, file = "sig.tsv", sep = "\t", col.names = NA, quote = FALSE)

# Split data set
sig_11111 <- sig[sig$MR_BS > 0 & sig$MR_DS > 0 & sig$MR_MH > 0 & sig$MR_MT > 0 & sig$MR_RB > 0, ]
sig_other <- sig[sig$MR_BS == 0 | sig$MR_DS == 0 | sig$MR_MH == 0 | sig$MR_MT == 0 | sig$MR_RB == 0, ]

# Perform T test
t.test(sig_11111$MeanDifference, sig_other$MeanDifference)

# Plot 
library(ggplot2)
sig_11111$Category <- "Found by All Tools"
sig_other$Category <- "Other"
all <- subset(sig_11111, select=c("MeanDifference", "Category"))
other <- subset(sig_other, select=c("MeanDifference", "Category"))
comb <- rbind(all, other)

p <- ggplot(comb, aes(x=Cat, y=MeanDifference, fill=Category)) + geom_boxplot()
p <- p+scale_color_manual(values=c("#4C76A3", "#E15759"))
ggsave('boxplot.pdf', plot=p, device="pdf")

