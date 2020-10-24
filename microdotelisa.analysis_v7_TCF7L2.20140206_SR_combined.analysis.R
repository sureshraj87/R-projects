## Aim: Analysis a TCF7L2 Oligos with Microdot.ELISA experiment.
## Author: Suresh Raju 
##Date: 22 July 2013 

## Idea:
## Start analysis seperately for the two individual measurements (June and july respectively)
## And check whether there is consideable deviation between two data measurements or not. 
## Normalize them separately with background corrected intensities
## Combine the June and July lot in the early final step

## Rationale:
## - Read in all data
## - Adjust oligo and column names, and split the data to sybr and protein intensities
## - Remove outliers (based on oligo-wise Z-score (=SD) cutoff (typically 3)
## - Harmonize NAs across corresponding protein and sybr intensities
## - Perform background correction (prot / relative sybr)
## - Normalize the intensity data (between arrays, using Cyclic loess method).
## - Convert from 'oligo x array' format to long table format. 
## - Subtract Scrambled oligo medians (or other quantile) from all prot intensities
## - Calculate final protein intensities as relative to the reference oligo


rm(list = ls())
date()
gc()


## Folder and library settings  

library(limma)
library(gplots)
library(Biostrings)
library(grid)
library(seqLogo)


out.data.marker <- "20140206"
infolder <- "/Users/sureshr/Project_files_scripts/TCF7L2_microdot.elisa/new_data_august_2013"

in.oligos <- "/Users/sureshr/Project_files_scripts/TCF7L2_microdot.elisa/microdotOligos_TCF7L2_all.long.txt"

setwd(infolder)
list.files()

infiles <- grep("E.txt$", list.files(), value = T)  ## possible to give as many as files you can.
infiles												## change the input as you wish.

int.names <- c("int.pos.E.sybr", "int.neg.E.sybr", "int.pos.E.tcf", "int.neg.E.tcf")
int.names.prot <- c("int.pos.E.tcf", "int.neg.E.tcf")
int.names.sybr <- c("int.pos.E.sybr", "int.neg.E.sybr")

reproduce.plots.part2 <- TRUE ## TRUE or FALSE
reproduce.plots.part3 <- TRUE ## TRUE or FALSE
reproduce.plots.part4 <- TRUE ## TRUE or FALSE
reproduce.plots.part5 <- TRUE ## TRUE or FALSE
resave.results.data <- TRUE ## TRUE or FALSE

z.cut.offs <- c(-2,3) ## abs cutoff for removing outliers based on z-scored data
max.na.per.array.pct <- 20 ## does not affect plotting, only the executed cutoff
scramble.quant.cutoff <- 0.9 ## 0.5 == median
ref.quant.cutoff <- 0.9


cols.oligo.groups <- c("black", "gray50", "brown", "green", "orange", "red") 
## in the order of ctrl, other, scrambled, ref, singlevariant and multivariant sets

## Define Scrambled and other oligo names
n.oligo.prefix.ref <- "7L2_2"
n.oligo.prefix.scr <- "7L2_scr"
n.oligo.prefix.single <- n.oligo.prefix.ref
n.oligo.prefix.mv <- "mv2" ## Only one accepted

n.oligos.other.in <- c("")
n.oligos.ctrl <- c("Control", "Mock")

n.scramble <- c(255, 371, 526, 598, 915)
n.ref <- "Ref2"
n.oligo.single.nb <- c(1:20)
n.oligo.mv.nb <- c(1:66)

n.oligo.ref.sep <- "."
n.oligo.scr.sep <- n.oligo.ref.sep
n.oligo.single.sep <- n.oligo.ref.sep
n.oligo.other.sep <- n.oligo.ref.sep
n.oligo.mv.sep <- "_"

n.oligo.ref <- paste(n.oligo.prefix.ref, n.ref, sep=".")
n.oligo.scr <- paste(n.oligo.prefix.scr, n.scramble, sep=".")
n.oligo.single.all.pos <- paste(paste(n.oligo.prefix.single, rep(n.oligo.single.nb, each=4), sep=n.oligo.single.sep), c("A", "C", "G", "T"), sep="")
n.oligo.mv.all.pos <- paste(n.oligo.prefix.mv, n.oligo.mv.nb, sep=n.oligo.mv.sep)
n.oligos.other <- gsub(" ", n.oligo.other.sep, n.oligos.other.in)

## To analyse the 2 individual measurement, make it very clear that how your data to be anaylsed. 
## In our case, Experients had been done on June and july, 2013 months.

june.measurement <- c("Lot.050613A", "Lot.050613B")
june.sybr.data <- paste(june.measurement, 488, sep = "_")
june.prtn.data <- paste(june.measurement, 543, sep = "_")

july.measurement <- c("Lot.190713A", "Lot.190713D")
july.sybr.data <- paste(july.measurement, 488, sep = "_")
july.prtn.data <- paste(july.measurement, 543, sep = "_")

# len.jun.aug.ms <- length(june.measurement) * length(june.measurement) / 2

## Construct a table for Oligo names, groups and plotting colors

oligos.ordered.all.pos <- data.frame(Oligo=c(n.oligos.ctrl, n.oligos.other, n.oligo.scr, n.oligo.ref, n.oligo.single.all.pos, n.oligo.mv.all.pos), stringsAsFactors=F)
oligos.ordered.all.pos[["gr.order"]] <- c(rep(1, length(n.oligos.ctrl)), rep(2, length(n.oligos.other)), rep(3, length(n.oligo.scr)), rep(4, length(n.oligo.ref)), rep(5, length(n.oligo.single.all.pos)), rep(6, length(n.oligo.mv.all.pos)))
oligos.ordered.all.pos[["gr.name"]] <- c(rep("Ctrl", length(n.oligos.ctrl)), rep("Other", length(n.oligos.other)), rep("Scrambled", length(n.oligo.scr)), rep("Ref", length(n.oligo.ref)), rep("Single", length(n.oligo.single.all.pos)), rep("Multi", length(n.oligo.mv.all.pos)))
oligos.ordered.all.pos[["o.order"]] <- c(1:length(n.oligos.ctrl), 1:length(n.oligos.other), 1:length(n.oligo.scr), 1:length(n.oligo.ref), 1:length(n.oligo.single.all.pos), 1:length(n.oligo.mv.all.pos))
oligos.ordered.all.pos[["o.col"]] <- NA
oligos.ordered.all.pos[["o.col"]][oligos.ordered.all.pos[["Oligo"]] %in% n.oligos.ctrl] <- cols.oligo.groups[1]
oligos.ordered.all.pos[["o.col"]][oligos.ordered.all.pos[["Oligo"]] %in% n.oligos.other] <- cols.oligo.groups[2]
oligos.ordered.all.pos[["o.col"]][oligos.ordered.all.pos[["Oligo"]] %in% n.oligo.scr] <- cols.oligo.groups[3]
oligos.ordered.all.pos[["o.col"]][oligos.ordered.all.pos[["Oligo"]] %in% n.oligo.ref] <- cols.oligo.groups[4]
oligos.ordered.all.pos[["o.col"]][oligos.ordered.all.pos[["Oligo"]] %in% n.oligo.single.all.pos] <- cols.oligo.groups[5]
oligos.ordered.all.pos[["o.col"]][oligos.ordered.all.pos[["Oligo"]] %in% n.oligo.mv.all.pos] <- cols.oligo.groups[6]
oligos.ordered.all.pos <- oligos.ordered.all.pos[oligos.ordered.all.pos[["Oligo"]] != "",]
rownames(oligos.ordered.all.pos) <- 1:nrow(oligos.ordered.all.pos)
oligos.ordered.all.pos


## Program has 5 parts.For compatibility Program is made into 5 section
## Part 1 does reading the intensity data and rename the oligos  
## Part 2 does outlier detection and removal of them. (also creates comparative maps for each array in the intensity data)
## Part 3 does removing all NAs by comparing both protein and SYBR intensity data. (More deeper view of outlier removal)
## Part 4.1 does the Normalization for the de-outliered intensity files and create informative maps 
## Part 
## Part 5 does make the final data for post processing

###########################################################
######### Part - 1 starts... ##########
###########################################################

measured.comb.int <- TRUE  ## TRUE or FALSE


## change the plot folder corresponding to the measurement folder


plot.folder <- paste(infolder, "/plots/new_measurement/", sep = "")
ouput.folder <- paste(infolder, "/output/new_measurement/", sep = "")

imp.column <- c("Oligo", "avg", "SD", "CV", "orig.order", "class", "class.color")

cat("\n", "Intensity files are reading and doing calc...", "\n")

## Reading the intensity data 
setwd(infolder)

int.ini <- list()
june.int <- list()
july.int <- list()
par(mfrow = c(1,2))

i = 1
for (i in 1:length(infiles)) {
	int.ini[[i]] <- read.table(infiles[i], header = T, sep = "\t") 
	int.ar <- grep("^Lot.", names(int.ini[[i]]), value = T)
	names(int.ini)[i] <- int.names[i]
	print(names(int.ini[[i]]))
	
	int.ini[[i]] <- cbind(int.ini[[i]][c("Spot.no.", "Oligo")], int.ini[[i]][names(int.ini[[i]]) %in% int.ar])
	print(head(int.ini[[i]]))
	print(nrow(int.ini[[i]]))
	print(names(int.ini[[i]]))
	names(int.ini[[i]])[names(int.ini[[i]]) %in% c("oligo", "Oligo")] <- "Oligo"
	int.ini[[i]][["Oligo"]] <- as.character(int.ini[[i]][["Oligo"]])

	## Oligo name changes in a more straightforward way
	int.ini[[i]][["Oligo"]][int.ini[[i]][["Oligo"]] == n.ref] <- n.oligo.ref
	int.ini[[i]][["Oligo"]][int.ini[[i]][["Oligo"]] %in% n.scramble] <- paste(n.oligo.prefix.scr, int.ini[[i]][["Oligo"]][int.ini[[i]][["Oligo"]] %in% n.scramble], sep=".")
	int.ini[[i]][["Oligo"]][int.ini[[i]][["Oligo"]] %in% grep("^[0-9]", grep("[ACGT]$", int.ini[[i]][["Oligo"]], value=T), value=T)] <- paste(n.oligo.prefix.single, grep("[ACGT]$", int.ini[[i]][["Oligo"]], value=T), sep=".")
	int.ini[[i]][["Oligo"]][int.ini[[i]][["Oligo"]] %in% grep("^[0-9]*$", int.ini[[i]][["Oligo"]], value=T)] <- paste(n.oligo.prefix.mv, int.ini[[i]][["Oligo"]][int.ini[[i]][["Oligo"]] %in% grep("^[0-9]*$", int.ini[[i]][["Oligo"]], value=T)], sep="_")
	int.ini[[i]][["Oligo"]][int.ini[[i]][["Oligo"]] %in% n.oligos.other.in] <- gsub(" ", n.oligo.other.sep, int.ini[[i]][["Oligo"]][int.ini[[i]][["Oligo"]] %in% n.oligos.other.in])
	print(int.ini[[i]][["Oligo"]])
		
		
	## Recalculate the stats. Note that apparently the original SD in the data txt file was
	## calculated using the STDEV.P Excel function (or similar) that assumes that the data is
	## for the entire population (here: of oligo intensities). True, this is the entire
	## MEASURED population, but still the N is fairly small and there may also be other
	## arrays to be measured, either in fact or at the very least in theory, so it is rather
	## a sample of the entire population. The latter is what the R function 'sd' calculates,
	## and in any case seems to be the correct one.
	
	int.ini[[i]][, int.ar][int.ini[[i]][, int.ar] == 0] <- NA
	int.ini[[i]][["avg"]] <- apply(int.ini[[i]][, int.ar], 1, mean, na.rm = T)
	int.ini[[i]][["SD"]] <- apply(int.ini[[i]][, int.ar], 1, sd, na.rm = T)
	int.ini[[i]][["CV"]] <- 100 * int.ini[[i]][["SD"]]/int.ini[[i]][["avg"]]
	int.ini[[i]][["orig.order"]] <- 1:nrow(int.ini[[i]])

	int.ini[[i]][["class"]] <- NA
	int.ini[[i]][["class"]][grep("^7L2_2", int.ini[[i]][["Oligo"]])] <- "Single"
	int.ini[[i]][["class"]][grep("Ref2$", int.ini[[i]][["Oligo"]])] <- "Ref2"
	int.ini[[i]][["class"]][grep("^mv2", int.ini[[i]][["Oligo"]])] <- "Multi"
	int.ini[[i]][["class"]][grep("_scr\\.", int.ini[[i]][["Oligo"]])] <- "Scramble"
	int.ini[[i]][["class"]][int.ini[[i]][["Oligo"]] %in% c("Control", "Mock")] <- "No oligo"
	int.ini[[i]][["class"]] <- as.factor(int.ini[[i]][["class"]])
	table(int.ini[[i]][["class"]])
	
	int.ini[[i]][["class.color"]] <- NA
			int.ini[[i]][["class.color"]][int.ini[[i]][["Oligo"]] %in% n.oligo.ref] <- cols.oligo.groups[4]
			int.ini[[i]][["class.color"]][int.ini[[i]][["Oligo"]] %in% n.oligo.scr] <- cols.oligo.groups[3]
			int.ini[[i]][["class.color"]][int.ini[[i]][["Oligo"]] %in% n.oligo.single.all.pos] <- cols.oligo.groups[5]
			int.ini[[i]][["class.color"]][int.ini[[i]][["Oligo"]] %in% n.oligo.mv.all.pos] <- cols.oligo.groups[6]
			int.ini[[i]][["class.color"]][int.ini[[i]][["Oligo"]] %in% n.oligos.ctrl] <- cols.oligo.groups[1]
			int.ini[[i]][["class.color"]][is.na(int.ini[[i]][["class"]])] <- cols.oligo.groups[2]
	table(int.ini[[i]][["class.color"]])

	plot(int.ini[[i]][["avg"]], las = 2, type = "o", col = "red", ylab = "Mean Intensity", main = paste("Mean Intensity of", int.names[i], "from direct measurement", sep = "_"))
	
	## Take the individual measurement (JUNE or july) or all arrays for further analysis
	
	june.ar <- int.ar[gsub("_.*", "", int.ar) %in% june.measurement]
	june.int[[i]] <- int.ini[[i]][names(int.ini[[i]]) %in% c(imp.column, june.ar)]
	# names(june.int[[i]])[names(june.int[[i]]) %in% imp.column] <- paste(imp.column, "june", sep = ".")
	names(june.int)[i] <- paste(int.names[i], "june", sep = "_")
	print(head(june.int[[i]]))
	
	
	july.ar <- int.ar[gsub("_.*", "", int.ar) %in% july.measurement]	
	july.int[[i]] <- int.ini[[i]][names(int.ini[[i]]) %in% c(imp.column, july.ar)]
	names(july.int)[i] <- paste(int.names[i], "july", sep = "_")
	print(head(july.int[[i]]))	
	
}

par(mfrow = c(1,1))

print(lapply(june.int, head))
print(lapply(june.int, nrow))
print(lapply(july.int, head))
print(lapply(july.int, nrow))

###########################################################
######## Part - 2 starts... #########
###########################################################

de.out.int.ini <- list()
de.out.june.int <- list()
de.out.july.int <- list()

all.int.data <- list()
all.int.data[1:4] <- june.int
all.int.data[5:8] <- july.int
names(all.int.data) <- names(c(june.int, july.int))

int.names <- names(all.int.data)

lapply(all.int.data, head)
lapply(all.int.data, nrow)


## important columns for all data frame

cat("\n", "Outlier detection starts...", "\n")

setwd(plot.folder)
i = 1
for (i in 1:length(all.int.data)) {
	
	int.ar <- grep("^Lot.", names(all.int.data[[i]]), value = T)
	tmp.int.arr <- all.int.data[[i]][, int.ar]
	tmp.int.arr[, int.ar][tmp.int.arr[, int.ar] == 0] <- NA
	tmp.int.arr <- t(tmp.int.arr[, int.ar])
	colnames(tmp.int.arr) <- all.int.data[[i]][["Oligo"]]
	
	## creating own color panel
	colors <- c(seq(-6, -2, length = 100), seq(-2, 0.5, length = 100), seq(0.5, 6, length = 200))
	# colors <- seq(-5,5, 0.1)
	mycol <- colorRampPalette(c("blue", "yellow", "brown"))(n = 299)

	
	## Z score calculation, Per oligo, convert all values to Z-scores using mean and SD.
	## Convert for combined, June and july measurements
	
	ini.ol <- all.int.data[[i]]
	n.oligo.u <- unique(ini.ol[["Oligo"]])
	k = 1
	for (k in 1:length(n.oligo.u)) {
		tmp.data <- as.matrix(ini.ol[which(ini.ol[["Oligo"]] == n.oligo.u[k]), ][, int.ar])
		print(k)
		print(mean(as.vector(tmp.data), na.rm = T))
		print(median(as.vector(tmp.data), na.rm = T))
		ini.ol[which(ini.ol[["Oligo"]] == n.oligo.u[k]), ][, int.ar] <- (tmp.data - mean(as.vector(tmp.data), na.rm = T))/sd(as.vector(tmp.data), na.rm = T)
		# print(as.matrix(ini.ol[which(ini.ol[["Oligo"]] == n.oligo.u[k]),][, int.ar]))
		cat("\n\n")
		
		}
	print(head(ini.ol))
	print(nrow(ini.ol))
	
	## Outlier detection and removal 
	## Identify oligos that deviate by more than 3 SD.			
	
	tmp.out <- ini.ol
	a <- ceiling(max(tmp.out[, int.ar], na.rm = T))
	b <- floor(min(tmp.out[, int.ar], na.rm = T))
	tmp.all <- t(tmp.out[, int.ar])
	colnames(tmp.all) <- tmp.out[["Oligo"]]

	## Outlier removal
	
	tmp.out[, int.ar][tmp.out[, int.ar] > z.cut.offs[2]] <- 0
	tmp.out[, int.ar][tmp.out[, int.ar] < z.cut.offs[1]] <- 0
	print(head(tmp.out))
	de.out <- t(tmp.out[, int.ar])
	colnames(de.out) <- ini.ol[["Oligo"]]
	
	##################
	# Removing outliers in the original intensity data:
	# Decided to remove the outliers wihch are existing more than the sd value 3
	
	de.out.int <- NULL
	de.out.int <- all.int.data[[i]]
	print(nrow(de.out.int))
	print(nrow(tmp.out))

	de.out.int <- data.matrix(de.out.int[, int.ar])
	de.out.int[is.na(de.out.int)] <- 0
	de.out.int[which(de.out.int == de.out.int) %in% which(0 == tmp.out[, int.ar])] <- 0
	de.out.int <- as.data.frame(de.out.int, stringsAsFactors = F)
	de.out.int[, int.ar][de.out.int[, int.ar] == 0] <- NA
	de.out.int <- cbind(all.int.data[[i]][imp.column], de.out.int)
	str(de.out.int)

	## storing the de outliered data into a list
	de.out.int.ini[[i]] <- de.out.int
	names(de.out.int.ini)[i] <- names(all.int.data)[i]
	de.out.arr <- t(de.out.int[, int.ar])
	colnames(de.out.arr) <- all.int.data[[i]][["Oligo"]]

	
	##### outlier removal in the calculated z score data#####
	heat.test <- NULL
	heat.test <- ini.ol
	print(nrow(heat.test))

	heat.test <- data.matrix(heat.test[, int.ar])
	heat.test[is.na(heat.test)] <- 0
	heat.test[which(heat.test == heat.test) %in% which(0 == tmp.out[, int.ar])] <- 0
	heat.test[heat.test == 0] <- NA
	heat.test <- data.frame(heat.test, stringsAsFactors = F)
	heat.test <- cbind(ini.ol[imp.column], heat.test)

	
	tmp.heat <- t(heat.test[, int.ar])
	# tmp.heat <- type.convert(tmp.heat, na.strings="NA", as.is=F, dec = ".")
	colnames(tmp.heat) <- heat.test[["Oligo"]]

		########
	
	# Replicate counting, without the loop
	df.rc <- data.frame(Oligo = unique(ini.ol[["Oligo"]]), oligo.repl = NA)
	df.rc[["oligo.repl"]] <- apply(df.rc, 1, function(x, y = ini.ol) length(which(y[["Oligo"]] %in% x[["Oligo"]])))
	ini.ol <- merge(ini.ol, df.rc, by = "Oligo")
	ini.ol <- ini.ol[order(ini.ol[["orig.order"]]), ]
	ini.ol[["Oligo_r"]] <- paste(ini.ol[["Oligo"]], ini.ol[["oligo.repl"]], sep = "_")
	uni.oid <- unique(ini.ol[["Oligo"]])
	
	## All the neccessary plots
	
	if (reproduce.plots.part2) {
		png(filename = paste(paste("Fig.1_Density.plot_original_data", int.names[i], sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		plot(density(tmp.int.arr, na.rm = T), xlab = "intensity values", main = paste("Density.plot_original_data", int.names[i], sep = "_"))
		dev.off()
	
		png(filename = paste(paste("Fig.2_Array.comparison_original_data", int.names[i], sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		heatmap.2(tmp.int.arr, trace = "none", dendrogram = "none", col = mycol, margin = c(5,10), Colv = F, Rowv = F, scale = "column", na.rm = T, keysize = 1, main = paste("Array.comparison_original_data", int.names[i], sep = "_"))
		dev.off()
	
		png(filename = paste(paste("Fig.3_Oligo.data_with_all.z.score", int.names[i], sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		boxplot(tmp.all, ylim = c(-a, a), las = 2, col = ini.ol[["class.color"]], ylab = "Z-score", main = paste(paste("Oligo.data_with_all.calculated.z.score", sep = "-"), int.names[i], sep = "_"), na.rm = T)
		abline(h = c(-3:-1, -1.5, 1:3, 1.5), col = c(rep("blue", 4), rep("red", 4)))
		dev.off()

		png(filename = paste(paste("Fig.4_Outlier.removed", int.names[i], "Oligo.data_with_z-score", "-2 to 3", sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		boxplot(de.out, ylim = c(-4, 4), las = 2, col = ini.ol[["class.color"]], ylab = "Z-score", main = paste(paste("Outlier.removed_oligo.data_z.score.of", "-2 to 3", sep = " "), int.names[i], sep = "_"), na.rm = T)
		abline(h = c(-3:-1, -1.5, 1:3, 1.5), col = c(rep("blue", 4), rep("red", 4)))
		abline(h = c(-2, 3), lwd = 2.5, col = c("blue", "red"))

		dev.off()

		png(filename = paste(paste("Fig.5_Density.plot_after_outliers_removal", int.names[i],	sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		plot(density(de.out.arr, na.rm = T), xlab = "intensity values", main = paste("Density.plot_after_outliers_removal", int.names[i], sep = "_"))
		dev.off()

		png(filename = paste(paste("Fig.6_Array_comparison_after_outliers.removal_original_data", int.names[i], sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		heatmap.2(de.out.arr, trace = "none", dendrogram = "none", col = mycol, margin = c(5,10), Colv = F, Rowv = F, scale = "column", na.rm = T, keysize = 1, main = paste("Array_comparison_after_outliers.removal_original_data", int.names[i], sep = "_"))
		dev.off()

		png(filename = paste(paste("Fig.7_Array_comparison_after_outliers.removal_calculated.zscore_data", int.names[i], sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		heatmap.2(tmp.heat, trace = "none", dendrogram = "none", col = mycol, margin = c(5,10), Colv = F, Rowv = F, scale = "none", na.rm = T, keysize = 1, main = paste("Array comparison for calculated z score data of", int.names[i], sep = "_"))
		dev.off()
	}
}

print(lapply(de.out.int.ini, head))
print(lapply(de.out.int.ini, nrow))
print(lapply(de.out.int.ini, names))

###################


###########################################################
######### Part - 3 starts... ##########
###########################################################


## More narrow and deeper view of outliers detection in both SYBR and protein data, 
## Idea:compare the outliers in original intensity data of both sybr and
## protein,then it can be removed vice versa and calculate the % of "NA"s.
## Generate relative intensity values for each oligos	

cat("\n", "Removal of missing values starts...", "\n")
int.names <- names(de.out.int.ini)
prtn.int <- list()
sybr.int <- list()
tmp.prtn <- list()
tmp.sybr <- list()

sybr.name <- grep("sybr", int.names, value = T)
prtn.name <- grep("tcf", int.names, value = T)

sybr.int[1:length(sybr.name)] <- de.out.int.ini[sybr.name]
names(sybr.int) <- sybr.name
prtn.int[1:length(prtn.name)] <- de.out.int.ini[prtn.name]
names(prtn.int) <- prtn.name

print(lapply(prtn.int, head))
print(lapply(prtn.int, nrow))
print(lapply(sybr.int, head))
print(lapply(sybr.int, nrow))

setwd(plot.folder)
h = 1
for (h in 1:length(prtn.int)) {
	prtn.arr <- grep("^Lot", names(prtn.int[[h]]), value = T)
	sybr.arr <- grep("^Lot", names(sybr.int[[h]]), value = T)
	table(is.na(prtn.int[[h]][, prtn.arr]))
	table(is.na(sybr.int[[h]][, sybr.arr]))
	prtn.int[[h]][, prtn.arr][is.na(prtn.int[[h]][, prtn.arr])] <- 0
	sybr.int[[h]][, sybr.arr][is.na(sybr.int[[h]][, sybr.arr])] <- 0
	prtn.data <- data.matrix(prtn.int[[h]][, prtn.arr])
	sybr.data <- data.matrix(sybr.int[[h]][, sybr.arr])

	## Removing missing values in both SYBR and Protein data by comparing each other
	
	prtn.data[sybr.data == 0] <- 0
	sybr.data[prtn.data == 0] <- 0

	## Then, convert the "new" missing values (ie zeros) back to NAs
	prtn.data[prtn.data == 0] <- NA
	sybr.data[sybr.data == 0] <- NA
	prtn.data <- data.frame(prtn.data, stringsAsFactors = F)
	prtn.data <- cbind(prtn.int[[h]][c("Oligo", "avg", "SD", "CV", "orig.order", "class", "class.color")], prtn.data)
	sybr.data <- data.frame(sybr.data, stringsAsFactors = F)
	sybr.data <- cbind(sybr.int[[h]][c("Oligo", "avg", "SD", "CV", "orig.order", "class", "class.color")], sybr.data)
	# print(head(prtn.data))
	# print(head(sybr.data))

	## The tabling done only on the measured intensities to enable verification
	## that the NA harmonization has worked as excpected.
	
	table(is.na(prtn.data[grep("^Lot", names(prtn.data))]))
	table(is.na(sybr.data[grep("^Lot", names(sybr.data))]))
	all(which(is.na(as.matrix(prtn.data[grep("^Lot", names(prtn.data))]))) == which(is.na(as.matrix(sybr.data[grep("^Lot", names(sybr.data))])))) 
	
	## Should be TRUE
	## Removal of the technical negative control spots that did not have any 
	## oligo on the arrays.
	prtn.data <- prtn.data[prtn.data[["class"]] %in% grep("^No", prtn.data[["class"]], invert = T, value = T), ]
	sybr.data <- sybr.data[sybr.data[["class"]] %in% grep("^No", sybr.data[["class"]], invert = T, value = T), ]


	## RESTRUCTURED THE NA PERCENTAGE CALCULATIONS AND CORRECTED A NUMBER
	## OF MISTAKES.
	## The task is to remove ARRAYs with lots of NAs, not OLIGOs.
	## The script proceeds as follows:
	## - the NA percentages are calculated for OLIGOs and corrected for those
	##   oligos that have data on more than 1 row
	## - plot the above; suffice to plot only either prtn or sybr, since the
	##   the NAs have already been harmonized
	## - the NA percentages are calculated for ARRAYs, and plotted
	## - identify ARRAYS with more than X% (X=20) NAs, and remove
	## - the NA percentages are re-calculated for OLIGOs as above, and plotted
	
	## Note that the all is done using the prtn.data and sybr.data objects
	## and only at the very last step assigned to tmp.prtn[[h]] etc.
	
	## NA percentage calculation per OLIGO (initial, before ARRAYS with
	## lots of NAs have been removed, to enable plotting already at this stage;
	## similar plots are generated also after the removal of ARRAYs with lots of NAs).

	prtn.data[["percent_NA"]] <- apply(prtn.data[, prtn.arr], 1, function(x) sum(is.na(x))/length(x)) * 100
	sybr.data[["percent_NA"]] <- apply(sybr.data[, sybr.arr], 1, function(x) sum(is.na(x))/length(x)) * 100
	
	## The original code removed the OLIGOS that have > 20% of missing values, when
	## the task was to remove ARRAYS that have lots of missing values.
	## There are oligos that have lots of missing values in either the June or July
	## array lots, but mostly valid values in the other lot, so removing those is an
	## error. Moreover, the calculations for the 'cor.percent.NA' were wrong,
	## effectively dividing all NA percentages with the same value when the
	## point was to "correct" the row-wise percentages only for such oligos that
	## needed the correction, ie those that have data on more than one row.
	## Furthermore, a "bad" array is worse than a "bad" oligo, so the "bad" arrays
	## should be removed first, before inspecting whether to remove any oligos.

	n <- by(prtn.data[["percent_NA"]], prtn.data[["Oligo"]], mean)
	prtn.data[["cor.percent.NA"]] <- as.vector(n[match(prtn.data[["Oligo"]], rownames(n))])
	n <- by(sybr.data[["percent_NA"]], sybr.data[["Oligo"]], mean)
	sybr.data[["cor.percent.NA"]] <- as.vector(n[match(sybr.data[["Oligo"]], rownames(n))])
	
	## Remove oligos which has more than 90% NAs. the oligos which has more NAs are introducing 
	## outliers values, mean and SD values cannot be calculated properly. 
	
	rem.oli.per <- 90
	if(max(prtn.data[["percent_NA"]], na.rm = T) >= rem.oli.per) {
		prtn.data <- subset(prtn.data, !percent_NA >= rem.oli.per)
	}
	if(max(sybr.data[["percent_NA"]], na.rm = T) > rem.oli.per) {
		sybr.data <- subset(sybr.data, !percent_NA > rem.oli.per)
	}
		
	if (reproduce.plots.part3) {
		
		## Plotting Oligo's percentage	
		png(filename = paste(paste("Fig.8.Oligo_NAs_percentage_initial", sybr.name[h], sep = "_"), "png", sep = "."), width = 12, height = 12, res = 200, unit = "in")
		par(mar = c(7, 4, 4, 2))
		par(mfrow = c(2, 1))
		barplot(sybr.data[["percent_NA"]], names = sybr.data[["Oligo"]], ylim = c(0, 100), ylab = "percentage of NAs", col = as.vector(sybr.data[["class.color"]]), main = paste("Oligo_NAs_percentage", sybr.name[h], sep = "_"), las = 2)
		abline(h = 20, lty = 2, col = "gray")
		barplot(sybr.data[["cor.percent.NA"]], names = sybr.data[["Oligo"]], ylim = c(0, 100), ylab = "percentage of NAs", col = as.vector(sybr.data[["class.color"]]), main = paste("Oligo_NAs_percentage", sybr.name[h], sep = "_"), las = 2)
		abline(h = 20, lty = 2, col = "gray")
		dev.off()
	}


	## NA percentage calculation per ARRAY
	
	percent.NA.pr.arr <- (apply(prtn.data[, prtn.arr], 2, function(x) sum(is.na(x)))/nrow(prtn.data)) * 100
	percent.NA.sy.arr <- (apply(sybr.data[, sybr.arr], 2, function(x) sum(is.na(x)))/nrow(sybr.data)) * 100
	names(percent.NA.sy.arr) <- NULL
	names(percent.NA.pr.arr) <- NULL

	## Below, the color scheme still had 25% as one of the cutoff points instead of
	## the 20% that we agreed to use. This has been corrected.

	## Assigning all calculated percentage values in one data frame
	percent.calc <- NULL
	percent.calc[["prtn.arrays"]] <- prtn.arr
	percent.calc[["sybr.arrays"]] <- sybr.arr
	percent.calc[["prtn.arr.percent"]] <- percent.NA.pr.arr
	percent.calc[["sybr.arr.percent"]] <- percent.NA.sy.arr
	percent.calc <- as.data.frame(percent.calc)
	percent.calc[["color.scheme"]][percent.calc[["prtn.arr.percent"]] <= 10] <- "green"
	# percent.calc[["color.scheme"]][percent.calc[["prtn.arr.percent"]] > 10 & percent.calc[["prtn.arr.percent"]] <= 25] <- "yellow"
	# percent.calc[["color.scheme"]][percent.calc[["prtn.arr.percent"]] > 25 & percent.calc[["prtn.arr.percent"]] <= 50] <- "red"
percent.calc[["color.scheme"]][percent.calc[["prtn.arr.percent"]] > 10 & percent.calc[["prtn.arr.percent"]] <= 20] <- "yellow"
	percent.calc[["color.scheme"]][percent.calc[["prtn.arr.percent"]] > 20 & percent.calc[["prtn.arr.percent"]] <= 50] <- "red"
	percent.calc[["color.scheme"]][percent.calc[["prtn.arr.percent"]] > 50] <- "black"

	print(head(percent.calc))
	print(nrow(percent.calc))

	## Plots the initial NA percentages for the ARRAYS. 
	## Note that the percentage ranges assigned above had the equal sign at the
	## upper limit whereas the plot legend did not have them except for the lowest,
	## 10% limit. This is now corrected.
	if (reproduce.plots.part3) {
		png(filename = paste(paste("Fig.9.Array_NAs_percentage_initial", prtn.name[h], sep = "_"),		"png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(10, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		barplot(percent.calc[["prtn.arr.percent"]], ylim = c(0, 100), ylab = "percentage of NAs", col = percent.calc[["color.scheme"]], names = percent.calc[["prtn.arrays"]], 	main = paste("Array_NAs_percentage", prtn.name[h], sep = "_"), las = 2)
		abline(h = 20, lty = 2, col = "gray")
		legend("topright", legend = c("<= 10%", "> 10% & <= 20%", "> 20% & <= 50%", "> 50%"), fill = c("green", "yellow", "red", "black"))
		dev.off()
	}


	## Get array names with too many NAs
	arr.id.many.nas.p <- as.character(percent.calc[["prtn.arrays"]][percent.calc[["prtn.arr.percent"]] > max.na.per.array.pct])
	arr.id.many.nas.s <- as.character(percent.calc[["sybr.arrays"]][percent.calc[["sybr.arr.percent"]] > max.na.per.array.pct])


	## Remove ARRAYS with lots of NAs and update the array name vectors.
	
	if (length(arr.id.many.nas.p) > 0) {
		prtn.data <- prtn.data[, -which(names(prtn.data) %in% arr.id.many.nas.p)]
		prtn.arr <- prtn.arr[-which(prtn.arr %in% arr.id.many.nas.p)]
	}
	if (length(arr.id.many.nas.s) > 0) {
		sybr.data <- sybr.data[, -which(names(sybr.data) %in% arr.id.many.nas.s)]
		sybr.arr <- sybr.arr[-which(sybr.arr %in% arr.id.many.nas.s)]
	}


	## Recalculate the NA percentages for OLIGOs as above (overwrites the 
	## pre-existing columns)
	
	prtn.data[["percent_NA"]] <- apply(prtn.data[, prtn.arr], 1, function(x) sum(is.na(x))/length(x)) * 100
	sybr.data[["percent_NA"]] <- apply(sybr.data[, sybr.arr], 1, function(x) sum(is.na(x))/length(x)) * 100

	n <- by(prtn.data[["percent_NA"]], prtn.data[["Oligo"]], mean)
	prtn.data[["cor.percent.NA"]] <- as.vector(n[match(prtn.data[["Oligo"]], rownames(n))])
	n <- by(sybr.data[["percent_NA"]], sybr.data[["Oligo"]], mean)
	sybr.data[["cor.percent.NA"]] <- as.vector(n[match(sybr.data[["Oligo"]], rownames(n))])

	if (reproduce.plots.part3) {
		## Plotting Oligo's percentage	
		png(filename = paste(paste("Fig.10.Oligo_NAs_percentage_after_handling_arrays", sybr.name[h], sep = "_"), "png", sep = "."), width = 12, height = 12, res = 200, unit = "in")
		par(mar = c(7, 4, 4, 2))
		par(mfrow = c(2, 1))
		barplot(sybr.data[["percent_NA"]], names = sybr.data[["Oligo"]], ylim = c(0, 100), ylab = "percentage of NAs", col = as.vector(sybr.data[["class.color"]]), main = paste("Oligo_NAs_percentage", sybr.name[h], sep = "_"), las = 2)
		abline(h = 20, lty = 2, col = "gray")
		barplot(sybr.data[["cor.percent.NA"]], names = sybr.data[["Oligo"]], ylim = c(0, 100), ylab = "percentage of NAs", col = as.vector(sybr.data[["class.color"]]), main = paste("Oligo_NAs_percentage", sybr.name[h], sep = "_"), las = 2)
		abline(h = 20, lty = 2, col = "gray")
		dev.off()
	}


	########################
	
	## Only the heatmaps remain to be plotted here.
	
	if (reproduce.plots.part3) {

		sybr.out.heat <- t(sybr.data[, sybr.arr])
		colnames(sybr.out.heat) <- sybr.data[["Oligo"]]
		png(filename = paste(paste("Fig.11.Array_comparison_after.removal.of.all_NAs", sybr.name[h], sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		heatmap.2(sybr.out.heat, trace = "none", dendrogram = "none", col = mycol, margin = c(5, 10), Colv = F, Rowv = F, scale = "column", na.rm = T, keysize = 1, main = paste("Array_comparison_after.removal.of.all_NAs", sybr.name[h], sep = "_"))
		dev.off()

		prtn.out.heat <- t(prtn.data[, prtn.arr])
		colnames(prtn.out.heat) <- prtn.data[["Oligo"]]
		png(filename = paste(paste("Fig.12.Array_comparison_after.removal.of.all_NAs", prtn.name[h], sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		heatmap.2(prtn.out.heat, trace = "none", dendrogram = "none", col = mycol, margin = c(5, 10), Colv = F, Rowv = F, scale = "column", na.rm = T, keysize = 1, main = paste("Array_comparison_after.removal.of.all_NAs", prtn.name[h], sep = "_"))
		dev.off()
	}

	## Now, finally assign the objects to the tmp.prtn and tmp.sybr lists
	tmp.prtn[[h]] <- prtn.data
	tmp.sybr[[h]] <- sybr.data

	### Background correction
	

	################################
	## SYBR relative intensity data can be calculated by total mean value of SYBR intensity 
	## In order to correct the background intensity in protein data, relative intensity 
	## data will be calculated by relative SYBR intensity.

	total.prtn.mean <- mean(as.matrix(tmp.prtn[[h]][, prtn.arr]), na.rm = T)
	total.sybr.mean <- mean(as.matrix(tmp.sybr[[h]][, sybr.arr]), na.rm = T)

	tmp.sybr[[h]][, paste("rel", sybr.arr, sep = "_")] <- tmp.sybr[[h]][, sybr.arr]/total.sybr.mean

	## Using the original '^rel' to 'grep' for the column names worked correctly
	## when calculating the mean but not for SD since that included also the
	## freshly calculated means, ie also 'rel.int_mean' gets grepped using '^rel'.
	## Therefore, it is grep search was changed to '^rel_' for both mean and SD. 

	tmp.sybr[[h]][["rel.int_mean"]] <- apply(tmp.sybr[[h]][, grep("^rel_", names(tmp.sybr[[h]]),value = T), ], 1, function(x) mean(x, na.rm = T))
	tmp.sybr[[h]][["rel.int_sd"]] <- apply(tmp.sybr[[h]][, grep("^rel_", names(tmp.sybr[[h]]), value = T), ], 1, function(x) sd(x, na.rm = T))

	print(lapply(tmp.sybr, head))
	print(lapply(tmp.sybr, nrow))

	## The bg correction was done using the original, not relative SYBR values that
	## were calculated above. The point is to correct the protein intensities for
	## RELATIVE technical differences due to especially variations in the dotted 
	## oligo amounts. Using RELATIVE differences also has the benefit of preserving
	## the protein intensities at near the original range.

	s = 1
	for (s in 1:length(prtn.arr)) {
		tmp.prtn[[h]][, paste("rel", prtn.arr[s], sep = "_")] <- tmp.prtn[[h]][, prtn.arr][s]/tmp.sybr[[h]][, paste("rel", sybr.arr[s], sep = "_")]
	}

	## Same grep mistake in the calc for SD as above; similarly corrected 	
	
	tmp.prtn[[h]][["rel.int_mean"]] <- apply(tmp.prtn[[h]][, grep("^rel_", names(tmp.prtn[[h]]), value = T), ], 1, function(x) mean(x, na.rm = T))
	tmp.prtn[[h]][["rel.int_sd"]] <- apply(tmp.prtn[[h]][, grep("^rel_", names(tmp.prtn[[h]]), value = T), ], 1, function(x) sd(x, na.rm = T))

	print(lapply(tmp.prtn, head))
	print(lapply(tmp.prtn, nrow))

	##################################				
	#### heatmap for background corrected samples

	## This heatmap is also now plotted only for protein data since the sybr data
	## was not corrected for anything and is therefore the same as above.

	if (reproduce.plots.part3) {
		# heat.bc.prtn <- t(tmp.prtn[[h]][, prtn.arr])
		heat.bc.prtn <- t(tmp.prtn[[h]][, paste("rel", prtn.arr, sep = "_")])
		colnames(heat.bc.prtn) <- tmp.prtn[[h]][["Oligo"]]
		png(filename = paste(paste("Fig.13.Array_comparison_background corrected", prtn.name[h], sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		heatmap.2(heat.bc.prtn, trace = "none", dendrogram = "none", col = mycol, margin = c(5, 10), Colv = F, Rowv = F, scale = "column", na.rm = T, keysize = 1, main = paste("Array_comparison_background corrected", prtn.name[h], sep = "_"))
		dev.off()

		}
}


## NAs removed data is assigned into a new object
## Only needed for protein data, because sybr data has already surved its purpose
## at the background removal stage above.
## Normalize the June and July data separately and do the background correction 
## (subtraction of scrambled oligos mean) and relative to the reference

tmp.na.rm.int.ini <- list()
tmp.na.rm.int.ini[1:length(tmp.prtn)] <- tmp.prtn
names(tmp.na.rm.int.ini) <- prtn.name

tmp.names <- c("int.pos.E.tcf_june", "int.pos.E.tcf_july")

na.rm.int.ini <- list()
na.rm.int.ini <- tmp.na.rm.int.ini[names(tmp.na.rm.int.ini) %in% tmp.names]
names(na.rm.int.ini) <- tmp.names

lapply(na.rm.int.ini, head)
lapply(na.rm.int.ini, nrow)
lapply(na.rm.int.ini, names)


###########################################################
			   ###### Part - 4.1 Starts.. #######
###########################################################

## RESTRUCTURED AS FOLLOWS:
## - get the de-outliered and NA-handled in intensities before and after
##   background-correction. Could be expanded to cover also the data
##   at steps of the processing
## - get the final, background-corrected data
## - normalize between the arrays, using all data points (i.e. without
##   averaging over those with data on multiple rows)
## - only then take into account the Ref2 oligo

cat("\n", "Normalization starts...", "\n")
pre.fin.int.data <- list()

setwd(plot.folder)
i = 1
for (i in 1:length(na.rm.int.ini)) {
	
	int.ar <- grep("^Lot\\.", names(na.rm.int.ini[[i]]), value = T)
	rel.ar <- grep("^rel_Lot\\.", names(na.rm.int.ini[[i]]), value = T)

	## For pre-bg-corrected values, normalize and calc mean and sd
	## The column callings are now more straightforwards and
	## less prone to misordering.
	
	na.rm.int.ini[[i]][, paste("norm", int.ar, sep = "_")] <- normalizeBetweenArrays(as.matrix(na.rm.int.ini[[i]][, int.ar]), method = "cyclicloess")
	na.rm.int.ini[[i]][["norm_int.mean"]] <- apply(na.rm.int.ini[[i]][, paste("norm", int.ar, sep = "_")], 1, mean, na.rm = T)
	na.rm.int.ini[[i]][["norm_int.sd"]] <- apply(na.rm.int.ini[[i]][, paste("norm", int.ar, sep = "_")], 1, sd, na.rm = T)
	int.torefnorm <- paste("norm", int.ar, sep = "_")
	print(int.torefnorm)

	na.rm.int.ini[[i]][, paste("norm", rel.ar, sep = "_")] <- normalizeBetweenArrays(as.matrix(na.rm.int.ini[[i]][, rel.ar]), method = "cyclicloess")
	na.rm.int.ini[[i]][["norm_rel.mean"]] <- apply(na.rm.int.ini[[i]][, paste("norm", rel.ar, sep = "_")], 1, mean, na.rm = T)
	na.rm.int.ini[[i]][["norm_rel.sd"]] <- apply(na.rm.int.ini[[i]][, paste("norm", rel.ar, sep = "_")], 1, sd, na.rm = T)
	rel.torefnorm <- paste("norm", rel.ar, sep = "_")

	setwd(plot.folder)
	if (reproduce.plots.part4) {
		png(filename = paste("Fig.14.Cyclic.Loess.Method_background.corrected_for", tmp.names[i], "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(5, 4, 4, 2) + 0.1)
		par(mfrow = c(2, 2))
		
		tmp.plot <- na.rm.int.ini[[i]]
		cn.plot <- grep("^rel_Lot", names(na.rm.int.ini[[i]]), value = T)
		col.plot <- rainbow(length(cn.plot))
		plot(tmp.plot[["rel.int_mean"]], tmp.plot[["rel.int_mean"]], xlim = c(0, 1.1 * max(tmp.plot[["rel.int_mean"]], na.rm = T)), ylim = c(0, max(as.matrix(na.rm.int.ini[[i]][,cn.plot]), na.rm = T)), type = "n", xlab = "mean intensity", ylab = "mean intensity", main = "Non-normalized_background.corrected_data")
		r = 1
		for (r in 1:length(cn.plot)) {
			# cat(max(tmp.plot[[cn.plot[r]]], na.rm=T), "\n")
			points(tmp.plot[["rel.int_mean"]], tmp.plot[[cn.plot[r]]], col = col.plot[r], pch = 20, cex = 1)
			lines(loess.smooth(tmp.plot[["rel.int_mean"]], tmp.plot[[cn.plot[r]]]), col = col.plot[r], lwd = 2)
		}

		cn.plot <- grep("^norm_rel_Lot", names(na.rm.int.ini[[i]]), value = T)
		plot(tmp.plot[["norm_rel.mean"]], tmp.plot[["norm_rel.mean"]], xlim = c(0, 1.1 * max(tmp.plot[["norm_rel.mean"]], na.rm = T)), ylim = c(0, max(as.matrix(na.rm.int.ini[[i]][,cn.plot]), na.rm = T)), type = "n", xlab = "mean intensity", ylab = "normalized mean intensity", main = "Normalized_Cyclic.Loess.method_background.corrected_data")
		r = 1
		for (r in 1:length(cn.plot)) {
			# cat(max(tmp.plot[[cn.plot[r]]], na.rm=T), "\n")
			points(tmp.plot[["rel.int_mean"]], tmp.plot[[cn.plot[r]]], col = col.plot[r], pch = 20, cex = 1)
			lines(loess.smooth(tmp.plot[["rel.int_mean"]], tmp.plot[[cn.plot[r]]]), col = col.plot[r], lwd = 2)
		}


		boxplot(tmp.plot[, grep("^rel_Lot", names(na.rm.int.ini[[i]]), value = T)], main = "Non-normalized", ylab = "test intensity", las = 2)
		boxplot(tmp.plot[, grep("^norm_rel_Lot", names(na.rm.int.ini[[i]]), value = T)], main = "Normalized using the Cyclic Loess method", ylab = "normalized test intensity", las = 2)
		par(mfrow = c(1, 1))
		dev.off()
		
		}
		
	## check for outliers in the normalized intensities using zscore calculation as above (in part-2)
	## Though the outliers have been discarded in the early part of the script, 
	## we may re consider of doing outlier detection over here as well
	
	int.out <- na.rm.int.ini[[i]]
	uni.oligo.o <- unique(int.out[["Oligo"]])
	norm.int.ar <- grep("^norm_rel_Lot", names(na.rm.int.ini[[i]]), value = T)
	k = 1
	for (k in 1:length(uni.oligo.o)) {
		tmp.out.data <- as.matrix(int.out[which(int.out[["Oligo"]] == uni.oligo.o[k]), ][, norm.int.ar])
		print(k)
		print(mean(as.vector(tmp.data), na.rm = T))
		print(median(as.vector(tmp.data), na.rm = T))
		int.out[which(int.out[["Oligo"]] == uni.oligo.o[k]), ][, norm.int.ar] <- (tmp.out.data - mean(as.vector(tmp.out.data), na.rm = T))/sd(as.vector(tmp.out.data), na.rm = T)
		# print(as.matrix(ini.ol[which(ini.ol[["Oligo"]] == n.oligo.u[k]),][, int.ar]))
		cat("\n\n")
		
		}
	# int.out <- cbind(int.out [,imp.column], int.out[, norm.int.ar])
	int.out[["norm_rel.mean"]] <- apply(int.out[, norm.int.ar], 1, function(x) mean(x, na.rm = T))
	int.out[["norm_rel.sd"]] <- apply(int.out[, norm.int.ar], 1, function(x) sd(x, na.rm = T))
	print(head(int.out))
	print(nrow(int.out))
	
	cat("\n", "maximum zscore value for the normalized intensities", max(int.out[, norm.int.ar], na.rm = T), "\n")
	cat("\n", "minimum zscore value for the normalized intensities", min(int.out[, norm.int.ar], na.rm = T), "\n")
	
	min.ylim <- min(int.out[, norm.int.ar], na.rm = T)
	max.ylim <- max(int.out[, norm.int.ar], na.rm = T)
	
	if(reproduce.plots.part4){
	png(filename = paste("Fig.15.Outlier detection and removal in normalized int data", tmp.names[i], "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
	par(mar = c(5, 4, 4, 2) + 0.1)
	par(mfrow=c(1,2))
	boxplot(t(int.out[, norm.int.ar]), las = 2, col = int.out[["class.color"]], ylab = "Z-score", ylim = c(min.ylim, max.ylim), main = paste(paste("Normalized.int.data with outliers", sep = "-"), int.names[i], sep = "_"), names = int.out[["Oligo"]])
	abline(h = c(-3:-1, -1.5, 1:3, 1.5), col = c(rep("blue", 4), rep("red", 4)))
	
	
	int.out[, norm.int.ar][int.out[, norm.int.ar] > z.cut.offs[2]] <- 0
	int.out[, norm.int.ar][int.out[, norm.int.ar] < z.cut.offs[1]] <- 0
	print(head(int.out))
	de.int.out <- t(int.out[, norm.int.ar])
	colnames(de.int.out) <- int.out[["Oligo"]]
	
	
	boxplot(de.int.out, las = 2, col = int.out[["class.color"]], ylab = "Z-score", ylim = c(min.ylim, max.ylim), main = paste(paste("De-outliered normalized.int.data with the z.score.of", "-2 to 3", sep = " "), int.names[i], sep = "_"))
	abline(h = c(-3:-1, -1.5, 1:3, 1.5), col = c(rep("blue", 4), rep("red", 4)))
	abline(h = c(-2, 3), lwd = 2.5, col = c("blue", "red"))
	
	par(mfrow=c(1,1))
	}
	dev.off()
	
	## Remove the outliers in the original normalized intensities
	int.out[is.na(int.out)] <- 0
	
	tmp.de.int.out <- NULL
	tmp.de.int.out <- na.rm.int.ini[[i]][, norm.int.ar]
	tmp.de.int.out <- data.matrix(tmp.de.int.out)
	tmp.de.int.out[is.na(tmp.de.int.out)] <- 0
	tmp.de.int.out[which(tmp.de.int.out == tmp.de.int.out) %in% which(0 == int.out[, norm.int.ar])] <- 0
	tmp.de.int.out <- as.data.frame(tmp.de.int.out, stringsAsFactors = F)
	tmp.de.int.out[, norm.int.ar][tmp.de.int.out[, norm.int.ar] == 0] <- NA
	tmp.de.int.out[["norm_rel.mean"]] <- apply(tmp.de.int.out[, norm.int.ar], 1, function(x) mean(x, na.rm = T))
	tmp.de.int.out[["norm_rel.sd"]] <- apply(tmp.de.int.out[, norm.int.ar], 1, function(x) sd(x, na.rm = T))
	tmp.de.int.out <- cbind(na.rm.int.ini[[i]][c("Oligo", "orig.order", "class", "class.color")], tmp.de.int.out)
	print(head(tmp.de.int.out))
	print(nrow(tmp.de.int.out))
	
	## view of normalized intensity data (like Cyclic Loess method) and observe the outliers
	## For the sake of comparision, previous scatter plot (normalized intensity data with outliers) have added  in this plot 

##################################################################
	if(reproduce.plots.part4){
	png(filename = paste("Fig.16.Comparision between normalized and outlierd removed normalized int data", tmp.names[i], "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
	par(mar = c(5, 4, 4, 2) + 0.1)
	par(mfrow = c(1, 2))
	cn.plot <- grep("^norm_rel_Lot", names(na.rm.int.ini[[i]]), value = T)
		plot(tmp.plot[["norm_rel.mean"]], tmp.plot[["norm_rel.mean"]], xlim = c(0, 1.1 * max(tmp.plot[["norm_rel.mean"]], na.rm = T)), ylim = c(0, max(as.matrix(na.rm.int.ini[[i]][,cn.plot]), na.rm = T)), type = "n", xlab = "mean intensity", ylab = "normalized mean intensity", main = "Normalized_Cyclic.Loess.method_background.corrected_data")
		r = 1
		for (r in 1:length(cn.plot)) {
			# cat(max(tmp.plot[[cn.plot[r]]], na.rm=T), "\n")
			points(tmp.plot[["rel.int_mean"]], tmp.plot[[cn.plot[r]]], col = col.plot[r], pch = 20, cex = 1)
			lines(loess.smooth(tmp.plot[["rel.int_mean"]], tmp.plot[[cn.plot[r]]]), col = col.plot[r], lwd = 2)
		}
	
	
	plot(tmp.de.int.out[["norm_rel.mean"]], tmp.de.int.out[["norm_rel.mean"]], xlim = c(0, 1.1 * max(tmp.de.int.out[["norm_rel.mean"]], na.rm = T)), ylim = c(0, max(as.matrix(na.rm.int.ini[[i]][,norm.int.ar]), na.rm = T)), type = "n", xlab = "mean intensity", ylab = "normalized mean intensity", main = "Outlier removed_Normalized.data")
		r = 1
		for (r in 1:length(norm.int.ar)) {
			points(tmp.de.int.out[["norm_rel.mean"]], tmp.de.int.out[[norm.int.ar[r]]], col = col.plot[r], pch = 20, cex = 1)
			lines(loess.smooth(tmp.de.int.out[["norm_rel.mean"]], tmp.de.int.out[[norm.int.ar[r]]]), col = col.plot[r], lwd = 2)
		}
	}
	dev.off()
	par(mfrow = c(1,1))
##################################################################
		
	## Additional box plot view of relative intensity means
	# norm.int <- NULL
	# norm.int <- na.rm.int.ini[[i]][, c("Oligo", paste("norm", rel.ar, sep = "_"))]
	# norm.int[["norm_rel_mean"]] <- apply(norm.int[, paste("norm", rel.ar, sep = "_")], 1, mean, na.rm = T)
	# print(head(norm.int))
	
	## Take the outlier removed background corrected and normlaized intensity data for further anaysis 
	## factor level re-ordering
	tmp.de.int.out[["Oligo"]] <- as.factor(tmp.de.int.out[["Oligo"]])
	tmp.de.int.out[["Oligo"]] <- reorder(tmp.de.int.out[["Oligo"]], na.rm.int.ini[[i]][["Oligo"]])
	
	## reducing the background intensities of scrambled oligos.
	
	tmp.norm.int <- tmp.de.int.out

	## Taking mean of scrambled oligos.
	# scr.mean <- mean(as.matrix(subset(tmp.norm.int, Oligo %in% grep("scr", tmp.norm.int[["Oligo"]], value=T))[paste("norm", rel.ar, sep="_")]), na.rm=T)
	# tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")] <- tmp.norm.int[paste("norm", rel.ar, sep = "_")] - scr.mean
	
	## Taking 90% probability of scrambled oligos.
	scr.quan.90 <- quantile(as.matrix(subset(tmp.norm.int, Oligo %in% grep("scr", tmp.norm.int[["Oligo"]], value=T))[paste("norm", rel.ar, sep="_")]), na.rm = T, probs = scramble.quant.cutoff)
	tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")] <- tmp.norm.int[paste("norm", rel.ar, sep = "_")] - scr.quan.90
	tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")][tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")] < 0] <- 0
	# # # tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")][tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")] == 0] <- NA

	
	## calculates reltive intensities to the reference oligo
	
	# tmp.norm.int[is.na(tmp.norm.int)] <- 0	
	# ref.mean <- mean(as.matrix(subset(tmp.norm.int, Oligo %in% grep("Ref", tmp.norm.int[["Oligo"]], value=T))[paste("b.c.norm", rel.ar, sep="_")]), na.rm=T)
	# tmp.norm.int[paste("rel.norm", rel.ar, sep = "_")] <- tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")] / ref.mean

	ref.quan.90 <- mean(as.matrix(subset(tmp.norm.int, Oligo %in% grep("Ref", tmp.norm.int[["Oligo"]], value=T))[paste("b.c.norm", rel.ar, sep="_")]), na.rm=T, probs = ref.quant.cutoff)
	tmp.norm.int[paste("rel.norm", rel.ar, sep = "_")] <- tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")] / ref.quan.90
		tmp.norm.int[paste("rel.norm", rel.ar, sep = "_")][tmp.norm.int[paste("rel.norm", rel.ar, sep = "_")] < 0] <- 0	
	# # # tmp.norm.int[paste("rel.norm", rel.ar, sep = "_")][tmp.norm.int[paste("rel.norm", rel.ar, sep = "_")] == 0] <- NA
	tmp.norm.int[["b.c.mean"]] <- apply(tmp.norm.int[paste("b.c.norm", rel.ar, sep = "_")], 1, function(x) mean(x, na.rm=T))
	tmp.norm.int[["rel.mean"]] <- apply(tmp.norm.int[paste("rel.norm", rel.ar, sep = "_")], 1, function(x) mean(x, na.rm=T))
	print(head(tmp.norm.int))
	
	pre.fin.int.data[[i]] <- tmp.norm.int
			
}
	
	
names(pre.fin.int.data) <- tmp.names
lapply(pre.fin.int.data, head)
lapply(pre.fin.int.data, nrow)

###########################################################
			###### Part - 4.2 Starts.. #######			
###########################################################
	
## Normalized and background corrected intensities (June and July lots) are combined together

cat("\n", "Combine the June and July lots together and Making the final data", "\n")

tmp.name <- "int.pos.E.tcf"
comb.int.data <- NULL
int.data.june <- pre.fin.int.data[["int.pos.E.tcf_june"]]
print(head(int.data.june))
print(nrow(int.data.june))

int.data.july <- pre.fin.int.data[["int.pos.E.tcf_july"]]
print(head(int.data.july))
print(nrow(int.data.july))

all(int.data.june[["Oligo"]] %in% int.data.july[["Oligo"]])

## finding missing oligos in both june and July slot before cinbine them togehter

miss.oli.june <- int.data.july[["Oligo"]][!int.data.july[["Oligo"]] %in% int.data.june[["Oligo"]]]
miss.oli.june <- as.data.frame(miss.oli.june)
miss.oli.june[2:length(int.data.june)] <- NA
names(miss.oli.june) <- names(int.data.june)
print(miss.oli.june)
int.data.june <- rbind(int.data.june, miss.oli.june)
names(int.data.june) <- paste(names(int.data.june), "june", sep = ".")

miss.oli.july <- int.data.june[["Oligo.june"]][!int.data.june[["Oligo.june"]] %in% int.data.july[["Oligo"]]]
miss.oli.july <- as.data.frame(miss.oli.july)
miss.oli.july[2:length(int.data.july)] <- NA
names(miss.oli.july) <- names(int.data.july)
print(miss.oli.july)
int.data.july <- rbind(int.data.july, miss.oli.july)
names(int.data.july) <- paste(names(int.data.july), "july", sep = ".")

## ordering data accoding to the Oligo information
int.data.june[["Oligo.june"]] <- as.character(int.data.june[["Oligo.june"]])
int.data.june <- int.data.june[order(int.data.june[["Oligo.june"]]), ]
int.data.july[["Oligo.july"]] <- as.character(int.data.july[["Oligo.july"]])
int.data.july <- int.data.july[order(int.data.july[["Oligo.july"]]), ]

comb.int.data <- cbind(int.data.june, int.data.july)
comb.int.data <- comb.int.data[!names(comb.int.data) %in% "Oligo.july"]
names(comb.int.data)[names(comb.int.data) %in% "Oligo.june"] <- "Oligo"

print(head(comb.int.data))
print(nrow(comb.int.data))


fin.int.data <- NULL
# i = 1
# for(i in 1:length(pre.fin.int.data)) {		
		
		norm.rel.ar <- grep("^norm_rel_Lot", names(comb.int.data), value = T)
		
		norm.test <- comb.int.data[, c("Oligo", norm.rel.ar)]
		# # # norm.test[norm.rel.ar][norm.test[norm.rel.ar] == 0] <- NA
		colnames(norm.test)[colnames(norm.test) %in% norm.rel.ar] <- gsub(".[01].*", "", norm.rel.ar)
		
		print(head(norm.test))
		
		norm.int.plot <- NULL
		for(l in 2:length(norm.test)){
			tmp.norm.int.plot <- norm.test[c(1, l)]
			norm.int.plot <- rbind(norm.int.plot, tmp.norm.int.plot)
		}
		
		print(head(norm.int.plot))
		print(nrow(norm.int.plot))
		
		## Making a data frame for backround corrected intensities. Mean of scramled Oligos
		## substracted from the normalized intensities.
		
		b.c.rel.ar <- grep("^b.c.norm_rel", names(comb.int.data), value = T)
		b.c.norm.test <- comb.int.data[, c("Oligo", b.c.rel.ar)]
		# # # b.c.norm.test[b.c.rel.ar][b.c.norm.test[b.c.rel.ar] == 0] <- NA
		colnames(b.c.norm.test)[colnames(b.c.norm.test) %in% b.c.rel.ar] <- gsub(".[01].*", "", b.c.rel.ar)
		print(head(b.c.norm.test))
		
		b.c.norm.int.plot <- NULL
		for(m in 2:length(b.c.norm.test)){
			tmp.b.c.norm.int.plot <- b.c.norm.test[c(1, m)]
			b.c.norm.int.plot <- rbind(b.c.norm.int.plot, tmp.b.c.norm.int.plot)
		}
		print(head(b.c.norm.int.plot))
		print(nrow(b.c.norm.int.plot))
		
		## Making a data frame for relative intensities. Mean of reference Oligos
		## divided from the background corrected intensities.
		
		rel.norm.rel.ar <- grep("^rel.norm_rel_Lot", names(comb.int.data), value = T)
		rel.norm.test <- comb.int.data[, c("Oligo", rel.norm.rel.ar)]
		# # # rel.norm.test[rel.norm.rel.ar][rel.norm.test[rel.norm.rel.ar] == 0] <- NA
		colnames(rel.norm.test)[colnames(rel.norm.test) %in% rel.norm.rel.ar] <- gsub(".[01].*", "", rel.norm.rel.ar)
		print(head(rel.norm.test))
				
		rel.norm.int.plot <- NULL
		for(m in 2:length(rel.norm.test)){
			tmp.rel.norm.int.plot <- rel.norm.test[c(1, m)]
			rel.norm.int.plot <- rbind(rel.norm.int.plot, tmp.rel.norm.int.plot)
		}
		print(head(rel.norm.int.plot))
		print(nrow(rel.norm.int.plot))	
		                        
		## making a long data frame which contains all intensities. for plotting purpose 
		
		fin.int.data <- cbind(norm.int.plot, b.c.norm.int.plot[2], rel.norm.int.plot[2])
		fin.int.data[["Oligo"]] <- as.factor(fin.int.data[["Oligo"]])
		fin.int.data[["Oligo"]] <- reorder(fin.int.data[["Oligo"]])
		
		## data frame for oligo color
		ordered.oligo.names <- NULL
		ordered.oligo.names[["oligo"]] <- levels(unique(fin.int.data[["Oligo"]]))
		ordered.oligo.names <- as.data.frame(ordered.oligo.names)
		ordered.oligo.names[["oligo.col"]] <- NA
		ordered.oligo.names[["oligo.col"]][ordered.oligo.names[["oligo"]] %in% grep("7L2_2.Ref2", ordered.oligo.names[["oligo"]], value = T)] <- "green"
		ordered.oligo.names[["oligo.col"]][ordered.oligo.names[["oligo"]] %in% grep("[ACGT]$", ordered.oligo.names[["oligo"]], value = T)] <- "red"
		ordered.oligo.names[["oligo.col"]][ordered.oligo.names[["oligo"]] %in% grep("scr", ordered.oligo.names[["oligo"]], value = T)] <- "grey"
		ordered.oligo.names[["oligo.col"]][ordered.oligo.names[["oligo"]] %in% grep("^mv", ordered.oligo.names[["oligo"]], value = T)] <- "orange"
		
		
		print(head(fin.int.data))
		print(nrow(fin.int.data))
		
		## boxplot for normalized intensities
		setwd(plot.folder)		
		png(filename = paste("Fig.17.boxplot_nomalized_intensities_quantile.90", tmp.name, "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(6, 5, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		boxplot(fin.int.data[["norm_rel_Lot"]] ~ fin.int.data[["Oligo"]] , las = 2, col = ordered.oligo.names[["oligo.col"]], ylab = "Mean Value")
		par(mfrow = c(1, 1))
		dev.off()
		
		png(filename = paste("Fig.18.boxplot_b.c.nomalized_intensities_quantile.90", tmp.name, "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(6, 5, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		boxplot(fin.int.data[["b.c.norm_rel_Lot"]] ~ fin.int.data[["Oligo"]] , las = 2, ylab = "Mean Value", col = ordered.oligo.names[["oligo.col"]])
		par(mfrow = c(1, 1))
		dev.off()
		
		png(filename = paste("Fig.19.boxplot_rel_normalized_intensities_with_reference_intensities_quantile.90", tmp.name, "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(6, 5, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		boxplot(fin.int.data[["rel.norm_rel_Lot"]] ~ fin.int.data[["Oligo"]] , las = 2, ylab = "Relative Mean Value", col = ordered.oligo.names[["oligo.col"]])
		abline(h = c(0, 0.5, 1), lty=2)
		par(mfrow = c(1, 1))
		dev.off()
		
		setwd(ouput.folder)
		## save the final processed intensity data from this part.	
		write.table(fin.int.data, file = paste(paste("pre.screened.intensity.data_all_quantile.90", tmp.name, sep = "_"), "txt", sep = "."), quote = F, sep = "\t", row.names = F)
		write.table(fin.int.data, file = paste(paste("pre.screened.intensity.data_all_quantile.90", tmp.name, sep = "_"), "txt", sep = "."), quote = F, sep = "\t", row.names = F)
	
	# }
# }

print(head(fin.int.data))
print(nrow(fin.int.data))

###########################################################
######### Part - 5 starts... ##########
###########################################################

setwd(infolder)
oligo.seq <- read.table("microdotOligos_TCF7L2_all.long.txt", header = T, sep = "\t", stringsAsFactors = F)
print(head(oligo.seq))
print(nrow(oligo.seq))


## Making final data for post processing.This data should contain all neccessary information 
## calculating the occurence of each oligo 

cat("\n", "Post processing and producing final data", "\n")

setwd(plot.folder)
post.sv.int.data <- NULL
post.mv.int.data <- NULL
# i=1
# for (i in 1:length(fin.int.data)) {
	
	## finding the missing oligos from the single varaint oligo set
	
	sv.oligo <- unique(grep("[ACGT]$", fin.int.data[["Oligo"]], value = T))
	mv.oligo <- unique(grep("^mv", fin.int.data[["Oligo"]], value = T))
	all.sv.oligo <- paste(paste("7L2_2", rep(1:20, each = 4), sep = "."), c("A", "G", "C", "T"), sep = "")
	missing.sv.oligos <- all.sv.oligo[!all.sv.oligo %in% sv.oligo]

	sv.oli.int.1 <- subset(fin.int.data, Oligo %in% sv.oligo)
	sv.oli.int.1 <- rbind(sv.oli.int.1, subset(fin.int.data, Oligo %in% "7L2_2.Ref2"))
	sv.oli.int <- aggregate(sv.oli.int.1[["norm_rel_Lot"]], by = list(sv.oli.int.1[["Oligo"]]), mean, na.rm = T)
	sv.oli.int <- merge(sv.oli.int, aggregate(sv.oli.int.1[["b.c.norm_rel_Lot"]], by = list(sv.oli.int.1[["Oligo"]]), 	mean, na.rm = T), by = "Group.1")
	sv.oli.int <- merge(sv.oli.int, aggregate(sv.oli.int.1[["rel.norm_rel_Lot"]], by = list(sv.oli.int.1[["Oligo"]]), 	mean, na.rm = T), by = "Group.1")
	names(sv.oli.int) <- c("Oligo", "norm_rel_Lot", "b.c.norm_rel_Lot", "rel.norm_rel_Lot")
	
	## adding missing oligos or reference oligos in the single oligo variant set
	 
	sv.oli.int <- rbind(sv.oli.int, data.frame(cbind(Oligo = missing.sv.oligos, norm_rel_Lot = 1, b.c.norm_rel_Lot = 1, rel.norm_rel_Lot = 1), stringsAsFactors = F))
	sv.oli.int[["Oligo"]] <- as.character(sv.oli.int[["Oligo"]])
	sv.oli.int[["norm_rel_Lot"]] <- as.numeric(sv.oli.int[["norm_rel_Lot"]])
	sv.oli.int[["b.c.norm_rel_Lot"]] <- as.numeric(sv.oli.int[["b.c.norm_rel_Lot"]])
	sv.oli.int[["rel.norm_rel_Lot"]] <- as.numeric(sv.oli.int[["rel.norm_rel_Lot"]])
	sv.oli.int[is.na(sv.oli.int)] <- 0

	print(head(sv.oli.int))
	print(nrow(sv.oli.int))	
	
	## calculating the oocurence of each oligo in the single variant set
	
	fin.sv.int <- NULL
	sv.plot <- data.frame(Oligo = c("A", "C", "G", "T"))
	p = 1
	for (p in 1:20) {
		tmp.sv.oli.int <- subset(sv.oli.int, Oligo %in% paste(paste("7L2_2", p, sep = "."), c("A", "G", "C", "T"), sep = ""))
		tmp.ref.int <- subset(sv.oli.int, Oligo %in% "7L2_2.Ref2")
		for (q in 1:nrow(tmp.sv.oli.int)) {
			tmp.sv.oli.int[["occ.oligo"]][q] <- tmp.sv.oli.int[["rel.norm_rel_Lot"]][q]/(sum(tmp.sv.oli.int[["rel.norm_rel_Lot"]], na.rm = T))
		}
		fin.sv.int <- rbind(fin.sv.int, tmp.sv.oli.int)
		tmp.sv.plot <- tmp.sv.oli.int[, c("Oligo", "occ.oligo")]
		tmp.sv.plot <- tmp.sv.plot[order(tmp.sv.plot[["Oligo"]]),]
		names(tmp.sv.plot) <- c("Oligo", p)
		tmp.sv.plot[["Oligo"]] <- gsub(".*[0-9]", "", tmp.sv.plot[["Oligo"]])
		sv.plot <- merge(sv.plot, tmp.sv.plot, by = "Oligo")
		
	}
	fin.sv.int[["variant.pos"]] <- gsub("^7L2_2.", "", fin.sv.int[["Oligo"]])
	print(head(fin.sv.int))
	print(nrow(fin.sv.int))

	## sequence Logo for single variant set. 
	
	print(sv.plot)
		
	m.sv.plot <- sv.plot[! names(sv.plot) %in% "Oligo"]
	
	pwm <- makePWM(m.sv.plot)
	
	## Multi variant oligos do not have the mutation information. In order to find the mutation 
	## and their position info, obtain information from the oligo seq file. The occurence of  
	## multi variant (MV) oligos are identified based on single variant set.  
	## binding sequences are also obtained from the seq file for MV set. 
	
	mv.oli.int.1 <- subset(fin.int.data, Oligo %in% mv.oligo)
	mv.oli.int.1 <- rbind(mv.oli.int.1, subset(fin.int.data, Oligo %in% "7L2_2.Ref2"))
	mv.oli.int <- aggregate(mv.oli.int.1[["norm_rel_Lot"]], by = list(mv.oli.int.1[["Oligo"]]), mean, na.rm = T)
	mv.oli.int <- merge(mv.oli.int, aggregate(mv.oli.int.1[["b.c.norm_rel_Lot"]], by = list(mv.oli.int.1[["Oligo"]]), 	mean, na.rm = T), by = "Group.1")
	mv.oli.int <- merge(mv.oli.int, aggregate(mv.oli.int.1[["rel.norm_rel_Lot"]], by = list(mv.oli.int.1[["Oligo"]]), 	mean, na.rm = T), by = "Group.1")
	names(mv.oli.int) <- c("Oligo", "norm_rel_Lot", "b.c.norm_rel_Lot", "rel.norm_rel_Lot")
	
	mv.oligo.seq <- subset(oligo.seq, name %in% grep("mv2", oligo.seq[["name"]], value = T))
	ref.seq <- subset(oligo.seq, name %in% "7L2_2.Ref2")
	mv.oligo.seq <- rbind(ref.seq, mv.oligo.seq)
	mv.oligo.seq[["bind.seq"]] <- apply(mv.oligo.seq, 1, function(x) substr(x[["seq"]], 12, 31))
	ref.bind.seq <- subset(mv.oligo.seq, name %in% "7L2_2.Ref2")[["bind.seq"]]
	
	## This loop gets the variant and their postional information from the single variant set
	## for the multi variant set. 
	
	var.info <- NULL
	for (j in 1:nrow(mv.oligo.seq)) {
		tmp.seq <- mv.oligo.seq[["bind.seq"]][j]
		test.mv.pos <- NULL
		for (k in 1:nchar(ref.bind.seq)) {
			if (!unlist(strsplit(tmp.seq, ""))[k] %in% unlist(strsplit(ref.bind.seq, ""))[k]){
				tmp.mv.seq.pos <- paste(k, unlist(strsplit(tmp.seq, ""))[k][! unlist(strsplit(tmp.seq, ""))[k] %in% unlist(strsplit(ref.bind.seq, ""))[k]], sep = "")
				test.mv.pos <- paste(test.mv.pos, tmp.mv.seq.pos, sep = ".")
			}
		}
		tmp.var.info <- paste(mv.oligo.seq[["name"]][j], test.mv.pos, sep = "")
		var.info <- c(var.info, tmp.var.info)
	}
	
	mv.oligo.seq[["variant.info"]] <- var.info
	names(mv.oligo.seq)[names(mv.oligo.seq) %in% "name"] <- "Oligo"
	print(head(mv.oligo.seq))
	print(nrow(mv.oligo.seq))
	
	mv.oli.int <- merge(mv.oli.int, mv.oligo.seq, by = "Oligo")
	mv.oli.int <- mv.oli.int[c("Oligo", "norm_rel_Lot", "b.c.norm_rel_Lot", "rel.norm_rel_Lot", "core", "bind.seq", "variant.info")]
	mv.oli.int <- subset(mv.oli.int, ! Oligo %in% "7L2_2.Ref2")
	print(head(mv.oli.int))
	print(nrow(mv.oli.int))
	
	## collecting the effect of the multi variant set on single variation set
	
	tmp.mv <- mv.oli.int[c("Oligo", "variant.info")]
	tmp.mv[["var.pos"]] <- apply(tmp.mv, 1, function(x) unlist(strsplit(x[["variant.info"]], "[.]")))
	tmp.mv[["var.pos"]] <- apply(tmp.mv, 1, function(x) x[["var.pos"]][x[["var.pos"]] %in% paste(rep(1:20, each = 4), c("A", "G", "C", "T"), sep = "")])
				
	variants <- c("A", "G", "C", "T")
	y = 1
	z = 1
	tmp.sv.mv.int <- NULL
	for (y in 1:nrow(tmp.mv)) {
		tmp.mv.int <- NULL
		for (z in 1:length(variants)) {
		 	tmp.var.int <- NULL
		 	tmp.var.int <- subset(fin.sv.int, variant.pos %in% tmp.mv[["var.pos"]][[y]])["rel.norm_rel_Lot"][z, ]
		 	tmp.mv.int <- c(tmp.mv.int, tmp.var.int)
		 }
		 tmp.sv.mv.int <- rbind(tmp.sv.mv.int, tmp.mv.int)
	}
	
	tmp.mv[paste("effect.sv", 1:length(variants), sep = ".")] <- tmp.sv.mv.int
	tmp.mv <- tmp.mv[! names(tmp.mv) %in% "var.pos"]
	print(head(tmp.mv))
	print(nrow(tmp.mv))
	
	## merging MV set sequence info and effect of SV set on MV set information together
	mv.oli.int <- mv.oli.int[! names(mv.oli.int) %in% "variant.info"]
	fin.mv.int <- merge(mv.oli.int, tmp.mv, by = "Oligo")
	
	print(head(fin.mv.int))
	print(nrow(fin.mv.int))
	
	## single varaint's effects are extracted for each multi variant oligo
	## In order to evaluate the influences of MV set, we further need to go for 
	## stepwise calculation (classifier converter)
	## The classifier converter allows the relative intensity values of SV set for each MV set
		
	class.con <- data.frame(cbind(min = c(0, 1, 0.7, 0.4, 0, 0.6, 0.3, 0), max = c(Inf, Inf, 1.0, 0.7, 0.4, 1, 0.6, 0.3), target = c(Inf, 1.0, 0.8, 0.5, 0.25, 0.75, 0.45, 0)), row.names = c("no.change", "class.v1", "class.v2.q1", "class.v2.q2", "class.v2.q3", "class.v3.q1", "class.v3.q2", "class.v3.q3"))
	
	## Prediction model.1 can be made by with out changing any intensity values of SV effect 
	## on MV set. 
	
	fin.mv.int[paste("no.ch", 1:4, sep =".")] <- fin.mv.int[paste("effect.sv", 1:4, sep =".")]
	
	## Similiarly, Prediction model.2 has relative intensities of SV set which have a effect 
	## on MV set and intesnity value which is more than 1 has been turned as 1
	
	fin.mv.int[paste("class.v1", 1:4, sep =".")] <- fin.mv.int[paste("no.ch", 1:4, sep =".")]
	fin.mv.int[paste("class.v1", 1:4, sep =".")][fin.mv.int[paste("class.v1", 1:4, sep =".")] > class.con["class.v1",][["min"]]] <- class.con["class.v1",][["target"]]
	
	## Prediction model.3 has modified relative intensities, intensities which are in the range 
	## of 0.7-1, 0.4-0.7 & 0.1-0.4 are modified into 0.8, 0.50 & 0.25 respectively
	fin.mv.int[paste("class.v2", 1:4, sep =".")] <- fin.mv.int[paste("class.v1", 1:4, sep =".")]
	fin.mv.int[paste("class.v2", 1:4, sep =".")][fin.mv.int[paste("class.v2", 1:4, sep =".")] > class.con["class.v2.q1",][["min"]] & fin.mv.int[paste("class.v2", 1:4, sep =".")] < class.con["class.v2.q1",][["max"]]] <- class.con["class.v2.q1",][["target"]]
	fin.mv.int[paste("class.v2", 1:4, sep =".")][fin.mv.int[paste("class.v2", 1:4, sep =".")] > class.con["class.v2.q2",][["min"]] & fin.mv.int[paste("class.v2", 1:4, sep =".")] < class.con["class.v2.q2",][["max"]]] <- class.con["class.v2.q2",][["target"]]
	fin.mv.int[paste("class.v2", 1:4, sep =".")][fin.mv.int[paste("class.v2", 1:4, sep =".")] > class.con["class.v2.q3",][["min"]] & fin.mv.int[paste("class.v2", 1:4, sep =".")] < class.con["class.v2.q3",][["max"]]] <- class.con["class.v2.q3",][["target"]]	
	
	## Prediction model.4 has modified relative intensities, intensities which are in the range 
	## of 0.6-0.9, 0.3-0.6 & 0.0 - 0.3 are modified into 0.75, 0.45 & 0 respectively
	fin.mv.int[paste("class.v3", 1:4, sep =".")] <- fin.mv.int[paste("class.v1", 1:4, sep =".")]
	fin.mv.int[paste("class.v3", 1:4, sep =".")][fin.mv.int[paste("class.v3", 1:4, sep =".")] > class.con["class.v3.q1",][["min"]] & fin.mv.int[paste("class.v3", 1:4, sep =".")] < class.con["class.v3.q1",][["max"]]] <- class.con["class.v3.q1",][["target"]]
	fin.mv.int[paste("class.v3", 1:4, sep =".")][fin.mv.int[paste("class.v3", 1:4, sep =".")] > class.con["class.v3.q2",][["min"]] & fin.mv.int[paste("class.v3", 1:4, sep =".")] < class.con["class.v3.q2",][["max"]]] <- class.con["class.v3.q2",][["target"]]
	fin.mv.int[paste("class.v3", 1:4, sep =".")][fin.mv.int[paste("class.v3", 1:4, sep =".")] > class.con["class.v3.q3",][["min"]] & fin.mv.int[paste("class.v3", 1:4, sep =".")] < class.con["class.v3.q3",][["max"]]] <- class.con["class.v3.q3",][["target"]]
	
	
	x = 1
	for (x in 1:nrow(fin.mv.int)) {
		## making prediction model for no change model
		# fin.mv.int[paste("no.ch", 1:4, sep = ".")][x,] <- sort(fin.mv.int[paste("no.ch", 1:4, sep = ".")][x,])
		# fin.mv.int[paste("no.ch", 1:4, sep =".")][x,][is.na(fin.mv.int[paste("no.ch", 1:4, sep =".")])[x,]] <- 1
		fin.mv.int[["pred.mod.no.ch"]][x] <- prod(fin.mv.int[paste("no.ch", 1:4, sep = ".")][x,], na.rm = T)
		
		## making prediction for classifier version 1		
		fin.mv.int[["pred.mod.cl.v1"]][x] <- prod(fin.mv.int[paste("class.v1", 1:4, sep = ".")][x,], na.rm = T)
				
		## making prediction for classifier version 2		
		fin.mv.int[["pred.mod.cl.v2"]][x] <- prod(fin.mv.int[paste("class.v2", 1:4, sep = ".")][x,], na.rm = T)
				
		## making prediction for classifier version 3				
		fin.mv.int[["pred.mod.cl.v3"]][x] <- prod(fin.mv.int[paste("class.v3", 1:4, sep = ".")][x,], na.rm = T)
				
	}
	print(head(fin.mv.int))
	print(nrow(fin.mv.int))
	
	## creating post processing data. Basically, post processing or final data contains all
	## neccesary information such as nomalized, background corrected and relatve intensities 
	## along with the occrrence of each oligo's info. 
	
	post.sv.int.data <- fin.sv.int
	post.sv.int.data[["Oligo"]] <- as.factor(post.sv.int.data[["Oligo"]])
	post.sv.int.data[["Oligo"]] <- reorder(post.sv.int.data[["Oligo"]])
	post.mv.int.data <- fin.mv.int
	post.mv.int.data[["Oligo"]] <- as.factor(post.mv.int.data[["Oligo"]])
	post.mv.int.data[["Oligo"]] <- reorder(post.mv.int.data[["Oligo"]])

	
	# names(post.sv.int.data)[i] <- tmp.names[i]
	print(head(post.sv.int.data))
	print(tail(post.sv.int.data))
	print(nrow(post.sv.int.data))
	# names(post.mv.int.data)[i] <- tmp.names[i]
	print(head(post.mv.int.data))
	print(tail(post.mv.int.data))
	print(nrow(post.mv.int.data))
	
	## plots for single variant and all oligo sets
	setwd(plot.folder)
	if (reproduce.plots.part5) {
		
		png(filename = paste(paste("Fig.20.Sequence.logo_single.variant.occurrence_ic.scale.F_quantile.90", tmp.name, sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(8, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		seqLogo(pwm)
		dev.off()
	
		png(filename = paste(paste("Fig.21.boxplot.for.occurence_of_each_oligo_quantile.90", tmp.name, sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mar = c(8, 4, 4, 2) + 0.1)
		par(mfrow = c(1, 1))
		boxplot(post.sv.int.data[["occ.oligo"]] ~ post.sv.int.data[["Oligo"]], , ylim = c(0, 1), las = 2)
		par(mfrow = c(1,1))
		dev.off()	
		
		png(filename = paste(paste("Fig.22.Measured relative intensity vs predicted model_for.MV.set_quantile.90", tmp.name, sep = "_"), "png", sep = "."), width = 12, height = 9, res = 200, unit = "in")
		par(mfrow = c(2,2))
		plot(post.mv.int.data[["rel.norm_rel_Lot"]], post.mv.int.data[["pred.mod.no.ch"]], col = "red", cex = 0.7, pch = 19, xlab = "Measured Relative intensities", ylab = "Predicted Model", main = paste("measured relative intensity vs predicted model_no.ch for", tmp.name, sep = "_"), xlim = c(0, max(post.mv.int.data[["pred.mod.no.ch"]], na.rm = T)), ylim = c(0, max(post.mv.int.data[["pred.mod.no.ch"]], na.rm = T)))
		abline(a = 0, b = 1)
		plot(post.mv.int.data[["rel.norm_rel_Lot"]], post.mv.int.data[["pred.mod.cl.v1"]], col = "red", cex = 0.7, pch = 19, xlab = "Measured Relative intensities", ylab = "Predicted Model", main = paste("measured relative intensity vs predicted model_v1 for", tmp.name, sep = "_"), xlim = c(0, max(post.mv.int.data[["pred.mod.no.ch"]], na.rm = T)), ylim = c(0, max(post.mv.int.data[["pred.mod.no.ch"]], na.rm = T)))
		abline(a = 0, b = 1)
		plot(post.mv.int.data[["rel.norm_rel_Lot"]], post.mv.int.data[["pred.mod.cl.v2"]], col = "red", cex = 0.7, pch = 19, xlab = "Measured Relative intensities", ylab = "Predicted Model", main = paste("measured relative intensity vs predicted model_v2 for", tmp.names[i], sep = "_"), xlim = c(0, max(post.mv.int.data[["pred.mod.no.ch"]], na.rm = T)), ylim = c(0, max(post.mv.int.data[["pred.mod.no.ch"]], na.rm = T)))
		abline(a = 0, b = 1)
		plot(post.mv.int.data[["rel.norm_rel_Lot"]], post.mv.int.data[["pred.mod.cl.v3"]], col = "red", cex = 0.7, pch = 19, xlab = "Measured Relative intensities", ylab = "Predicted Model", main = paste("measured relative intensity vs predicted model_v3 for", tmp.names[i], sep = "_"), xlim = c(0, max(post.mv.int.data[["pred.mod.no.ch"]], na.rm = T)), ylim = c(0, max(post.mv.int.data[["pred.mod.no.ch"]], na.rm = T)))
		abline(a = 0, b = 1)
		par(mfrow = c(1,1))
		dev.off()
	}
	
	
	setwd(ouput.folder)
	## saving the post proceesing or final data.
	write.table(post.sv.int.data, file = paste(paste("post.processing.SV.set.intensity.data_quantile.90", tmp.name, sep = "_"), "txt", sep = "."), quote = F, sep = "\t", row.names = F)
	write.table(post.mv.int.data, file = paste(paste("post.processing.MV.set.intensity.data_quantile.90", tmp.name, sep = "_"), "txt", sep = "."), quote = F, sep = "\t", row.names = F)

# }

head(post.sv.int.data)
head(post.mv.int.data)
names(post.sv.int.data)
names(post.mv.int.data)


###################################################################

list.files()
date()
gc()

###################################################################

