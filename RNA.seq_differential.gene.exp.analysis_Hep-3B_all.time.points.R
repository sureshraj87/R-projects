## Author : Suresh Raju
## Date : 12 May 2016


## Objective
## Differential gene expression analysis

rm(list=ls())
gc()


## Folder and initial settings. 

filefolder <- "/Volumes/SureshData/RNAseq_Hep3B.siST.timecourse_Raju_2015_hg38_all_042016/results"
plotfolder <- "/Volumes/SureshData/RNAseq_Hep3B.siST.timecourse_Raju_2015_hg38_all_042016/plots"

library(gplots)
library(amap)

version.marker <- gsub(" ","",gsub("..:..:..", "", date()))

tp.ids <- c("03h","06h","09h","12h","15h","18h","24h","48h","mock","untreated")
tp.col <- rev(c("red", "darkorange", "purple", "dodgerblue", "blue", "darkgreen", "green", "grey50", "orange4", "plum2"))

######################################################################

## Read all the necessary files

setwd(filefolder)

load("genes.all.merged_all.time.points__g.all.Rdata")
load("genes.all.cuffdiff__g.diff.Rdata")
load("genes.siginificant.cuffdiff__tmp.sig.g.Rdata")

## Generating PCA plots for all samples

head(g.all);nrow(g.all)
pca.gene.fpkm <- g.all[,grep("^f.", names(g.all), value = T)]
names(pca.gene.fpkm) <- gsub("f.", "", names(pca.gene.fpkm))
rownames(pca.gene.fpkm) <- g.all[["tracking_id"]]
head(pca.gene.fpkm);nrow(pca.gene.fpkm)
m.pca.gene.fpkm <- t(pca.gene.fpkm)
m.pca.gene.fpkm <- prcomp(m.pca.gene.fpkm, cor = T)
summary(m.pca.gene.fpkm)
predict(m.pca.gene.fpkm)
# library(ggbiplot)
plot(m.pca.gene.fpkm, type = "lines")

setwd(plotfolder)
pdf(paste(paste("PCA for raw data", sep="_"), "pdf", sep="."))
par(mar = c(6, 5, 3, 1))
par(mfrow = c(1,1))
biplot(m.pca.gene.fpkm, expand = 0)
dev.off()


## Cluster dendrogram for all samples

rownames(pca.gene.fpkm) <- NULL
tmp.clusters <- t(pca.gene.fpkm)
clusters.gene.fpkm <- dist(tmp.clusters)
clusters.gene.fpkm <- hclust(clusters.gene.fpkm)


pdf(paste(paste("Cluster dendrogram for raw data", sep="_"), "pdf", sep="."))
par(mar = c(6, 5, 3, 1))
par(mfrow = c(1,1))
plot(clusters.gene.fpkm, las = 2, xlab = "Time points", labels = NULL, hang = 0.15, ann = T)
dev.off()

########################################################################

## Find out the outliers , calculate z score 

rel.gene.fpkm <- g.all[, c("tracking_id", "gene_short_name", "tss_id", "locus", grep("f.", names(g.all), value = T))]
tmp.rel.gene.fpkm <- rel.gene.fpkm
fpkm.all.tp <- grep("^f", names(tmp.rel.gene.fpkm), value = T)


j = 1
for (j in 1:length(rel.gene.fpkm[["tracking_id"]])) {
	tmp.data <- as.matrix(tmp.rel.gene.fpkm[which(tmp.rel.gene.fpkm[["tracking_id"]] == rel.gene.fpkm[["tracking_id"]][j]), ][grep("^f", names(tmp.rel.gene.fpkm), value = T)])
	tmp.data[tmp.data == 0] <- NA
	print(j)
	print(mean(as.vector(tmp.data), na.rm = T))
	print(median(as.vector(tmp.data), na.rm = T))
	tmp.rel.gene.fpkm[which(tmp.rel.gene.fpkm[["tracking_id"]] == rel.gene.fpkm[["tracking_id"]][j]), ][grep("^f", names(tmp.rel.gene.fpkm), value = T)] <- (tmp.data - mean(as.vector(tmp.data), na.rm = T))/sd(as.vector(tmp.data), na.rm = T)
	cat("\n\n")
	
}

head(tmp.rel.gene.fpkm);nrow(tmp.rel.gene.fpkm)

## Outlier detection and removal 
## Identify the gene that deviate by more than 3 SD.	

z.cut.offs <- c(-2,3)


tmp.rel.gene.fpkm[, fpkm.all.tp][tmp.rel.gene.fpkm[, fpkm.all.tp] < z.cut.offs[1]] <- 0
tmp.rel.gene.fpkm[, fpkm.all.tp][tmp.rel.gene.fpkm[, fpkm.all.tp] > z.cut.offs[2]] <- 0
head(tmp.rel.gene.fpkm);nrow(tmp.rel.gene.fpkm)

## Removing outliers in the original FPKM file 

out.rm.rel.gene.fpkm <- rel.gene.fpkm
out.rm.rel.gene.fpkm <- data.matrix(out.rm.rel.gene.fpkm[, fpkm.all.tp])
out.rm.rel.gene.fpkm[is.na(out.rm.rel.gene.fpkm)] <- 0
out.rm.rel.gene.fpkm[which(out.rm.rel.gene.fpkm == out.rm.rel.gene.fpkm) %in% which(0 == tmp.rel.gene.fpkm[, fpkm.all.tp])] <- 0
out.rm.rel.gene.fpkm <- as.data.frame(out.rm.rel.gene.fpkm, stringsAsFactors = F)
out.rm.rel.gene.fpkm[out.rm.rel.gene.fpkm == 0] <- NA
out.rm.rel.gene.fpkm <- cbind(rel.gene.fpkm[c("tracking_id", "gene_short_name", "tss_id", "locus")], out.rm.rel.gene.fpkm)

head(out.rm.rel.gene.fpkm);nrow(out.rm.rel.gene.fpkm)

save(out.rm.rel.gene.fpkm, file="genes.all.outlier.detected.all.time.points__ out.tmp.rel.gene.fpkm.Rdata")



######## Heat maps and further plots for raw fpkm values for all genes ########

library(gplots)
library(amap)

## Heat map for raw data for all genes

heat.raw <- g.all
heat.raw <- heat.raw[!duplicated(heat.raw[["gene_short_name"]]), ]
heat.raw <- heat.raw[, c("gene_short_name", fpkm.all.tp)]
names(heat.raw)[1] <- "gene"
head(heat.raw);nrow(heat.raw)
heat.raw.plot <- heat.raw[, fpkm.all.tp]
names(heat.raw.plot) <- gsub("f.", "", names(heat.raw.plot))
rownames(heat.raw.plot) <- heat.raw[["gene"]]
heat.raw.plot <- head(heat.raw.plot, 5000)
names(heat.raw.plot)
heat.side.col <- c(rep(tp.col[1:8], each = 3), rep(tp.col[1:8], each = 3), rep(tp.col[9:10], each = 3))
heat.raw.plot <- as.matrix(heat.raw.plot)
str(heat.raw.plot)

qtile.cut <- quantile(heat.raw.plot, 0.99)
heat.raw.plot[heat.raw.plot > qtile.cut] <- qtile.cut

heat.breaks.raw <- seq(0, max(heat.raw.plot), max(heat.raw.plot)/200)

mycol.raw <-colorRampPalette(c("lightblue", "dodgerblue", "lightgreen", "palevioletred3", "darkorange", "red"), bias=3)(n=length(heat.breaks.raw)-1)


setwd(plotfolder)
pdf(paste(paste("Heatmap for 5000 genes_fpkm values_non scaled", sep="_"), "pdf", sep="."))
par(mar = c(4, 7, 2, 1))
par(mfrow = c(1,2))
heatmap.2(heat.raw.plot, dendrogram="column", scale="none", col= mycol.raw, breaks = heat.breaks.raw, lhei = c(2,8), trace="none", margins=c(6, 13), na.rm = T, main="Top 5000 genes_scale.none_raw data", labRow = "", ColSideColors = heat.side.col)
dev.off()


pdf(paste(paste("Heatmap for 5000 genes_fpkm values_scaled", sep="_"), "pdf", sep="."))
par(mar = c(4, 7, 2, 1))
par(mfrow = c(2,1))
heatmap.2(heat.raw.plot, distfun = function(x) dist(x, method = "correlation"), dendrogram = "none", scale="row", col= bluered(100), Colv = F, Rowv = F, trace="none", margins=c(6, 13), lhei=c(2, 8), main="Top 5000.genes_Correlation_scaled", labRow = "", ColSideColors = heat.side.col)
par(mfrow = c(1,1))
dev.off()

################################################################


### calculate the the relative mean for siTCF7L2
## Replace the Inf and NaN values 
nrow(rel.gene.fpkm);names(rel.gene.fpkm)

repl <- 0:2
i = 1
j = 1
for (i in 1:length(tp.ids[1:8])){
	for (j in 1:length(repl)){
		rel.gene.fpkm[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")] <- log2(rel.gene.fpkm[paste(paste("f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]/rel.gene.fpkm[paste("m.f.siS", tp.ids[i], sep = "_")])
		## Remove NaN with NAs
		rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% NaN] <- NA
		## Remove Inf and -Inf values
		rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% Inf] <- 0
		a <- max(rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]], na.rm = T)
		rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% 0] <- a		
		
		rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% -Inf] <- 0
		b <- min(rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]], na.rm = T)
		rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% 0] <- b

		## Make changes in relative siT
		rel.gene.fpkm[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")] <- log2(rel.gene.fpkm[paste(paste("f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]/rel.gene.fpkm[paste("m.f.siS", tp.ids[i], sep = "_")])		
		rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% NaN] <- NA	

		rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% Inf] <- 0
		c <- max(rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]], na.rm = T)
		rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% 0] <- c
		
		rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% -Inf] <- 0
		d <- min(rel.gene.fpkm[[paste(paste("rel.f.siT", tp.ids[i], sep = "_"), repl[j], sep = "_")]], na.rm = T)	
		rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]][rel.gene.fpkm[[paste(paste("rel.f.siS", tp.ids[i], sep = "_"), repl[j], sep = "_")]] %in% 0] <- d	
	
		
			
	}
	rel.gene.fpkm[paste("m.rel.f.siS", tp.ids[i], sep = "_")] <- apply(rel.gene.fpkm[grep(paste("^rel.f.siS", tp.ids[i], sep = "_"), names(rel.gene.fpkm), value = T)], 1, mean, na.rm = T)
	rel.gene.fpkm[paste("m.rel.f.siT", tp.ids[i], sep = "_")] <- apply(rel.gene.fpkm[grep(paste("^rel.f.siT", tp.ids[i], sep = "_"), names(rel.gene.fpkm), value = T)], 1, mean, na.rm = T)
	rel.gene.fpkm[paste("sd.rel.f.siT", tp.ids[i], sep = "_")] <- apply(rel.gene.fpkm[grep(paste("^rel.f.siT", tp.ids[i], sep = "_"), names(rel.gene.fpkm), value = T)], 1, sd, na.rm = T)
	rel.gene.fpkm[paste("sd.rel.f.siS", tp.ids[i], sep = "_")] <- apply(rel.gene.fpkm[grep(paste("^rel.f.siS", tp.ids[i], sep = "_"), names(rel.gene.fpkm), value = T)], 1, sd, na.rm = T)
	
}

head(rel.gene.fpkm);nrow(rel.gene.fpkm);names(rel.gene.fpkm)

######################################################################

## Take significant genes which are only present in real time points (03h to 48h)

head(g.diff);nrow(g.diff)
g.diff[["sig.col"]] <- NA
g.diff[["sig.col"]][g.diff[["p_value"]] <= 0.05] <- "red"
g.diff[["sig.col"]][g.diff[["p_value"]] > 0.05] <- "blue"
g.diff[["log2.fold_change."]][g.diff[["log2.fold_change."]] == Inf] <- 0
g.diff[["log2.fold_change."]][g.diff[["log2.fold_change."]] == -Inf] <- 0
g.diff[["log2.fold_change."]][g.diff[["log2.fold_change."]] == Inf] <- max(g.diff[["log2.fold_change."]])
g.diff[["log2.fold_change."]][g.diff[["log2.fold_change."]] == -Inf] <- min(g.diff[["log2.fold_change."]])

pdf(paste("P value vs Log2 fold change for all genes", "pdf", sep="."))
plot(g.diff[["log2.fold_change."]], -log10(g.diff[["p_value"]]), cex = 0.5, pch = 19, ylab = "-log10 p values", xlab = "Log fold change", col = g.diff[["sig.col"]])
legend("topright", legend = c("Up-regulated", "Down-regulated"), pch = 19, col = c("red", "blue"), cex = 0.8, inset = 0.01)
dev.off()

head(tmp.sig.g);nrow(tmp.sig.g)
gene.sig.tp.only <- subset(tmp.sig.g, sample_2 %in% grep("^siTCF7L2", tmp.sig.g[["sample_2"]], value = T))
gene.sig.tp.only <- subset(gene.sig.tp.only, significant %in% "yes")
head(gene.sig.tp.only); nrow(gene.sig.tp.only)

## Some exploratory plots

# pdf(paste("P value vs Log2 fold change for all significant genes", "pdf", sep="."))
# plot(gene.sig.tp.only[["log2.fold_c hange."]], -log10(gene.sig.tp.only[["p_value"]]), cex = 0.5, pch = 19, ylab = "p values", xlab = "Log fold change", col = gene.sig.tp.only[["t.col"]], xlim = c(min(g.diff[["log2.fold_change."]]), max(g.diff[["log2.fold_change."]])), ylim = c(min(-log10(g.diff[["p_value"]])), max(-log10(g.diff[["p_value"]]))))
# legend("topright", legend = unique(gene.sig.tp.only[["t.point"]]), pch = 19, col = unique(gene.sig.tp.only[["t.col"]]), cex = 0.8)
# dev.off()

pdf(paste("Log2 fold values for significant genes", "pdf", sep="."))
plot(gene.sig.tp.only[["log2.fold_change."]], pch = 19, cex = 0.6, col = gene.sig.tp.only[["t.col"]], ylab = "Log 2 Fold change", xlab = "Genes")
abline(a = 0, b = 0, col = "red")
legend("topright", legend = unique(gene.sig.tp.only[["t.point"]]), pch = 19, col = unique(gene.sig.tp.only[["t.col"]]), cex = 0.8)
dev.off()


## Collect the gene sets according to their time points. 

# # pdf(paste("Log2 fold values for all significant genes", "pdf", sep="."))

# i = 1
# tp.genes.sig <- list()
# plot(gene.sig.tp.only[["log2.fold_change."]], type = "n", xlim = c(0, length(unique(gene.sig.tp.only[["tracking_id"]]))), ylab = "Log2 fold change", xlab = "Genes")
# for (i in 1:length(orig.tp.ids)){
	# tmp <- subset(gene.sig.tp.only, t.point %in% orig.tp.ids[i])
	# tmp <- tmp[rev(order(tmp[["gene"]])),]
	# tp.genes.sig[[i]] <- tmp
	# points(tp.genes.sig[[i]][["log2.fold_change."]], pch = 19, col = tp.genes.sig[[i]][["t.col"]], cex = 0.5)
	
# }
# abline(a = 0, b = 0, col = "red")
# legend("topright", legend = unique(gene.sig.tp.only[["t.point"]]), pch = 19, col = unique(gene.sig.tp.only[["t.col"]]), cex = 0.8)
# dev.off()
# lapply(tp.genes.sig, head); lapply(tp.genes.sig, nrow)


## save significant gene list
setwd(outputfolder)

write.table(gene.sig.tp.only, file = "siginificant.gene.list_in.real.time.points.only.txt", sep = "", row.names = F, col.names = F, quote = F)

######################################################################

## Bar plot for top significant genes.
orig.tp.ids <- grep("[0-9]", tp.ids, value = T)
tmp.sig.genes <- unique(gene.sig.tp.only[["gene"]])
tmp.sig.genes <- tmp.sig.genes[! tmp.sig.genes %in% grep("[,.]", tmp.sig.genes, value =T)]
tmp.sig.genes <- tmp.sig.genes[! tmp.sig.genes %in% grep("Metazoa_SRP", tmp.sig.genes, value =T)]
sig.genes <- subset(gene.sig.tp.only, gene %in% tmp.sig.genes)

sig.genes.tp <- subset(rel.gene.fpkm, gene_short_name %in% tmp.sig.genes)
names(sig.genes.tp)[2] <- "gene"

top.sig.genes <- sig.genes.tp
head(top.sig.genes);nrow(top.sig.genes)

p.x <- 1:length(orig.tp.ids)

## Plot the fold change (relative to the scrambled) in all time points for top listed genes

setwd(plotfolder)
pdf(paste("Fold change for top listed genes", "pdf", sep="."))
i = 1
# par (mfrow = c(5,5))
for(i in 1:nrow(top.sig.genes)) {
	tmp.fpkm <- top.sig.genes[i,]
	siS.tmp.fpkm <- tmp.fpkm[names(tmp.fpkm) %in% grep("^m.rel.f.siS", names(tmp.fpkm), value = T)]
	siT.tmp.fpkm <- tmp.fpkm[names(tmp.fpkm) %in% grep("^m.rel.f.siT", names(tmp.fpkm), value = T)]
	if(all(siS.tmp.fpkm > -0.5)){		
		plot(t(siS.tmp.fpkm), type = "o", ylim = c(-5,5), pch = 19, col = "grey", ylab = "Fold change", xlab = "Time points", main = tmp.fpkm[["gene"]], xaxt = "n", cex = 0.7)
		s.m <- tmp.fpkm[grep("^m.rel.f.siS", names(tmp.fpkm), value = T)]
		s.sd <- tmp.fpkm[grep("^sd.rel.f.siS", names(tmp.fpkm), value = T)]
		sd.bar <- 0.02
		segments(p.x, t(s.m)-t(s.sd), p.x, t(s.m)+t(s.sd), col = "grey")
		segments(p.x-sd.bar, t(s.m)-t(s.sd), p.x+sd.bar, t(s.m)-t(s.sd))
		segments(p.x-sd.bar, t(s.m)+t(s.sd), p.x+sd.bar, t(s.m)+t(s.sd))		
		points(t(siT.tmp.fpkm), type = "o", pch = 19, col = "red", cex = 0.7)
		t.m <- tmp.fpkm[grep("^m.rel.f.siT", names(tmp.fpkm), value = T)]
		t.sd <- tmp.fpkm[grep("^sd.rel.f.siT", names(tmp.fpkm), value = T)]
		sd.bar <- 0.02
		segments(p.x, t(t.m)-t(t.sd), p.x, t(t.m)+t(t.sd), col = "red")
		segments(p.x-sd.bar, t(t.m)-t(t.sd), p.x+sd.bar, t(t.m)-t(t.sd))
		segments(p.x-sd.bar, t(t.m)+t(t.sd), p.x+sd.bar, t(t.m)+t(t.sd))		
		mtext(text = orig.tp.ids, side = 1, las = 1, line = 0.5, at = seq(1, length(orig.tp.ids), by = 1), cex = 0.6)
		legend("topright", legend = c("siScr", "siTCF"), pch = 19, col = c("grey", "red"), cex = 0.7, lty = 1)		
	}	
	
}
dev.off()
par(mfrow = c(1,1))


pdf(paste(paste("barplot.for.DEGs_log2.value", version.marker, sep="_"), "pdf", sep="."))
par(mar = c(4, 7, 2, 1))
par(mfrow = c(1, 1))
barplot(sig.genes[["log2.fold_change."]], col = sig.genes[["t.col"]], horiz = T, xlab = "Fold Change", space = 1, names.arg = sig.genes[["gene"]], las = 2, cex = 0.2, cex.lab = 0.6, cex.axis = 0.6, border = F, main = "Fold change")
legend("topright", legend = unique(gene.sig.tp.only[["t.point"]]), pch = 19, col = unique(gene.sig.tp.only[["t.col"]]), cex = 0.8)
dev.off()



## Heat map with significant gene list.
## combine both cuffdiff and cuffnorm 

tmp.heat <- sig.genes.tp[grep("^f.*", names(sig.genes.tp), value = T)]
tmp.heat <- cbind(sig.genes.tp[["gene"]], tmp.heat)
names(tmp.heat)[1] <- "gene"
tmp.heat <- tmp.heat[!duplicated(tmp.heat[["gene"]]),]
head(tmp.heat);nrow(tmp.heat)
tmp.heat1 <- tmp.heat[grep("^f.*", names(tmp.heat), value = T)]
names(tmp.heat1) <- gsub("f.", "", names(tmp.heat1))
rownames(tmp.heat1) <- tmp.heat[["gene"]]
tmp.heat1 <- as.matrix(tmp.heat1)
str(tmp.heat1)

qtile.cut <- quantile(tmp.heat1, 0.99)
tmp.heat1[tmp.heat1 > qtile.cut] <- qtile.cut

heat.breaks <- seq(0, max(tmp.heat1), max(tmp.heat1)/200)

setwd(plotfolder)
pdf(paste(paste("Heatmap for significant genes fpkm values_non scaled", sep="_"), "pdf", sep="."))

par(mar = c(4, 7, 2, 1))
par(mfrow = c(1,1))
mycol <-colorRampPalette(c("lightblue", "dodgerblue", "lightgreen", "palevioletred3", "darkorange", "red"), bias=3)(n=length(heat.breaks)-1)

heatmap.2(tmp.heat1, dendrogram="both", scale="none", col= mycol, breaks=heat.breaks, lhei = c(2,8), trace="none", margins=c(6, 13), na.rm = T, main="DEGs_scale.none", labRow = "", ColSideColors = heat.side.col)
dev.off()

pdf(paste(paste("Heatmap for significant genes fpkm values_scaled", sep="_"), "pdf", sep="."))
par(mar = c(4, 7, 2, 1))
par(mfrow = c(1,1))

heatmap.2(tmp.heat1, distfun = function(x) dist(x, method = "correlation"), dendrogram = "both", scale="row", col= bluered(100), Colv = F, Rowv = F, trace="none", margins=c(6, 13), lhei=c(2, 8), main="DEGs_Correlation_scaled", labRow = "", ColSideColors = heat.side.col)

par(mfrow = c(1,1))
dev.off()


## Heat map only for siTCF7L2, DEGs

tmp.heat2 <- sig.genes.tp[grep("^m.f.siT", names(sig.genes.tp), value = T)]
tmp.heat2 <- cbind(sig.genes.tp[["gene"]], tmp.heat2)
names(tmp.heat2)[1] <- "gene"
tmp.heat2 <- tmp.heat2[!duplicated(tmp.heat2[["gene"]]),]
head(tmp.heat2);nrow(tmp.heat2)
tmp.heat2 <- tmp.heat2[grep("^m.f.siT", names(tmp.heat2), value = T)]
names(tmp.heat2) <- gsub("m.f.siT_", "", names(tmp.heat2))
rownames(tmp.heat2) <- tmp.heat2[["gene"]]
# tmp.heat2 <- tmp.heat2[order(tmp.heat2[["03h"]]),]
tmp.heat2 <- as.matrix(tmp.heat2)
str(tmp.heat2)

qtile.cut <- quantile(tmp.heat2, 0.99)
tmp.heat2[tmp.heat2 > qtile.cut] <- qtile.cut

heat.breaks <- seq(0, max(tmp.heat2), max(tmp.heat2)/200)

pdf(paste(paste("Heatmap for significant genes_siTCF7L2.only fpkm values_scaled", sep="_"), "pdf", sep="."))
par(mar = c(4, 7, 2, 1))
par(mfrow = c(1,1))

heatmap.2(tmp.heat2, distfun = function(x) dist(x, method = "correlation"), dendrogram = "none", scale="row", col= bluered(100), Colv = F, Rowv = F, trace="none", margins=c(6, 13), lhei=c(2, 8), main="Significant.genes_Correlation_scaled_siTCF7L2", labRow = "")

dev.off()


## Name genes-of-interest to plot below
## NOTE: Could be further developed by reading these from (even multiple) cuffdiff output files,
## and by supplying the plots with qvalues at all time points on top of siT boxes, whether the
## gene was a DEG at that time point or not.



## Select the most responsing genes to prior treatment

imp.sig.genes <- tmp.sig.g[!duplicated(tmp.sig.g[["gene"]]), ]
imp.sig.genes <- imp.sig.genes[order(imp.sig.genes[["q_value"]]),]

tmp.g.all <- g.all[order(-g.all[["m.f.siS_03h"]]),]

# genes.oi <- imp.sig.genes[["gene"]]

genes.oi <- tmp.sig.genes

setwd(plotfolder)
pdf(paste(paste("boxplots.for.genes.of.interest", version.marker, sep="_"), "pdf", sep="."))
i=1
for(i in 1:length(genes.oi)){
	
	cat("Processing: ", genes.oi[i], "\n")
	tmp <- subset(g.all, gene_short_name == genes.oi[i])
	tmp <- tmp[,c("gene_short_name", grep("^f\\.", names(tmp), value=T))]
	if(nrow(tmp) == 1){
		rownames(tmp) <- paste(tmp[,c("gene_short_name")], collapse=".")
		tmp.plot <- as.data.frame(t(tmp[,grep("^f\\.", names(tmp))]))
		tmp.plot[["ini.grp"]] <- gsub("h_[[:digit:]]$", "", gsub("^f\\.", "", rownames(tmp.plot)))
		tmp.plot[["grp"]] <- paste(gsub("^.*_", "", tmp.plot[["ini.grp"]]), "h_", gsub("_.*$", "", tmp.plot[["ini.grp"]]), sep="")
		tmp.plot[["grp"]][tmp.plot[["grp"]] %in% c("24h_mock", "24h_untreated")] <- c("mock", "untreated")
		tmp.plot[["grp"]] <- as.factor(tmp.plot[["grp"]])
		tmp.plot[["ini.grp"]] <- NULL
		print(tmp.plot)
		plot(tmp.plot[,c(2,1)], main=colnames(tmp.plot)[1], ylab="FPKM", col=c(rep(c("dodgerblue", "darkorange"),8), "yellow", "green"), las =2, xlab = "")
		points(as.numeric(tmp.plot[,2]), tmp.plot[,1])
		legend("topright", inset=c(0.01,0.01), fill=c("dodgerblue", "darkorange", "yellow", "green"), legend=c("siS", "siT", "mock", "untreated"))
		
	} else {
		warning("Duplicates present for this gene!")
 	}
	
	cat("\n\n\n")
	
}
dev.off()

