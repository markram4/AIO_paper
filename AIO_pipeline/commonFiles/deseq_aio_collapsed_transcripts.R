library(methods)
library("DESeq2")
library("gplots")
library("RColorBrewer")
library(tidyverse)
library(ggrepel)
library("ggsci")
library(ggpubr)

CountData <- read.table("/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/deseq/merged_input_deseq.txt",header = T) %>% column_to_rownames("Gene")

labels <- colnames(CountData) 
condition <- sub("\\..*", "", labels)

smpInfo = data.frame(labels, condition)


# Build DESeq data set
dds <- DESeqDataSetFromMatrix(countData=CountData, colData=smpInfo, design= ~condition)
dds$condition <- relevel(dds$condition, ref="wt_red")
dds <- DESeq(dds)


## wt_red vs dcl_red
outPrefix1 <- "/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/deseq/DESeq_AIO_wt_vs_dcl_red.rnd1.collapsed"

normTable <- counts(dds, normalized = TRUE)
write.table(normTable, file=paste(outPrefix1,".normalizedRawValues.txt",sep=''), quote=F, sep="\t", row.names=T)
rawTable <- counts(dds, normalized = F)
write.table(rawTable, file=paste(outPrefix1,".RawValues.txt",sep=''), quote=F, sep="\t", row.names=T)

res.dds <- results(dds, contrast = c("condition", "dcl_red","wt_red"))
res.dds <- subset(res.dds, !is.na(res.dds$padj))
resSig.dds <- res.dds[ abs(res.dds$padj) < .05, ]
resSig20.dds <- res.dds[ abs(res.dds$padj) < .20, ]
resAll.dds <- res.dds[ abs(res.dds$padj) < 2, ]

write.table(resSig.dds, file=paste(outPrefix1,".FDR05.txt",sep=''), quote=F, sep="\t", row.names=T)
write.table(resSig20.dds, file=paste(outPrefix1,".FDR20.txt",sep=''), quote=F, sep="\t", row.names=T)
write.table(resAll.dds, file=paste(outPrefix1,".all.txt",sep=''), quote=F, sep="\t", row.names=T)

pdf(paste(outPrefix1,"pdf",sep="."), title="DESeq plots")
plotMA(res.dds,main = "dcl_vs_wt",ylim=c(-7,8),alpha = 0.05)
hist(res.dds$pvalue,xlab="p value",main = "dcl_vs_wt", breaks=100, col="lightblue")
dev.off()


#---------------
## wt_red vs wt.fullRed
outPrefix2 <- "/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/deseq/DESeq_AIO_red_vs_fullRed.rnd1.collapsed"

res.dds <- results(dds, contrast = c("condition", "wt_fullRed","wt_red"))
res.dds <- subset(res.dds, !is.na(res.dds$padj))
resSig.dds <- res.dds[ abs(res.dds$padj) < .05, ]
resSig20.dds <- res.dds[ abs(res.dds$padj) < .20, ]
resAll.dds <- res.dds[ abs(res.dds$padj) < 2, ]

write.table(resSig.dds, file=paste(outPrefix2,".FDR05.txt",sep=''), quote=F, sep="\t", row.names=T)
write.table(resSig20.dds, file=paste(outPrefix2,".FDR20.txt",sep=''), quote=F, sep="\t", row.names=T)
write.table(resAll.dds, file=paste(outPrefix2,".all.txt",sep=''), quote=F, sep="\t", row.names=T)
write.table(normTable, file=paste(outPrefix2,".normalizedRawValues.txt",sep=''), quote=F, sep="\t", row.names=T)

pdf(paste(outPrefix2,"pdf",sep="."), title="DESeq plots")
plotMA(res.dds,main = "fullRed_vs_red",ylim=c(-8,8),alpha = 0.05)
hist(res.dds$pvalue,xlab="p value",main = "fullRed_vs_red", breaks=100, col="lightblue")
dev.off()

## wt_fullRed vs dcl_red
outPrefix3 <- "/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/deseq/DESeq_AIO_wt_fullRed_vs_dcl_red.rnd1.collapsed"

res.dds <- results(dds, contrast = c("condition", "dcl_red","wt_fullRed"))
res.dds <- subset(res.dds, !is.na(res.dds$padj))
resSig.dds <- res.dds[ abs(res.dds$padj) < .05, ]
resSig20.dds <- res.dds[ abs(res.dds$padj) < .20, ]
resAll.dds <- res.dds[ abs(res.dds$padj) < 2, ]

write.table(resSig.dds, file=paste(outPrefix3,".FDR05.txt",sep=''), quote=F, sep="\t", row.names=T)
write.table(resSig20.dds, file=paste(outPrefix3,".FDR20.txt",sep=''), quote=F, sep="\t", row.names=T)
write.table(resAll.dds, file=paste(outPrefix3,".all.txt",sep=''), quote=F, sep="\t", row.names=T)
write.table(normTable, file=paste(outPrefix3,".normalizedRawValues.txt",sep=''), quote=F, sep="\t", row.names=T)


pdf(paste(outPrefix3,"pdf",sep="."), title="DESeq plots")
plotMA(res.dds,main = "dcl_vs_wt_fullRed",ylim=c(-8,8),alpha = 0.05)
hist(res.dds$pvalue,xlab="p value",main = "dcl_vs_wt_fullRed", breaks=100, col="lightblue")
dev.off()

##--------------
## Clustering
vsd <- varianceStabilizingTransformation(dds)
vsdMat <- assay(vsd)
colnames(vsdMat) <- with(colData(dds), paste(labels, sep=" : "))

hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsVS = dist(t(vsdMat))
matVS = as.matrix(distsVS)
rownames(matVS) <- colnames(matVS) <- with(colData(dds), paste(labels, sep=" : "))

pdf("/cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round1/data/6_minimap/2_map_to_capture/output/mapped_to_targets/bedFiles/deseq/DESeq_AIO.clustering_plots.pdf", title="DESeq plots")
heatmap.2(matVS, trace="none",col=rev(hmcol),margin=c(13,13))
plotPCA(vsd, intgroup="condition")

correlation <- cor(vsdMat)
corrstr <- apply(correlation, c(1,2), function(x) sprintf("%.3g", x))
colors2 <- colorRampPalette(c("white","darkblue"))(150)
heatmap.2(correlation, col=colors2, scale="none", trace="none", main="correlation between samples", cellnote=corrstr, notecex=0.4, notecol="black")
dev.off()

#----
# Plot using ggplot2
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Add sample labels to the PCA data
pcaData$Sample <- rownames(pcaData)

pcaPlot.gg <-ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = Sample,shape=condition)) +
  geom_point(size = 2) + scale_color_manual(values=c("#B03060","#853061","#C97795"))+
  #geom_text_repel(seed=123,size=2,force=0.5) +
  xlab(paste0("PC1: ", round(100 * attr(pcaData, "percentVar")[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaData, "percentVar")[2], 1), "% variance")) +
  theme_bw() + themes+
  theme(legend.position = "top")
