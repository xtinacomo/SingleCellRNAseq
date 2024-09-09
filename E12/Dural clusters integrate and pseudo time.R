library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(cowplot)


Idents(fibroblasts.final2) <- "Cluster"
Idents(mm10.merged.fibroblasts.E14) <- "Cluster"

fibroblasts.final2@meta.data$Cluster <- paste("E12.5", fibroblasts.final2@meta.data$Cluster, sep = "_")
mm10.merged.fibroblasts.E14@meta.data$Cluster <- paste("E14", mm10.merged.fibroblasts.E14@meta.data$Cluster, sep = "_")

#Normalize
norm.fibroblastsE12.5 <- NormalizeData(fibroblasts.final2)
norm.fibroblastsE14 <- NormalizeData(mm10.merged.fibroblasts.E14)

#Find variable features
E12.5fibro <-  FindVariableFeatures(norm.fibroblastsE12.5, selection.method = "vst",
                                    nfeatures = 2000, verbose = FALSE)
E14fibro <-  FindVariableFeatures(norm.fibroblastsE14, selection.method = "vst",
                                  nfeatures = 2000, verbose = FALSE)

#Find Anchors
reference_list <- c(E12.5fibro, E14fibro)
fibroblast_anchors <- FindIntegrationAnchors(object.list = reference_list, dims = 1:30)

#Integrate
fibroblast_integrated <- IntegrateData(anchorset = fibroblast_anchors, dims = 1:30)
DefaultAssay(fibroblast_integrated) <- "integrated"
saveRDS(fibroblast_integrated, file = "fibroblast_integrated.RDS")


#Scale integrated seurat object
fibroblast_integrated <- readRDS("fibroblast_integrated.RDS")
DefaultAssay(fibroblast_integrated) <- "integrated"
fibroblast_integrated <- ScaleData(fibroblast_integrated, verbose = FALSE)

#Run PCA
fibroblast_integrated <- RunPCA(fibroblast_integrated, npcs = 30, verbose = FALSE)

#Run UMAP and plot 
#May need to adjust dims parameter 
#grouped by Cluster since at this point it is equivalent to a sample ID, I just want to see if my samples overlap/cluster together
fibroblast_integrated <- RunUMAP(fibroblast_integrated, reduction = "pca", dims = 1:30)
DimPlot(fibroblast_integrated, reduction = "umap", group.by = "seurat_clusters")
DimPlot(fibroblast_integrated, reduction = "umap", pt.size = 0.1, split.by = "orig.ident", ncol = 2)
DimPlot(fibroblast_integrated, reduction = "umap", group.by = "Cluster")

#Find neighbors
Idents(fibroblast_integrated) <- "Sample"
fibroblast_integrated <- FindNeighbors(fibroblast_integrated, dims = 1:15)
fibroblast_integrated <- FindClusters(fibroblast_integrated, resolution = .3)
#5 clusters 

#should noodle more with the dims parameter 
fibroblast_integrated <- RunUMAP(fibroblast_integrated, reduction = "pca", dims = 1:30)
DimPlot(fibroblast_integrated, reduction = "umap", group.by = "seurat_clusters")
DimPlot(fibroblast_integrated, reduction = "umap", label = TRUE, pt.size = 0.2)
DimPlot(fibroblast_integrated, reduction = "umap", group.by = "Cluster")
DimPlot(fibroblast_integrated, reduction = "umap", pt.size = 0.5, split.by = "orig.ident", ncol = 2)

Fig2_pia_markers <- c("Ngfr","Bgn","S100a6","Rdh10","Cxcl12","Mfap2")
Fig3_arachnoid_markers <- c("Slc6a13","Crabp2","Angptl2","Aebp1","Wnt6","Aldh1a2")
Fig4_dura_markers <- c("Fxyd5","Foxp1","Nov","Smoc2","Ndrg1","Kctd12")
# 
FeaturePlot(fibroblast_integrated, features = "Vtn", max.cutoff = 4)
#clusters 2 and 1 are arachnoid, cluster 5 is dura and pia is cluster 0 

#further subcluster - cluster 5 for dural cells ___________________________________________________________________________-
Idents(fibroblast_integrated) <- "seurat_clusters"
fibroblast_integrated_8 <- subset(fibroblast_integrated,idents = c('5'))
DimPlot(fibroblast_integrated_8, reduction = 'umap',label = TRUE, pt.size = 0.2) + xlim(-10,5) + ylim(-8,4)

fibroblast_integrated_8 <- FindVariableFeatures(fibroblast_integrated_8, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(fibroblast_integrated_8)
VariableFeaturePlot(fibroblast_integrated_8) + ylim(0,10)

fibroblast_integrated_8 <- ScaleData(fibroblast_integrated_8, verbose = F)
fibroblast_integrated_8 <- RunPCA(fibroblast_integrated_8, npcs = 40, verbose = F)
ElbowPlot(fibroblast_integrated_8, ndims = 40)
fibroblast_integrated_8 <- FindNeighbors(fibroblast_integrated_8, reduction = "pca", dims = 1:30)
fibroblast_integrated_8 <- FindClusters(fibroblast_integrated_8, resolution = .2, algorithm = 1)
fibroblast_integrated_8 <- RunUMAP(fibroblast_integrated_8, reduction = "pca", dims = 1:40)
DimPlot(fibroblast_integrated_8, reduction = "umap", label = TRUE, pt.size = 1)
DimPlot(fibroblast_integrated_8, reduction = "umap", pt.size = 0.1, split.by = "orig.ident", ncol = 2)

dc0_1 <- FindMarkers(fibroblast_integrated_8, ident.1 = "0", ident.2 = "1")
write.csv(dc0_1, 'dc0_1.csv')
FeaturePlot(fibroblast_integrated_8, features = c("Birc5"))


#Slingshot
fibroblast_integrated_8.sce <- as.SingleCellExperiment(fibroblast_integrated_8)
fibroblast_integrated_8_slingshot <- slingshot(fibroblast_integrated_8.sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')
summary(fibroblast_integrated_8_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(fibroblast_integrated_8_slingshot$slingPseudotime_1, breaks=100)]

#Lineage structure estimation
plot(reducedDims(fibroblast_integrated_8_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_8_slingshot), lwd=2, col='black')

plot(reducedDims(fibroblast_integrated_8_slingshot)$PCA, col = brewer.pal(9,'Set1')[fibroblast_integrated_8_slingshot$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_8_slingshot), lwd=2, type = 'lineages', col = 'black')

# fit negative binomial GAM; Identify temporally expressed genes 

require(gam)
t <- fibroblast_integrated_8_slingshot$slingPseudotime_1
# for time, only look at the 100 most variable genes
Y <- log1p(assays(fibroblast_integrated_8_slingshot)$logcounts)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]

Y <- Y[var100,]
Y
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

#Cdh1 Cells 
summary(fibroblast_integrated_8_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(fibroblast_integrated_8_slingshot$slingPseudotime_1, breaks=100)]
plot(reducedDims(fibroblast_integrated_8_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_8_slingshot), lwd=2, col='black')
SlingshotDataSet(fibroblast_integrated_8_slingshot)

#Heatmap
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(fibroblast_integrated_8_slingshot)$norm[topgenes, order(t, na.last = NA)]
heatclus <- fibroblast_integrated_8_slingshot$seurat_clusters[order(t, na.last = NA)]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])


#further subcluster - cluster  for arachnoid cells ___________________________________________________________________________-
Idents(fibroblast_integrated) <- "seurat_clusters"
fibroblast_integrated_2_1 <- subset(fibroblast_integrated,idents = c('2'))
DimPlot(fibroblast_integrated_2_1, reduction = 'umap',label = TRUE, pt.size = 0.2) + xlim(-10,5) + ylim(-8,4)

fibroblast_integrated_2_1 <- FindVariableFeatures(fibroblast_integrated_2_1, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(fibroblast_integrated_2_1)
VariableFeaturePlot(fibroblast_integrated_2_1) + ylim(0,10)

fibroblast_integrated_2_1 <- ScaleData(fibroblast_integrated_2_1, verbose = F)
fibroblast_integrated_2_1 <- RunPCA(fibroblast_integrated_2_1, npcs = 40, verbose = F)
ElbowPlot(fibroblast_integrated_2_1, ndims = 40)
fibroblast_integrated_2_1 <- FindNeighbors(fibroblast_integrated_2_1, reduction = "pca", dims = 1:30)
fibroblast_integrated_2_1 <- FindClusters(fibroblast_integrated_2_1, resolution = .4, algorithm = 1)
fibroblast_integrated_2_1 <- RunUMAP(fibroblast_integrated_2_1, reduction = "pca", dims = 1:40)
DimPlot(fibroblast_integrated_2_1, reduction = "umap", label = TRUE, pt.size = 1)
DimPlot(fibroblast_integrated_2_1, reduction = "umap", pt.size = 0.1, split.by = "orig.ident", ncol = 2)

arach_pia <- FindMarkers(fibroblast_integrated_2_1, ident.1 = "5", ident.2 = "3")
write.csv(arach_pia, 'arach_pia_DEG.csv')



#Slingshot
fibroblast_integrated_5.sce <- as.SingleCellExperiment(fibroblast_integrated_2_1)
fibroblast_integrated_5_slingshot <- slingshot(fibroblast_integrated_5.sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')
summary(fibroblast_integrated_5_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(fibroblast_integrated_5_slingshot$slingPseudotime_1, breaks=100)]

#Lineage structure estimation
plot(reducedDims(fibroblast_integrated_5_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_5_slingshot), lwd=2, col='black')

plot(reducedDims(fibroblast_integrated_5_slingshot)$PCA, col = brewer.pal(9,'Set1')[fibroblast_integrated_5_slingshot$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_5_slingshot), lwd=2, type = 'lineages', col = 'black')

# fit negative binomial GAM; Identify temporally expressed genes 

require(gam)
t <- fibroblast_integrated_5_slingshot$slingPseudotime_1
# for time, only look at the 100 most variable genes
Y <- log1p(assays(fibroblast_integrated_5_slingshot)$logcounts)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]

Y <- Y[var100,]
Y
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

#Cdh1 Cells 
summary(fibroblast_integrated_5_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(fibroblast_integrated_5_slingshot$slingPseudotime_1, breaks=100)]
plot(reducedDims(fibroblast_integrated_5_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_5_slingshot), lwd=2, col='black')
SlingshotDataSet(fibroblast_integrated_5_slingshot)

#Heatmap
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(fibroblast_integrated_2_1_slingshot)$logcounts[topgenes, order(t, na.last = NA)]
heatclus <- fibroblast_integrated_2_1_slingshot$seurat_clusters[order(t, na.last = NA)]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])


apc_52 <- FindMarkers(fibroblast_integrated_2_1, ident.1 = "5", ident.2 = "2")
write.csv(apc_52, 'apc52_DEG.csv')


#___________________________________________________________________________________________________________________
#pia cluster 0 

Idents(fibroblast_integrated) <- "seurat_clusters"
fibroblast_integrated_0 <- subset(fibroblast_integrated,idents = c('0'))
DimPlot(fibroblast_integrated_0, reduction = 'umap',label = TRUE, pt.size = 0.2) + xlim(-10,5) + ylim(-8,4)

fibroblast_integrated_0 <- FindVariableFeatures(fibroblast_integrated_0, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(fibroblast_integrated_0)
VariableFeaturePlot(fibroblast_integrated_0) + ylim(0,10)

fibroblast_integrated_0 <- ScaleData(fibroblast_integrated_0, verbose = F)
fibroblast_integrated_0 <- RunPCA(fibroblast_integrated_0, npcs = 40, verbose = F)
ElbowPlot(fibroblast_integrated_0, ndims = 40)
fibroblast_integrated_0 <- FindNeighbors(fibroblast_integrated_0, reduction = "pca", dims = 1:30)
fibroblast_integrated_0 <- FindClusters(fibroblast_integrated_0, resolution = .4, algorithm = 1)
fibroblast_integrated_0 <- RunUMAP(fibroblast_integrated_0, reduction = "pca", dims = 1:40)
DimPlot(fibroblast_integrated_0, reduction = "umap", label = TRUE, pt.size = 1)
DimPlot(fibroblast_integrated_0, reduction = "umap", pt.size = 0.1, split.by = "orig.ident", ncol = 2)

pia_only_DEG <- FindMarkers(fibroblast_integrated_0, ident.1 = "3", ident.2 = "0")
write.csv(pia_only_DEG, 'pia_only_DEG.csv')



#Slingshot
fibroblast_integrated_0.sce <- as.SingleCellExperiment(fibroblast_integrated_0)
fibroblast_integrated_0_slingshot <- slingshot(fibroblast_integrated_0.sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')
summary(fibroblast_integrated_0_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(fibroblast_integrated_0_slingshot$slingPseudotime_1, breaks=100)]

#Lineage structure estimation
plot(reducedDims(fibroblast_integrated_0_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_0_slingshot), lwd=2, col='black')

plot(reducedDims(fibroblast_integrated_0_slingshot)$PCA, col = brewer.pal(9,'Set1')[fibroblast_integrated_0_slingshot$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_0_slingshot), lwd=2, type = 'lineages', col = 'black')

# fit negative binomial GAM; Identify temporally expressed genes 

require(gam)
t <- fibroblast_integrated_0_slingshot$slingPseudotime_1
# for time, only look at the 100 most variable genes
Y <- log1p(assays(fibroblast_integrated_0_slingshot)$logcounts)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]

Y <- Y[var100,]
Y
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

#Cdh1 Cells 
summary(fibroblast_integrated_0_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(fibroblast_integrated_0_slingshot$slingPseudotime_1, breaks=100)]
plot(reducedDims(fibroblast_integrated_0_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_0_slingshot), lwd=2, col='black')
SlingshotDataSet(fibroblast_integrated_0_slingshot)

#Heatmap
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(fibroblast_integrated_2_1_slingshot)$logcounts[topgenes, order(t, na.last = NA)]
heatclus <- fibroblast_integrated_2_1_slingshot$seurat_clusters[order(t, na.last = NA)]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])


apc_52 <- FindMarkers(fibroblast_integrated_2_1, ident.1 = "5", ident.2 = "2")
write.csv(apc_52, 'apc52_DEG.csv')


