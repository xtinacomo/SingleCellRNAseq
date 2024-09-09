#further subcluster
Idents(fibroblasts.E12.5.further) <- "seurat_clusters"
fibroblasts.E12.5.further2 <- subset(fibroblasts.E12.5.further,idents = c('1', '2'))
nonfibroblasts.E12.5.further2 <- subset(fibroblasts.E12.5, idents = c('0'))
DimPlot(fibroblasts.E12.5.further2, reduction = 'umap',label = TRUE, pt.size = 0.2) + xlim(-10,5) + ylim(-8,4)

fibroblasts.E12.5.further2 <- FindVariableFeatures(fibroblasts.E12.5.further2, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(fibroblasts.E12.5.further2)
VariableFeaturePlot(fibroblasts.E12.5.further2) + ylim(0,10)

fibroblasts.E12.5.further2 <- ScaleData(fibroblasts.E12.5.further2, verbose = F)
fibroblasts.E12.5.further2 <- RunPCA(fibroblasts.E12.5.further2, npcs = 40, verbose = F)
ElbowPlot(fibroblasts.E12.5.further2, ndims = 40)
fibroblasts.E12.5.further2 <- FindNeighbors(fibroblasts.E12.5.further2, reduction = "pca", dims = 1:30)
fibroblasts.E12.5.further2 <- FindClusters(fibroblasts.E12.5.further2, resolution = 0.4, algorithm = 1)
fibroblasts.E12.5.further2 <- RunUMAP(fibroblasts.E12.5.further2, reduction = "pca", dims = 1:40)
DimPlot(fibroblasts.E12.5.further2, reduction = "umap", label = TRUE, pt.size = 1)
DimPlot(fibroblasts.E12.5.further2, reduction = "umap", pt.size = 0.1, split.by = "orig.ident", ncol = 2)

FeaturePlot(fibroblasts.E12.5.further2, features = Fig2_pia_markers, max.cutoff = 5)
FeaturePlot(fibroblasts.E12.5.further2, features = Fig3_arachnoid_markers, max.cutoff = 5)
FeaturePlot(fibroblasts.E12.5.further2, features = Fig4_dura_markers, max.cutoff = 5)
FeaturePlot(fibroblasts.E12.5.further2, features = c("Tbx18"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Cxcl12"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Wnt6"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Rdh10"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Crabp2"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Aldh1a2"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Aebp1"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Cdh1"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Cldn11"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Smoc2"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Ndrg1"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Slc6a13"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("S100a6"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Ngfr"))
FeaturePlot(fibroblasts.E12.5.further2, features = c("Angptl2"))



