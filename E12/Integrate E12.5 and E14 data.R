library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(cowplot)


Idents(fibroblasts.E12.5.further) <- "Cluster"
Idents(mm10.merged.fibroblasts.E14) <- "Cluster"

fibroblasts.E12.5.further@meta.data$Cluster <- paste("E12.5", fibroblasts.E12.5.further@meta.data$Cluster, sep = "_")
mm10.merged.fibroblasts.E14@meta.data$Cluster <- paste("E14", mm10.merged.fibroblasts.E14@meta.data$Cluster, sep = "_")

#Normalize
norm.fibroblastsE12.5.further <- NormalizeData(fibroblasts.E12.5.further)
norm.fibroblastsE14 <- NormalizeData(mm10.merged.fibroblasts.E14)

#Find variable features
E12.5fibro.further <-  FindVariableFeatures(norm.fibroblastsE12.5.further, selection.method = "vst",
                                   nfeatures = 2000, verbose = FALSE)
E14fibro <-  FindVariableFeatures(norm.fibroblastsE14, selection.method = "vst",
                                  nfeatures = 2000, verbose = FALSE)

#Find Anchors
reference_list <- c(E12.5fibro.further, E14fibro)
fibroblast_anchors <- FindIntegrationAnchors(object.list = reference_list, dims = 1:30)

#Integrate
fibroblast_integrated <- IntegrateData(anchorset = fibroblast_anchors, dims = 1:30)
DefaultAssay(fibroblast_integrated) <- "integrated"
saveRDS(fibroblast_integrated, file = "fibroblast_integrated.RDS")


#Scale integrated seurat object
fibroblast_integrated <- readRDS("fibroblast_integrated.RDS")
DefaultAssay(fibroblast_integrated) <- "integrated"
fibroblast_integrated <- ScaleData(fibroblast_integrated, verbose = FALSE)

#Count all of the "Fibroblast-like" Cells
CountCells <- table(fibroblast_integrated$manualCluster)
CountCells

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
fibroblast_integrated <- FindClusters(fibroblast_integrated, resolution = 0.5)
#5 clusters 

#should noodle more with the dims parameter 
fibroblast_integrated <- RunUMAP(fibroblast_integrated, reduction = "pca", dims = 1:30)
DimPlot(fibroblast_integrated, reduction = "umap", group.by = "seurat_clusters")
DimPlot(fibroblast_integrated, reduction = "umap", label = TRUE, pt.size = 0.2)
DimPlot(fibroblast_integrated, reduction = "umap", group.by = "Cluster")
DimPlot(fibroblast_integrated, reduction = "umap", pt.size = 0.1, split.by = "orig.ident", ncol = 2)

cluster0.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '0') #compares 7 to all other clusters 
cluster1.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '1') #compares 7 to all other clusters 
cluster2.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '2') #compares 7 to all other clusters 
cluster3.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '3') #compares 7 to all other clusters 
cluster4.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '4') #compares 7 to all other clusters 
cluster5.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '5') #compares 7 to all other clusters 
cluster6.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '6') #compares 7 to all other clusters 
cluster7.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '7') #compares 7 to all other clusters 
cluster8.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '8') #compares 7 to all other clusters 
cluster9.markers <- FindMarkers(fibroblast_integrated, only.pos = T, ident.1 = '9') #compares 7 to all other clusters 


#Feature Plots:
ClusterMarkers <- c("Camk2n1","Postn","Vegfa","Ntm","S100a6",
                    "Pltp","Islr","Itih5","Wnt4","Enpp2","Aldh1a2","Ngfr",
                    "Crabp2","Ogn","Cdkn1c","Alcam","Ptgds","Wnt6",
                    "Fxyd5","Dkk2","Mgp","1500015O10Rik","Anxa3",
                    "Hist1h1b","Hist1h2ae","2810417H13Rik","Top2a","Hmgb2",
                    "Cenpa","Cenpf","Tubb4b","Birc5","Ube2c")

FeaturePlot(fibroblast_integrated, features = ClusterMarkers[1:5], max.cutoff = 5)
FeaturePlot(fibroblast_integrated, features = ClusterMarkers[6:12], max.cutoff = 5)
FeaturePlot(fibroblast_integrated, features = ClusterMarkers[13:18])
FeaturePlot(fibroblast_integrated, features = ClusterMarkers[19:23], max.cutoff = 5)
NewMarkers <- c("S100a6", "Pdgfra", "Aldh1a1", "Camk2n1")

pia_markers <- c("Ngfr","Bgn","S100a6","Rdh10","Cxcl12","Mfap2")
arachnoid_markers <- c("Slc6a13","Crabp2","Angptl2","Aebp1","Wnt6","Aldh1a2")
dura_markers <- c("Fxyd5","Foxp1","Nov","Smoc2","Ndrg1","Kctd12")

FeaturePlot(fibroblast_integrated, features = NewMarkers, max.cutoff = 5)
FeaturePlot(fibroblast_integrated, features = arachnoid_markers, max.cutoff = 5)
FeaturePlot(fibroblast_integrated, features = dura_markers, max.cutoff = 5)
FeaturePlot(fibroblast_integrated, features = c("Tbx18"))
FeaturePlot(fibroblast_integrated, features = c("Rdh10"))
FeaturePlot(fibroblast_integrated, features = c("Kcnj8"))
#Remove clusters that are not fibroblasts and plot what is 
Idents(fibroblast_integrated) <- "seurat_clusters"
fibroblasts <- subset(fibroblast_integrated, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8"))

Idents(fibroblasts) <- "Sample"
num <- table(Idents(fibroblasts))
num

#To fetch the old data cluster IDs
p1 <- DimPlot(mm10.merged.fibroblasts.E14, reduction = "umap", group.by = "seurat_clusters")
p2 <- DimPlot(fibroblast_integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

FeaturePlot(mm10.merged.fibroblasts.E12.5, features = c("Rdh10"))
FeaturePlot(fibroblast_integrated, features = c("Ptgfr"))
FeaturePlot(fibroblast_integrated, features = c("Crabp2"))
FeaturePlot(fibroblast_integrated, features = c("Rdh10"))
FeaturePlot(fibroblast_integrated, features = c("Fxyd5"))
FeaturePlot(fibroblast_integrated, features = c("S100a6"))
FeaturePlot(fibroblast_integrated, features = c("Sox9"))
FeaturePlot(fibroblast_integrated, features = c("Tmem119"))
FeaturePlot(fibroblast_integrated, features = c("Crym"))
FeaturePlot(fibroblast_integrated, features = c("Cldn5"))
FeaturePlot(fibroblast_integrated, features = c("Tbx18"))
FeaturePlot(fibroblast_integrated, features = c("Cd68"))
FeaturePlot(mm10.merged.fibroblasts.E12.5, features = c("Des"))
FeaturePlot(fibroblast_integrated, features = c("Cldn1"))
FeaturePlot(fibroblast_integrated, features = c("Cldn11"))
FeaturePlot(fibroblast_integrated, features = c("Zeb1"))
FeaturePlot(fibroblast_integrated, features = c("Tsc2"))
FeaturePlot(fibroblast_integrated, features = c("Pdgfra"))
FeaturePlot(fibroblast_integrated, features = c("Ngfr"))
FeaturePlot(fibroblast_integrated, features = c("Aldh1a2"))


#further subcluster - cluster 3
Idents(fibroblast_integrated) <- "seurat_clusters"
fibroblast_integrated_3 <- subset(fibroblast_integrated,idents = c('3'))
DimPlot(fibroblast_integrated_3, reduction = 'umap',label = TRUE, pt.size = 0.2) + xlim(-10,5) + ylim(-8,4)

fibroblast_integrated_3 <- FindVariableFeatures(fibroblast_integrated_3, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(fibroblast_integrated_3)
VariableFeaturePlot(fibroblast_integrated_3) + ylim(0,10)

fibroblast_integrated_3 <- ScaleData(fibroblast_integrated_3, verbose = F)
fibroblast_integrated_3 <- RunPCA(fibroblast_integrated_3, npcs = 40, verbose = F)
ElbowPlot(fibroblast_integrated_3, ndims = 40)
fibroblast_integrated_3 <- FindNeighbors(fibroblast_integrated_3, reduction = "pca", dims = 1:30)
fibroblast_integrated_3 <- FindClusters(fibroblast_integrated_3, resolution = 0.4, algorithm = 1)
fibroblast_integrated_3 <- RunUMAP(fibroblast_integrated_3, reduction = "pca", dims = 1:40)
DimPlot(fibroblast_integrated_3, reduction = "umap", label = TRUE, pt.size = 1)
DimPlot(fibroblast_integrated_3, reduction = "umap", pt.size = 0.1, split.by = "orig.ident", ncol = 2)

FeaturePlot(fibroblast_integrated_3, features = c("Rdh10"))
FeaturePlot(fibroblast_integrated_3, features = c("Ptgfr"))
FeaturePlot(fibroblast_integrated_3, features = c("Crabp2"))
FeaturePlot(fibroblast_integrated_3, features = c("Rdh10"))
FeaturePlot(fibroblast_integrated_3, features = c("Fxyd5"))
FeaturePlot(fibroblast_integrated_3, features = c("S100a6"))
FeaturePlot(fibroblast_integrated_3, features = c("Sox9"))
FeaturePlot(fibroblast_integrated_3, features = c("Tmem119"))
FeaturePlot(fibroblast_integrated_3, features = c("Crym"))
FeaturePlot(fibroblast_integrated_3, features = c("Cldn5"))
FeaturePlot(fibroblast_integrated_3, features = c("Tbx18"))
FeaturePlot(fibroblast_integrated, features = c("Cd68"))
FeaturePlot(fibroblast_integrated_3, features = c("Des"))
FeaturePlot(fibroblast_integrated_3, features = c("Cldn1"))
FeaturePlot(fibroblast_integrated_3, features = c("Cldn11"))
FeaturePlot(fibroblast_integrated_3, features = c("Zeb1"))
FeaturePlot(fibroblast_integrated_3, features = c("Tsc2"))
FeaturePlot(fibroblast_integrated_3, features = c("Pdgfra"))
FeaturePlot(fibroblast_integrated_3, features = c("Ngfr"))
FeaturePlot(fibroblast_integrated_3, features = c("Aldh1a2"))
FeaturePlot(fibroblast_integrated_3, features = c("Cldn11"))
FeaturePlot(fibroblast_integrated_3, features = c("Smoc2"))
FeaturePlot(fibroblast_integrated_3, features = c("Ndrg1"))
FeaturePlot(fibroblast_integrated_3, features = c("Slc6a13"))
FeaturePlot(fibroblast_integrated_3, features = c("S100a6"))
FeaturePlot(fibroblast_integrated_3, features = c("Ngfr"))
FeaturePlot(fibroblast_integrated_3, features = c("Atp1a2"))

cluster2.markers2 <- FindMarkers(fibroblast_integrated_3, only.pos = T, ident.1 = '2') #compares 7 to all other clusters 
Wnt6pos <- subset(fibroblast_integrated, subset = Wnt6 >1)
table(Idents(Wnt6pos))
