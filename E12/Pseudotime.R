

#Load Libraries from tutorial 

knitr::opts_chunk$set(fig.align="center", cache=TRUE,error=FALSE, #stop on error
                      fig.width=5, fig.height=5, autodep=TRUE, out.width="600px", out.height="600px",
                      results="markup", echo=TRUE, eval=TRUE)
#knitr::opts_knit$set(stop_on_error = 2L) #really make it stop
#knitr::dep_auto()
options(getClass.msg=FALSE)
graphics:::par(pch = 16, las = 1)
set.seed(12345) ## for reproducibility
library(SingleCellExperiment)
library(slingshot, quietly = TRUE)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(mclust)
library(gam)


# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")



read.csv('E12.5_script.R')
read.csv('Seurat.merge.filter_nFeatures_mt_FibroblastOnly.R')
read.csv('E12.5 and E14 merged.R')
fibroblast_integrated.sce <- as.SingleCellExperiment(fibroblast_integrated)


# generate synthetic count data representing a single lineage
means <- rbind(
  # non-DE genes
  matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
         ncol = 300, byrow = TRUE),
  # early deactivation
  matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # late deactivation
  matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # early activation
  matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # late activation
  matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # transient
  matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
         ncol = 300, byrow = TRUE)
)
counts <- apply(means,2,function(cell_means){
  total <- rnbinom(1, mu = 7500, size = 4)
  rmultinom(1, total, cell_means)
})
rownames(counts) <- paste0('G',1:750)
colnames(counts) <- paste0('c',1:300)
fibroblast_integrated.sce <- SingleCellExperiment(assays = List(counts = counts))

# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(fibroblast_integrated.sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
fibroblast_integrated.sce <- fibroblast_integrated.sce[geneFilter, ]

#Normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(fibroblast_integrated.sce)$norm <- FQnorm(assays(fibroblast_integrated.sce)$counts)

#Dimensionality Reduction
pca <- prcomp(t(log1p(assays(fibroblast_integrated.sce)$counts)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)


rd2 <- umap(t(log1p(assays(fibroblast_integrated.sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(fibroblast_integrated.sce) <- SimpleList(PCA = rd1, UMAP = rd2)

#Mcluster
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(fibroblast_integrated.sce)$seurat_clusters <- cl1

#K cluster
library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 6)$cluster
colData(fibroblast_integrated.sce)$kmeans <- cl2
plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)
####__________________________________________________________________________________________________

fibroblast_integrated.sce <- as.SingleCellExperiment(fibroblast_integrated)

#Slingshot
fibroblast_integrated_slingshot <- slingshot(fibroblast_integrated.sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')
summary(fibroblast_integrated_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(fibroblast_integrated_slingshot$slingPseudotime_1, breaks=100)]

#Lineage structure estimation
plot(reducedDims(fibroblast_integrated_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_slingshot), lwd=2, col='black')

plot(reducedDims(fibroblast_integrated_slingshot)$PCA, col = brewer.pal(9,'Set1')[fibroblast_integrated_slingshot$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_slingshot), lwd=2, type = 'lineages', col = 'black')

# fit negative binomial GAM; Identify temporally expressed genes 

require(gam)
t <- fibroblast_integrated$slingPseudotime_1
# for time, only look at the 100 most variable genes
Y <- log1p(assays(fibroblast_integrated_slingshot)$logcounts)
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
summary(fibroblast_integrated_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(fibroblast_integrated_slingshot$slingPseudotime_1, breaks=100)]
plot(reducedDims(fibroblast_integrated_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(fibroblast_integrated_slingshot), lwd=2, col='black')
SlingshotDataSet(fibroblast_integrated_slingshot)

#Heatmap
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(fibroblast_integrated_slingshot)$norm[topgenes, order(t, na.last = NA)]
heatclus <- fibroblast_integrated_slingshot$seurat_clusters[order(t, na.last = NA)]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])


#Get lineages
lin1 <- getLineages(rd1, cl1, start.clus = '0')
lin1

plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 2, pch = 16)
lines(lin1, lwd = 3, col = 'black')

#start 6 end 1
lin5.1 <- getLineages(rd1, cl1, start.clus= '3', end.clus = '0')
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin5.1, lwd = 3, col = 'black', show.constraints = TRUE)

#start 6 end 2
lin6.5 <- getLineages(rd1, cl1, start.clus= '6', end.clus = '5')
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin6.5, lwd = 3, col = 'black', show.constraints = TRUE)

#start 5 end 3
lin5.3 <- getLineages(rd1, cl1, start.clus= '5', end.clus = '3')
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin5.3, lwd = 3, col = 'black', show.constraints = TRUE)

#start 5 end 4
lin5.4 <- getLineages(rd1, cl1, start.clus= '5', end.clus = '4')
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin5.4, lwd = 3, col = 'black', show.constraints = TRUE)

#start 5 end 6
lin5.6 <- getLineages(rd1, cl1, start.clus= '5', end.clus = '6')
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin5.6, lwd = 3, col = 'black', show.constraints = TRUE)
#Curve
crv5 <- getCurves(lin5)
crv5

plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(crv5, lwd = 3, col = 'black')
