#Slingshot data on E12.5 data - Cluster starting point 3

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
library(Seurat)
install.packages("umap")
library(umap)
#load data
read.csv('E12.5_script.R')

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")
#Run Slingshot #Tutorial 
#sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')

#____________________________________________________________--

E12.5_subclusters.sce <- as.SingleCellExperiment(fibroblasts.E12.5.further)

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
E12.5_subclusters.sce <- SingleCellExperiment(assays = List(counts = counts))

# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(E12.5_subclusters.sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
E12.5_subclusters.sce <- E12.5_subclusters.sce[geneFilter, ]


#Normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(E12.5_subclusters.sce)$norm <- FQnorm(assays(E12.5_subclusters.sce)$counts)

#Reduce dimensionality
pca <- prcomp(t(log1p(assays(E12.5_subclusters.sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

rd2 <- umap(t(log1p(assays(E12.5_subclusters.sce)$norm)))
colnames(rd2) <- c('UMAP1', "UMAP2")

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)


reducedDims(E12.5_subclusters.sce) <- SimpleList(PCA = rd1, UMAP = rd2)
cl1 <- Mclust(rd1)$classification
colData(E12.5_subclusters.sce)$GMM <- cl1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(E12.5_subclusters.sce)$kmeans <- cl2
plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

#Run Slingshot #My data
E12.5_subclusters_slingshot <- slingshot(E12.5_subclusters.sce, clusterLabels = "GMM", reducedDim = 'PCA')
E12.5_subclusters_slingshot$slingPseudotime_1

#E12.5_subclusters_slingshot
summary(E12.5_subclusters_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(E12.5_subclusters_slingshot$slingPseudotime_1, breaks=100)]
plot(reducedDims(E12.5_subclusters_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(E12.5_subclusters_slingshot), lwd=2, col='black')
SlingshotDataSet(E12.5_subclusters_slingshot)

plot(reducedDims(E12.5_subclusters_slingshot)$PCA, col = brewer.pal(9,'Set1')[E12.5_subclusters_slingshot$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(E12.5_subclusters_slingshot), lwd=2, type = 'lineages', col = 'black')

## fit negative binomial GAM
E12.5_subclusters_slingshot <- fitGAM(E12.5_subclusters_slingshot)

# test for dynamic expression
ATres <- associationTest(E12.5_subclusters_slingshot)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(E12.5_subclusters_slingshot$slingPseudotime_1, na.last = NA)
heatdata <- assays(E12.5_subclusters_slingshot)$counts[topgenes, pst.ord]
heatclus <- E12.5_subclusters_slingshot$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

topgenes



E12.5_subclusters_slingshot <- slingshot(E12.5_subclusters_slingshot, clusterLabels = 'GMM', reducedDim = 'PCA',
                  approx_points = 5)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(E12.5_subclusters_slingshot$slingPseudotime_1, breaks=100)]

plot(reducedDims(sim5)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim5), lwd=2, col='black')


#We can also see how the lineage structure was intially estimated by the
#cluster-based minimum spanning tree by using the `type` argument.

#{r plot_curve_2}
plot(reducedDims(E12.5_subclusters_slingshot)$PCA, col = brewer.pal(9,'Set1')[E12.5_subclusters_slingshot$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(E12.5_subclusters_slingshot), lwd=2, col='black')
SlingshotDataSet(E12.5_subclusters_slingshot)

library(tradeSeq)

# fit negative binomial GAM
E12.5_subclusters_slingshot <- fitGAM(E12.5_subclusters_slingshot)

# test for dynamic expression
ATres <- associationTest(E12.5_subclusters_slingshot)

# test for dynamic expression
require(gam)
t <- E12.5_subclusters_slingshot$slingPseudotime_1
# for time, only look at the 100 most variable genes
head(E12.5_subclusters_slingshot)
Y <- log1p(assays(E12.5_subclusters_slingshot)$logcounts)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})


#Heatmap of top genes 
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(E12.5_subclusters.sce)$norm[topgenes, order(t, na.last = NA)]
heatclus <- E12.5_subclusters.sce$seurat_clusters[order(t, na.last = NA)]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])
topgenes

lin1 <- getLineages(rd, cl, start.clus = '1')
lin1
