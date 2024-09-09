E12.5_clusters.sce <- as.SingleCellExperiment(fibroblasts.final)

# # generate synthetic count data representing a single lineage
# means <- rbind(
#   # non-DE genes
#   matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
#          ncol = 300, byrow = TRUE),
#   # early deactivation
#   matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
#   # late deactivation
#   matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
#   # early activation
#   matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
#   # late activation
#   matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
#   # transient
#   matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
#          ncol = 300, byrow = TRUE)
# )
# counts <- apply(means,2,function(cell_means){
#   total <- rnbinom(1, mu = 7500, size = 4)
#   rmultinom(1, total, cell_means)
# })
# rownames(counts) <- paste0('G',1:750)
# colnames(counts) <- paste0('c',1:300)
# E12.5_clusters.sce <- SingleCellExperiment(assays = List(counts = counts))
# 
# E12.5_clusters.sce <- as.SingleCellExperiment(fibroblasts.final)
# 
# 
# # filter genes down to potential cell-type markers
# # at least M (15) reads in at least N (15) cells
# geneFilter <- apply(assays(E12.5_clusters.sce)$counts,1,function(x){
#   sum(x >= 3) >= 10
# })
# E12.5_clusters.sce <- E12.5_clusters.sce[geneFilter, ]
# 
# #Normalization
# FQnorm <- function(counts){
#   rk <- apply(counts,2,rank,ties.method='min')
#   counts.sort <- apply(counts,2,sort)
#   refdist <- apply(counts.sort,1,median)
#   norm <- apply(rk,2,function(r){ refdist[r] })
#   rownames(norm) <- rownames(counts)
#   return(norm)
# }
# assays(E12.5_clusters.sce)$norm <- FQnorm(assays(E12.5_clusters.sce)$counts)
# 
# #Dimensionality Reduction
# pca <- prcomp(t(log1p(assays(E12.5_clusters.sce)$norm)), scale. = FALSE)
# rd1 <- pca$x[,1:2]
# 
# plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
# 
# 
# rd2 <- umap(t(log1p(assays(E12.5_clusters.sce)$counts)))
# colnames(rd2) <- c('UMAP1', 'UMAP2')
# 
# #plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
# #reducedDims(E12.5_clusters.sce) <- SimpleList(PCA = rd1, UMAP = rd2)
# 
# #Mcluster
# library(mclust, quietly = TRUE)
# cl1 <- Mclust(rd1)$classification
# colData(E12.5_clusters.sce)$seurat_clusters <- cl1
# 
# #K cluster
# library(RColorBrewer)
# plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
# 
# #cl2 <- kmeans(rd1, centers = 6)$cluster
# c#olData(E12.5_clusters.sce)$kmeans <- cl2
# #plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)


#Run Slingshot #My data
timepoint_slingshot <- slingshot(E12.5_clusters.sce, clusterLabels = "seurat_clusters", reducedDim = 'PCA')
timepoint_slingshot$slingPseudotime_1

#Cdh1 Cells 
summary(timepoint_slingshot$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(timepoint_slingshot$slingPseudotime_1, breaks=100)]
plot(reducedDims(timepoint_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(timepoint_slingshot), lwd=2, col='black')
SlingshotDataSet(timepoint_slingshot)

#Lineage structure estimation
plot(reducedDims(timepoint_slingshot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(timepoint_slingshot), lwd=2, col='black')

plot(reducedDims(timepoint_slingshot)$PCA, col = brewer.pal(9,'Set1')[timepoint_slingshot$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(timepoint_slingshot), lwd=2, type = 'lineages', col = 'black')

# fit negative binomial GAM; Identify temporally expressed genes 

require(gam)
t <- timepoint_slingshot$slingPseudotime_1
# for time, only look at the 100 most variable genes
Y <- log1p(assays(timepoint_slingshot)$counts)
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

#Heatmap
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(timepoint_slingshot)$counts[topgenes, order(t, na.last = NA)]
heatclus <- timepoint_slingshot$seurat_clusters[order(t, na.last = NA)]
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

#Get lineages
lin1 <- getLineages(rd1, cl1, start.clus = '0')
lin1

plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 3, col = 'black')

#start 6 end 1
lin5.1 <- getLineages(rd2, cl2, start.clus= '1', end.clus = '4')
plot(rd2, col = brewer.pal(9,"Set1")[cl2], asp = 1, pch = 16)
lines(lin5.1, lwd = 3, col = 'black', show.constraints = TRUE)

s#start 6 end 2
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

