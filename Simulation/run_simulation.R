source("Simulation/Simulation.R")
set.seed(2020)
library(umap)






## generate dataset: X, AS, A0, W(ground truth), true.label
gamma <- simulate_gamma(iteration, K, D)
A0 <- gamma$gamma
W <- simulate_w(D,K,A0, p.0 = p.0, q.0=q.0, p.neg1= p.neg1, q.neg1 = q.neg1, 
                big_w_mean=big_w_mean, big_tau_w = big_tau_w, 
                small_w_mean=small_w_mean, small_tau_w = small_tau_w)
AS <- simulate_AS(D,K,W,corr=corr)
plot(as.vector(W),as.vector(AS), xlab = 'w', ylab = 'AS')
cor(as.vector(W),as.vector(AS))
cor(as.vector(W),as.vector(AS),method='spearman')
label.output <- simulate_label(D,N,K,prob_k)
X <- simulate_x(D,N,W,label.output$label, mean_var_ratio = mean_var_ratio)
true.H <- label.output$H


## umap visualization
X.umap = umap(t(X))
plot(X.umap$layout,col=label.output$label)




# seur <- NormalizeData(seur) 
# seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
# seur <- ScaleData(seur, slot='counts')
# seur <- RunPCA(seur,slot='counts', verbose = TRUE,features = rownames(seur))




# Idents(seur) <- 'label'
# RidgePlot(seur, slot = 'counts',features = rownames(seur))


## get louvain cluster
rownames(X) <- c(1:D)
colnames(X) <- c(1:N)
seur <- CreateSeuratObject(X)
seur@meta.data$label <- label.output$label[1,]
seur <- NormalizeData(seur) 
seur@assays$RNA@data = seur@assays$RNA@counts
seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
seur <- ScaleData(seur)
seur <- RunPCA(seur,verbose = TRUE,features = rownames(seur))
seur <- FindNeighbors(seur, dims = 1:6)
seur <- FindClusters(seur)
## Senario 1,2,3
# seur <- FindClusters(seur, resolution = .6)
## Senario 4
# seur <- FindClusters(seur, resolution = 1)
## Senario 5
# seur <- FindClusters(seur, resolution = 1.2)

seur <- RunUMAP(seur, dims = 1:5)
DimPlot(seur, reduction='umap', group.by = 'label')

DimPlot(seur, reduction='umap',group.by='seurat_clusters')

