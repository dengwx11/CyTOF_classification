source("Simulation.R")
set.seed(2020)



gamma <- simulate_gamma(10, K, D)
A0 <- gamma$gamma
W <- simulate_w(D,K,A0, p.0 = 0.15, q.0=.3, p.neg1= 0.1, q.neg1 = 0.15)
AS <- simulate_AS(D,K,W,corr=0.3)
plot(as.vector(W),as.vector(AS), xlab = 'w', ylab = 'AS')
cor(as.vector(W),as.vector(AS))
cor(as.vector(W),as.vector(AS),method='spearman')
label.output <- simulate_label(D,N,K, prob_k = c(3,1,2,3,3,1,2))
X <- simulate_x(D,N,W,label.output$label, mean_var_ratio = 3)
true.H <- label.output$H

rownames(X) <- c(1:D)
colnames(X) <- c(1:N)
seur <- CreateSeuratObject(X)
seur@meta.data$label <- label.output$label[1,]

# seur <- NormalizeData(seur) 
# seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
seur <- ScaleData(seur, slot='counts')
seur <- RunPCA(seur,slot='counts', verbose = T,features = rownames(seur))

seur <- FindNeighbors(seur, dims = 1:5)
seur <- FindClusters(seur, resolution = 0.8)
seur <- RunUMAP(seur, dims = 1:5)
DimPlot(seur, reduction='umap', group.by = 'label')

Idents(seur) <- 'label'
RidgePlot(seur, slot = 'counts',features = rownames(seur))
