source("Simulation.R")
set.seed(2020)


## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 10
small_w_mean = 0.5
small_tau_w = 10
p.0 = 0.1
q.0 = 0.2
p.neg1 = 0.05
q.neg1 = 0.1
mean_var_ratio = 0.5
corr = 0.5
prob_k = c(1,2,2,3,3)
#prob_k = c(3,1,2,3,3,1,2,1)




gamma <- simulate_gamma(iteration, K, D)
A0 <- gamma$gamma
W <- simulate_w(D,K,A0, p.0 = p.0, q.0=q.0, p.neg1= p.neg1, q.neg1 = q.neg1)
AS <- simulate_AS(D,K,W,corr=corr)
plot(as.vector(W),as.vector(AS), xlab = 'w', ylab = 'AS')
cor(as.vector(W),as.vector(AS))
cor(as.vector(W),as.vector(AS),method='spearman')
label.output <- simulate_label(D,N,K, prob_k = prob_k)
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

sum=0
for(i in 1:10) sum = sum + i^2