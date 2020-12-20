source("Simulation.R")
set.seed(2020)

## senario 1
set.seed(2020)
K = 5 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

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
mean_var_ratio = 10
corr = 0.5
prob_k = c(1,2,2,3,3)
#### rst<-run(X,0.4,.5,3.5,2,AS,A0,D,K,N, epsilon = 10^(-4))
#### rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=150)




## senario 2
set.seed(2020)
K = 8 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 8
small_w_mean = 1
small_tau_w = 8
p.0 = 0.15
q.0 = 0.3
p.neg1 = 0.1
q.neg1 = 0.15
mean_var_ratio = 5
corr = 0.2
#prob_k = c(1,2,2,3,3)
prob_k = c(3,1,2,3,3,1,2,1)
### rst <- run(X,1,10,30,10,AS,A0,D,K,N, epsilon = 10^(-3))
### rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=150)

## senario 3
set.seed(2020)
K = 8 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 8
small_w_mean = 1
small_tau_w = 8
p.0 = 0.15
q.0 = 0.3
p.neg1 = 0.1
q.neg1 = 0.15
mean_var_ratio = 3
corr = 0.3
#prob_k = c(1,2,2,3,3)
prob_k = c(3,1,2,3,3,1,2,1)
### rst<-run(X,1,15,40,10,AS,A0,D,K,N, epsilon = 10^(-3))
### rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=150)

## senario 4
set.seed(2020)
K = 8 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 7
small_w_mean = 1
small_tau_w = 7
p.0 = 0.2
q.0 = 0.4
p.neg1 = 0.1
q.neg1 = 0.2
mean_var_ratio = 3
corr = 0.2
#prob_k = c(1,2,2,3,3)
prob_k = c(3,1,2,3,3,1,2,1)
### rst<-run(X,3,45,60,10,AS,A0,D,K,N, epsilon = 10^(-3))


## senario 5
set.seed(2020)
K = 8 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 5
small_w_mean = 1
small_tau_w = 5
p.0 = 0.2
q.0 = 0.4
p.neg1 = 0.1
q.neg1 = 0.2
mean_var_ratio = 2
corr = 0.2
#prob_k = c(1,2,2,3,3)
prob_k = c(3,1,2,3,3,1,2,1)
### rst<-run(X,3,45,60,10,AS,A0,D,K,N, epsilon = 10^(-3))





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
seur <- FindClusters(seur, resolution = 1)
seur <- RunUMAP(seur, dims = 1:5)
DimPlot(seur, reduction='umap', group.by = 'label')

DimPlot(seur, reduction='umap',group.by='seurat_clusters')