# /gpfs/loomis/project/zhao/bz234/Results/Tomo

set.seed(2020)
source('run_para.R')

## Input
## inhouse
A <- as.matrix(read.csv('./write/matrix_A/inhouse_w_hvg_covid_s_35_a.csv',sep=' '))
S <- readRDS('./write/matrix_S/hvgs_inhouse_w_covid_s.rds')
X <- as.matrix(read.table('./write/matrix_X/X_inhouse.txt'))
A0 <- read.csv('./write/matrix_A0/A0_inhouse.csv')
# truth_ct = read.table('./write/true_label/celltype_inhouse.txt', sep='\t',header=T)
metadata = read.table('./write/true_label/metadata_inhouse.txt', sep='\t',header=T)
truth_ct = metadata$celltype
colnames(A0)=c("Gene","CD4+T cell","CD8+T cell", "naïve B cells","NK cell", "classical monocytes", "non-classical monocytes")
rownames(A0) = A0$Gene
A0 <- as.matrix(A0[-7,-1])
S <- as.matrix(S[,colnames(A0)])
AS <- t(A) %*% S
# X <- X[-7,truth_ct$cellID]
X <- X[-7,]

## BCR 5 celltypes
S <- readRDS('write/matrix_S/hvgs_bcr_w_tomo32_s.rds')
# S <-readRDS('write/matrix_S/hvgs_bcr_w_tomo35_s.rds')
A0 <- read.csv('./write/matrix_A0/A0_BCR.csv')
W <- readRDS('write/matrix_W/bcr_w.rds')
A <- read.table('write/matrix_A/bcr_w_hvg_tomo_s_32_a.csv')
X <- as.matrix(read.table('write/matrix_X/X_BCR.txt'))

## get louvain cluster
rownames(X) <- c(1:D)
colnames(X) <- c(1:N)
seur <- CreateSeuratObject(X)
seur@meta.data$label <- truth_ct
seur <- NormalizeData(seur) 
seur@assays$RNA@data = seur@assays$RNA@counts
seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
seur <- ScaleData(seur)
seur <- RunPCA(seur,verbose = TRUE,features = rownames(seur))
seur <- FindNeighbors(seur, dims = 1:6)
seur <- FindClusters(seur, resolution = 0.3)
seur <- RunUMAP(seur, dims = 1:5)
DimPlot(seur, reduction='umap', group.by = 'label')

DimPlot(seur, reduction='umap',group.by='seurat_clusters')

## ground truth
Idents(seur) <- "label"
W <- AverageExpression(seur,slot="counts")$RNA
W <- as.matrix(W[,colnames(A0)])


## Other para
# lambda1 = 0.5
# lambda2 = 0.5
# mu = 1
D = nrow(X)
K = ncol(AS)
N = ncol(X)


## 
rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=250,depth=4)
rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=0)
rst<-run(X,1,5,5,5,AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=1000)


##
H_est <- data.frame(rst$H)
H_est <- cbind(matrix(colnames(A0),ncol=1),H_est)


get_truth <- function(true_ct){
    return(which(as.character(H_est[,1]) == true_ct))
}
infer_max <- function(truth,h){
    return(1*(h[truth]==max(h)))
}
predict_realdata<-function(h){
    loc = which(h == max(h))
    return(colnames(A0)[loc])
}

celltype_pred <- apply(rst$H, 2, predict_realdata)
ARI.louvain <- adjustedRandIndex(celltype_pred, seur$seurat_clusters)
print(ARI.louvain)
ARI.truth <- adjustedRandIndex(celltype_pred, seur$label)
print(ARI.truth)
plot(as.vector(W),as.vector(rst$W))
cor(as.vector(W),as.vector(rst$W))

#truth = sapply(as.character(truth_ct), get_truth)
cnt_max = 0
for(i in 1:length(truth_ct)){
    cnt_max = cnt_max + 1*(truth_ct[i]==celltype_pred[i])
}  
cnt_max
print(cnt_max/N)

seur$pred = celltype_pred
DimPlot(seur, reduction='umap', group.by = 'pred')
DimPlot(seur, reduction='umap', group.by = 'label')