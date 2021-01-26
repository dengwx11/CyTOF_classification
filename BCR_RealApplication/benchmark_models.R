options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

sample <- as.character(args[1])
ct.cnt <- as.integer(args[2])
print(c(sample,ct.cnt))

set.seed(2020)
source('Opt/run_para.R')
source('Opt/preparation.R')
suppressPackageStartupMessages(library(SummarizedExperiment))
library(Seurat)

## BCR
S <- readRDS('write/matrix_S/bcr_w_7celltype_tomo32_s.rds')
# S <-readRDS('write/matrix_S/hvgs_bcr_w_tomo35_s.rds')
A0 <- read.csv('./write/matrix_A0/A0_BCR.csv')
W <- readRDS('write/matrix_W/bcr_w.rds')
A <- read.table('write/matrix_A/bcr_w_7celltype_hvg_tomo_s_32_a.csv')
#X <- as.matrix(read.table('write/matrix_X/X_BCR.txt'))
BCR <- readRDS('/gpfs/loomis/project/zhao/bz234//Data/BCR/BCR_celltype_annotated.rds')
X <- assay(BCR[, colData(BCR)$marker_class == "type"])
X <- t(X)
cofactor <- 5
X <- asinh(X / cofactor)
colnames(X) <- c(1:ncol(X))
metadata.BCR<-rowData(BCR)
metadata.BCR <- data.frame(metadata.BCR@listData)
rownames(metadata.BCR) <- colnames(X)

## clean the input
rownames(A0) <- A0[,1]
A0 <- A0[,-1]
rownames(S) <- rownames(A)
colnames(S) <- c("memory.B.cells","naïve.B.cells","CD4.T.cells","DC","monocytes","NK.cells","CD8.T.cells")

S <- S[,colnames(A0)]
colnames(W) <- c('CD4.T.cells','CD8.T.cells','DC','memory.B.cells','monocytes','naïve.B.cells','NK.cells')


W <- W[,colnames(A0)]

A0<-as.matrix(A0)
if(ct.cnt == 5){
S <- S[,-c(1,5)]
A0 <- A0[,-c(1,5)]
W <- W[,-c(1,5)]
}else if(ct.cnt == 6){
S <- S[,-c(5)]
A0 <- A0[,-c(5)]
W <- W[,-c(5)]
}
AS <- t(A)%*%S

if(ct.cnt == 5){
    idx.DC = which(metadata.BCR$population_id == 'DC')
    idx.mB = which(metadata.BCR$population_id == 'memory B cells')
    metadata.BCR <- metadata.BCR[-c(idx.mB,idx.DC),]
    X <- X[,-c(idx.mB,idx.DC)]
}else if(ct.cnt== 6){
    idx = which(metadata.BCR$population_id == 'DC')
    metadata.BCR <- metadata.BCR[-idx,]
    X <- X[,-idx]
}


metadata.BCR$population_id <- droplevels(metadata.BCR$population_id)


#sample.id.list <- unique(metadata.BCR$sample_id)
#W.true <- list()
#for(i in c(1:length(sample.id.list))){


    X.sample <- X[,metadata.BCR$sample_id %in% sample]
    metadata.sample <- metadata.BCR[metadata.BCR$sample_id %in% sample,]
    seur <- CreateSeuratObject( X.sample)
    seur@meta.data$label <- metadata.sample$population_id
    seur <- NormalizeData(seur) 
    seur@assays$RNA@data = seur@assays$RNA@counts
    seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
    seur <- ScaleData(seur)
    seur <- RunPCA(seur,verbose = TRUE,features = rownames(seur))
    seur <- FindNeighbors(seur, dims = 1:6)
    seur <- FindClusters(seur, resolution = 0.3)
    seur <- RunUMAP(seur, dims = 1:5)

    ## ground truth
    Idents(seur) <- "label"
    W <- AverageExpression(seur,slot="counts")$RNA
if(ct.cnt == 5){
    colnames(W) <- c('CD4.T.cells','CD8.T.cells','monocytes','naïve.B.cells','NK.cells')
}else if(ct.cnt== 6){
    colnames(W) <- c('CD4.T.cells','CD8.T.cells','memory.B.cells','monocytes','naïve.B.cells','NK.cells')
}else{
    colnames(W) <- c('CD4.T.cells','CD8.T.cells','DC','memory.B.cells','monocytes','naïve.B.cells','NK.cells')
}
    W <- as.matrix(W[,colnames(A0)])
 #   W.true[[i]] <- W

    D = nrow(X.sample)
    K = ncol(AS)
    N = ncol(X.sample)

print(dim(X.sample))
print(dim(AS))
print(dim(A0))
print(c(D,K,N))

rst.para<-runOptimalPara(X.sample,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=3,lambda1.on=T,lambda2.on=T)
rst<-run(X.sample,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
rst.para.woas<-runOptimalPara(X.sample,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=3,lambda1.on=F,lambda2.on=T)
rst.woas<-run(X.sample,rst.para.woas$para$lambda1,rst.para.woas$para$lambda2,rst.para.woas$para$mu,rst.para.woas$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
rst.para.woa0<-runOptimalPara(X.sample,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=3,lambda1.on=T,lambda2.on=F)
rst.woa0<-run(X.sample,rst.para.woa0$para$lambda1,rst.para.woa0$para$lambda2,rst.para.woa0$para$mu,rst.para.woa0$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
rst.para.woa0as<-runOptimalPara(X.sample,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=3,lambda1.on=F,lambda2.on=F)           
rst.woa0as<-run(X.sample,rst.para.woa0as$para$lambda1,rst.para.woa0as$para$lambda1,rst.para.woa0as$para$mu,rst.para.woa0as$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)

truth_ct = seur$label
H_est <- data.frame(rst$H)
H_est <- cbind(matrix(colnames(A0),ncol=1),H_est)
H_est_woas <- data.frame(rst.woas$H)
H_est_woas <- cbind(matrix(colnames(A0),ncol=1),H_est_woas)
H_est_woa0 <- data.frame(rst.woa0$H)
H_est_woa0 <- cbind(matrix(colnames(A0),ncol=1),H_est_woa0)
H_est_woa0as <- data.frame(rst.woa0as$H)
H_est_woa0as <- cbind(matrix(colnames(A0),ncol=1),H_est_woa0as)
print("1")
if(ct.cnt == 5){
    H_est[,1] <- c("naïve B cells","CD4 T-cells","CD8 T-cells","monocytes","NK cells")
    H_est_woas[,1] <- c("naïve B cells","CD4 T-cells","CD8 T-cells","monocytes","NK cells")
    H_est_woa0[,1] <- c("naïve B cells","CD4 T-cells","CD8 T-cells","monocytes","NK cells")
    H_est_woa0as[,1] <- c("naïve B cells","CD4 T-cells","CD8 T-cells","monocytes","NK cells")
}else if(ct.cnt== 6){
#	print(H_est[,1:4])
    H_est[,1] <- c("memory B cells","naïve B cells","CD4 T-cells","CD8 T-cells","monocytes","NK cells")
    H_est_woas[,1] <- c("memory B cells","naïve B cells","CD4 T-cells","CD8 T-cells","monocytes","NK cells")
    H_est_woa0[,1] <- c("memory B cells","naïve B cells","CD4 T-cells","CD8 T-cells","monocytes","NK cells")
    H_est_woa0as[,1] <- c("memory B cells","naïve B cells","CD4 T-cells","CD8 T-cells","monocytes","NK cells")
}else{
    H_est[,1] <- c("memory B cells","naïve B cells","CD4 T-cells","CD8 T-cells","DC","monocytes","NK cells")
    H_est_woas[,1] <- c("memory B cells","naïve B cells","CD4 T-cells","CD8 T-cells","DC","monocytes","NK cells")
    H_est_woa0[,1] <- c("memory B cells","naïve B cells","CD4 T-cells","CD8 T-cells","DC","monocytes","NK cells")
    H_est_woa0as[,1] <- c("memory B cells","naïve B cells","CD4 T-cells","CD8 T-cells","DC","monocytes","NK cells")
}
print("2")
celltype_pred <- apply(rst$H, 2, predict_realdata, H = H_est)
celltype_pred_woas <- apply(rst.woas$H, 2, predict_realdata, H = H_est_woas)
celltype_pred_woa0 <- apply(rst.woa0$H, 2, predict_realdata, H = H_est_woa0)
celltype_pred_woa0as <- apply(rst.woa0as$H, 2, predict_realdata, H = H_est_woa0as)
print("3")
print(celltype_pred)
print(celltype_pred_woas)

seur$pred = celltype_pred
seur$pred_woas = celltype_pred_woas
seur$pred_woa0 = celltype_pred_woa0
seur$pred_woa0as = celltype_pred_woa0as

rst$sample.id = sample
rst.woas$sample.id = sample
rst.woa0$sample.id = sample
rst.woa0as$sample.id = sample

seur$sample.id = sample
rst$H = H_est
rst.woas$H = H_est_woas
rst.woa0$H = H_est_woa0
rst.woa0as$H = H_est_woa0as

saveRDS(rst,paste0('/home/bz234/project/Results/CyTOF/BCR_Cytof/rst_models/',ct.cnt,'celltype/rst_',sample,".rds"))
saveRDS(rst.woas,paste0('/home/bz234/project/Results/CyTOF/BCR_Cytof/rst_models/',ct.cnt,'celltype/rst_woas_',sample,".rds"))
saveRDS(rst.woa0,paste0('/home/bz234/project/Results/CyTOF/BCR_Cytof/rst_models/',ct.cnt,'celltype/rst_woa0_',sample,".rds"))
saveRDS(rst.woa0as,paste0('/home/bz234/project/Results/CyTOF/BCR_Cytof/rst_models/',ct.cnt,'celltype/rst_woa0as_',sample,".rds"))

saveRDS(seur,paste0('/home/bz234/project/Results/CyTOF/BCR_Cytof/rst_models/',ct.cnt,'celltype/seur_',sample,".rds"))
