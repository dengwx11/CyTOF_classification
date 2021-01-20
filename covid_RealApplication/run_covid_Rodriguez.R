options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

expr.dir <- as.character(args[1])
A0.dir <- as.character(args[2])
sample.id <- as.character(args[3])


set.seed(2021)
setwd('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/CyTOF_classification')
source('run_para.R')
source('preparation.R')
library(Seurat)

# load input
expr <- readRDS(expr.dir)
A0 <- read.csv(A0.dir,header =T, row.names=1)
AS <- matrix(rep(0, nrow(A0)*ncol(A0)),nrow = nrow(A0), ncol=ncol(A0))
lineage_channels = rownames(A0)
X <- expr[,lineage_channels]
X <- t(as.matrix(X))
A0 <- as.matrix(A0)

    seur <- CreateSeuratObject( X)
    seur <- NormalizeData(seur) 
    seur@assays$RNA@data = seur@assays$RNA@counts
    seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
    seur <- ScaleData(seur)
    seur <- RunPCA(seur,verbose = TRUE,features = rownames(seur))
    seur <- FindNeighbors(seur, dims = 1:6)
    seur <- FindClusters(seur, resolution = 1)
    seur <- RunUMAP(seur, dims = 1:5)


    D = nrow(X)
    K = ncol(AS)
    N = ncol(X)

## run
rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=3,lambda1.on=F,lambda2.on=T)
rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)

H_est <- data.frame(rst$H)
H_est <- cbind(matrix(colnames(A0),ncol=1),H_est)
celltype_pred <- apply(rst$H, 2, predict_realdata)

saveRDS(rst,paste0('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/covid_Rodriguez/only_A0/rst_',sample.id,".rds"))
saveRDS(seur,paste0('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/covid_Rodriguez/only_A0/seur_',sample.id,".rds"))