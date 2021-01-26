options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

expr.dir <- as.character(args[1])
sample.id <- as.character(args[2])


set.seed(2021)
setwd('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/CyTOF_classification')
source('Opt/run_para.R')
source('Opt/preparation.R')
library(Seurat)

# load input
expr <- readRDS(expr.dir)
A0 <- read.csv('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/CyTOF_classification/write/matrix_A0/A0_covid_Rodrigeuz_in_Rami.csv',header =T, row.names=1)
#AS <- matrix(rep(0, nrow(A0)*ncol(A0)),nrow = nrow(A0), ncol=ncol(A0))
AS <- readRDS("~/zhao-data/new_cytof/write/covid_Rode_in_covid_rami_as.rds")
lineage_channels = rownames(A0)
X <- expr[,lineage_channels]
X <- t(as.matrix(X))
A0 <- as.matrix(A0)
colnames(A0) <- c('Classical monocytes','NC & IM monocytes','effector CD4+T cell',
                 'Memory CD4 T','Naive CD4 T','effector CD8+T cell','Memory CD8 T','Naive CD8 T','Tregs',
                 'Naive B','Plasma cells','Memory B','NK cell','Myeloid DC','Plasmacytoid DC')
AS <- AS[,colnames(A0)]

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

saveRDS(rst,paste0('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/covid_Rodriguez/covid_rami_A_S/rst_',sample.id,".rds"))
saveRDS(seur,paste0('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/covid_Rodriguez/covid_rami_A_S/seur_',sample.id,".rds"))