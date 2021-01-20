options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

seurat_directory <- as.character(args[1])
saverX_directory <- as.character(args[2])
geneList_directory <- as.character(args[3])
ADTList_directory <- as.character(args[4]) 
output_directory <- as.character(args[5])

library(Seurat)
library(glmnet)

## load data
CITE <- readRDS(seurat_directory)
RNA.saverx <- readRDS(saverX_directory)
geneList <- read.csv(geneList_directory, row.names = 1)[,1]
ADTList <- as.character(read.table(ADTList_directory, header = T)[,2])

## Antibody expression profile matrix
E_ADT <- GetAssayData(CITE[['ADT']], slot='counts')
E_ADT <- E_ADT[match(ADTList,rownames(E_ADT)),]
cofactor <- 5
E_ADT <- asinh(E_ADT / cofactor)

## saverX RNA expression profile matrix
if(class(RNA.saverx) == 'matrix'){
        E_RNA <- RNA.saverx[match(geneList, rownames(RNA.saverx)), colnames(E_ADT)]
} else{
        E_RNA <- RNA.saverx$estimate[match(geneList, rownames(RNA.saverx$estimate)), colnames(E_ADT)]
}

## Lasso function
elnet <- function(j){
    set.seed(2020)
    #sample.idx <- sample(c(1:ncol(E_RNA)),15000)
    sample.idx <- c(1:ncol(E_RNA))
    x.train <- t(E_RNA)[sample.idx,]
    y.train <- E_ADT[j,sample.idx]
#   x.test <- t(E_RNA)[-sample.idx,]
#    y.test <- E_ADT[j,-sample.idx]

#    fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5)
    fit5 <- cv.glmnet(x.train, y.train, type.measure="mse", family="gaussian", alpha=.5)
    return(as.matrix(coef(fit5)))
}

rst <- lapply(c(1:nrow(E_ADT)), elnet)
A <- Reduce(cbind, rst)
colnames(A) <- rownames(E_ADT)

## save A
write.table(A, output_directory, quote=F)
