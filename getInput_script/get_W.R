options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

ADT_directory <- as.character(args[1])
datatype <- as.character(args[2]) ## CyTOF/CITE
output_directory <- as.character(args[3])

library(Seurat)
library(HDCytoData)

## load data and get ADT expression profile matrix
ADT <- readRDS(ADT_directory)
if(datatype == 'CyTOF'){

    E_ADT <- assay(ADT[, colData(ADT)$marker_class == "type"])
    E_ADT <- t(E_ADT)

    colnames(E_ADT) <- c(1:ncol(E_ADT))
    metadata<-rowData(ADT)
    metadata <- data.frame(metadata@listData)
    rownames(metadata) <- colnames(E_ADT)

    celltype <- metadata$population_id

}else{
    E_ADT <- GetAssayData(ADT[['ADT']], slot='counts')
    E_ADT <- as.matrix(E_ADT)

    celltype <- Idents(ADT) # inhouse
    #celltype <- ADT$celltype_clean # covid
}

## normalization 
cofactor <- 5
E_ADT <- asinh(E_ADT / cofactor)

## get W
df <- data.frame(celltype = celltype, t(E_ADT))
W <- aggregate(. ~ celltype, df, mean)
rownames(W) <-W$celltype
W <- W[,-1]
W <- t(as.matrix(W))


W <- as.matrix(W)
saveRDS(W, output_directory)


