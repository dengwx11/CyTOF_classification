options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

seurat_directory <- as.character(args[1])
saverX_directory <- as.character(args[2])
geneList_directory <- as.character(args[3])
ADT_data <- as.character(args[4]) ## BCR/inhouse/covid
output_geneList_directory  <- as.character(args[5])
output_directory <- as.character(args[6])

library(Seurat)
library(dplyr)

## load data
seur <- readRDS(seurat_directory)
RNA.saverx <- readRDS(saverX_directory)
geneList <- as.character(read.table(geneList_directory)[,1])
all_RNAs <- as.character(read.table('/gpfs/loomis/project/zhao/bz234/Results/CyTOF/Tomo/geneList/intersect_RNAs.txt')[[1]])
DefaultAssay(seur) <- 'RNA'

## Get Top 15 DEGs for each cell type
if(ADT_data == 'BCR'){
    Idents(seur) <- seur$BCR_celltype
}else if(ADT_data == 'inhouse'){
    Idents(seur) <- seur$inhouse_celltype
}else if(ADT_data == 'covid'){
    Idents(seur) <- seur$celltype_clean
}
seur.markers <- FindAllMarkers(seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seur.markers <- seur.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) 
geneList <- c(seur.markers$gene, geneList)

## Identify the 50 most highly variable genes
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 50)
geneList <- c(VariableFeatures(seur), geneList)
geneList <- unique(geneList)

## Not in RNAs
geneList <- all_RNAs[match(geneList, all_RNAs)]
geneList <- geneList[!is.na(geneList)]
write.csv(geneList, output_geneList_directory)

if(class(RNA.saverx) == 'matrix'){
	E_RNA <- RNA.saverx[as.character(geneList), colnames(seur)]
} else{
	E_RNA <- RNA.saverx$estimate[as.character(geneList), colnames(seur)]
}

## RNA signature matrix
if(ADT_data == 'BCR'){
    df <- data.frame(celltype = seur$BCR_celltype, t(E_RNA))
}else if(ADT_data == 'inhouse'){
    df <- data.frame(celltype = seur$inhouse_celltype, t(E_RNA))
}else if(ADT_data == 'covid'){
    df <- data.frame(celltype = seur$celltype_clean, t(E_RNA))
}

S <- aggregate(. ~ celltype, df, mean)
rownames(S) <-S$celltype
S <- S[,-1]
S <- t(as.matrix(S))

S <- rbind(1, S)
#rownames(S)[1] <- rownames(A)[1]

saveRDS(S, output_directory)
