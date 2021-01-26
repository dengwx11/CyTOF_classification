options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

sample <- as.character(args[1])

library(Seurat)
library(ggplot2)
library(mltools)
library(data.table)
library(lsa)
library(aricode)
library(kBET)
library(cowplot)
options(repr.plot.width=14, repr.plot.height=10)
library(viridis)
library(mclust)
library("hrbrthemes")
source('Opt/preparation.R')

seur <- readRDS(paste0('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype/seur_', sample, '.rds'))
rst <- readRDS(paste0('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype/rst_', sample, '.rds'))
truth <- seur$label

celltype_pred <- seur$pred
rownames(rst$H) <- rst$H[,1]
H_est_subsetting <- data.frame(rst$H[,-1])
H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))
celltypes <- c('CD4 T-cells', 'CD8 T-cells', 'monocytes', 'naïve B cells', 'NK cells')
H_est_subsetting <- H_est_subsetting[celltypes,]

cutoff.list <- c(1:100)/100                          
df <- data.frame(cutoff = cutoff.list)
df$acc <- df$N <- df$cos <- df$nmi <- df$ari <- df$sil <- 0              
                          
for(j in c(1:length(cutoff.list))){
    ## cnt_max
    cnt_max = 0
    N_subsetting = 0
    idx.subsetting = c()
        for(i in 1:length(truth)){
            if(max(H_est_subsetting[,i])>cutoff.list[j]){
                idx.subsetting = c(idx.subsetting,i)
                N_subsetting = N_subsetting+1
                cnt_max = cnt_max + 1*(truth[i]==celltype_pred[i])
            }

        }
    if(N_subsetting == 0) {
        df$acc[j] <- 0
    } else {
       df$acc[j]<-cnt_max/N_subsetting
    } 
    if(length(idx.subsetting) > 0) {
        df$N[j] <- N_subsetting
    
        truth_sub <- truth[idx.subsetting]
        truth_fact_sub <- as.factor(truth_sub)

        
        celltype_pred_sub <- seur$pred[idx.subsetting]
        celltype_pred_fact_sub <- as.factor(celltype_pred_sub)
        celltype_pred_fact_sub <- factor(celltype_pred_fact_sub, levels = levels(truth_fact_sub))

        ## ARI
        ari <- adjustedRandIndex(celltype_pred_sub, truth_sub)
        df$ari[j] <- ari
        
        ## Cosine similarity
        truth_onehot_sub <- as.data.frame(t(one_hot(as.data.table(truth_fact_sub))))
        cos_sim <- mean(mapply(cosine, truth_onehot_sub, as.data.frame(H_est_subsetting[,idx.subsetting])))
        df$cos[j] <- cos_sim

        ## NMI
        pred_onehot_sub <- as.data.frame(t(one_hot(as.data.table(celltype_pred_fact_sub))))
        nmi <- mean(mapply(NMI, truth_onehot_sub, as.data.frame(pred_onehot_sub)))
        df$nmi[j] <- nmi

        ## Silhouette
        if(length(unique(celltype_pred_sub)) == 1) {
            df$sil[j] <- 0
        } else{
            pca.data <- list()
            pca.data$x <- seur@reductions$pca@cell.embeddings[idx.subsetting,]
            sil <- batch_sil(pca.data, as.numeric(as.factor(celltype_pred_sub)))
            df$sil[j] <- sil
        }   
    }
}
write.table(df, paste0('/home/bz234/project/Results/CyTOF/BCR_Cytof/subset_metrics/', sample, '.txt'))