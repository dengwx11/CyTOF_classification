library(mltools)
library(data.table)
library(lsa)
library(aricode)
library(kBET)
library(cowplot)
library(ggplot2)
library(mclust)

seur_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype', pattern='*seur*', full=T, ignore.case = TRUE)
rst_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype', pattern='*rst*', full=T, ignore.case = TRUE)
BCR_XL_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype', pattern=glob2rx('*BCR-XL.rds'), full=T, ignore.case = TRUE)
Reference_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype', pattern=glob2rx('*Reference.rds'), full=T, ignore.case = TRUE)
rst_BCR_files <- intersect(rst_files,BCR_XL_files)
rst_Reference_files <- intersect(rst_files,Reference_files)
seur_BCR_files <- intersect(seur_files,BCR_XL_files)
seur_Reference_files <- intersect(seur_files,Reference_files)

## BCR
acc_raw <- c()
acc_subsetting <- c()
cos_raw <- c()
cos_subsetting <- c()
nmi_raw <- c()
nmi_subsetting <- c()
ari_raw <- c()
ari_subsetting <- c()
sil_raw <- c()
sil_subsetting <- c()
for(sample.idx in 1:length(rst_BCR_files)){
    rst <- readRDS(rst_BCR_files[sample.idx])
    seur <- readRDS(seur_BCR_files[sample.idx])
    truth <- seur$label
    pred <- seur$pred
    H <- rst$H
    rownames(H) <- H[,1]
    H <- H[,-1]
    H <- H[c('CD4 T-cells', 'CD8 T-cells', 'monocytes', 'naïve B cells', 'NK cells'),]
    truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))
#     celltype_pred_fact <- as.factor(pred)
#     celltype_pred_fact <- factor(celltype_pred_fact, levels = as.character(1:dim(table(truth))))
    pred_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(pred)))))
    
    ## raw 
    # Accuracy
    cnt_max = 0
    N = ncol(seur)
    for(i in 1:length(truth)){
        cnt_max = cnt_max + 1*(truth[i]==pred[i])
    }
    print(paste0(rst_BCR_files[sample.idx]," : ", cnt_max))
    print(cnt_max/N)
    acc_raw <- append(acc_raw,cnt_max/N)
    
    # Cosine similarity
    cos_raw <- append(cos_raw, mean(mapply(cosine, truth_onehot, as.data.frame(H))))
    
    # NMI
    nmi_raw <- append(nmi_raw, mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot))))
    
    # ARI
    ari_raw <- append(ari_raw, adjustedRandIndex(pred, truth))
    
    # Silhouette
    pca.data <- list()
    pca.data$x <- seur@reductions$pca@cell.embeddings
    sil_raw <- append(sil_raw, batch_sil(pca.data, as.numeric(as.factor(pred))))
    
    ## subsetting
    # Accuracy
    H_est_subsetting <- apply(H,2,function(x) x/sum(x))
    cnt_max = 0
    N_subsetting = 0
    idx.subsetting = c()
    for(i in 1:length(truth)){
        if(max(H_est_subsetting[,i])>0.5){
            idx.subsetting = append(idx.subsetting,i)
            N_subsetting = N_subsetting+1
            cnt_max = cnt_max + 1*(truth[i]==pred[i])
        }

    }
    cnt_max
    print(cnt_max/N_subsetting)
    acc_subsetting <- append(acc_subsetting,cnt_max/N_subsetting)
                              
    # Cosine similarity
    truth_sub <- truth[idx.subsetting]
    truth_fact_sub <- as.factor(truth_sub)
#     truth_fact_sub <- factor(truth_fact_sub, levels = as.character(1:dim(table(truth_sub))))
    truth_onehot_sub <- as.data.frame(t(one_hot(as.data.table(truth_fact_sub))))
    cos_subsetting <- append(cos_subsetting, mean(mapply(cosine, truth_onehot_sub, 
                                                         as.data.frame(H_est_subsetting[,idx.subsetting]))))
    # NMI
    celltype_pred_sub <- pred[idx.subsetting]
    celltype_pred_sub <- as.factor(celltype_pred_sub)
    celltype_pred_sub <- factor(celltype_pred_sub, levels = levels(truth_fact_sub))
    pred_onehot_sub <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred_sub)))))
    nmi_subsetting <- append(nmi_subsetting, mean(mapply(NMI, truth_onehot_sub, pred_onehot_sub)))
                              
    # ARI
    ari_subsetting <- append(ari_subsetting, adjustedRandIndex(celltype_pred_sub, truth_sub))
                              
    # Silhouette
    pca.data <- list()
    pca.data$x <- seur@reductions$pca@cell.embeddings[idx.subsetting,]
    sil_subsetting <- append(sil_subsetting, batch_sil(pca.data, as.numeric(as.factor(celltype_pred_sub))))
}

BCR_raw <- data.frame(
    accuracy = acc_raw,
    cos = cos_raw,
    nmi = nmi_raw,
    ari = ari_raw,
    sil = sil_raw
)
BCR_sub <- data.frame(
    accuracy = acc_subsetting,
    cos = cos_subsetting,
    nmi = nmi_subsetting,
    ari = ari_subsetting,
    sil = sil_subsetting
)
write.table(BCR_raw, '/home/bz234/project/Results/CyTOF/BCR_Cytof/benchmark_raw_subset/BCR_raw_5type.txt')
write.table(BCR_sub, '/home/bz234/project/Results/CyTOF/BCR_Cytof/benchmark_raw_subset/BCR_sub_5type.txt')

## Reference
acc_raw <- c()
acc_subsetting <- c()
cos_raw <- c()
cos_subsetting <- c()
nmi_raw <- c()
nmi_subsetting <- c()
ari_raw <- c()
ari_subsetting <- c()
sil_raw <- c()
sil_subsetting <- c()
for(sample.idx in 1:length(rst_BCR_files)){
    rst <- readRDS(rst_Reference_files[sample.idx])
    seur <- readRDS(seur_Reference_files[sample.idx])
    truth <- seur$label
    pred <- seur$pred
    H <- rst$H
    rownames(H) <- H[,1]
    H <- H[,-1]
    H <- H[c('CD4 T-cells', 'CD8 T-cells', 'monocytes', 'naïve B cells', 'NK cells'),]
    truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))
#     celltype_pred_fact <- as.factor(pred)
#     celltype_pred_fact <- factor(celltype_pred_fact, levels = as.character(1:dim(table(truth))))
    pred_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(pred)))))
    
    ## raw 
    # Accuracy
    cnt_max = 0
    N = ncol(seur)
    for(i in 1:length(truth)){
        cnt_max = cnt_max + 1*(truth[i]==pred[i])
    }
    print(paste0(rst_BCR_files[sample.idx]," : ", cnt_max))
    print(cnt_max/N)
    acc_raw <- append(acc_raw,cnt_max/N)
    
    # Cosine similarity
    cos_raw <- append(cos_raw, mean(mapply(cosine, truth_onehot, as.data.frame(H))))
    
    # NMI
    nmi_raw <- append(nmi_raw, mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot))))
    
    # ARI
    ari_raw <- append(ari_raw, adjustedRandIndex(pred, truth))
    
    # Silhouette
    pca.data <- list()
    pca.data$x <- seur@reductions$pca@cell.embeddings
    sil_raw <- append(sil_raw, batch_sil(pca.data, as.numeric(as.factor(pred))))
    
    ## subsetting
    # Accuracy
    H_est_subsetting <- apply(H,2,function(x) x/sum(x))
    cnt_max = 0
    N_subsetting = 0
    idx.subsetting = c()
    for(i in 1:length(truth)){
        if(max(H_est_subsetting[,i])>0.5){
            idx.subsetting = append(idx.subsetting,i)
            N_subsetting = N_subsetting+1
            cnt_max = cnt_max + 1*(truth[i]==pred[i])
        }

    }
    cnt_max
    print(cnt_max/N_subsetting)
    acc_subsetting <- append(acc_subsetting,cnt_max/N_subsetting)
                              
    # Cosine similarity
    truth_sub <- truth[idx.subsetting]
    truth_fact_sub <- as.factor(truth_sub)
#     truth_fact_sub <- factor(truth_fact_sub, levels = as.character(1:dim(table(truth_sub))))
    truth_onehot_sub <- as.data.frame(t(one_hot(as.data.table(truth_fact_sub))))
    cos_subsetting <- append(cos_subsetting, mean(mapply(cosine, truth_onehot_sub, 
                                                         as.data.frame(H_est_subsetting[,idx.subsetting]))))
    # NMI
    celltype_pred_sub <- pred[idx.subsetting]
    celltype_pred_sub <- as.factor(celltype_pred_sub)
    celltype_pred_sub <- factor(celltype_pred_sub, levels = levels(truth_fact_sub))
    pred_onehot_sub <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred_sub)))))
    nmi_subsetting <- append(nmi_subsetting, mean(mapply(NMI, truth_onehot_sub, pred_onehot_sub)))
                              
    # ARI
    ari_subsetting <- append(ari_subsetting, adjustedRandIndex(celltype_pred_sub, truth_sub))
                              
    # Silhouette
    pca.data <- list()
    pca.data$x <- seur@reductions$pca@cell.embeddings[idx.subsetting,]
    sil_subsetting <- append(sil_subsetting, batch_sil(pca.data, as.numeric(as.factor(celltype_pred_sub))))
}

ref_raw <- data.frame(
    accuracy = acc_raw,
    cos = cos_raw,
    nmi = nmi_raw,
    ari = ari_raw,
    sil = sil_raw
)
ref_sub <- data.frame(
    accuracy = acc_subsetting,
    cos = cos_subsetting,
    nmi = nmi_subsetting,
    ari = ari_subsetting,
    sil = sil_subsetting
)
write.table(ref_raw, '/home/bz234/project/Results/CyTOF/BCR_Cytof/benchmark_raw_subset/ref_raw_5type.txt')
write.table(ref_sub, '/home/bz234/project/Results/CyTOF/BCR_Cytof/benchmark_raw_subset/ref_sub_5type.txt')

BCR_raw_plot <- data.frame(
    Metrics = c('ARI', 'Accuracy', 'NMI', 'Cosine Similarity', 'ASW'),
    Values = colMeans(BCR_raw)[c('ari', 'accuracy', 'nmi', 'cos', 'sil')]
)
BCR_raw_plot$Metrics <- factor(BCR_raw_plot$Metrics,levels = c('ASW', 'Cosine Similarity', 'NMI', 'Accuracy', 'ARI'))
BCR_sub_plot <- data.frame(
    Metrics = c('ARI', 'Accuracy', 'NMI', 'Cosine Similarity', 'ASW'),
    Values = colMeans(BCR_sub)[c('ari', 'accuracy', 'nmi', 'cos', 'sil')]
)
BCR_sub_plot$Metrics <- factor(BCR_sub_plot$Metrics,levels = c('ASW', 'Cosine Similarity', 'NMI', 'Accuracy', 'ARI'))
ref_raw_plot <- data.frame(
    Metrics = c('ARI', 'Accuracy', 'NMI', 'Cosine Similarity', 'ASW'),
    Values = colMeans(ref_raw)[c('ari', 'accuracy', 'nmi', 'cos', 'sil')]
)
ref_raw_plot$Metrics <- factor(ref_raw_plot$Metrics,levels = c('ASW', 'Cosine Similarity', 'NMI', 'Accuracy', 'ARI'))
ref_sub_plot <- data.frame(
    Metrics = c('ARI', 'Accuracy', 'NMI', 'Cosine Similarity', 'ASW'),
    Values = colMeans(ref_sub)[c('ari', 'accuracy', 'nmi', 'cos', 'sil')]
)
ref_sub_plot$Metrics <- factor(ref_sub_plot$Metrics,levels = c('ASW', 'Cosine Similarity', 'NMI', 'Accuracy', 'ARI'))

dat_lst <- list(BCR_raw_plot, BCR_sub_plot, ref_raw_plot, ref_sub_plot)
plotlist <- list()
title_vec <- c('BCR Raw 5 Cell types', 'BCR Subsetting 5 Cell types', 'Reference Raw 5 Cell types', 
               'Reference Subsetting 5 Cell types')
for(i in 1:length(dat_lst)){
  plotlist[[i]] <- ggplot(data = dat_lst[[i]], aes(x = Metrics, y = Values, fill = Metrics, width=.7)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(y = 'Values', title = title_vec[i]) +
  coord_flip() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) + 
  ylim(0, .8)
}
pll <- plot_grid(plotlist = plotlist, ncol=4)
ggsave('/home/bz234/project/Results/CyTOF/BCR_Cytof/benchmark_raw_subset/bar_plots.png', pll, width = 30, height = 5)
