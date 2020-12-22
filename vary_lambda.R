source("Simulation.R")
source('run_para.R')
library(mltools)
library(data.table)
library(lsa)
library(aricode)
library(kBET)

gamma <- simulate_gamma(iteration, K, D)
A0 <- gamma$gamma
W <- simulate_w(D,K,A0, p.0 = p.0, q.0=q.0, p.neg1= p.neg1, q.neg1 = q.neg1, 
                big_w_mean=big_w_mean, big_tau_w = big_tau_w, 
                small_w_mean=small_w_mean, small_tau_w = small_tau_w)
AS <- simulate_AS(D,K,W,corr=corr)
label.output <- simulate_label(D,N,K,prob_k)
X <- simulate_x(D,N,W,label.output$label, mean_var_ratio = mean_var_ratio)
true.H <- label.output$H

truth <- label.output$label
truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))

rownames(X) <- c(1:D)
colnames(X) <- c(1:N)
seur <- CreateSeuratObject(X)
seur@meta.data$label <- label.output$label[1,]
seur <- NormalizeData(seur) 
seur@assays$RNA@data = seur@assays$RNA@counts
seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
seur <- ScaleData(seur)
seur <- RunPCA(seur,verbose = TRUE,features = rownames(seur), seed.use = NULL)
seur <- FindNeighbors(seur, dims = 1:6)
seur <- FindClusters(seur, resolution = .6)

rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=2,lambda1.on=T,lambda2.on=T)
rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)

accu_vec <- c()
cos_vec <- c()
ari_vec <- c()
nmi_vec <- c()
sil_vec <- c()

## Vary lambda2
for(lambda2 in seq(0, 60, 5)){
    rst<-run(X,rst.para$para$lambda1,lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
    celltype_pred <- apply(rst$H, 2, predict)
    pred_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred)))))
    cnt_max = 0
    for(i in 1:length(truth)){
        cnt_max = cnt_max + infer_max(truth[i], rst$H[,i])
    }  

    ## Accuracy
    accu <- cnt_max / N
    accu_vec <- c(accu_vec, accu)

    ## Cosine similarity
    cos_sim <- mean(mapply(cosine, truth_onehot, as.data.frame(rst$H)))
    cos_vec <- c(cos_vec, cos_sim)

    ## ARI
    ari <- adjustedRandIndex(celltype_pred, truth)
    ari_vec <- c(ari_vec, ari)

    ## NMI
    nmi <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot)))
    nmi_vec <- c(nmi_vec, nmi)

    ## Silhouette
    pca.data <- list()
    pca.data$x <- seur@reductions$pca@cell.embeddings
    sil <- batch_sil(pca.data, celltype_pred)
    sil_vec <- c(sil_vec, sil)
}

dat <- matrix(,length(accu_vec),ncol=5)
dat[,1] <- accu_vec
dat[,2] <- cos_vec
dat[,3] <- ari_vec
dat[,4] <- nmi_vec
dat[,5] <- sil_vec

png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda2.png")
matplot(seq(0, 60, 5), dat, type = c("b"),pch=1,col = 1:5, xlab = 'lambda 2', ylab = 'Value', main = 'Varying lambda 2 (Senario 1)') #plot
legend('right', legend = c('Accuracy', 'Cosine', 'ARI', 'NMI', 'ASW'), col=1:5, pch=1)
dev.off()

png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda2_ex0.png")
matplot(seq(5, 60, 5), dat[2:length(accu_vec),1:4], type = c("b"),pch=1,col = 1:4, xlab = 'lambda 2', ylab = 'Value', main = 'Varying lambda 2 (Senario 1)') #plot
legend('right', legend = c('Accuracy', 'Cosine', 'ARI', 'NMI'), col=1:4, pch=1)
dev.off()

png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda2_ASW.png")
plot(seq(0, 60, 5), sil_vec, xlab = 'lambda 2', ylab = 'Value', main = 'Varying lambda 2 ASW (Senario 1)')
lines(seq(0, 60, 5), sil_vec)
dev.off()

accu_vec_l1 <- c()
cos_vec_l1 <- c()
ari_vec_l1 <- c()
nmi_vec_l1 <- c()
sil_vec_l1 <- c()

## Vary lambda1
for(lambda1 in seq(0, 5, 0.4)){
    rst<-run(X,lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
    celltype_pred <- apply(rst$H, 2, predict)
    pred_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred)))))
    cnt_max = 0
    for(i in 1:length(truth)){
        cnt_max = cnt_max + infer_max(truth[i], rst$H[,i])
    }  

    ## Accuracy
    accu <- cnt_max / N
    accu_vec_l1 <- c(accu_vec_l1, accu)

    ## Cosine similarity
    cos_sim <- mean(mapply(cosine, truth_onehot, as.data.frame(rst$H)))
    cos_vec_l1 <- c(cos_vec_l1, cos_sim)

    ## ARI
    ari <- adjustedRandIndex(celltype_pred, truth)
    ari_vec_l1 <- c(ari_vec_l1, ari)

    ## NMI
    nmi <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot)))
    nmi_vec_l1 <- c(nmi_vec_l1, nmi)

    ## Silhouette
    pca.data <- list()
    pca.data$x <- seur@reductions$pca@cell.embeddings
    sil <- batch_sil(pca.data, celltype_pred)
    sil_vec_l1 <- c(sil_vec_l1, sil)
}

dat_l1 <- matrix(,length(accu_vec_l1),ncol=5)
dat_l1[,1] <- accu_vec_l1
dat_l1[,2] <- cos_vec_l1
dat_l1[,3] <- ari_vec_l1
dat_l1[,4] <- nmi_vec_l1
dat_l1[,5] <- sil_vec_l1

png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda1.png")
matplot(seq(0, 5, 0.4), dat_l1, type = c("b"),pch=1,col = 1:5, xlab = 'lambda 1', ylab = 'Value', main = 'Varying lambda 1 (Senario 1)') #plot
legend('right', legend = c('Accuracy', 'Cosine', 'ARI', 'NMI', 'ASW'), col=1:5, pch=1)
dev.off()

png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda1_ex_asw.png")
matplot(seq(0, 5, 0.4), dat_l1[,1:4], type = c("b"),pch=1,col = 1:4, xlab = 'lambda 1', ylab = 'Value', main = 'Varying lambda 1 (Senario 1)') #plot
legend('right', legend = c('Accuracy', 'Cosine', 'ARI', 'NMI'), col=1:4, pch=1)
dev.off()

png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda1_ASW.png")
plot(seq(0, 5, 0.4), sil_vec_l1, xlab = 'lambda 1', ylab = 'Value', main = 'Varying lambda 1 ASW (Senario 1)')
lines(seq(0, 5, 0.4), sil_vec)
dev.off()