source("Opt/Simulation.R")
source('Opt/run_para.R')
library(mltools)
library(data.table)
library(lsa)
library(aricode)
library(kBET)

# gamma <- simulate_gamma(iteration, K, D)
# A0 <- gamma$gamma
# W <- simulate_w(D,K,A0, p.0 = p.0, q.0=q.0, p.neg1= p.neg1, q.neg1 = q.neg1, 
#                 big_w_mean=big_w_mean, big_tau_w = big_tau_w, 
#                 small_w_mean=small_w_mean, small_tau_w = small_tau_w)
# AS <- simulate_AS(D,K,W,corr=corr)
# label.output <- simulate_label(D,N,K,prob_k)
# X <- simulate_x(D,N,W,label.output$label, mean_var_ratio = mean_var_ratio)
# true.H <- label.output$H

# truth <- label.output$label
# truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))

# rownames(X) <- c(1:D)
# colnames(X) <- c(1:N)
# seur <- CreateSeuratObject(X)
# seur@meta.data$label <- label.output$label[1,]
# seur <- NormalizeData(seur) 
# seur@assays$RNA@data = seur@assays$RNA@counts
# seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
# seur <- ScaleData(seur)
# seur <- RunPCA(seur,verbose = TRUE,features = rownames(seur), seed.use = NULL)
# seur <- FindNeighbors(seur, dims = 1:6)
# seur <- FindClusters(seur, resolution = .6)

# rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=2,lambda1.on=T,lambda2.on=T)
# rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
#             AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)

accu_vec <- c()
cos_vec <- c()
ari_vec <- c()
nmi_vec <- c()
# sil_vec <- c()

## Vary lambda2

## Senario 1, 2
# for(lambda2 in seq(10, 50, 4)){
## Senario 3
# for(lambda2 in seq(25, 65, 4)){
## Senario 4
# for(lambda2 in seq(2, 7, 0.4)){
## Senario 5
for(lambda2 in seq(4, 10, 0.5)){
    rst<-run(X,rst.para$para$lambda1,lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
    celltype_pred <- apply(rst$H, 2, predict)
    truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))
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

#     ## Silhouette
#     pca.data <- list()
#     pca.data$x <- seur@reductions$pca@cell.embeddings
#     sil <- batch_sil(pca.data, celltype_pred)
#     sil_vec <- c(sil_vec, sil)
}

dat_l2_5 <- data.frame(
    lambda2 = seq(4, 10, 0.5),
    accu = accu_vec,
    cos = cos_vec, 
    ari = ari_vec,
    nmi = nmi_vec
)
dat_l2_5 <- melt(dat_l2_5, id.vars = 'lambda2', variable.name = 'Metrics')
saveRDS(dat_l2_5, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/dat_l2_5.rds')

l2_pll <- list()
l2_pll[[5]] <- ggplot(dat_l2_5, aes(lambda2,value)) + 
geom_point(aes(colour = Metrics), shape = 4) +
labs(y = 'Annotation Metric', x = bquote(lambda[2]), title = 'Scenario 5') + 
geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
# ylim(0.8, 1) +
theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
# geom_smooth(aes(x = lambda2, y = value, color = Metrics))
pll_2 <- plot_grid(plotlist = l2_pll, ncol=5, align = 'h')
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/Plots/vary_lambda2.png', pll_2, width = 25, height = 5)

# png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda2.png")
# matplot(seq(0, 60, 5), dat, type = c("b"),pch=1,col = 1:5, xlab = 'lambda 2', ylab = 'Value', main = 'Varying lambda 2 (Senario 1)') #plot
# legend('right', legend = c('Accuracy', 'Cosine', 'ARI', 'NMI', 'ASW'), col=1:5, pch=1)
# dev.off()

# png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda2_ex0.png")
# matplot(seq(5, 60, 5), dat[2:length(accu_vec),1:4], type = c("b"),pch=1,col = 1:4, xlab = 'lambda 2', ylab = 'Value', main = 'Varying lambda 2 (Senario 1)') #plot
# legend('right', legend = c('Accuracy', 'Cosine', 'ARI', 'NMI'), col=1:4, pch=1)
# dev.off()

# png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda2_ASW.png")
# plot(seq(0, 60, 5), sil_vec, xlab = 'lambda 2', ylab = 'Value', main = 'Varying lambda 2 ASW (Senario 1)')
# lines(seq(0, 60, 5), sil_vec)
# dev.off()

accu_vec_l1 <- c()
cos_vec_l1 <- c()
ari_vec_l1 <- c()
nmi_vec_l1 <- c()

## Vary lambda1
## Senario 1,2,3
# for(lambda1 in seq(0, 5, 0.4)){
## Senario 4
# for(lambda1 in seq(0, 0.04, 0.003)){
## Senario 5
for(lambda1 in seq(0, 0.4, 0.03)){    
    rst<-run(X,lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
    celltype_pred <- apply(rst$H, 2, predict)
    truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))
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

    # ## Silhouette
    # pca.data <- list()
    # pca.data$x <- seur@reductions$pca@cell.embeddings
    # sil <- batch_sil(pca.data, celltype_pred)
    # sil_vec_l1 <- c(sil_vec_l1, sil)
}

dat_l1_5 <- data.frame(
    lambda1 = seq(0, 0.4, 0.03),
    accu = accu_vec_l1,
    cos = cos_vec_l1, 
    ari = ari_vec_l1,
    nmi = nmi_vec_l1
)
dat_l1_5 <- melt(dat_l1_5, id.vars = 'lambda1', variable.name = 'Metrics')
saveRDS(dat_l1_5, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/dat_l1_5.rds')

l1_pll <- list()
l1_pll[[5]] <- ggplot(dat_l1_5, aes(lambda1,value)) + 
geom_point(aes(colour = Metrics), shape = 4) +
labs(y = 'Annotation Metric', x = bquote(lambda[1]), title = 'Scenario 5') + 
geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
# ylim(0.8, 1) +
theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
# geom_smooth(aes(x = lambda2, y = value, color = Metrics))
pll_1 <- plot_grid(plotlist = l1_pll, ncol=5, align = 'h')
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/Plots/vary_lambda1.png', pll_1, width = 25, height = 5)

# png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda1.png")
# matplot(seq(0, 5, 0.4), dat_l1, type = c("b"),pch=1,col = 1:5, xlab = 'lambda 1', ylab = 'Value', main = 'Varying lambda 1 (Senario 1)') #plot
# legend('right', legend = c('Accuracy', 'Cosine', 'ARI', 'NMI', 'ASW'), col=1:5, pch=1)
# dev.off()

# png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda1_ex_asw.png")
# matplot(seq(0, 5, 0.4), dat_l1[,1:4], type = c("b"),pch=1,col = 1:4, xlab = 'lambda 1', ylab = 'Value', main = 'Varying lambda 1 (Senario 1)') #plot
# legend('right', legend = c('Accuracy', 'Cosine', 'ARI', 'NMI'), col=1:4, pch=1)
# dev.off()

# png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_vary_lambda1_ASW.png")
# plot(seq(0, 5, 0.4), sil_vec_l1, xlab = 'lambda 1', ylab = 'Value', main = 'Varying lambda 1 ASW (Senario 1)')
# lines(seq(0, 5, 0.4), sil_vec)
# dev.off()