source("Opt/Simulation.R")
source('Opt/run_para.R')
library(ggplot2)
library(scales)
library(cowplot)

# loops <- 100
# corr_vec <- c()
# i <- 0
# while(i < loops) {
#     ## generate dataset: X, AS, A0, W(ground truth), true.label
#     gamma <- simulate_gamma(iteration, K, D)
#     A0 <- gamma$gamma
#     W <- simulate_w(D,K,A0, p.0 = p.0, q.0=q.0, p.neg1= p.neg1, q.neg1 = q.neg1, 
#                 big_w_mean=big_w_mean, big_tau_w = big_tau_w, 
#                 small_w_mean=small_w_mean, small_tau_w = small_tau_w)
#     AS <- simulate_AS(D,K,W,corr=corr)
#     label.output <- simulate_label(D,N,K,prob_k)
#     X <- simulate_x(D,N,W,label.output$label, mean_var_ratio = mean_var_ratio)
#     true.H <- label.output$H

#     rownames(X) <- c(1:D)
#     colnames(X) <- c(1:N)
#     seur <- CreateSeuratObject(X)
#     seur@meta.data$label <- label.output$label[1,]
#     seur <- NormalizeData(seur) 
#     seur@assays$RNA@data = seur@assays$RNA@counts
#     seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
#     seur <- ScaleData(seur)
#     seur <- RunPCA(seur,verbose = TRUE,features = rownames(seur), seed.use = NULL)
#     seur <- FindNeighbors(seur, dims = 1:6)
#     seur <- FindClusters(seur, resolution = .6)

#     rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=2,lambda1.on=T,lambda2.on=T)
#     rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
#             AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
#     corr_vec <- c(corr_vec, cor(as.vector(rst$W), as.vector(W)))
#     print(corr_vec)
#     print(i)
#     i <- i + 1
# }

# write.table(corr_vec, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/w_corr_vec.txt')

run_eachW <- function(X,lambda1,lambda2,mu, eta, AS, A0, D, K, N, epsilon = 10^(-4),fixed_loop = 0) {
    ## Initialization
    W_iter <- matrix(1,nrow = D, ncol = K)
    H <- matrix(1,nrow = K, ncol = N)
    corr_vec <- c()

    L.current <- compute_L(X,W,H,lambda1,lambda2,mu,eta)
    L.min <- L.current
    L.prev <- 100000
    L.save <- c(L.current)
    k=1
     while(abs(L.prev - L.current) > epsilon){
    #while(abs(L.prev - L.current) > epsilon && L.min >= L.current){
    #for(i in 1:1000){
        W_iter <- update_W(X, W_iter, H, mu, AS, A0, D, K, lambda1, lambda2)
        corr_vec <- c(corr_vec, cor(as.vector(W_iter), as.vector(W)))
        L.prev <- L.current
        L.current <- compute_L(X,W_iter,H,lambda1,lambda2,mu,eta)
        #print(paste0("update W, L change ", L.current-L.prev))
        H <- update_H(X, W_iter, H, K, N, eta)
        L.prev <- L.current
        L.current <- compute_L(X,W_iter,H,lambda1,lambda2,mu,eta)
        #print(paste0("update H, L change ", L.current-L.prev))
        L.save <- c(L.save, L.current)
        if(L.min>L.current) L.min=L.current
        #print(L.current)
        k = k+1
        if(fixed_loop >0 && k>fixed_loop) break
    }

    print(paste0("iteration = ",k))
    return(corr_vec)
}
corr_vec_4 <- run_eachW(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
saveRDS(corr_vec_4, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/corr_vec_4.rds')

corr_vec_1 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/corr_vec_1.rds')
corr_vec_2 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/corr_vec_2.rds')
corr_vec_3 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/corr_vec_3.rds')
corr_vec_4 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/corr_vec_4.rds')
corr_vec_5 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/corr_vec_5.rds')
corr_vec_6 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/corr_vec_6.rds')

corr_vec_lst <- list()
for(i in 1:6) {
    corr_vec_lst[[i]] <- data.frame(
    x = 1:300, 
    y = eval(parse(text=paste0('corr_vec_', i)))[1:300]
    )   
}
p_corr_lst <- list()
for(i in 1:length(corr_vec_lst)) {
    p_corr_lst[[i]] <- ggplot(data = corr_vec_lst[[i]], mapping = aes(x = x, y = y)) +
    geom_line() +
    labs(x = 'Epoch', y = 'Correlation', title = paste('Scenario', i)) +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14))
}
pll_corr <- plot_grid(plotlist = p_corr_lst, ncol=3)
title <- ggdraw() + draw_label("Correlation between true W and estimated W", fontface='bold', size = 16)
pll <- plot_grid(title, pll_corr, ncol=1, rel_heights=c(0.1, 1))
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/Plots/w_corr_v1.png', pll, width = 15, height = 10)
# png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_w_corr.png") 
# plot(corr_vec, type = 'l', xlab = 'Simulation', ylab = 'Correlation', main = 'Correlation between W and W^hat (Senario 1)')
# dev.off() 

## Loss
loss_4 <- rst$L.save
saveRDS(loss_4, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/loss_4.rds')

loss_1 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/loss_1.rds')
loss_2 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/loss_2.rds')
loss_3 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/loss_3.rds')
loss_4 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/loss_4.rds')
loss_5 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/loss_5.rds')
loss_6 <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/loss_6.rds')

loss_lst <- list()
for(i in 1:6) {
    loss_lst[[i]] <- data.frame(
        x = 1:100,
        y = eval(parse(text=paste0('loss_', i)))[1:100]
    )
}
plotlist <- list()
for(i in 1:length(loss_lst)) {
    plotlist[[i]] <- ggplot(data = loss_lst[[i]], mapping = aes(y = y, x = x)) + 
    geom_point() +
    geom_line() +
    scale_y_continuous(trans='log2', ,
    breaks = trans_breaks("log2", function(x) 2^x),
    labels = trans_format("log2", math_format(2^.x))) +
    labs(title=paste("Scenario", i), x ="Epoch", y = "Log loss") +
    theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
}
pll_loss <- plot_grid(plotlist = plotlist, ncol=3)
title_loss <- ggdraw() + draw_label("Log loss function", fontface='bold', size = 16)
pll_l <- plot_grid(title_loss, pll_loss, ncol=1, rel_heights=c(0.1, 1))
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/Plots/log_loss_v1.png', pll_l, width = 15, height = 10)

# png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_loss.png")
# plot(rst$L.save, xlab = 'Epoch', ylab = 'Loss', main = 'Loss (Senario 1)')
# dev.off() 