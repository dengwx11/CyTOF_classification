source("Simulation.R")
source('run_para.R')
library(ggplot2)

loops <- 100
corr_vec <- c()
i <- 0
while(i < loops) {
    ## generate dataset: X, AS, A0, W(ground truth), true.label
    gamma <- simulate_gamma(iteration, K, D)
    A0 <- gamma$gamma
    W <- simulate_w(D,K,A0, p.0 = p.0, q.0=q.0, p.neg1= p.neg1, q.neg1 = q.neg1, 
                big_w_mean=big_w_mean, big_tau_w = big_tau_w, 
                small_w_mean=small_w_mean, small_tau_w = small_tau_w)
    AS <- simulate_AS(D,K,W,corr=corr)
    label.output <- simulate_label(D,N,K,prob_k)
    X <- simulate_x(D,N,W,label.output$label, mean_var_ratio = mean_var_ratio)
    true.H <- label.output$H

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
    corr_vec <- c(corr_vec, cor(as.vector(rst$W), as.vector(W)))
    print(corr_vec)
    print(i)
    i <- i + 1
}

write.table(corr_vec, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/w_corr_vec.txt')

png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_w_corr.png") 
plot(corr_vec, type = 'l', xlab = 'Simulation', ylab = 'Correlation', main = 'Correlation between W and W^hat (Senario 1)')
dev.off() 

## Loss
png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_log_loss.png")
plot(rst$L.save, log="y", xlab = 'Epoch', ylab = 'Log loss', main = 'Loss (Senario 1 log scale)')
dev.off() 

png("/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_loss.png")
plot(rst$L.save, xlab = 'Epoch', ylab = 'Loss', main = 'Loss (Senario 1)')
dev.off() 