source("Simulation.R")
set.seed(2020)



gamma <- simulate_gamma(100, K, D)
A0 <- gamma$gamma
W <- simulate_w(D,K,A0)
AS <- simulate_AS(D,K,W,corr=0.9)
plot(as.vector(w),as.vector(AS), xlab = 'w', ylab = 'AS')
cor(as.vector(w),as.vector(AS))
cor(as.vector(w),as.vector(AS),method='spearman')
label.output <- simulate_label(D,N,K)
X <- simulate_x(D,N,w,label.output$label)
true.H <- label.output$H

rownames(X) <- c(1:D)
colnames(X) <- c(1:N)
seur <- CreateSeuratObject(X)
seur@meta.data$label <- label.output$label[1,]


Idents(seur) <- 'label'
RidgePlot(seur, slot = 'counts',features = rownames(seur))
