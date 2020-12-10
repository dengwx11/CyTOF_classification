set.seed(2020)
source("getPara.R")




## step 1 : get an initial eta by kkt condition
rst.eta = run(X,1,10,50,50,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=1500)
plot(as.vector(W),as.vector(rst.eta$W))
W.ica = rst.eta$W
H.ica = rst.eta$H
eta.kkt = getEta1(X,W.ica,H.ica)
print(eta.kkt)

## step 2 : get the initial lambdas by screening
lambda_1.sequence = c(10^(-1),10^(0),10^(1),10^(2))
lambda_2.sequence = c(2*10^(-1),2*10^(0),2*10^(1),2*10^(2))
ans.lambda1 = getLambda1(object = 1,lambda_1.sequence,10, 10,eta.kkt,fixed_loop=1000)
lambda_1.init1 = ans.lambda1$lambda_1.init1
lambda_2.init1 = ans.lambda1$lambda_2.init1
print(c(lambda_1.init1,lambda_2.init1,ans.lambda1$ARI))
ans.lambda2 = getLambda1(object = 2,lambda_2.sequence,lambda_1.init1, 10,eta.kkt,fixed_loop=1000)
lambda_1.init2 = ans.lambda2$lambda_1.init1
lambda_2.init2 = ans.lambda2$lambda_2.init1
print(c(lambda_1.init2,lambda_2.init2,ans.lambda2$ARI))

## step 3 : get mu by W reliability
mu.sequence = c(10^(-1),10^(0),10^(1),10^(2))
ans.mu =  getMu1(X, lambda_1.init2,lambda_2.init2,mu.sequence,eta.kkt,epsilon = .01,fixed_loop=1500)
mu.init = ans.mu$mu.init
idx = ans.mu$idx
print(c(mu.init,idx))


## step 4 : get an second initial eta by kkt condition
rst.eta = run(X,lambda_1.init2,lambda_2.init2,mu.init,eta.kkt,AS,A0,D,K,N, fixed_loop=1500)
plot(as.vector(W),as.vector(rst.eta$W))
W.ica = rst.eta$W
H.ica = rst.eta$H
eta.kkt = getEta1(X,W.ica,H.ica)
print(eta.kkt)

## step 5 : select the best lambda
lambda_2.init3 = getLambda2(object = 2,depth=3,lambda_2.init2,lambda_1.init2, mu.init,eta.kkt,fixed_loop=1500)
print(lambda_2.init3)
lambda_1.init3 = getLambda2(object = 1,depth=3,lambda_1.init2,lambda_2.init3, mu.init,eta.kkt,fixed_loop=1500)
print(lambda_1.init3)


## step 6 : select the best mu
mu.init2 =  getMu2(X, depth = 3, lambda_1.init3,lambda_2.init3,mu.init,eta.kkt,epsilon = .01,fixed_loop=1500)
print(c(mu.init2))

## step 7 : select the best lambda
rst.lambda2 = getLambda2(object = 2,depth=3,lambda_2.init3,lambda_1.init2, mu.init2,eta.kkt,fixed_loop=1500)
lambda_2.init4 = rst.lambda2$lambda
print(lambda_2.init4)
rst.lambda1 = getLambda2(object = 1,depth=1,lambda_1.init3,lambda_2.init4, mu.init2,eta.kkt,fixed_loop=1500)
lambda_1.init4 = rst.lambda1$lambda
print(lambda_1.init4)
rst = rst.lambda1$rst

##################################


rst<-run(X,lambda_1.init4,lambda_2.init4,mu.init2,eta.kkt,AS,A0,D,K,N, epsilon = .1,fixed_loop=1500)
celltype_pred <- apply(rst$H, 2, predict)*2+10
adjustedRandIndex(celltype_pred, seur$seurat_clusters)
cor(as.vector(W),as.vector(rst$W))


rst2<-run(X,1,lambda_2.init4,mu.init2,eta.kkt,AS,A0,D,K,N, epsilon = 10^(-3))
celltype_pred <- apply(rst2$H, 2, predict)
adjustedRandIndex(celltype_pred, seur$seurat_clusters)

# k=4
# rst = rst.seq[[k]]
# plot(as.vector(W),as.vector(rst$W))
# hist(as.vector(rst$W),30)
# max(rst$W)

######### 

# dist.matrix <- dist(x = seur@reductions$umap@cell.embeddings)
# clusters <- seur$seurat_clusters
# sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
# median(sil[, 3])
# hist(sil[, 3],30)