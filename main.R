set.seed(2020)
source('run_para.R')
library(ggplot2)
library(nnls)



## lambda1 and lambda2 could be found by screening on |X-WH|_F


# ## full penalization
# rst<-run(X,1,65,60,10,AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=1000)
# rst<-run(X,0,10,30,10,AS,A0,D,K,N, epsilon = 10^(-3))
# rst<-run(X,1,0,30,10,AS,A0,D,K,N, epsilon = 10^(-3))
# rst<-run(X,0,0,30,10,AS,A0,D,K,N, epsilon = 10^(-3))


rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=2,lambda1.on=T,lambda2.on=T)
saveRDS(rst.para, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/rst_para/rst_para_6.rds')
rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
rst.para.woas<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=2,lambda1.on=F,lambda2.on=T)
rst.woas<-run(X,rst.para.woas$para$lambda1,rst.para.woas$para$lambda2,rst.para.woas$para$mu,rst.para.woas$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
rst.para.woa0<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=2,lambda1.on=T,lambda2.on=F)
rst.woa0<-run(X,rst.para.woa0$para$lambda1,rst.para.woa0$para$lambda2,rst.para.woa0$para$mu,rst.para.woa0$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
rst.para.woa0as<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=2,lambda1.on=F,lambda2.on=F)           
rst.woa0as<-run(X,rst.para.woa0as$para$lambda1,rst.para.woa0as$para$lambda1,rst.para.woa0as$para$mu,rst.para.woa0as$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)

## NNLS
nnls_comp <- function(x, W) {
    return(nnls(W, x)$x)
}
H_nnls <- apply(X,2,nnls_comp, W = W)

plot(as.vector(W),as.vector(rst$W))
cor(as.vector(W),as.vector(rst$W))
#plot(as.vector(true.H),as.vector(rst$H))



## True postive rate


truth = label.output$label
# cnt = 0
# for(i in 1:length(truth)){
#     cnt = cnt + infer(truth[i], rst$H[,i])
# }

cnt_max = 0
for(i in 1:length(truth)){
    cnt_max = cnt_max + infer_max(truth[i], rst$H[,i])
}  
cnt_max_woas = 0
for(i in 1:length(truth)){
    cnt_max_woas = cnt_max_woas + infer_max(truth[i], rst.woas$H[,i])
} 
cnt_max_woa0 = 0
for(i in 1:length(truth)){
    cnt_max_woa0 = cnt_max_woa0 + infer_max(truth[i], rst.woa0$H[,i])
} 
cnt_max_woa0as = 0
for(i in 1:length(truth)){
    cnt_max_woa0as = cnt_max_woa0as + infer_max(truth[i], rst.woa0as$H[,i])
} 
cnt_max_nnls = 0
for(i in 1:length(truth)){
    cnt_max_nnls = cnt_max_nnls + infer_max(truth[i], H_nnls[,i])
} 

print(cnt_max)
print(cnt_max/N)


## prediction visualization
celltype_pred <- apply(rst$H, 2, predict)
celltype_pred_woas <- apply(rst.woas$H, 2, predict)
celltype_pred_woa0 <- apply(rst.woa0$H, 2, predict)
celltype_pred_woa0as <- apply(rst.woa0as$H, 2, predict)
celltype_pred_nnls <- apply(H_nnls, 2, predict)

#celltype_pred <- as.character(celltype_pred)
df.plot <- data.frame(x = X.umap$layout[,1],y=X.umap$layout[,2],
                            true_label = factor(label.output$label[1,]),
                            pred_label = factor(celltype_pred))
ggplot(df.plot, aes(x=x,y=y,color=true_label)) + geom_point() 

ggplot(df.plot, aes(x=x,y=y,color=pred_label)) + geom_point()                         
plot(X.umap$layout,col=label.output$label)
plot(X.umap$layout,col=celltype_pred)

seur$pred = celltype_pred
seur$pred_woas = celltype_pred_woas
seur$pred_woa0 = celltype_pred_woa0
seur$pred_woa0as = celltype_pred_woa0as
seur$pred_nnls = celltype_pred_nnls

DimPlot(seur, reduction='umap', group.by = 'pred')

adjustedRandIndex(celltype_pred, seur$seurat_clusters)