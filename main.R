set.seed(2020)
source('run_opt.R')
library(ggplot2)




## lambda1 and lambda2 could be found by screening on |X-WH|_F


## full penalization
rst<-run(X,2,100,60,70,AS,A0,D,K,N, epsilon = 10^(-3))
rst<-run(X,0,10,30,10,AS,A0,D,K,N, epsilon = 10^(-3))
rst<-run(X,1,0,30,10,AS,A0,D,K,N, epsilon = 10^(-3))
rst<-run(X,0,0,30,10,AS,A0,D,K,N, epsilon = 10^(-3))
plot(as.vector(W),as.vector(rst$W))
cor(as.vector(W),as.vector(rst$W))
plot(as.vector(true.H),as.vector(rst$H))





## True postive rate
get_truth <- function(true_ct){
    return(which(H_empirical$celltype == true_ct))
}
infer <- function(truth, h){
    return(1*(h[truth]>0))
}
infer_max <- function(truth,h){
    return(1*(h[truth]==max(h)))
}
predict<-function(h){
    return(which(h == max(h)))
}

truth = label.output$label
cnt = 0
for(i in 1:length(truth)){
    cnt = cnt + infer(truth[i], rst$H[,i])
}

cnt_max = 0
for(i in 1:length(truth)){
    cnt_max = cnt_max + infer_max(truth[i], rst$H[,i])
}  
print(c(cnt,cnt_max))
print(cnt_max/N)



celltype_pred <- apply(rst$H, 2, predict)
#celltype_pred <- as.character(celltype_pred)
df.plot <- data.frame(x = X.umap$layout[,1],y=X.umap$layout[,2],
                            true_label = factor(label.output$label[1,]),
                            pred_label = factor(celltype_pred))
ggplot(df.plot, aes(x=x,y=y,color=true_label)) + geom_point() 
ggplot(df.plot, aes(x=x,y=y,color=pred_label)) + geom_point()                         
plot(X.umap$layout,col=label.output$label)
plot(X.umap$layout,col=celltype_pred)