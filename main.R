set.seed(2020)
source('run_opt.R')

K = 5 # cell types number ## K could be larger
D = 10 # surface markers number
N = 500 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## lambda1 and lambda2 could be found by screening on |X-WH|_F

## full penalization
rst<-run(X,0.5,0.5,1,AS,A0,D,K,N)
plot(as.vector(W),as.vector(rst$W))
plot(as.vector(true.H),as.vector(rst$H))

## without AS
rst<-run(X,0,1,1,AS,A0,D,K,N)
plot(as.vector(W),as.vector(rst$W))
plot(as.vector(true.H),as.vector(rst$H))

## without A0
rst<-run(X,1,0,1,AS,A0,D,K,N)
plot(as.vector(W),as.vector(rst$W))
plot(as.vector(true.H),as.vector(rst$H))

## without any penalization
rst<-run(X,0,0,1,AS,A0,D,K,N)
plot(as.vector(W),as.vector(rst$W))
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


predict<-function(h){
    return(which(h == max(h)))
}
celltype_pred <- apply(rst$H, 2, predict)
celltype_pred <- as.character(celltype_pred)