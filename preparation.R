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
    return(which(h == max(h))[1])
}
predict_realdata<-function(h){
    loc = which(h == max(h))[1]
    return(H_est[,1][loc])
}