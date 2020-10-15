set.seed(2020)
source('run_opt.R')

rst<-run(X,1,1,1,AS,A0,D,K,N)
plot(as.vector(W),as.vector(rst$W))
plot(as.vector(true.H),as.vector(rst$H))