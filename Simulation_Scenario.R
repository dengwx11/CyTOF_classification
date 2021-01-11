## senario 1
set.seed(2020)
K = 5 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 10
small_w_mean = 0.5
small_tau_w = 10
p.0 = 0.1
q.0 = 0.2
p.neg1 = 0.05
q.neg1 = 0.1
mean_var_ratio = 10
corr = 0.5
prob_k = c(1,2,2,3,3)
#### rst<-run(X,0.4,.5,3.5,2,AS,A0,D,K,N, epsilon = 10^(-4))
#### rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=150)




## senario 2
set.seed(2020)
K = 8 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 8
small_w_mean = 1
small_tau_w = 8
p.0 = 0.15
q.0 = 0.3
p.neg1 = 0.1
q.neg1 = 0.15
mean_var_ratio = 5
corr = 0.2
#prob_k = c(1,2,2,3,3)
prob_k = c(3,1,2,3,3,1,2,1)
### rst <- run(X,1,10,30,10,AS,A0,D,K,N, epsilon = 10^(-3))
### rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=150)

## senario 3
set.seed(2020)
K = 8 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 8
small_w_mean = 1
small_tau_w = 8
p.0 = 0.15
q.0 = 0.3
p.neg1 = 0.1
q.neg1 = 0.15
mean_var_ratio = 3
corr = 0.3
#prob_k = c(1,2,2,3,3)
prob_k = c(3,1,2,3,3,1,2,1)
### rst<-run(X,1,15,40,10,AS,A0,D,K,N, epsilon = 10^(-3))
### rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=150)

## senario 4
set.seed(2020)
K = 8 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 7
small_w_mean = 1
small_tau_w = 7
p.0 = 0.2
q.0 = 0.4
p.neg1 = 0.1
q.neg1 = 0.2
mean_var_ratio = 3
corr = 0.2
#prob_k = c(1,2,2,3,3)
prob_k = c(3,1,2,3,3,1,2,1)
### rst<-run(X,3,45,60,10,AS,A0,D,K,N, epsilon = 10^(-3))


## senario 5
set.seed(2020)
K = 8 # cell types number ## K could be larger
D = 10 # surface markers number
N = 2000 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

## simulation parameters:
iteration = 100
big_w_mean = 2
big_tau_w = 5
small_w_mean = 1
small_tau_w = 5
p.0 = 0.2
q.0 = 0.4
p.neg1 = 0.1
q.neg1 = 0.2
mean_var_ratio = 2
corr = 0.2
#prob_k = c(1,2,2,3,3)
prob_k = c(3,1,2,3,3,1,2,1)
### rst<-run(X,3,45,60,10,AS,A0,D,K,N, epsilon = 10^(-3))
