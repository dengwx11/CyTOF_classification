set.seed(2020)
source('run_para.R')
library(ggplot2)
library(tidyverse)
theme_set(theme_bw(base_size=16))
library(nnls)
#library(dplyr)

rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=2,lambda1.on=T,lambda2.on=T)
rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)

compare1.true = 2*rst.para$para$lambda1/rst.para$para$mu 
compare2.true = 2*rst.para$para$lambda2/rst.para$para$mu
compare1.predicted <- W/AS
compare2.predicted <- W/A0    

hist(compare1.predicted,30)
hist(compare2.predicted,30)

compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))


summarize(compare1.df,median=median(predicted))

  ggplot(compare1.df,aes(x=predicted)) +
  geom_density(alpha=0.3,size=1)+ 
  scale_x_log10()+
  geom_vline(data = median(compare1.predicted), aes(xintercept = median), size=1.5)+
  geom_vline(data = compare1.true, aes(xintercept = compare1.true), size=1.5)+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and AS")+
  theme(legend.position="bottom")