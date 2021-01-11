set.seed(2020)
source('run_para.R')
library(ggplot2)
library(dplyr)
#library(tidyverse)
theme_set(theme_bw(base_size=16))

rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=3,lambda1.on=T,lambda2.on=T)
rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)

compare1.true = (2*rst.para$para$lambda1)/rst.para$para$mu
compare2.true = (2*rst.para$para$lambda2)/rst.para$para$mu
compare1.predicted <- W/AS
compare2.predicted <- W/A0    

#hist(compare1.predicted,30)
#hist(compare2.predicted,30)



##### ggplot2

### scenario 1
## 2*lambda1/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df <- summarize(compare1.df,median=median(predicted))
med_ratio_df$true <- compare1.true
print(med_ratio_df)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare1.density.s1.pdf')
  ggplot(compare1.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-11,13)+
  geom_vline(data = med_ratio_df, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and AS (Scenario 1)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df$median+3.7, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=0.3))+
  geom_text(aes(x=med_ratio_df$true-2, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=.3))
dev.off()
## 2*lambda2/mu
compare2.df <- data.frame(predicted = as.vector(compare2.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df.2 <- summarize(compare2.df,median=median(predicted))
med_ratio_df.2$true <- compare2.true
print(med_ratio_df.2)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare2.density.s1.pdf')
  ggplot(compare2.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-3,4.5)+
  geom_vline(data = med_ratio_df.2, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df.2, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and A0 (Scenario 1)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df.2$median+1, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df.2$true-1, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=.9))
dev.off()

### scenario 2
## 2*lambda1/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df <- summarize(compare1.df,median=median(predicted))
med_ratio_df$true <- compare1.true
print(med_ratio_df)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare1.density.s2.pdf')
  ggplot(compare1.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-4,7.5)+
  geom_vline(data = med_ratio_df, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and AS (Scenario 2)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df$median+2, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=1.4))+
  geom_text(aes(x=med_ratio_df$true-1.2, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=1.4))
dev.off()
## 2*lambda2/mu
compare2.df <- data.frame(predicted = as.vector(compare2.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df.2 <- summarize(compare2.df,median=median(predicted))
med_ratio_df.2$true <- compare2.true
print(med_ratio_df.2)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare2.density.s2.pdf')
  ggplot(compare2.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-3,4.5)+
  geom_vline(data = med_ratio_df.2, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df.2, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and A0 (Scenario 2)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df.2$median+1, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df.2$true-1, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=.9))
dev.off()


### scenario 3
## (2*lambda1)/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df <- summarize(compare1.df,median=median(predicted))
med_ratio_df$true <- compare1.true
print(med_ratio_df)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare1.density.s3.pdf')
  ggplot(compare1.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(min(compare1.predicted)-.5,max(compare1.predicted)+.5)+
  geom_vline(data = med_ratio_df, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and AS (Scenario 3)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df$median+4, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df$true-2.5, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=1))
dev.off()
## (2*lambda2)/mu
compare2.df <- data.frame(predicted = as.vector(compare2.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df.2 <- summarize(compare2.df,median=median(predicted))
med_ratio_df.2$true <- compare2.true
print(med_ratio_df.2)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare2.density.s3.pdf')
  ggplot(compare2.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-.5,4)+
  geom_vline(data = med_ratio_df.2, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df.2, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and A0 (Scenario 3)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df.2$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df.2$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df.2$true-.5, label=paste0("True ratio\n",round(med_ratio_df.2$true,3)), y=.9))
dev.off()

### scenario 4
## (2*lambda1)/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df <- summarize(compare1.df,median=median(predicted))
med_ratio_df$true <- compare1.true
print(med_ratio_df)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare1.density.s4.pdf')
  ggplot(compare1.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(min(compare1.predicted)-.5,max(compare1.predicted)+.5)+
  geom_vline(data = med_ratio_df, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and AS (Scenario 4)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df$median+3, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df$true-2, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=1))
dev.off()
## (2*lambda2)/mu
compare2.df <- data.frame(predicted = as.vector(compare2.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df.2 <- summarize(compare2.df,median=median(predicted))
med_ratio_df.2$true <- compare2.true
print(med_ratio_df.2)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare2.density.s4.pdf')
  ggplot(compare2.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-.5,4)+
  geom_vline(data = med_ratio_df.2, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df.2, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and A0 (Scenario 4)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df.2$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df.2$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df.2$true-.5, label=paste0("True ratio\n",round(med_ratio_df.2$true,3)), y=.9))
dev.off()

### scenario 5
## (2*lambda1)/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df <- summarize(compare1.df,median=median(predicted))
med_ratio_df$true <- compare1.true
print(med_ratio_df)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare1.density.s5.pdf')
  ggplot(compare1.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(min(compare1.predicted)-.5,max(compare1.predicted)+.5)+
  geom_vline(data = med_ratio_df, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and AS (Scenario 5)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df$median+3, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df$true-2, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=1))
dev.off()
## (2*lambda2)/mu
compare2.df <- data.frame(predicted = as.vector(compare2.predicted),
                            W = as.vector(W), AS = as.vector(AS))
med_ratio_df.2 <- summarize(compare2.df,median=median(predicted))
med_ratio_df.2$true <- compare2.true
print(med_ratio_df.2)
pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare2.density.s5.pdf')
  ggplot(compare2.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-.5,4)+
  geom_vline(data = med_ratio_df.2, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df.2, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and A0 (Scenario 5)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df.2$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df.2$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df.2$true-.5, label=paste0("True ratio\n",round(med_ratio_df.2$true,3)), y=.9))
dev.off()