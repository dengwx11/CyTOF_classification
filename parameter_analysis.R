set.seed(2020)
source('run_para.R')
library(ggplot2)
library(dplyr)
#library(tidyverse)
theme_set(theme_bw(base_size=16))
rst.para<-runOptimalPara(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=50,depth=3,lambda1.on=T,lambda2.on=T)
rst.para <- readRDS('/Users/mac/Desktop/Yale/Hongyu/CyTOF/rst_para/rst_para_6.rds')
rst<-run(X,rst.para$para$lambda1,rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)

compare1.true = (2*rst.para$para$lambda1)/rst.para$para$mu
write.table(compare1.true, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_true_s6.txt')

compare2.true = (2*rst.para$para$lambda2)/rst.para$para$mu
write.table(compare2.true, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_true_s6.txt')

compare1.predicted <- W/AS
compare2.predicted <- W/A0    

#hist(compare1.predicted,30)
#hist(compare2.predicted,30)



##### ggplot2

### scenario 1
## 2*lambda1/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare1.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s1.txt')

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
write.table(compare2.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s1.txt')
                            
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
  geom_text(aes(x=med_ratio_df.2$median+1, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=0.75))+
  geom_text(aes(x=med_ratio_df.2$true-1, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=.75))
dev.off()

### scenario 2
## 2*lambda1/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare1.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s2.txt')

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
write.table(compare2.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s2.txt')                            
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
  geom_text(aes(x=med_ratio_df.2$median+1, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=0.75))+
  geom_text(aes(x=med_ratio_df.2$true-1, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=.75))
dev.off()


### scenario 3
## (2*lambda1)/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare1.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s3.txt')

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
write.table(compare2.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s3.txt')                            

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

### scenario 3.5
## (2*lambda1)/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare1.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s4.txt')

med_ratio_df <- summarize(compare1.df,median=median(predicted))
med_ratio_df$true <- compare1.true
print(med_ratio_df)
# pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare1.density.s3.pdf')
  ggplot(compare1.df,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(min(compare1.predicted)-.5,max(compare1.predicted)+.5)+
  geom_vline(data = med_ratio_df, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Ratio between predicted W and AS (Scenario 4)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df$median+4, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df$true-2.5, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=1))
# dev.off()
## (2*lambda2)/mu
compare2.df <- data.frame(predicted = as.vector(compare2.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare2.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s4.txt')                            

med_ratio_df.2 <- summarize(compare2.df,median=median(predicted))
med_ratio_df.2$true <- compare2.true
print(med_ratio_df.2)
# pdf('/gpfs/ysm/pi/zhao-data/wd262/new_cytof/figure/compare2.density.s3.pdf')
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
# dev.off()

### scenario 4
## (2*lambda1)/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare1.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s5.txt')

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
       subtitle="Ratio between predicted W and AS (Scenario 5)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df$median+3, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df$true-2, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=1))
dev.off()
## (2*lambda2)/mu
compare2.df <- data.frame(predicted = as.vector(compare2.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare2.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s5.txt')

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
       subtitle="Ratio between predicted W and A0 (Scenario 5)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df.2$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df.2$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df.2$true-.5, label=paste0("True ratio\n",round(med_ratio_df.2$true,3)), y=.9))
dev.off()

### scenario 5
## (2*lambda1)/mu
compare1.df <- data.frame(predicted = as.vector(compare1.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare1.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s6.txt')

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
       subtitle="Ratio between predicted W and AS (Scenario 6)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df$median+3, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df$true-2, label=paste0("True ratio\n",round(med_ratio_df$true,3)), y=1))
dev.off()
## (2*lambda2)/mu
compare2.df <- data.frame(predicted = as.vector(compare2.predicted),
                            W = as.vector(W), AS = as.vector(AS))
write.table(compare2.df, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s6.txt')

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
       subtitle="Ratio between predicted W and A0 (Scenario 6)")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df.2$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df.2$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df.2$true-.5, label=paste0("True ratio\n",round(med_ratio_df.2$true,3)), y=.9))
dev.off()

### Grid plot
compare1_df_s1 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s1.txt')
compare2_df_s1 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s1.txt')
compare1_df_s2 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s2.txt')
compare2_df_s2 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s2.txt')
compare1_df_s3 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s3.txt')
compare2_df_s3 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s3.txt')
compare1_df_s4 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s4.txt')
compare2_df_s4 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s4.txt')
compare1_df_s5 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s5.txt')
compare2_df_s5 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s5.txt')
compare1_df_s6 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_df_s6.txt')
compare2_df_s6 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_df_s6.txt')

compare1_true_s1 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_true_s1.txt')
compare2_true_s1 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_true_s1.txt')
compare1_true_s2 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_true_s2.txt')
compare2_true_s2 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_true_s2.txt')
compare1_true_s3 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_true_s3.txt')
compare2_true_s3 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_true_s3.txt')
compare1_true_s4 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_true_s4.txt')
compare2_true_s4 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_true_s4.txt')
compare1_true_s5 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_true_s5.txt')
compare2_true_s5 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_true_s5.txt')
compare1_true_s6 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare1_true_s6.txt')
compare2_true_s6 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/parameter_analysis/compare2_true_s6.txt')

med_ratio_df1_s1 <- summarize(compare1_df_s1,median=median(predicted))
med_ratio_df1_s1$true <- unlist(compare1_true_s1)
med_ratio_df2_s1 <- summarize(compare2_df_s1,median=median(predicted))
med_ratio_df2_s1$true <- unlist(compare2_true_s1)
med_ratio_df1_s2 <- summarize(compare1_df_s2,median=median(predicted))
med_ratio_df1_s2$true <- unlist(compare1_true_s2)
med_ratio_df2_s2 <- summarize(compare2_df_s2,median=median(predicted))
med_ratio_df2_s2$true <- unlist(compare2_true_s2)
med_ratio_df1_s3 <- summarize(compare1_df_s3,median=median(predicted))
med_ratio_df1_s3$true <- unlist(compare1_true_s3)
med_ratio_df2_s3 <- summarize(compare2_df_s3,median=median(predicted))
med_ratio_df2_s3$true <- unlist(compare2_true_s3)
med_ratio_df1_s4 <- summarize(compare1_df_s4,median=median(predicted))
med_ratio_df1_s4$true <- unlist(compare1_true_s4)
med_ratio_df2_s4 <- summarize(compare2_df_s4,median=median(predicted))
med_ratio_df2_s4$true <- unlist(compare2_true_s4)
med_ratio_df1_s5 <- summarize(compare1_df_s5,median=median(predicted))
med_ratio_df1_s5$true <- unlist(compare1_true_s5)
med_ratio_df2_s5 <- summarize(compare2_df_s5,median=median(predicted))
med_ratio_df2_s5$true <- unlist(compare2_true_s5)
med_ratio_df1_s6 <- summarize(compare1_df_s6,median=median(predicted))
med_ratio_df1_s6$true <- unlist(compare1_true_s6)
med_ratio_df2_s6 <- summarize(compare2_df_s6,median=median(predicted))
med_ratio_df2_s6$true <- unlist(compare2_true_s6)

was_lst <- list()
was_lst[[1]] <- ggplot(compare1_df_s1,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-11,13)+
  geom_vline(data = med_ratio_df1_s1, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df1_s1, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 1")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df1_s1$median+5, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df1_s1$median,3)), y=.3))+
  geom_text(aes(x=med_ratio_df1_s1$true-3, label=paste0("True Ratio\n",round(med_ratio_df1_s1$true,3)), y=.3))
was_lst[[2]] <- ggplot(compare1_df_s2,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-4,7.5)+
  geom_vline(data = med_ratio_df1_s2, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df1_s2, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 2")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df1_s2$median+2, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df1_s2$median,3)), y=1.4))+
  geom_text(aes(x=med_ratio_df1_s2$true-1.2, label=paste0("True Ratio\n",round(med_ratio_df1_s2$true,3)), y=1.4))
was_lst[[3]] <- ggplot(compare1_df_s3,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(min(compare1_df_s3$predicted)-.5,max(compare1_df_s3$predicted)+.5)+
  geom_vline(data = med_ratio_df1_s3, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df1_s3, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 3")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df1_s3$median+6, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df1_s3$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df1_s3$true-4, label=paste0("True Ratio\n",round(med_ratio_df1_s3$true,3)), y=1))
was_lst[[4]] <- ggplot(compare1_df_s4,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(min(compare1_df_s4$predicted)-.5,max(compare1_df_s4$predicted)+.5)+
  geom_vline(data = med_ratio_df1_s4, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df1_s4, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 4")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df1_s4$median+4, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df1_s4$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df1_s4$true-2.5, label=paste0("True Ratio\n",round(med_ratio_df1_s4$true,3)), y=1))
was_lst[[5]] <- ggplot(compare1_df_s5,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(min(compare1_df_s5$predicted)-.5,max(compare1_df_s5$predicted)+.5)+
  geom_vline(data = med_ratio_df1_s5, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df1_s5, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 5")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df1_s5$median+3.5, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df1_s5$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df1_s5$true-2.5, label=paste0("True Ratio\n",round(med_ratio_df1_s5$true,3)), y=1))
was_lst[[6]] <- ggplot(compare1_df_s6,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(min(compare1_df_s6$predicted)-.5,max(compare1_df_s6$predicted)+.5)+
  geom_vline(data = med_ratio_df1_s6, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df1_s6, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 6")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df1_s6$median+3.5, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df1_s6$median,3)), y=1))+
  geom_text(aes(x=med_ratio_df1_s6$true-2.5, label=paste0("True Ratio\n",round(med_ratio_df1_s6$true,3)), y=1))
pll_was <- plot_grid(plotlist = was_lst, ncol=3)
title_was <- ggdraw() + draw_label("Ratio between predicted W and AS", fontface='bold', size = 16)
pll_was <- plot_grid(title_was, pll_was, ncol=1, rel_heights=c(0.1, 1))
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/Plots/para_w_as.png', pll_was, width = 15, height = 10)

wa0_lst <- list()
wa0_lst[[1]] <- ggplot(compare2_df_s1,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-3,4.5)+
  geom_vline(data = med_ratio_df2_s1, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df2_s1, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 1")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df2_s1$median+1.5, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df2_s1$median,3)), y=0.75))+
  geom_text(aes(x=med_ratio_df2_s1$true-1, label=paste0("True Ratio\n",round(med_ratio_df2_s1$true,3)), y=.75))
wa0_lst[[2]] <- ggplot(compare2_df_s2,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-3,4.5)+
  geom_vline(data = med_ratio_df2_s2, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df2_s2, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 2")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df2_s2$median+1.5, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df2_s2$median,3)), y=0.75))+
  geom_text(aes(x=med_ratio_df2_s2$true-1, label=paste0("True Ratio\n",round(med_ratio_df2_s2$true,3)), y=.75))
wa0_lst[[3]] <- ggplot(compare2_df_s3,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-.5,4)+
  geom_vline(data = med_ratio_df2_s3, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df2_s3, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 3")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df2_s3$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df2_s3$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df2_s3$true-.5, label=paste0("True Ratio\n",round(med_ratio_df2_s3$true,3)), y=.9))
wa0_lst[[4]] <- ggplot(compare2_df_s4,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-.5,4)+
  geom_vline(data = med_ratio_df2_s4, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df2_s4, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 4")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df2_s4$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df2_s4$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df2_s4$true-.5, label=paste0("True Ratio\n",round(med_ratio_df2_s4$true,3)), y=.9))
wa0_lst[[5]] <- ggplot(compare2_df_s5,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-.5,4)+
  geom_vline(data = med_ratio_df2_s5, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df2_s5, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 5")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df2_s5$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df2_s5$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df2_s5$true-.5, label=paste0("True Ratio\n",round(med_ratio_df2_s5$true,3)), y=.9))
wa0_lst[[6]] <- ggplot(compare2_df_s6,aes(x=predicted)) +
  geom_density(fill="dodgerblue", alpha=0.5)+ 
 # scale_x_log10()+
  xlim(-.5,4)+
  geom_vline(data = med_ratio_df2_s6, aes(xintercept = median), size=1.5,color="red")+
  geom_vline(data = med_ratio_df2_s6, aes(xintercept = true), size=1.5,color="orange")+
  labs(x= "Predicted Ratio",
       subtitle="Scenario 6")+
  theme(legend.position="bottom")+
  geom_text(aes(x=med_ratio_df2_s6$median+.8, label=paste0("Predicted Ratio\nMedian\n",round(med_ratio_df2_s6$median,3)), y=0.9))+
  geom_text(aes(x=med_ratio_df2_s6$true-.5, label=paste0("True Ratio\n",round(med_ratio_df2_s6$true,3)), y=.9))
pll_wa0 <- plot_grid(plotlist = wa0_lst, ncol=3)
title_wa0 <- ggdraw() + draw_label("Ratio between predicted W and A0", fontface='bold', size = 16)
pll_wa0 <- plot_grid(title_wa0, pll_wa0, ncol=1, rel_heights=c(0.1, 1))
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/Plots/para_w_a0.png', pll_wa0, width = 15, height = 10)
