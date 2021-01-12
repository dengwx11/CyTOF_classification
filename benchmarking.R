library(mltools)
library(data.table)
library(lsa)
library(aricode)
library(kBET)
library(cowplot)

truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))
pred_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred)))))
pred_onehot_woas <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred_woas)))))
pred_onehot_woa0 <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred_woa0)))))
## For senario 1
# celltype_pred_fact <- as.factor(celltype_pred_woa0as)
# celltype_pred_fact <- factor(celltype_pred_fact, levels = c("1", "2", "3", "4", "5"))
## For senario 3, 3.5
celltype_pred_fact <- as.factor(celltype_pred_woa0as)
celltype_pred_fact <- factor(celltype_pred_fact, levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
pred_onehot_woa0as <- as.data.frame(t(one_hot(as.data.table(celltype_pred_fact))))
pred_onehot_louvain <- as.data.frame(t(one_hot(as.data.table(seur$seurat_clusters))))
pred_onehot_nnls <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred_nnls)))))

## Cosine similarity
cos_sim <- mean(mapply(cosine, truth_onehot, as.data.frame(rst$H)))
cos_sim_woas <- mean(mapply(cosine, truth_onehot, as.data.frame(rst.woas$H)))
cos_sim_woa0 <- mean(mapply(cosine, truth_onehot, as.data.frame(rst.woa0$H)))
cos_sim_woa0as <- mean(mapply(cosine, truth_onehot, as.data.frame(rst.woa0as$H)))
cos_sim_louvain <- mean(mapply(cosine, truth_onehot, pred_onehot_louvain))
cos_sim_nnls <- mean(mapply(cosine, truth_onehot, pred_onehot_nnls))

## NMI
nmi <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot)))
nmi_woas <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_woas)))
nmi_woa0 <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_woa0)))
nmi_woa0as <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_woa0as)))
nmi_louvain<- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_louvain)))
nmi_nnls<- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_nnls)))

## ARI
ari <- adjustedRandIndex(celltype_pred, truth)
ari_woas <- adjustedRandIndex(celltype_pred_woas, truth)
ari_woa0 <- adjustedRandIndex(celltype_pred_woa0, truth)
ari_woa0as <- adjustedRandIndex(celltype_pred_woa0as, truth)
ari_louvain <- adjustedRandIndex(as.numeric(seur$seurat_clusters), truth)
ari_nnls <- adjustedRandIndex(celltype_pred_nnls, truth)

## Accuracy
accu <- cnt_max / N
accu_woas <- cnt_max_woas / N
accu_woa0 <- cnt_max_woa0 / N
accu_woa0as <- cnt_max_woa0as / N
accu_louvain <- mean(pred_onehot_louvain == truth)
accu_nnls <- cnt_max_nnls / N

## Silhouette
pca.data <- list()
pca.data$x <- seur@reductions$pca@cell.embeddings
sil <- batch_sil(pca.data, celltype_pred)
sil_woas <- batch_sil(pca.data, celltype_pred_woas)
sil_woa0 <- batch_sil(pca.data, celltype_pred_woa0)
sil_woa0as <- batch_sil(pca.data, celltype_pred_woa0as)
sil_louvain <- batch_sil(pca.data, as.numeric(seur$seurat_clusters))
sil_nnls <- batch_sil(pca.data, celltype_pred_nnls)

# pars <- data.frame(
#     names = rep(c('Accuracy', 'ARI', 'Cosine Similarity', 'NMI', 'ASW'), 5),
#     vars = c(accu, ari, cos_sim, nmi, sil, 
#              accu_woas, ari_woas, cos_sim_woas, nmi_woas, sil_woas,
#              accu_woa0, ari_woa0, cos_sim_woa0, nmi_woa0, sil_woa0,
#              accu_woa0as, ari_woa0as, cos_sim_woa0as, nmi_woa0as, sil_woa0as,
#              accu_louvain, ari_louvain, cos_sim_louvain, nmi_louvain, sil_louvain),
#     group = c(rep('Full', 5), rep('w/o AS', 5), rep('w/o A0', 5), rep('w/o A0 AS', 5), rep('Louvain', 5))
# )

pars_s4 <- data.frame(
    models = c(rep('Full', 5), rep('w/o AS', 5), rep('w/o A0', 5), rep('w/o A0 AS', 5), rep('Louvain', 5), rep('NNLS', 5)),
    vars = c(accu, ari, cos_sim, nmi, sil, 
             accu_woas, ari_woas, cos_sim_woas, nmi_woas, sil_woas,
             accu_woa0, ari_woa0, cos_sim_woa0, nmi_woa0, sil_woa0,
             accu_woa0as, ari_woa0as, cos_sim_woa0as, nmi_woa0as, sil_woa0as,
             accu_louvain, ari_louvain, cos_sim_louvain, nmi_louvain, sil_louvain,
             accu_nnls, ari_nnls, cos_sim_nnls, nmi_nnls, sil_nnls),
    metric = rep(c('Accuracy', 'ARI', 'Cosine Similarity', 'NMI', 'ASW'), 6)
)
pars_s6$models <- factor(pars_s6$models,levels = c("NNLS", "Louvain", "w/o A0 AS", "w/o A0", "w/o AS", "Full"))
pars_s6$metric <- factor(pars_s6$metric,levels = c('ASW', 'Cosine Similarity', 'NMI', 'Accuracy', 'ARI'))

write.table(pars_s1, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s1.txt')
write.table(pars_s2, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s2.txt')
write.table(pars_s3, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s3.txt')
write.table(pars_s4, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s4.txt')
write.table(pars_s5, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s5.txt')
write.table(pars_s6, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s6.txt')

pars_s1 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s1.txt')
pars_s2 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s2.txt')
pars_s3 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s3.txt')
pars_s5 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s5.txt')
pars_s6 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/bcmk_s6.txt')

pars_lst <- list(pars_s1, pars_s2, pars_s3, pars_s4, pars_s5, pars_s6)
plotlist <- list()
for(i in 1:length(pars_lst)){
  plotlist[[i]] <- ggplot(data = pars_lst[[i]], aes(x = models, y = vars, fill = metric, width=.7)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(y = 'Values', title = paste('Senario', i)) +
  coord_flip() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14))
}
pll <- plot_grid(plotlist = plotlist, ncol=3)
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/Plots/benchmark_v1.png', pll, width = 15, height = 6)

Plo


