library(mltools)
library(data.table)
library(lsa)
library(aricode)
library(kBET)

truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))
pred_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred)))))
pred_onehot_woas <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred_woas)))))
pred_onehot_woa0 <- as.data.frame(t(one_hot(as.data.table(as.factor(celltype_pred_woa0)))))
celltype_pred_fact <- as.factor(celltype_pred_woa0as)
celltype_pred_fact <- factor(celltype_pred_fact, levels = c("1", "2", "3", "4", "5"))
pred_onehot_woa0as <- as.data.frame(t(one_hot(as.data.table(celltype_pred_fact))))
pred_onehot_louvain <- as.data.frame(t(one_hot(as.data.table(seur$seurat_clusters))))

## Cosine similarity
cos_sim <- mean(mapply(cosine, truth_onehot, as.data.frame(rst$H)))
cos_sim_woas <- mean(mapply(cosine, truth_onehot, as.data.frame(rst.woas$H)))
cos_sim_woa0 <- mean(mapply(cosine, truth_onehot, as.data.frame(rst.woa0$H)))
cos_sim_woa0as <- mean(mapply(cosine, truth_onehot, as.data.frame(rst.woa0as$H)))
cos_sim_louvain <- mean(mapply(cosine, truth_onehot, pred_onehot_louvain))

## NMI
nmi <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot)))
nmi_woas <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_woas)))
nmi_woa0 <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_woa0)))
nmi_woa0as <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_woa0as)))
nmi_louvain<- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot_louvain)))

## ARI
ari <- adjustedRandIndex(celltype_pred, truth)
ari_woas <- adjustedRandIndex(celltype_pred_woas, truth)
ari_woa0 <- adjustedRandIndex(celltype_pred_woa0, truth)
ari_woa0as <- adjustedRandIndex(celltype_pred_woa0as, truth)
ari_louvain <- adjustedRandIndex(as.numeric(seur$seurat_clusters), truth)

## Accuracy
accu <- cnt_max / N
accu_woas <- cnt_max_woas / N
accu_woa0 <- cnt_max_woa0 / N
accu_woa0as <- cnt_max_woa0as / N
accu_louvain <- mean(pred_onehot_louvain == truth)

## Silhouette
pca.data <- list()
pca.data$x <- seur@reductions$pca@cell.embeddings
sil <- batch_sil(pca.data, celltype_pred)
sil_woas <- batch_sil(pca.data, celltype_pred_woas)
sil_woa0 <- batch_sil(pca.data, celltype_pred_woa0)
sil_woa0as <- batch_sil(pca.data, celltype_pred_woa0as)
sil_louvain <- batch_sil(pca.data, celltype_pred_woa0as)

pars <- data.frame(
    names = rep(c('Accuracy', 'ARI', 'Cosine Similarity', 'NMI', 'ASW'), 5),
    vars = c(accu, ari, cos_sim, nmi, sil, 
             accu_woas, ari_woas, cos_sim_woas, nmi_woas, sil_woas,
             accu_woa0, ari_woa0, cos_sim_woa0, nmi_woa0, sil_woa0,
             accu_woa0as, ari_woa0as, cos_sim_woa0as, nmi_woa0as, sil_woa0as,
             accu_louvain, ari_louvain, cos_sim_louvain, nmi_louvain, sil_louvain),
    group = c(rep('Full', 5), rep('w/o AS', 5), rep('w/o A0', 5), rep('w/o A0 AS', 5), rep('Louvain', 5))
)
p <- ggplot(data = pars, aes(x = names, y = vars, fill = group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle('Senario 1') + 
  coord_flip()
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/senario1_benchmark.png', p)



