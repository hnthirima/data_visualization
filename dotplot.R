library(ggplot2)
library(forcats)

GO_terms <- read.csv("/fh/fast/clurman_b/user/hthirima/ManuscriptFigures_Feb2021/Hct116_RNAseq_clustering/Clusters_GOterms_dotplot.csv", header = T)

#class(GO_terms)

GO_Terms=factor(GO_terms$GO_Term, levels = GO_terms$GO_Term[order(GO_terms$logAdjpvalue)])

ggplot(GO_terms) +
  geom_point(aes(x = -log(Adjpvalue), y = GO_Terms, Rowv=NULL, Colv=NULL, color = Cluster, size = -(logAdjpvalue))) +
  theme_classic() +
  theme(text = element_text(size=15), axis.text.y = element_text(colour = "black"))+
  xlab("-log10(Adj.pvalue)") +
  labs(size="-log10(Adj.pvalue)")+
  scale_y_discrete(limits=rev)


dev.off()
