library(tidyverse)
library(ggrepel)

# Read in Hct116 RNA-Seq data
Hct_RNAseq <- read.csv("/fh/fast/clurman_b/user/hthirima/RNAseq_Hct116_FW/Analysis/DGE/RNAseq_Hct_FW_proteincoding.csv")

## tab3 : marks genes that are differentially expressed in Del or R
tab3 <- Hct_RNAseq %>%
  mutate(sig = ifelse(Del_vs_WT.FDR < 0.05 & abs(Del_vs_WT.log2FC) > log2(1.5) |
                        R_vs_WT.FDR < 0.05 & abs(R_vs_WT.log2FC) > log2(1.5), 
                      'yes', 'no'))

tab3$genelabels <- ""
tab3$genelabels[tab3$gene_name == "KCNQ5"] <- T
tab3$genelabels[tab3$gene_name == "GNG11"] <- T
tab3$genelabels[tab3$gene_name == "MAML2"] <- T
tab3$genelabels[tab3$gene_name == "ITGA2"] <- T
tab3$genelabels[tab3$gene_name == "CIITA"] <- T
#tab3$genelabels[tab3$gene_name == "NEDD9"] <- T
tab3

# this scatter plot shows genes that are significanlty differentially expressed in Del OR R. (FDR < 0.05 and log2FC(1.5))

tab3 %>%
  ggplot(aes(x = R_vs_WT.log2FC, y = Del_vs_WT.log2FC,
             alpha = sig, fill = sig)) +
  geom_point(shape = 21, size = 1.5) +
  theme_classic() +
  scale_alpha_manual(values = c('yes' = 1, 'no' = 0.05)) +
  scale_fill_manual(values = c('yes' = 'grey60', 'no' = 'lightgrey')) +
  labs(x = 'log2FoldChange(Fbw7R/+ / WT)', y = 'log2FoldChange(Fbw7-/- / WT)') +
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), linetype = 'dashed') +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 'dashed') +
  theme(legend.position = 'none') +
  geom_point(data = tab3[tab3$genelabels == T,], color = "red") +
  geom_text_repel(aes(x = R_vs_WT.log2FC, y = Del_vs_WT.log2FC, label = ifelse(genelabels == T, tab3$gene_name, "")), color = "red", hjust=-20,vjust=-3) +
  xlim(-12,12) +
  ylim(-12,12)


 

