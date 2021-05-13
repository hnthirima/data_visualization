library(tidyverse)
library(dplyr)
library(ggplot2)

# read the data
trans_cts_1 <- read.delim("./RNAseq_CPM.txt")
#trans_cts <- trans_cts_1[-c(1)]
sample_info <- read.csv("./SampleInfo.csv", header = TRUE)
# Create a matrix from our table of counts. (exclude the gene column and coerce to a matrix)
pca_matrix <- trans_cts_1 %>%
  select(-Gene_name) %>%
  as.matrix()

# Assign row names to the matrix
rownames(pca_matrix) <- trans_cts_1$Gene_name

# Transpose the matrix so that rows = samples and columns = variables (genes)
pca_matrix <- t(pca_matrix)

# Perform the PCA
sample_pca <- prcomp(pca_matrix, scale. = TRUE)

head(sample_pca)

pc_scores <- sample_pca$x

variances <- (100*((sample_pca$sdev)^2))/sum((sample_pca$sdev)^2)
variances

#pc_scores %>%
  as_tibble(rownames = "sample") %>%
  ggplot(aes(x=PC1, y=PC2)) + 
  geom_point()

#pc_scores %>%
    as_tibble(rownames = "sample") %>%
    full_join(sample_info, by = "sample") %>%
    ggplot(aes(x=PC1, y=PC2, colour = factor(Condition), shape = factor(Replicate))) + 
    geom_point()
  
#pc_scores %>%
    as_tibble(rownames = "sample") %>%
    full_join(sample_info, by = "sample") %>%
    ggplot(aes(x=PC1, y=PC2, colour = factor(sample), shape = factor(Condition))) + 
    geom_point()
dev.off()

RNAseq_PCA.pdf <- pc_scores %>%
  as_tibble(rownames = "sample") %>%
  full_join(sample_info, by = "sample") %>%
  ggplot(aes(x=PC1, y=PC2, colour = factor(Condition))) + 
  geom_point(size = 3)
RNAseq_PCA.pdf

ggsave("RNAseq_PCA.pdf", device = NULL, dpi = 300)

samplepca_x <- sample_pca$x
write.csv(samplepca_x, file = "./RNAseq_PCA_samplepca_x.csv")
write.csv(samplepca_rotation, file = "./RNAseq_PCA_samplepca_rotation.csv")
