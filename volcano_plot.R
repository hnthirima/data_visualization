Hct_H3K27ac_DelWT <- read.csv("<PATH_TO_FILE>")
Hct_H3K27ac_DelWT <- Hct_H3K27ac_DelWT %>%
  mutate(threshold = factor(case_when(logFC > 0.6 & FDR  < 0.05 ~ "Up in mutant",
                                      logFC < -0.6 & FDR < 0.05 ~ "Down in mutant",
                                      TRUE ~ "Non-diff.")))
ggplot(Hct_H3K27ac_DelWT, aes(x = logFC, y=-log10(FDR))) + 
  geom_point(aes(color=Hct_H3K27ac_DelWT$threshold), size = 0.8) +
  scale_color_manual(name = " ", values = c("Up in mutant" = "red", "Down in mutant"="blue", "Non-diff."="grey"))  +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  theme(legend.text = element_text(size = 20), text = element_text(size = 20)) + 
  xlim(-6,6) +
  labs(title = "H3K27ac in Fbw7-/-")
