library(ggplot2)

Frequency <- read.csv("./Frequency_edited.csv")
Frequency <- data.frame(Frequency)
Frequency$Gene_Regions <- factor(Frequency$Gene_Regions, levels = Frequency$Gene_Regions)
My_Theme = theme(axis.title.x = element_text(size = 18),
                 axis.text.x = element_text(size = 16, angle = 90),
                 axis.title = element_text(size = 16))
Jun_WT <- ggplot(data = Frequency, aes(x= Gene_Regions, y=Jun_Delfcpos0.01_Freq)) +
  geom_bar(stat = "identity", fill="black" ) +
  theme_minimal() +
  My_Theme +
  ggtitle("peaks upregulated in Fbw7-/-") +
  xlab("Gene_Regions") +
  ylab("Percentage of peaks")  
Jun_WT
