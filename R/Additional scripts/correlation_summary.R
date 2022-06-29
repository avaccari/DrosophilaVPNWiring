library(corrplot)
library(psych)

LC4 <- read.csv("data/LC4_3-way-cor.csv")
LPLC2 <- read.csv("data/LPLC2_3-way-cor.csv")

# 009900 green
# 993399 purple

# Plot correlation between two parameters
ggplot(LPLC2, aes(x = Glomerulus, y = RF.Map)) +
  geom_point(size = 2, col = "#009900", fill = "#009900", alpha = 0.5) +
  # geom_smooth(method="loess",formula= 'y ~ x', col="red", se=TRUE)+
  # ggrepel::geom_label_repel(data = subset(LPLC2, row == "PVLP071" & column == "PVLP076"), aes(label = paste0(row, "::", column))) +
  theme_classic() +
  ylab("RF Map, um") +
  xlab("Synaptic Topography, um") +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 13, angle = 0),
    axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
    plot.title = element_text(face = "bold", color = "blue", size = 15),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    aspect.ratio = 1
  ) +
  # geom_vline(xintercept =0, colour="black", linetype = "dashed", lwd=1.5)+
  # annotate('text', size=5, col="red",
  #        label=toString(paste(" r =", round(pc, 2))),
  #       x=0.7 * max(merge$n.x.norm),
  #      y=0.7 * max(merge$n.y.norm),
  #     size=5)+
  xlim(0, 10) +
  ylim(0, 40)




# Visualize correlation between all 3 parameters

col2 <- colorRampPalette(c("#0033CC", "#CCCCCC", "#CC0033"))

y <- cor(LPLC2[3:5])
corrplot(y, type = "upper", method = "ellipse", order = "alpha", diag = FALSE, tl.srt = 0, addCoef.col = "white", col = col2(10), tl.col = "white", tl.cex = 1, tl.offset = 1, tl.pos = "td")

pairs.panels(LPLC2[3:5],
  smooth = FALSE, density = TRUE, lm = FALSE, method = "spearman", jiggle = TRUE, hist.col = "#CC0033", rug = TRUE, breaks = 10,
  cex.cor = 0.5, smoother = FALSE, stars = TRUE, ci = TRUE, ellipses = FALSE
)

corr <- cor.test(c(LC4$RF.Map), c(LC4$Connectivity), method = "spearman")
