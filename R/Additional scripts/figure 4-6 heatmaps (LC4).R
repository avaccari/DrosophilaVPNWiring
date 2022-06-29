library(tidyverse)
library(cowplot)

rm(list = ls())

# 1. LC4
load("data/LC4_matrices.Rdata")

## Order
clust0 <- hclust(dist(pcorr_c))
# order0 <- clust0$labels[clust0$order]
order0 <- c(
  "DNp11", "PLP219", "SAD064", "PVLP026", "PVLP137", "PVLP130", "PVLP024", "AVLP452", "AVLP259",
  "PVLP151", "PVLP136", "DNp03", "DNp04", "PVLP021", "PVLP122", "PLP220", "AMMC-A1", "Giant Fiber", "PVLP141", "PVLP010", "PVLP022",
  "LHAD1g1", "PVLP142", "DNp05", "DNp02"
)
labels0 <- cutree(clust0, k = 5)[order0]

## Theme
theme0 <- theme(
  axis.title = element_blank(),
  axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5)
)

## Plot 1
pdata1 <- pcorr_c %>%
  as.data.frame() %>%
  rownames_to_column("cell1") %>%
  gather(key = "cell2", value = "corr1", -cell1)

p1 <- ggplot(pdata1, aes(x = cell1, y = cell2, fill = corr1)) +
  geom_tile(col = "grey20") +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_x_discrete(limits = order0) +
  scale_y_discrete(limits = order0) +
  geom_hline(yintercept = c(5.5), size = 2) +
  geom_vline(xintercept = c(5.5), size = 2) +
  geom_hline(yintercept = c(19.5), size = 2) +
  geom_vline(xintercept = c(19.5), size = 2) +
  theme0

## Plot 2
pdata2 <- dist_mtrx_glom %>%
  as.data.frame() %>%
  rownames_to_column("cell1") %>%
  gather(key = "cell2", value = "dist1", -cell1) %>%
  mutate(dist1 = dist1 * 0.008)

p3 <- ggplot(pdata2, aes(x = cell1, y = cell2, fill = dist1)) +
  geom_tile(col = "grey20") +
  # scale_fill_distiller(palette = "Blues", direction = 1) +
  # scale_fill_gradient2(low = "white", high = "blue") +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_x_discrete(limits = order0) +
  scale_y_discrete(limits = order0) +
  geom_hline(yintercept = c(5.5), size = 2) +
  geom_vline(xintercept = c(5.5), size = 2) +
  geom_hline(yintercept = c(19.5), size = 2) +
  geom_vline(xintercept = c(19.5), size = 2) +
  theme0

## Plot 3
pdata3 <- dist_mtrx_lob %>%
  as.data.frame() %>%
  rownames_to_column("cell1") %>%
  gather(key = "cell2", value = "dist1", -cell1) %>%
  mutate(dist1 = dist1 * 0.008)

# ggplot(pdata0, aes(x = post.type, y = pre.type, fill = ifelse(n > 200, 200, n)))

p2 <- ggplot(pdata3, aes(x = cell1, y = cell2, fill = dist1)) +
  geom_tile(col = "grey20") +
  # scale_fill_distiller(palette = "Oranges", direction = -1) +
  # scale_fill_gradientn(low = "white", high = "blue") +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_x_discrete(limits = order0) +
  scale_y_discrete(limits = order0) +
  geom_hline(yintercept = c(5.5), size = 2) +
  geom_vline(xintercept = c(5.5), size = 2) +
  geom_hline(yintercept = c(19.5), size = 2) +
  geom_vline(xintercept = c(19.5), size = 2) +
  theme0

## Combine
plot_grid(p1, p2, p3, nrow = 1)
