# Original file to visualize VPN gradients in the lobula

library(tidyverse)
library(natverse)
library(RColorBrewer)

options(stringsAsFactors = F)


# 0. Load and Connect
con0 <- readRDS("data/hemibrain_con0.rds")
nlist0 <- readRDS("data/nlist1.rds")

# 1. Plot LC-gradients
grad0 <- con0 %>%
  filter(pre.type == "LC4", post.type == "DNp02") %>%
  group_by(pre.bodyID, pre.type, post.type) %>%
  count() %>%
  arrange(desc(n)) %>%
  pull(pre.bodyID)
pal0 <- colorRampPalette(brewer.pal(9, "RdBu"))(length(grad0))

plist0 <- nlist0[as.character(grad0)]
nclear3d(plotengine = "plotly")
plot3d(plist0, col = pal0, plotengine = "plotly")
