library(tidyverse)
library(magrittr)
library(neuprintr)
library(igraph)
library(ggplot2)
library(ggpubr)
library(plotly)
# library(scran)
library(factoextra)
library(RColorBrewer)
library(gplots)
library(cluster)
library(philentropy)
library(dendextend)
library(colorspace)
library(reshape2)
library(scales)
library(VennDiagram)

# Neuprint login info
conn <- neuprint_login(
  server = "https://neuprint.janelia.org/",
  # Use your personal token here
  token = # your token here #,
  dataset = "hemibrain:v1.2.1"
)

# More digits for the seconds in time
options(digits.secs = 6)

# Function for scaling edge arrows by edge width,
# from: https://github.com/jevansbio/igraphhack
source("R/VPN_cluster_analysis/3d_party/igraphplot2.R")

environment(plot.igraph2) <- asNamespace("igraph")
environment(igraph.Arrows2) <- asNamespace("igraph")

###### Misc #####
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

###### SI LC graphs  ######
# takes LC string and list of LC-output matrices
# outputs an igraph object with the LC nodes
# labeled by their three k means clusters and their outputs
get_LC_out_network_collapse <- function(lc, mat_list, syn_thresh) {
  lck <- readRDS(sprintf("output/LC_kmeans/%s_cluster_assignments_k%i.rds", lc, 3))
  # in denotes edges going from LCs (columns) to outputs (rows)
  g <- graph_from_incidence_matrix(mat_list[[lc]][-c(1:3)],
    weighted = TRUE,
    directed = TRUE,
    mode = "in"
  )

  # convert to edgelist
  el <- igraph::as_edgelist(g, names = TRUE) %>% as.data.frame()
  el$weight <- E(g)$weight
  el <- filter(el, weight >= syn_thresh)
  # all LC6->LC6 connections are <5 so they get filtered out
  # mat_list[[lc]][mat_list[[lc]]$name == lc,]

  # Rename LCs by kmeans
  # type.to = replace(type.to, !(is.na(DN.Type.to)), DN.Type.to[which( !(is.na(DN.Type.to)) )])
  el <- left_join(el, lck, by = c("V1" = "bodyid")) %>%
    dplyr::mutate(V1 = replace(
      V1,
      !(is.na(cluster)),
      cluster[which(!(is.na(cluster)))]
    )) %>%
    select(-cluster) %>%
    left_join(lck, by = c("V2" = "bodyid")) %>%
    dplyr::mutate(V2 = replace(
      V2,
      !(is.na(cluster)),
      cluster[which(!(is.na(cluster)))]
    )) %>%
    select(-cluster) %>%
    group_by(V1, V2) %>%
    summarize(weight = sum(weight))

  g <- graph_from_data_frame(el, directed = TRUE)

  # make bipartite
  V(g)[V(g)$name %in% c("1", "2", "3")]$type <- TRUE
  V(g)[!(V(g)$name %in% colnames(mat_list[[lc]]))]$type <- FALSE

  # assign shapes
  V(g)[V(g)$type == TRUE]$shape <- "square"
  V(g)[V(g)$type == FALSE]$shape <- "circle"

  # color by k means cluster
  V(g)$color <- "gray"

  V(g)[V(g)$name == "1"]$color <- "#7DBFA7"
  V(g)[V(g)$name == "2"]$color <- "#EC936A"
  V(g)[V(g)$name == "3"]$color <- "#919FC8"

  # size by degree, capping at 35 size
  V(g)$strength <- strength(g, V(g), mode = "all")
  V(g)[V(g)$strength > 100]$strength <- 100

  #<=2 for LC4, <=1 for LC6
  g <- delete_vertices(g, V(g)[degree(g) <= 1])
  # color edges by the start node
  edge_start <- ends(g, es = E(g), names = F)[, 1]
  edge_col <- V(g)$color[edge_start]

  E(g)[E(g)$weight > 1000]$weight <- 1000
  return(list("graph" = g, "edge_col" = edge_col))
}

# take graph and scaling factor
# returns fruchterman reingold layout matrix
# adjust nodes in the middle (x) to push them left or right
# also spreads out y direction
layout_LC <- function(g, m, adjustx, scaley, adjustlc) {
  set.seed(100000 * as.POSIXlt(Sys.time())$sec)
  LO <- layout_with_fr(g)
  LO <- LO * m # 2.3 for LC4/LPLC1, 1.6 for LC6, 3 for LPLC2
  # get x and y just for outputs
  x <- as.data.frame(LO[-(1:3), ])$V1
  y <- as.data.frame(LO[-(1:3), ])$V2
  # get indices of nodes in the middle and move them to the left or right
  left <- which(x <= mean(x) & x >= (mean(x) - 0.5 * sd(x)))
  right <- which(x >= mean(x) & x <= (mean(x) + 0.5 * sd(x)))
  center <- which(x >= (mean(x) - 0.1 * sd(x)) & x <= (mean(x) + 0.1 * sd(x)))

  x[left] <- x[left] - adjustx
  x[right] <- x[right] + adjustx
  y[center] <- y[center] * scaley
  lcs <- cbind(rep(mean(x), 3), c(LO[1, 2] + adjustlc, LO[2, 2], LO[3, 2] - adjustlc))
  newLO <- rbind(lcs, cbind(x, y))
  return(newLO)
}


### Get PCA object of LCs by their outputs ####
# takes lc, a matrix of LCs as columns and output neurons as rows
# (neuprint_simple_connectivity)
# takes DN, a df of DNs and bodyids
# filters either to just DNs or top outputs above a total synapse count
# rm_self is to filter out LC-LC connections, input lc type in str
# syn_count removes connections below a threshold
lc_getpca <- function(lc, DN = NULL, syn_count = NULL, rm_self = FALSE, str = NULL) {
  lc$output %<>% as.character()
  if (!is.null(DN)) {
    tmp <- data.frame(
      "output" = lc$output,
      "lc_syn" = rowSums(lc[-c(1:3)],
        na.rm = TRUE
      )
    )
    lc <- filter(lc, output %in% tmp[tmp$lc_syn > syn_count, "output"])

    lc <- inner_join(DN, lc, by = c("bodyid" = "output"))
    rownames(lc) <- lc$DN.Type

    if (rm_self) {
      lc <- filter(lc, !(is.na(type)) & type != str)
    }
    lc <- t(as.matrix(lc[-c(1:4)]))
    rownames(lc) <- separate(as.data.frame(rownames(lc)),
      col = "rownames(lc)",
      into = c("bodyid", NA)
    )$bodyid
    lc[is.na(lc)] <- 0
  } else {
    tmp <- data.frame(
      "output" = lc$output,
      "lc_syn" = rowSums(lc[-c(1:3)], na.rm = TRUE)
    )
    lc <- filter(lc, output %in% tmp[tmp$lc_syn > syn_count, "output"])
    if (rm_self) {
      lc <- filter(lc, type != str)
    }
    rownames(lc) <- lc$output
    lc <- t(as.matrix(lc[-c(1:3)]))
    rownames(lc) <- separate(as.data.frame(rownames(lc)),
      col = "rownames(lc)", into = c("bodyid", NA)
    )$bodyid
    lc[is.na(lc)] <- 0
  }

  # run PCA
  lc_pca <- prcomp(lc, scale = TRUE)
  p1 <- fviz_eig(lc_pca, main = str_c("Scree plot for ", str))

  # return(list("PCA" = lc_pca, "eigen_plot" = p1, "cluster_plot" = p2, 'pca_df' = lc))
  return(list("PCA" = lc_pca, "eigen_plot" = p1, "pca_df" = lc))
}

#### Cluster LCs by their components  #####
# takes PCA object and number of components and clusters, and performs k means
# clustering
# returns k means clustering and plot
lc_pca_clust <- function(lc_pca, num_comp, num_clust, str) {
  # Initialize the seed for the 'nstart' initial locations in kmeans
  set.seed(100000 * as.POSIXlt(Sys.time())$sec)

  # Cut the PCA to at most the desired number of PC
  comp <- lc_pca$x[, 1:min(num_comp, ncol(lc_pca$x))]

  # Evaluate the withinss for 2->5 clusters in kmeans
  wss <- (nrow(comp) - 1) * sum(apply(comp, 2, var)) # Withinss for k = 1
  for (i in 2:15) {
    wss[i] <- sum(kmeans(comp, i, nstart = 25, iter.max = 1000)$withinss)
  }

  #
  p1 <- ggplot(data.frame(1:15, wss), aes(x = X1.15, y = wss)) +
    geom_point() +
    geom_line() +
    xlab("Number of Clusters") +
    ylab("Within groups sum of squares") +
    ggtitle(lc)

  # Evaluate kmeans for the specified number of clusters
  k <- kmeans(comp, num_clust, nstart = 25, iter.max = 1000)
  p2 <- ggplot(
    as.data.frame(comp),
    aes(x = PC1, y = PC2, col = as.character(k$clust))
  ) +
    geom_point() +
    guides(col = guide_legend(title = "k means cluster")) +
    ggtitle(sprintf("Clusters of %s neurons", str))

  return(list("withinss" = p1, "kmeans" = k, "clusters" = p2))
}

#### Plot ROI outputs for a set of LC clusters #####
# determine if lc clusters have different patterns of output
# takes k the output of kmeans clustering on PCA of LC neurons
lc_clust_roi <- function(k, str) {
  clust <- names(sort(table(k$clust)))
  clust_ROIs <- data.frame()
  for (c in clust) {
    c1 <- names(k$cluster[k$cluster == c])
    c1_rois <- neuprint_get_roiInfo(c1)
    c1_rois <- data.frame(
      roiname = names(colSums(c1_rois[-1], na.rm = TRUE)),
      syn_count = colSums(c1_rois[-1], na.rm = TRUE),
      cluster = rep(sprintf("cluster_%s", c), ncol(c1_rois) - 1)
    )
    c1_rois <- separate(c1_rois,
      col = "roiname",
      into = c("roi", "prepost"),
      sep = "[.]",
      extra = "merge"
    )
    c1_rois <- filter(c1_rois, prepost == "pre" & syn_count > 10)
    # take average number of synapses per neuron
    c1_rois$syn_count <- c1_rois$syn_count / length(c1)
    clust_ROIs <- rbind(clust_ROIs, c1_rois)
  }
  p1 <- ggplot(data = clust_ROIs, aes(x = roi, y = syn_count, fill = cluster)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    scale_fill_brewer(palette = "Spectral") +
    labs(x = "ROI", y = "Average Number of Presynaptic Terminals Per neuron") +
    ggtitle(sprintf("Projection to ROIs of different %s clusters", str)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(list("ROI_df" = clust_ROIs, "ROI_plot" = p1))
}


###### K means Spatial Centers #####
# takes k means object of LC clusters and returns the spatial locations of
# presynaptic terminals of each cluster
# takes median of presynaptic terminals
lc_clust_syn <- function(k, str) {
  clust <- names(sort(table(k$clust)))
  clust_syn <- data.frame()
  for (c in clust) {
    c1 <- names(k$cluster[k$cluster == c])
    df <- neuprint_get_synapses(c1)
    syn <- filter(df, prepost == 0)
    tmp <- data.frame(
      "xyz" = c("x", "y", "z"),
      "syn_loc" = c(median(syn$x), median(syn$y), median(syn$z)),
      "cluster" = rep(sprintf("cluster_%s", c), 3)
    )
    clust_syn <- rbind(clust_syn, tmp)
  }

  p1 <- ggplot(data = clust_syn, aes(x = xyz, y = syn_loc, fill = cluster)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    scale_fill_brewer(palette = "Spectral") +
    labs(x = "Location", y = "Location of presynaptic terminal") +
    ggtitle(sprintf("Spatial location of pre-synaptic terminals for different %s clusters", str)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(list("syn_df" = clust_syn, "syn_plot" = p1))
}

#### Makes opponent connectivity gradient plot for LC ####
# x axis is LC ordered from fewest to highest synapses for 1 DN, overlaid with
# another DN
# y axis is number of synapses with DN
## LC_outputs is list of LC outputs from simple_connectivity
## LC is string of LC neuron type, DN_interest is list of strings of DN names
# makes either a double bar plot or a scatterplot
opponent_gradient_plot <- function(LC_outputs, DN, LC, DN_interest, plot_type = "double_bar") {
  plots <- list()
  lc <- LC_outputs[[LC]]
  pairs <- combn(DN_interest, 2)
  for (i in 1:ncol(pairs)) {
    dn1 <- pairs[1, i]
    dn2 <- pairs[2, i]
    dn_id1 <- filter(DN, DN.Type == dn1) %>% pull(bodyid)
    dn_id2 <- filter(DN, DN.Type == dn2) %>% pull(bodyid)

    # make data frame of LC4s and number synapses to DN
    syn_dist <- data.frame(
      as.numeric(lc[lc$output == dn_id1, -c(1:3)]),
      as.numeric(lc[lc$output == dn_id2, -c(1:3)]),
      colnames(lc)[-c(1:3)]
    ) %>%
      setNames(c("syn1", "syn2", "lc_bodyid"))
    syn_dist[is.na(syn_dist$syn1), "syn1"] <- 0
    syn_dist[is.na(syn_dist$syn2), "syn2"] <- 0
    # order by number synapses of first DN
    syn_dist <- syn_dist[order(syn_dist$syn1), ]
    syn_dist["ID"] <- 1:nrow(syn_dist)
    t <- cor.test(syn_dist$syn1, syn_dist$syn2)
    if (plot_type == "double_bar") {
      # barplot
      plots[[i]] <- ggplot(data = syn_dist, aes(width = 0.9)) +
        geom_bar(aes(x = ID, y = syn1, fill = "b"),
          stat = "identity",
          position = "identity", colour = "black", alpha = 0.5
        ) +
        geom_bar(aes(x = ID, y = syn2, fill = "r"),
          stat = "identity",
          position = "identity", colour = "black", alpha = 0.5
        ) +
        scale_fill_manual(
          name = "DN",
          values = c("b" = "blue", "r" = "red"),
          labels = c(dn1, dn2)
        ) +
        ggtitle(sprintf("%s connectivity to %s and %s", LC, dn1, dn2)) +
        labs(x = sprintf("LC ID (ordered by %s)", dn1), y = "Number Synapses ") +
        theme_classic() +
        annotate("text",
          x = max(syn_dist$syn1) - median(syn_dist$syn1) / 2, y = max(syn_dist$syn2),
          label = sprintf("\u03C1 = %.2f; p-value = %.1e", t$estimate, t$p.value)
        )
    }
    if (plot_type == "scatter") {
      plots[[i]] <- ggplot(data = syn_dist, aes(x = syn1, y = syn2)) +
        geom_point(size = 1) +
        geom_smooth(method = lm, se = FALSE, linetype = "dashed", col = "gray", size = 0.5) +
        ggtitle(sprintf("%s connectivity to %s and %s", LC, dn1, dn2)) +
        labs(
          x = sprintf("Number of synapses to %s", dn1),
          y = sprintf("Number of synapses to %s", dn2)
        ) +
        theme_classic() +
        annotate("text",
          x = max(syn_dist$syn1) - median(syn_dist$syn1) / 2, y = max(syn_dist$syn2),
          label = sprintf("\u03C1 = %.2f; p-value = %.1e", t$estimate, t$p.value)
        )
    }
  }
  return(plots)
}

#### Get correlation heatmap of outputs ####
# takes list of output matrices (LC cols, output rows), name of LC,
# cutoff for rho and number of clusters for dendrogram display
cor_map <- function(mat_list, lc_str, rho, c_group) {
  cor_mat <- cor(t(mat_list[[lc_str]][-c(1:3)]), method = "spearman")
  # outputs with >1 corr >|rho|
  cor_mat <- cor_mat[
    cor_mat[1:ncol(cor_mat)] > rho | cor_mat[1:ncol(cor_mat)] < -rho,
    cor_mat[1:ncol(cor_mat)] > rho | cor_mat[1:ncol(cor_mat)] < -rho
  ]

  dend1 <- as.dendrogram(hclust(dist(cor_mat)))
  dend1 <- color_branches(dend1, k = c_group, col = rainbow_hcl) # add color to the lines
  dend1 <- color_labels(dend1, k = c_group, col = rainbow_hcl) # add color to the labels

  # reorder the dendrogram, must incl. `agglo.FUN = mean`
  dend1 <- reorder(dend1, rowMeans(cor_mat, na.rm = T), agglo.FUN = mean)

  # get the color of the leaves (labels) for `heatmap.2`
  col_labels <- get_leaves_branches_col(dend1)
  col_labels <- col_labels[order(order.dendrogram(dend1))]
  h1 <- heatmap.2(cor_mat,
    scale = "none", col = bluered(100), Rowv = dend1,
    trace = "none", density.info = "none",
    breaks = seq(-1, 1, length.out = 101), labRow = "", labCol = "",
    RowSideColors = col_labels,
    Colv = "Rowv", dendrogram = "none", symm = TRUE, key = FALSE
  )
  return(h1)
}
