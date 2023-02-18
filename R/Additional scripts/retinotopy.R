# Created on Wed Aug  4 13:25:59 2021
# Name:    retinotopy.R
# Purpose: Evaluate if there is a relationship between the distances of the
#          end-points centroids in the lobula  and in the glomerulus
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - One 3D plot showing the used endpoints and the evaluated centroids both in
#   the glomerulus and the lobula
# - Two heatmaps showing the distances between neurons' centroids in the lobula
#   and in the glomerulus
# - Two plots exploring the potential correlation between the distances in the
#   lobula and in the glomerulus.
#
#
# Copyright (c) 2021 Andrea Vaccari
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#
# TODO:
# - Can be optimized by plotting all the endpoints and the synapses without
#   iterating, but the current version allows to select specific neurons to
#   visualize.
# - The computation of distances can be optimized.
# - project on the usual plane for the lobula and try those distances
# - try to do separately DV and AP using just the AP and DV coordinates
# - add kendall tau or spearman for the sorted distances
# - check if it is worth ordering the pairs in the glomerulus correlation plot
#   as those in the lobula
# - check what is going on with LC9 and LC10


# Import required libraries
library(tidyverse)
library(natverse)
library(zeallot)
library(pracma)
library(ks)
library(matrixStats)
library(egg)
library(corrplot)
library(aricode)



# Clean everything up ----
# Except the connection and skeleton files if they are already loaded.
items <- ls()
items <- items[items != "con"] # Comment this line to reload the connection file
items <- items[items != "nlist"] # Comment this line to reload the skeleton file
rm(list = items)

# Close any open plotting windows
while (dev.cur() > 1) {
  dev.off()
}
while (rgl.cur() > 0) {
  rgl.close()
}

# Load datasets
if (!exists("con")) {
  con <- readRDS("data/hemibrain_con0.rds")
}
if (!exists("nlist")) {
  nlist <- readRDS("data/nlist1.rds")
}

# Source local files
source("R/aux_functions.R")








###############################################################################
# Define items to analyze here
pre_type <- "LPLC4"

# Get the correct planes
lob_planes <- get.plane(pre_type, "lobula")
glo_planes <- get.plane(pre_type, "glomerulus")

# Plot window size
win_siz <- 1000

###############################################################################


# Identify body all posts and body IDs based on pre ----
# Extract all pre.type synapses excluding with pre
pre <- con %>%
  filter(pre.type == pre_type, post.type != pre_type)

# Find all the pre neurons
pre_IDs <- as.character(unique(pre$pre.bodyID))
pre_neu <- nlist[pre_IDs]
# pre_neu <- nlist[c("1907587934", "1838257401")]  # A-P
# pre_neu <- nlist[c("1838257401", "1281303666")]  # P-P
# Display all the neurons defined by pre
open3d()
par3d("windowRect" = c(100, 100, win_siz, win_siz))

# Add neuron skeletons (to visualize the somas)
plot3d(pre_neu, soma = TRUE, col = "light gray", add = TRUE)

# Show all endpoints
for (i in 1:length(pre_neu)) {
  neu <- pre_neu[[i]]
  ep <- neu$d[neu$EndPoints, ] %>%
    select(X, Y, Z)
  plot3d(ep, col = "black", size = 1, add = TRUE)
}

# and the lobula planes
planes3d(a = lob_planes[, 1], b = lob_planes[, 2], c = lob_planes[, 3], d = lob_planes[, 4], alpha = 0.2, add = TRUE)

# Select preserved end points for each neuron
for (i in 1:length(pre_neu)) {
  neu <- pre_neu[[i]]
  ep <- neu$d[neu$EndPoints, ]
  pts <- ep
  for (p in 1:nrow(lob_planes)) {
    pts <- pts %>%
      mutate(LO = lob_planes[p, 5] * (lob_planes[p, 1] * X + lob_planes[p, 2] * Y + lob_planes[p, 3] * Z + lob_planes[p, 4])) %>%
      filter(LO < 0) %>%
      select(X, Y, Z)
  }
  if (i == 1) {
    end_pts <- pts
  } else {
    end_pts <- rbind(end_pts, pts)
  }
}

# Show the preserved points
plot3d(end_pts, col = "green", size = 3, add = TRUE)

# Evaluates the lobula centroids for each body ID ----
# Calculate centroid of lobula end points
for (i in 1:length(pre_neu)) {
  neu <- pre_neu[[i]]
  pts <- neu$d[neu$EndPoints, ]
  for (p in 1:nrow(lob_planes)) {
    pts <- pts %>%
      mutate(LO = lob_planes[p, 5] * (lob_planes[p, 1] * X + lob_planes[p, 2] * Y + lob_planes[p, 3] * Z + lob_planes[p, 4])) %>%
      filter(LO < 0) %>%
      select(X, Y, Z)
  }
  pts <- pts %>%
    colMeans()
  if (i == 1) {
    lob_ctrs <- pts
  } else {
    lob_ctrs <- rbind(lob_ctrs, pts)
  }
}

# SHow in red on the 3d plot
plot3d(lob_ctrs, col = "red", size = 5, add = TRUE)

# Create a dataframe of lobula centroids
lob_ctrs_df <- data.frame(lob_ctrs) %>%
  rename(
    lob.X = X,
    lob.Y = Y,
    lob.Z = Z
  )
lob_ctrs_df$bodyid <- pre_IDs

# Add the synapses
for (n in pre_neu) {
  pre_coors <- pre %>%
    filter(pre.bodyID == n$bodyid) %>%
    pull(pre.coors) %>%
    xyzmatrix()
  plot3d(pre_coors, size = 1, add = TRUE)
}

# and the glomerulus planes
planes3d(a = glo_planes[, 1], b = glo_planes[, 2], c = glo_planes[, 3], d = glo_planes[, 4], alpha = 0.2, add = TRUE)

# Filter for glomerulus synapses using all the planes
for (n in pre_neu) {
  pre.glo <- pre %>%
    filter(pre.bodyID == n$bodyid)
  for (i in 1:nrow(glo_planes)) {
    pre.glo <- pre.glo %>%
      filter(pre.glo %>%
        pull(pre.coors) %>%
        xyzmatrix() %>%
        as_tibble() %>%
        mutate(GLO = glo_planes[i, 5] * (glo_planes[i, 1] * X + glo_planes[i, 2] * Y + glo_planes[i, 3] * Z + glo_planes[i, 4]) > 0) %>%
        select(GLO))
  }
  # Show filtered synapses
  pre_kept <- pre.glo %>%
    pull(pre.coors) %>%
    xyzmatrix()
  plot3d(pre_kept, col = "green", size = 3, add = TRUE)
}

# Evaluates the glomerulus centroids for each body ID ----
# Calculate centroid of glomerulus synapses
for (i in 1:length(pre_neu)) {
  neu <- pre_neu[[i]]
  syn <- pre %>% filter(pre.bodyID == neu$bodyid)
  for (p in 1:nrow(glo_planes)) {
    syn <- syn %>%
      filter(syn %>%
        pull(pre.coors) %>%
        xyzmatrix() %>%
        as_tibble() %>%
        mutate(GLO = glo_planes[p, 5] * (glo_planes[p, 1] * X + glo_planes[p, 2] * Y + glo_planes[p, 3] * Z + glo_planes[p, 4]) > 0) %>%
        select(GLO))
  }
  pts <- syn %>%
    pull(pre.coors) %>%
    xyzmatrix() %>%
    colMeans()
  if (i == 1) {
    glo_ctrs <- pts
  } else {
    glo_ctrs <- rbind(glo_ctrs, pts)
  }
}

# SHow in red on the 3d plot
plot3d(glo_ctrs, col = "red", size = 5, add = TRUE)

# Create a dataframe of glomerulus centroids
glo_ctrs_df <- data.frame(glo_ctrs) %>%
  rename(
    glo.X = X,
    glo.Y = Y,
    glo.Z = Z
  )
glo_ctrs_df$bodyid <- pre_IDs

# Merge the dataframes and drop NAs
ctrs_df <- merge(lob_ctrs_df, glo_ctrs_df, by = "bodyid") %>%
  na.omit()

# Create a dataframe with one entry for all possible combinations of bodyIDs
# and corresponding distances in glomerulus and lobula

# Evaluate all possible combinations of body IDs
com_full <- combn(ctrs_df$bodyid, 2)

# Create the dataframe
dist_df <- data.frame(t(com_full)) %>%
  rename(
    Id1 = X1,
    Id2 = X2
  )

# Add empty columns to the dataframe
dist_df$lob.dist <- 0
dist_df$glo.dist <- 0

# Create a couple of array to hold the distances saparately for glomerulus and
# lobula
glo_dist_mtrx <- matrix(0, nrow = nrow(ctrs_df), ncol = nrow(ctrs_df))
rownames(glo_dist_mtrx) <- ctrs_df$bodyid
colnames(glo_dist_mtrx) <- ctrs_df$bodyid
lob_dist_mtrx <- matrix(0, nrow = nrow(ctrs_df), ncol = nrow(ctrs_df))
rownames(lob_dist_mtrx) <- ctrs_df$bodyid
colnames(lob_dist_mtrx) <- ctrs_df$bodyid


# Calculate distances in lobula and glomerulus for each pair and add to the
# dataframe
for (r in 1:nrow(dist_df)) {
  id1 <- dist_df[r, 1]
  id2 <- dist_df[r, 2]
  ctrs1 <- ctrs_df[ctrs_df$bodyid == id1, ]
  ctrs2 <- ctrs_df[ctrs_df$bodyid == id2, ]
  lob1 <- ctrs1[2:4]
  lob2 <- ctrs2[2:4]
  glo1 <- ctrs1[5:7]
  glo2 <- ctrs2[5:7]
  lob.diff <- lob1 - lob2
  lob.dist <- sqrt(sum(lob.diff * lob.diff))
  glo.diff <- glo1 - glo2
  glo.dist <- sqrt(sum(glo.diff * glo.diff))
  dist_df[r, ]$lob.dist <- lob.dist
  dist_df[r, ]$glo.dist <- glo.dist
  lob_dist_mtrx[id1, id2] <- lob.dist
  lob_dist_mtrx[id2, id1] <- lob.dist
  glo_dist_mtrx[id1, id2] <- glo.dist
  glo_dist_mtrx[id2, id1] <- glo.dist
}

# Normalize to 1 for easier comparison
lob_dist_mtrx_n <- lob_dist_mtrx / max(lob_dist_mtrx)
glo_dist_mtrx_n <- glo_dist_mtrx / max(glo_dist_mtrx)
dist_df$lob.dist.n <- dist_df$lob.dist / max(dist_df$lob.dist)
dist_df$glo.dist.n <- dist_df$glo.dist / max(dist_df$glo.dist)

# Evaluate the correlation between the two
corr <- cor(c(dist_df$lob.dist), c(dist_df$glo.dist), method = "pearson", use = "complete.obs")

# Evaluate the normalized mutual information
nmi <- NMI(round(dist_df$lob.dist), round(dist_df$glo.dist))

# Evaluate the Spearman's coefficient
spear <- cor(c(dist_df$lob.dist), c(dist_df$glo.dist), method = "spearman", use = "complete.obs")

# Evaluate the Kendall's coefficient
tau <- cor(c(dist_df$lob.dist), c(dist_df$glo.dist), method = "kendall", use = "complete.obs")

# Print coefficients
cat(paste(pre_type, ":\n"))
cat(" Correlation coefficient (Pearson):", corr, "\n")
cat(" Normalized Mutual Information:", nmi, "\n")
cat(" Correlation coefficient (Spearman):", spear, "\n")
cat(" Correlation coefficient (Kendall):", tau, "\n")

# Plot the correlation matrices
corrplot(lob_dist_mtrx_n,
  is.corr = FALSE, # It is not a correlation matrix
  method = "color", # Color the background
  # type='upper',  # Only upper diagonal
  order = "hclust", # Hierarchical
  hclust.method = "median",
  # addCoef.col='black',  # Add values in black
  # number.cex=0.6,  # values size
  # diag=FALSE,  # No diagonal
  tl.col = "black"
) # Color of the labels in black

corrplot(glo_dist_mtrx_n,
  is.corr = FALSE, # It is not a correlation matrix
  method = "color", # Color the background
  # type='upper',  # Only upper diagonal
  order = "hclust", # Hierarchical
  hclust.method = "median",
  # addCoef.col='black',  # Add values in black
  # number.cex=0.6,  # values size
  # diag=FALSE,  # No diagonal
  tl.col = "black"
) # Color of the labels in black


# Plot the correlation
ggplot() +
  ylab("Distance between centroids in the lobula, mcm") +
  xlab("Distance between centroids in the glomerulus, mcm") +
  theme_classic() +
  geom_point(
    data = dist_df,
    aes(x = 0.008 * glo.dist, y = 0.008 * lob.dist),
    size = 1,
    col = "steelblue"
  ) +
  geom_smooth(
    data = dist_df,
    aes(x = 0.008 * glo.dist, y = 0.008 * lob.dist),
    method = "lm",
    formula = "y ~ x",
    col = "red",
    se = TRUE
  ) +
  annotate("text",
    col = "red", size = 8,
    label = toString(paste("r =", round(corr, 2))),
    x = 0.8 * 0.008 * max(dist_df$glo.dist),
    y = 0.8 * 0.008 * max(dist_df$lob.dist),
    size = 5
  ) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
    axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
    plot.title = element_text(face = "bold", color = "blue", size = 15),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  )

# Using the normalized data
ggplot() +
  ylab("Distance between centroids in the lobula (normalized to max)") +
  xlab("Distance between centroids in the glomerulus (normalized to max)") +
  theme_classic() +
  geom_point(
    data = dist_df,
    aes(x = glo.dist.n, y = lob.dist.n),
    size = 1,
    col = "steelblue"
  ) +
  geom_smooth(
    data = dist_df,
    aes(x = glo.dist.n, y = lob.dist.n),
    method = "lm",
    formula = "y ~ x",
    col = "red",
    se = TRUE
  ) +
  annotate("text",
    col = "red", size = 8,
    label = toString(paste("r =", round(corr, 2))),
    x = 0.8,
    y = 0.8,
    size = 5
  ) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
    axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
    plot.title = element_text(face = "bold", color = "blue", size = 15),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  )
