# Created on Sun Jun 27 21:50:13 2021
# Name:    labula.R
# Purpose: Evaluates distribution of connections of the two post in the lobula
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - A 3d plot showing all the selected endpoints considered in the lobula, the
#   selection plane, the centroid of each pre neuron calculated using only
#   the selected end points, the projection of the end points and the center of
#   mass on the selection plane and the axes used to build the 2d image
# - A 2d plot showing the projection of the centroids for the given VPN
# - A 2d plot comparing the projections of the centroids for the two
#   posts within the projection of the lobula on the selection plane. The plot
#   also shows the weighted median (and weighted box plots) for the centers of
#   mass (weighted by the number of synapses with pre)
# - A 2d plot showing the distance between the weighted medians (two versions).
# - A 2d plot showing the gradient of synapses sorted using the position of the
#   centroids with respect to either the anterior-posterior or the
#   dorsal-ventral as specified by the user
#
#
# Copyright (c) 2022 Andrea Vaccari
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
# - duplicate the script and allow users to specify three groups of posts that
#   should be shown in different colors on the loblula. 3 for LC4 and 5 for
#   LPLC2. (maybe it could be done in such a way that users can select how many
#   different clusters they need)


# Import required libraries
library(tidyverse)
library(natverse)
library(cetcolor)
library(zeallot)
library(pracma)
library(ks)
library(matrixStats)
library(egg)
library(ggpmisc)


# Clean everything up ----
# Except the connection and skeleton files if they are already loaded.
items <- ls()
items <- items[items != "con"] # Comment to reload the connection file
items <- items[items != "nlist"] # Comment to reload the skeleton file
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
pre_type <- "LC24"
post_type1 <- "mALC2" # Red
post_type2 <- "CL246" # Blue

# Grab the correct plane
# NOTE: the FIRST plane specified is also use as projection plane
planes <- get.plane(pre_type, "lobula")

# Using perceptually uniform color map
# More options here:
# https://cran.r-project.org/web/packages/cetcolor/vignettes/cet_color_schemes.html
c_map_single <- "r2" # Map to use for individual post plot
c_size_single <- 100 # Number of colors in the map
c_map_both <- "d9" # Map to use for both posts together
c_size_both <- 100 # Number of colors in the map

# Projection line a and b
# The a and b of the projection line to use to generate the contrcution figure.
# The coefficients for a specific VPN can be found running the
# lobula_ranking.R script.
# If all zeros, the figure will not be generated.
proj_v <- c(983.6229896, -0.5838685) # Top25 for LC4
# proj_v <- c(-476.2505,  0.2322217)  # Top25 for LPLC2
# proj_v <- c(0, 0)

# Artificial line offset to apply to the projection line for visualization purposes
line_offset <- 2000 # For LC4
# line_offset <- 0  # For LPLC2

# When generating the gradients, sort the data by dorsal-ventral. If FALSE, the
# data will be sorted by anterior-posterior
sort_by_dv <- FALSE

# Plot window size
win_siz <- 1000
###############################################################################







# Normalize the vector normal to the plane
w_norm <- planes[1, 1:3]
w0 <- planes[1, 4]

# Pre ----
# Extract all pre.type synapses
pre <- con %>% filter(pre.type == pre_type)

# Post ----
# Extract bodyIDs of pre synapses with post.type1
post1 <- pre %>%
  filter(post.type == post_type1) %>%
  group_by(pre.bodyID) %>%
  dplyr::count()

# Extract bodyIDs of pre synapses with post.type2
post2 <- pre %>%
  filter(post.type == post_type2) %>%
  group_by(pre.bodyID) %>%
  dplyr::count()

# Merge the two posts with the counts
post <- merge(post1,
  post2,
  by = "pre.bodyID",
  all = TRUE,
  suffix = c(".post1", ".post2")
)

# Clear potential NaN
post[is.na(post)] <- 0

# Get the neuron list from the posts bodyIDs
n_list <- nlist[as.character(post$pre.bodyID)]

# Drop potential NAs
n_list <- n_list[!is.na(names(n_list))]
dropped <- setdiff(as.character(post$pre.bodyID), names(n_list))
if (length(dropped) != 0) {
  post <- subset(post, post$pre.bodyID != dropped)
}

# 3D analysis ----
nopen3d()
par3d("windowRect" = c(100, 100, win_siz, win_siz))

# Add neuron skeletons (to visualize the somas)
plot3d(n_list, soma = TRUE, col = "light gray", add = TRUE)

# Show all endpoints
for (i in 1:length(n_list)) {
  neu <- n_list[[i]]
  ep <- neu$d[neu$EndPoints, ] %>%
    select(X, Y, Z)
  plot3d(ep, col = "black", size = 1, add = TRUE)
}


# and the planes
planes3d(
  a = planes[1, 1],
  b = planes[1, 2],
  c = planes[1, 3],
  d = planes[1, 4],
  col = "blue",
  alpha = 0.2,
  add = TRUE
)
planes3d(
  a = planes[-1, 1],
  b = planes[-1, 2],
  c = planes[-1, 3],
  d = planes[-1, 4],
  alpha = 0.2,
  add = TRUE
)

# Select preserved end points for each neuron
for (i in 1:length(n_list)) {
  neu <- n_list[[i]]
  ep <- neu$d[neu$EndPoints, ]
  pts <- ep
  for (p in 1:nrow(planes)) {
    pts <- pts %>%
      mutate(LO = planes[p, 5] * (planes[p, 1] * X + planes[p, 2] * Y + planes[p, 3] * Z + planes[p, 4])) %>%
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








# Project the preserved end points on the plane in yellow.
end_pts.proj <- plane.proj(end_pts, w_norm, -w0)
plot3d(end_pts.proj, col = "yellow", size = 1, add = TRUE)

# Find the centroid of the end points projections and plot in white
# This is the origin of the coordinate system
center <- colMeans(end_pts.proj)
plot3d(center[1], center[2], center[3], col = "white", size = 10, add = TRUE)

# Calculate centroid of lobula end points and plot in red
for (i in 1:length(n_list)) {
  neu <- n_list[[i]]
  pts <- neu$d[neu$EndPoints, ]
  for (p in 1:nrow(planes)) {
    pts <- pts %>%
      mutate(LO = planes[p, 5] * (planes[p, 1] * X + planes[p, 2] * Y + planes[p, 3] * Z + planes[p, 4])) %>%
      filter(LO < 0) %>%
      select(X, Y, Z)
  }
  pts <- pts %>%
    colMeans()
  if (i == 1) {
    ctrs <- pts
  } else {
    ctrs <- rbind(ctrs, pts)
  }
}

# Set the row names to be the body id
rownames(ctrs) <- post$pre.bodyID

# Clear rows with NaN
ctrs <- ctrs[!is.na(ctrs[, 1]), ]

plot3d(ctrs, col = "red", size = 5, add = TRUE)

# Project the centroid of the end points on the surface in blue
ctrs.proj <- plane.proj(ctrs, w_norm, -w0)
plot3d(ctrs.proj, col = "blue", size = 5, add = TRUE)

# 2D analysis ----
# Define coordinate system in the plane
# Use the projection of the z-axis on the plane to determine one of the new 2D
# axes
# If more than one plane are defined to select the lobula region, the first one
# is the plane used for the projection.
pointY <- c(plane.proj(t(center + c(0, 0, 2000)), w_norm, -w0))
pointY <- pointY - center
unitY <- -pointY / sqrt(sum(pointY * pointY))
arrow3d(center, center + 2000 * unitY, thikness = 1, col = "white")

pointX <- cross(w_norm, pointY)
unitX <- -pointX / sqrt(sum(pointX * pointX))
arrow3d(center, center + 2000 * unitX, thikness = 1, col = "white")

# Project data on 2D plane ----
# All considered endpoints
Xep <- as.matrix(sweep(end_pts.proj, 2, center)) %*% unitX
Yep <- as.matrix(sweep(end_pts.proj, 2, center)) %*% unitY
end_pts.plane <- data.frame(cbind(Xep, Yep))

# Centroids
Xcm <- as.matrix(sweep(ctrs.proj, 2, center)) %*% unitX
Ycm <- as.matrix(sweep(ctrs.proj, 2, center)) %*% unitY
ctrs.plane <- data.frame(cbind(Xcm, Ycm))

# Join the counts of synapses with the centroid dataset
ctrs.plane["pre.bodyID"] <- row.names(ctrs.plane)
ctrs.plane <- merge(ctrs.plane, post, by = "pre.bodyID", all = TRUE) %>% na.omit()

# 2d maps of synapses-weighted centroid projections in the lobula ----
# Create the color maps
col_single <- cet_pal(c_size_single, name = c_map_single)

# Evaluate parameters for the box plots
dist_box <- 0.1 * abs(max(min(end_pts.plane$X1), min(end_pts.plane$X2)))
width_box <- 0.75 * dist_box
X1_rng <- abs(max(end_pts.plane$X1) - min(end_pts.plane$X1))
X2_rng <- abs(max(end_pts.plane$X2) - min(end_pts.plane$X2))
seg_len <- 0.05 * min(X1_rng, X2_rng)

# Calculate weighted medians
p1_x1_wm <- weightedMedian(ctrs.plane$X1, w = ctrs.plane$n.post1, interpolate = FALSE)
p1_x2_wm <- weightedMedian(ctrs.plane$X2, w = ctrs.plane$n.post1, interpolate = FALSE)
p2_x1_wm <- weightedMedian(ctrs.plane$X1, w = ctrs.plane$n.post2, interpolate = FALSE)
p2_x2_wm <- weightedMedian(ctrs.plane$X2, w = ctrs.plane$n.post2, interpolate = FALSE)

# Distance between weighted medians
dist_wm <- sqrt((p1_x1_wm - p2_x1_wm)^2 + (p1_x2_wm - p2_x2_wm)^2)

# Calculate the convex hull
# For consistency use all the neuron in the pre so that it doesn't change based
# on the chosen posts, but just on the VPN.
hull <- pre %>%
  filter(post.type == pre_type) %>%
  group_by(pre.bodyID) %>%
  pull(pre.bodyID) %>%
  unique()

# Get the neuron list from the posts bodyIDs
n_list_hull <- nlist[as.character(hull)]

# Drop potential NAs
n_list_hull <- n_list_hull[!is.na(names(n_list_hull))]
# dropped <- setdiff(as.character(post$pre.bodyID), names(n_list_hull))
# post <- subset(post, post$pre.bodyID != dropped)


# Select preserved end points for each neuron
for (i in 1:length(n_list_hull)) {
  neu <- n_list_hull[[i]]
  ep <- neu$d[neu$EndPoints, ]
  pts <- ep
  for (p in 1:nrow(planes)) {
    pts <- pts %>%
      mutate(LO = planes[p, 5] * (planes[p, 1] * X + planes[p, 2] * Y + planes[p, 3] * Z + planes[p, 4])) %>%
      filter(LO < 0) %>%
      select(X, Y, Z)
  }
  if (i == 1) {
    end_pts_hull <- pts
  } else {
    end_pts_hull <- rbind(end_pts_hull, pts)
  }
}

# Project the hull end points on the 3D plane.
end_pts_hull.proj <- plane.proj(end_pts_hull, w_norm, -w0)


# Project the hull endpoints on the 2D plane
Xep_hull <- as.matrix(sweep(end_pts_hull.proj, 2, center)) %*% unitX
Yep_hull <- as.matrix(sweep(end_pts_hull.proj, 2, center)) %*% unitY
end_pts_hull.plane <- data.frame(cbind(Xep_hull, Yep_hull))

# Evaluate convex hull
lo <- chull(end_pts_hull.plane)

# Sort the data to generate the  gradient plots
ctrs.plane.sorted <- ctrs.plane %>%
  arrange(X1)
if (sort_by_dv == TRUE) {
  ctrs.plane.sorted <- ctrs.plane %>%
    arrange(X2)
}
# Add index as column
ctrs.plane.sorted$idx <- as.numeric(rownames(ctrs.plane.sorted))

# Evaluate the limits for the plots
limits1 <- c(0, max(ctrs.plane$n.post1))
breaks1 <- round(c(0, .25, .5, .75, 1) * max(ctrs.plane$n.post1))
limits2 <- c(0, max(ctrs.plane$n.post2))
breaks2 <- round(c(0, .25, .5, .75, 1) * max(ctrs.plane$n.post2))

# Plot projection of centroids
ggplot() +
  coord_fixed() +
  xlab("A-P axis, um") +
  ylab("D-V axis, um") +
  geom_polygon(
    data = end_pts_hull.plane[lo, ],
    aes(x = 0.008 * X1, y = 0.008 * X2),
    alpha = 0.3
  ) +
  geom_point(
    data = ctrs.plane,
    aes(x = 0.008 * X1, y = 0.008 * X2)
  ) +
  scale_color_gradientn(
    name = "synapses\ncount",
    colours = col_single,
    breaks = breaks1,
    limits = limits1
  ) +
  scale_size_continuous(
    name = "synapses\ncount",
    limits = limits1,
    breaks = breaks1,
    range = c(0, 7)
  ) +
  guides(colour = guide_legend(), size = guide_legend()) +
  ggtitle(paste(pre_type, "Dendritic map")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
    axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
    plot.title = element_text(face = "bold", color = "black", size = 15),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  )


# Generate the plot for post1
proj1 <- ggplot() +
  coord_fixed() +
  xlab("A-P axis, um") +
  ylab("D-V axis, um") +
  geom_polygon(
    data = end_pts_hull.plane[lo, ],
    aes(x = 0.008 * X1, y = 0.008 * X2),
    alpha = 0.3
  ) +
  # geom_point(data=end_pts.plane,
  #            aes(x=0.008 * X1, y=0.008 * X2),
  #            size=0.1,
  #            alpha=0.2) +
  geom_point(
    data = ctrs.plane,
    aes(x = 0.008 * X1, y = 0.008 * X2, colour = n.post1, size = n.post1)
  ) +
  scale_color_gradientn(
    name = "synapses\ncount",
    colours = col_single,
    breaks = breaks1,
    limits = limits1
  ) +
  scale_size_continuous(
    name = "synapses\ncount",
    limits = limits1,
    breaks = breaks1,
    range = c(0, 7)
  ) +
  guides(colour = guide_legend(), size = guide_legend()) +
  # geom_boxplot(data=ctrs.plane,
  #              aes(x=0.008 * X1,
  #                  y=0.008 * (min(end_pts.plane$X2) - dist_box),
  #                  weight=n.post1),
  #              width=0.008 * width_box,
  #              notch=TRUE) +
  # geom_boxplot(data=ctrs.plane,
  #              aes(y=0.008 * X2,
  #                  x=0.008 * (min(end_pts.plane$X1) - dist_box),
  #                  weight=n.post1),
  #              width=0.008 * width_box,
  #              notch=TRUE) +
  # geom_segment(aes(
  #   x = 0.008 * (p1_x1_wm - seg_len),
  #   xend = 0.008 * (p1_x1_wm + seg_len),
  #   y = 0.008 * p1_x2_wm,
  #   yend = 0.008 * p1_x2_wm
  # ), col = "black", size = 2) +
  # geom_segment(aes(
  #   y = 0.008 * (p1_x2_wm - seg_len),
  #   yend = 0.008 * (p1_x2_wm + seg_len),
  #   x = 0.008 * p1_x1_wm,
  #   xend = 0.008 * p1_x1_wm
  # ), col = "black", size = 2) +
  ggtitle(paste(pre_type, ">", post_type1, "(lobula projection)")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
    axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
    plot.title = element_text(face = "bold", color = "black", size = 15),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  )

# Generate the plot for post2
proj2 <- ggplot() +
  coord_fixed() +
  xlab("A-P axis, um") +
  ylab("D-V axis, um") +
  geom_polygon(
    data = end_pts_hull.plane[lo, ],
    aes(x = 0.008 * X1, y = 0.008 * X2),
    alpha = 0.3
  ) +
  # geom_point(data=end_pts.plane,
  #            aes(x=0.008 * X1, y=0.008 * X2),
  #            size=0.1,
  #            alpha=0.2) +
  geom_point(
    data = ctrs.plane,
    aes(x = 0.008 * X1, y = 0.008 * X2, colour = n.post2, size = n.post2)
  ) +
  scale_color_gradientn(
    name = "synapses\ncount",
    colours = col_single,
    breaks = breaks2,
    limits = limits2
  ) +
  scale_size_continuous(
    name = "synapses\ncount",
    limits = limits2,
    breaks = breaks2,
    range = c(0, 7)
  ) +
  guides(colour = guide_legend(), size = guide_legend()) +
  # geom_boxplot(data=ctrs.plane,
  #              aes(x=0.008 * X1,
  #                  y=0.008 * (min(end_pts.plane$X2) - dist_box),
  #                  weight=n.post2),
  #              width=0.008 * width_box,
  #              notch=TRUE) +
  # geom_boxplot(data=ctrs.plane,
  #              aes(y=0.008 * X2,
  #                  x=0.008 * (min(end_pts.plane$X1) - dist_box),
  #                  weight=n.post2),
  #              width=0.008 * width_box,
  #              notch=TRUE) +
  # geom_segment(aes(
  #   x = 0.008 * (p2_x1_wm - seg_len),
  #   xend = 0.008 * (p2_x1_wm + seg_len),
  #   y = 0.008 * p2_x2_wm,
  #   yend = 0.008 * p2_x2_wm
  # ), col = "black", size = 2) +
  # geom_segment(aes(
  #   y = 0.008 * (p2_x2_wm - seg_len),
  #   yend = 0.008 * (p2_x2_wm + seg_len),
  #   x = 0.008 * p2_x1_wm,
  #   xend = 0.008 * p2_x1_wm
  # ), col = "black", size = 2) +
  ggtitle(paste(pre_type, ">", post_type2, "(lobula projection)")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
    axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
    plot.title = element_text(face = "bold", color = "black", size = 15),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  )

# Combine and show the plots
grid.arrange(proj1, proj2, ncol = 2)

# Generate a single plot that shows the location of the two centroids
ggplot() +
  coord_fixed() +
  xlab("A-P axis, um") +
  ylab("D-V axis, um") +
  geom_polygon(
    data = end_pts_hull.plane[lo, ],
    aes(x = 0.008 * X1, y = 0.008 * X2),
    alpha = 0.3
  ) +
  # geom_point(data=end_pts.plane,
  #            aes(x=0.008 * X1, y=0.008 * X2),
  #            size=0.5,
  #            alpha=0.1) +
  # geom_point(data=ctrs.plane,
  #            aes(x=0.008 * X1, y=0.008 * X2, size=n.post1),
  #            colour='red',
  #            shape=1,
  #            stroke=1.5) +
  # geom_point(data=ctrs.plane,
  #            aes(x=0.008 * X1, y=0.008 * X2, size=n.post2),
  #            colour='blue',
  #            shape=1,
  #            stroke=1.5) +
  geom_boxplot(
    data = ctrs.plane,
    aes(
      x = 0.008 * X1,
      y = 0.008 * (min(end_pts.plane$X2) - dist_box),
      weight = n.post1
    ),
    width = 0.008 * width_box,
    notch = TRUE,
    col = "red", size = 1
  ) +
  geom_boxplot(
    data = ctrs.plane,
    aes(
      y = 0.008 * X2,
      x = 0.008 * (min(end_pts.plane$X1) - dist_box),
      weight = n.post1
    ),
    width = 0.008 * width_box,
    notch = TRUE,
    col = "red", size = 1
  ) +
  geom_segment(aes(
    x = 0.008 * (p1_x1_wm - seg_len),
    xend = 0.008 * (p1_x1_wm + seg_len),
    y = 0.008 * p1_x2_wm,
    yend = 0.008 * p1_x2_wm
  ),
  col = "red", size = 3
  ) +
  geom_segment(aes(
    y = 0.008 * (p1_x2_wm - seg_len),
    yend = 0.008 * (p1_x2_wm + seg_len),
    x = 0.008 * p1_x1_wm,
    xend = 0.008 * p1_x1_wm
  ),
  col = "red", size = 3
  ) +
  geom_boxplot(
    data = ctrs.plane,
    aes(
      x = 0.008 * X1,
      y = 0.008 * (min(end_pts.plane$X2) - 2 * dist_box),
      weight = n.post2
    ),
    width = 0.008 * width_box,
    notch = TRUE,
    col = "blue", size = 1
  ) +
  geom_boxplot(
    data = ctrs.plane,
    aes(
      y = 0.008 * X2,
      x = 0.008 * (min(end_pts.plane$X1) - 2 * dist_box),
      weight = n.post2
    ),
    width = 0.008 * width_box,
    notch = TRUE,
    col = "blue", size = 1
  ) +
  geom_segment(aes(
    x = 0.008 * (p2_x1_wm - seg_len),
    xend = 0.008 * (p2_x1_wm + seg_len),
    y = 0.008 * p2_x2_wm,
    yend = 0.008 * p2_x2_wm
  ),
  col = "blue", size = 3
  ) +
  geom_segment(aes(
    y = 0.008 * (p2_x2_wm - seg_len),
    yend = 0.008 * (p2_x2_wm + seg_len),
    x = 0.008 * p2_x1_wm,
    xend = 0.008 * p2_x1_wm
  ),
  col = "blue", size = 3
  ) +
  geom_segment(aes(
    x = 0.008 * p1_x1_wm,
    xend = 0.008 * p2_x1_wm,
    y = 0.008 * p1_x2_wm,
    yend = 0.008 * p2_x2_wm
  ),
  col = "yellow", size = 1
  ) +
  annotate("text",
    label = toString(paste(round(0.008 * dist_wm, 2), "um")),
    x = 0.008 * max(p1_x1_wm, p2_x1_wm),
    y = 0.008 * max(p1_x2_wm, p2_x2_wm),
    size = 5, col = "black"
  ) +
  ggtitle(paste("Weighted median distance", post_type1, "(red) <=>", post_type2, "(blue)")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
    axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
    plot.title = element_text(face = "bold", color = "black", size = 12),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  )


# Generate the gradient plot based on the position of the centroids along the axis specified
# by the user

# Max number of synapses
max_syn <- max(ctrs.plane.sorted$n.post1, ctrs.plane.sorted$n.post2)

if (sort_by_dv == FALSE) {
  grad1 <- ggplot(
    data = ctrs.plane.sorted,
    aes(x = 0.008 * X1, y = n.post1, colour = n.post1)
  ) +
    coord_fixed(ratio = (0.008 * diff(range(ctrs.plane.sorted$X1)) / max_syn)) +
    xlab("A-P axis, um") +
    ylab("Number of synapses") +
    ylim(0, max_syn) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    scale_color_gradientn(colours = cet_pal(nrow(ctrs.plane.sorted))) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, color = "black") +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(..rr.label..)),
      parse = TRUE
    ) +
    ggtitle(paste(pre_type, ">", post_type1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
      axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
      plot.title = element_text(face = "bold", color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15)
    )

  grad2 <- ggplot(
    data = ctrs.plane.sorted,
    aes(x = 0.008 * X1, y = n.post2, colour = n.post2)
  ) +
    coord_fixed(ratio = (0.008 * diff(range(ctrs.plane.sorted$X1)) / max_syn)) +
    xlab("A-P axis, um") +
    ylab("Number of synapses") +
    ylim(0, max_syn) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    scale_color_gradientn(colours = cet_pal(nrow(ctrs.plane.sorted))) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, color = "black") +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(..rr.label..)),
      parse = TRUE
    ) +
    ggtitle(paste(pre_type, ">", post_type2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
      axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
      plot.title = element_text(face = "bold", color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15)
    )

  # Calculate R^2 for the plots
  rsq1 <- summary(lm(n.post1 ~ X1, data = ctrs.plane.sorted))$r.squared
  rsq2 <- summary(lm(n.post2 ~ X1, data = ctrs.plane.sorted))$r.squared
  cat("R-squared for", post_type1, ":", rsq1)
  cat("R-squared for", post_type2, ":", rsq2)
} else {
  grad1 <- ggplot(
    data = ctrs.plane.sorted,
    aes(x = 0.008 * X2, y = n.post1, colour = n.post1)
  ) +
    coord_fixed(ratio = (0.008 * diff(range(ctrs.plane.sorted$X2)) / max_syn)) +
    xlab("D-V axis, um") +
    ylab("Number of synapses") +
    ylim(0, max_syn) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    scale_color_gradientn(colours = cet_pal(nrow(ctrs.plane.sorted))) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, color = "black") +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(..rr.label..)),
      parse = TRUE
    ) +
    ggtitle(paste(pre_type, ">", post_type1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
      axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
      plot.title = element_text(face = "bold", color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15)
    )

  grad2 <- ggplot(
    data = ctrs.plane.sorted,
    aes(x = 0.008 * X2, y = n.post2, colour = n.post2)
  ) +
    coord_fixed(ratio = (0.008 * diff(range(ctrs.plane.sorted$X2)) / max_syn)) +
    xlab("D-V axis, um") +
    ylab("Number of synapses") +
    ylim(0, max_syn) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    scale_color_gradientn(colours = cet_pal(nrow(ctrs.plane.sorted))) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, color = "black") +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(..rr.label..)),
      parse = TRUE
    ) +
    ggtitle(paste(pre_type, ">", post_type2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
      axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
      plot.title = element_text(face = "bold", color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15)
    )

  # Calculate R^2 for the plots
  rsq1 <- summary(lm(n.post1 ~ X2, data = ctrs.plane.sorted))$r.squared
  rsq2 <- summary(lm(n.post2 ~ X2, data = ctrs.plane.sorted))$r.squared
  cat("R-squared for", post_type1, ":", rsq1)
  cat("R-squared for", post_type2, ":", rsq2)
}
# Combine and show the plots
grid.arrange(grad1, grad2, ncol = 2)





# If a projection line is specified, generate a figure showing the projection
# on the user-specified line
if (any(proj_v != 0)) {

  # Add an artificial offset to better show the process
  proj_v <- proj_v - c(line_offset, 0)

  # Find the unit vector
  p0 <- c(0, t(as.matrix(c(1, 0))) %*% as.matrix(proj_v))
  p1 <- c(1, t(as.matrix(c(1, 1))) %*% as.matrix(proj_v))
  uv <- p1 - p0
  uv_mod <- sqrt(sum(uv * uv))
  uv_norm <- uv / uv_mod

  # Define the projection matrix
  p_matrix <- (uv_norm %*% t(uv_norm)) / c(t(uv_norm) %*% uv_norm)

  # Project the centers of mass on the lines
  p1_prj <- p_matrix %*% as.matrix(c(p1_x1_wm, p1_x2_wm) - p0) + p0
  p2_prj <- p_matrix %*% as.matrix(c(p2_x1_wm, p2_x2_wm) - p0) + p0

  # Calculate the projected distance
  dist_wm_prj <- sqrt((p1_prj[1] - p2_prj[1])^2 + (p1_prj[2] - p2_prj[2])^2)

  # Generate the image showing how the distance is evaluated
  ggplot() +
    coord_fixed() +
    xlab("A-P axis, um") +
    ylab("D-V axis, um") +
    geom_polygon(
      data = end_pts_hull.plane[lo, ],
      aes(x = 0.008 * X1, y = 0.008 * X2),
      alpha = 0.3
    ) +
    geom_segment(aes(
      x = 0.008 * (p1_x1_wm - seg_len),
      xend = 0.008 * (p1_x1_wm + seg_len),
      y = 0.008 * p1_x2_wm,
      yend = 0.008 * p1_x2_wm
    ),
    col = "red", size = 3
    ) +
    geom_segment(aes(
      y = 0.008 * (p1_x2_wm - seg_len),
      yend = 0.008 * (p1_x2_wm + seg_len),
      x = 0.008 * p1_x1_wm,
      xend = 0.008 * p1_x1_wm
    ),
    col = "red", size = 3
    ) +
    geom_segment(aes(
      x = 0.008 * (p2_x1_wm - seg_len),
      xend = 0.008 * (p2_x1_wm + seg_len),
      y = 0.008 * p2_x2_wm,
      yend = 0.008 * p2_x2_wm
    ),
    col = "blue", size = 3
    ) +
    geom_segment(aes(
      y = 0.008 * (p2_x2_wm - seg_len),
      yend = 0.008 * (p2_x2_wm + seg_len),
      x = 0.008 * p2_x1_wm,
      xend = 0.008 * p2_x1_wm
    ),
    col = "blue", size = 3
    ) +
    geom_abline(aes(
      slope = proj_v[2],
      intercept = 0.008 * proj_v[1]
    ),
    col = "blue", size = 2
    ) +
    geom_segment(aes(
      x = 0.008 * p1_prj[1],
      xend = 0.008 * p1_x1_wm,
      y = 0.008 * p1_prj[2],
      yend = 0.008 * p1_x2_wm
    ),
    col = "red",
    size = 1,
    linetype = "dotted"
    ) +
    geom_segment(aes(
      x = 0.008 * p2_prj[1],
      xend = 0.008 * p2_x1_wm,
      y = 0.008 * p2_prj[2],
      yend = 0.008 * p2_x2_wm
    ),
    col = "blue",
    size = 1,
    linetype = "dotted"
    ) +
    geom_segment(aes(
      x = 0.008 * p1_prj[1],
      xend = 0.008 * p2_prj[1],
      y = 0.008 * p1_prj[2],
      yend = 0.008 * p2_prj[2]
    ),
    col = "yellow", size = 1
    ) +
    annotate("text",
      label = toString(paste(round(0.008 * dist_wm_prj, 2), "um")),
      x = 0.008 * min(p1_x1_wm, p2_x1_wm),
      y = 0.008 * (proj_v[1] - 1000),
      size = 5, col = "black"
    ) +
    ggtitle(paste("Weighted median distance", post_type1, "(red) <=>", post_type2, "(blue)")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 0),
      axis.text.y = element_text(face = "bold", color = "black", size = 15, angle = 0),
      plot.title = element_text(face = "bold", color = "black", size = 12),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15)
    )
}
