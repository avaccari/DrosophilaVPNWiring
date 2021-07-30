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
# - A 2d plot comparing the projections of the centroids for the two
#   posts within the projection of the lobula on the selection plane. The plot
#   also shows the weighted median (and weighted box plots) for the centers of
#   mass (weighted by the number of synapses with pre)
# - A 2d plot showing the distance between the weighted medians.
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


# Import required libraries
library(tidyverse)
library(natverse)
library(cetcolor)
library(zeallot)
library(pracma)
library(ks)
library(matrixStats)
library(egg)


# Clean everything up ----
# Except the connection and skeleton files if they are already loaded.
items <- ls()
items <- items[items != 'con']  # Comment this line to reload the connection file
items <- items[items != 'nlist']  # Comment this line to reload the skeleton file
rm(list=items)

# Close any open plotting windows
while (dev.cur() > 1) { dev.off() }
while (rgl.cur() > 0) { rgl.close() }

# Load datasets
if (!exists('con')) {
  con <- readRDS("hemibrain_con0.rds")
}
if (!exists('nlist')) {
  nlist <- readRDS("nlist1.rds")
}

# Source local files
source('aux_functions.R')








###############################################################################
# Define items to analyze here
pre_type <- 'LC4'
post_type1 <- 'DNp02'  # Red
post_type2 <- 'DNp11'  # Blue

# Define a plane separating lobula from glomerulus (it will be used to only
# consider the lobulat in the analysis)
c(a, b, c, d) %<-% c(0.5, 0.37, -0.19, -5200)  # Original -5200; Above somas: -7000

# Using perceptually uniform color map
# More options here:
# https://cran.r-project.org/web/packages/cetcolor/vignettes/cet_color_schemes.html
c_map_single <- 'r2'  # Map to use for individual post plot
c_size_single <- 100  # Number of colors in the map
c_map_both <- 'd9'  # Map to use for both posts together
c_size_both <- 100  # Number of colors in the map

# Plot window size
win_siz <- 500
###############################################################################







# Normalize the vector normal to the plane
w <- c(a, b, c)
w_mod <- sqrt(sum(w *w))
w_norm <- w / w_mod
w0 <- d

# Pre ----
# Extract all pre.type synapses
pre <- con %>% filter(pre.type==pre_type)

# Post ----
# Extract bodyIDs of pre synapses with post.type1
post1 <- pre %>%
         filter(post.type==post_type1) %>%
         group_by(pre.bodyID) %>%
         dplyr::count()

# Extract bodyIDs of pre synapses with post.type2
post2 <- pre %>%
         filter(post.type==post_type2) %>%
         group_by(pre.bodyID) %>%
         dplyr::count()

# Merge the two posts with the counts
post <- merge(post1,
              post2,
              by='pre.bodyID',
              all=TRUE,
              suffix=c('.post1', '.post2'))

# Clear potential NaN
post[is.na(post)] <- 0

# Get the neuron list from the posts bodyIDs
n_list <- nlist[as.character(post$pre.bodyID)]

# 3D analysis ----
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Grab the coordinates of the end points of the neurons which are on the lobula
# side of the plane
# Plot the exluded points in gray
for (i in 1:length(n_list)) {
  neu <- n_list[[i]]
  ep <- neu$d[neu$EndPoints, ]
  pts <- ep %>%
         mutate(LO = w_norm[1]*X + w_norm[2]*Y + w_norm[3]*Z + w0) %>%
         filter(LO >= 0) %>%
         select(X, Y, Z)
  if (i == 1) {
    end_pts <- pts
  } else {
    end_pts <- rbind(end_pts, pts)
  }
}
plot3d(end_pts, col='gray', size=1)

# Plot the preserved points in red
for (i in 1:length(n_list)) {
  neu <- n_list[[i]]
  ep <- neu$d[neu$EndPoints, ]
  pts <- ep %>%
         mutate(LO = w_norm[1]*X + w_norm[2]*Y + w_norm[3]*Z + w0) %>%
         filter(LO < 0) %>%
         select(X, Y, Z)
  if (i == 1) {
    end_pts <- pts
  } else {
    end_pts <- rbind(end_pts, pts)
  }
}
plot3d(end_pts, col='red', size=1, add=TRUE)

# Plot the separating plane
planes3d(a=w_norm[1], b=w_norm[2], c=w_norm[3], d=w0, add=TRUE)

# Project the preserved end points on the plane in yellow. The convex hull of
# these points define the lobula we are analyzing.
end_pts.proj <- plane.proj(end_pts, w_norm, -w0)
plot3d(end_pts.proj, col='yellow', size=1, add=TRUE)

# Find the centroid of the end points projections and plot in white
center <- colMeans(end_pts.proj)
plot3d(center[1], center[2], center[3], col='white', size=10, add=TRUE)

# Calculate centroid of lobula end points and plot in blue
for (i in 1:length(n_list)) {
  neu <- n_list[[i]]
  pts <- neu$d[neu$EndPoints, ] %>%
         mutate(LO = w_norm[1]*X + w_norm[2]*Y + w_norm[3]*Z + w0) %>%
         filter(LO < 0) %>%
         select(X, Y, Z) %>%
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

plot3d(ctrs, col='blue', size=10, add=TRUE)

# Project the centroid of the end points on the surface in green
ctrs.proj <- plane.proj(ctrs, w_norm, -w0)
plot3d(ctrs.proj, col='green', size=10, add=TRUE)

# 2D analysis ----
# Define coordinate system in the plane
# Use the projection of the z-axis on the plane to determine one of the new 2D
# axes
pointY <- c(plane.proj(t(center + c(0, 0, 2000)), w_norm, -w0))
pointY <- pointY - center
unitY <- - pointY / sqrt(sum(pointY * pointY))
arrow3d(center, center + 2000 * unitY, thikness=1, col='white')

pointX <- cross(w_norm, pointY)
unitX <- - pointX / sqrt(sum(pointX * pointX))
arrow3d(center, center + 2000 * unitX, thikness=1, col='white')

# Project data on 2D plane ----
# All considered endpoints
Xep <- as.matrix(sweep(end_pts.proj, 2, center)) %*% unitX
Yep <- as.matrix(sweep(end_pts.proj, 2, center)) %*% unitY
end_pts.plane <- data.frame(cbind(Xep, Yep))

# Calculate and plot the convex hull
lo <- chull(end_pts.plane)

# Centroids
Xcm <- as.matrix(sweep(ctrs.proj, 2, center)) %*% unitX
Ycm <- as.matrix(sweep(ctrs.proj, 2, center)) %*% unitY
ctrs.plane <- data.frame(cbind(Xcm, Ycm))

# Join the counts of synapses with the centroid dataset
ctrs.plane['pre.bodyID'] <- row.names(ctrs.plane)
ctrs.plane <- merge(ctrs.plane, post, by='pre.bodyID', all=TRUE) %>% na.omit()

# 2d maps of synapses-weighted centroid projections in the lobula ----
# Create the color maps
col_single <- cet_pal(c_size_single, name=c_map_single)

# Evaluate parameters for the box plots
dist_box <- 0.1 * abs(max(min(end_pts.plane$X1), min(end_pts.plane$X2)))
width_box <- 0.75 * dist_box
X1_rng <- abs(max(end_pts.plane$X1) - min(end_pts.plane$X1))
X2_rng <- abs(max(end_pts.plane$X2) - min(end_pts.plane$X2))
seg_len <- 0.05 * min(X1_rng, X2_rng)

# Calculate weighted medians
p1_x1_wm <- weightedMedian(ctrs.plane$X1, w=ctrs.plane$n.post1, interpolate=FALSE)
p1_x2_wm <- weightedMedian(ctrs.plane$X2, w=ctrs.plane$n.post1, interpolate=FALSE)
p2_x1_wm <- weightedMedian(ctrs.plane$X1, w=ctrs.plane$n.post2, interpolate=FALSE)
p2_x2_wm <- weightedMedian(ctrs.plane$X2, w=ctrs.plane$n.post2, interpolate=FALSE)

# Distance between weighted medians
dist_wm = sqrt((p1_x1_wm - p2_x1_wm)^2 + (p1_x2_wm - p2_x2_wm)^2)

# Generate the plot for post1
proj1 <- ggplot() +
         coord_fixed() +
         xlab('A-P axis, um') +
         ylab('D-V axis, um') +
         geom_polygon(data=end_pts.plane[lo, ],
                      aes(x=0.008 * X1, y=0.008 * X2),
                      alpha=0.3) +
         #geom_point(data=end_pts.plane,
          #          aes(x=0.008 * X1, y=0.008 * X2),
           #         size=0.1,
            #        alpha=0.2) +
         geom_point(data=ctrs.plane,
                    aes(x=0.008 * X1, y=0.008 * X2, colour=n.post1, size=n.post1)) +
         scale_color_gradientn(colours=col_single) +
         #scale_size_continuous(range=c(0, 7)) +
         geom_boxplot(data=ctrs.plane,
                      aes(x=0.008 * X1,
                          y=0.008 * (min(end_pts.plane$X2) - dist_box),
                          weight=n.post1),
                      width=0.008 * width_box,
                      notch=TRUE) +
         geom_boxplot(data=ctrs.plane,
                      aes(y=0.008 * X2,
                          x=0.008 * (min(end_pts.plane$X1) - dist_box),
                          weight=n.post1),
                      width=0.008 * width_box,
                      notch=TRUE) +
         geom_segment(aes(x=0.008 * (p1_x1_wm - seg_len),
                          xend=0.008 * (p1_x1_wm + seg_len),
                          y=0.008 * p1_x2_wm,
                          yend=0.008 * p1_x2_wm), col="black", size=2) +
         geom_segment(aes(y=0.008 * (p1_x2_wm - seg_len),
                          yend=0.008 * (p1_x2_wm + seg_len),
                          x=0.008 * p1_x1_wm,
                          xend=0.008 * p1_x1_wm), col="black", size=2) +
         ggtitle(paste(pre_type, '>', post_type1, '(lobula projection)')) +
         theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.x = element_text(face="bold", color="black", size=15, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=15, angle=0),
        plot.title = element_text(face="bold", color="black", size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))

# Generate the plot for post2
proj2 <- ggplot() +
         coord_fixed() +
         xlab('A-P axis, um') +
         ylab('D-V axis, um') +
         geom_polygon(data=end_pts.plane[lo, ],
                      aes(x=0.008 * X1, y=0.008 * X2),
                      alpha=0.3) +
         #geom_point(data=end_pts.plane,
          #          aes(x=0.008 * X1, y=0.008 * X2),
           #         size=0.1,
            #        alpha=0.2) +
         geom_point(data=ctrs.plane,
                    aes(x=0.008 * X1, y=0.008 * X2, colour=n.post2, size=n.post2)) +
         scale_color_gradientn(colours=col_single) +
         #scale_size_continuous(range=c(0, 7)) +
         geom_boxplot(data=ctrs.plane,
                      aes(x=0.008 * X1,
                          y=0.008 * (min(end_pts.plane$X2) - dist_box),
                          weight=n.post2),
                      width=0.008 * width_box,
                      notch=TRUE) +
         geom_boxplot(data=ctrs.plane,
                      aes(y=0.008 * X2,
                          x=0.008 * (min(end_pts.plane$X1) - dist_box),
                          weight=n.post2),
                      width=0.008 * width_box,
                      notch=TRUE) +
        geom_segment(aes(x=0.008 * (p2_x1_wm - seg_len),
                          xend=0.008 * (p2_x1_wm + seg_len),
                          y=0.008 * p2_x2_wm,
                          yend=0.008 * p2_x2_wm), col="black", size=2) +
         geom_segment(aes(y=0.008 * (p2_x2_wm - seg_len),
                          yend=0.008 * (p2_x2_wm + seg_len),
                          x=0.008 * p2_x1_wm,
                          xend=0.008 * p2_x1_wm), col="black", size=2) +
         ggtitle(paste(pre_type, '>', post_type2, '(lobula projection)')) +
         theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.x = element_text(face="bold", color="black", size=15, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=15, angle=0),
        plot.title = element_text(face="bold", color="black", size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))

# Combine and show the plots
grid.arrange(proj1, proj2, ncol=2)

# Generate a single plot that shows the location of the two centroids
ggplot() +
  coord_fixed() +
  xlab('A-P axis, um') +
  ylab('D-V axis, um') +
  geom_polygon(data=end_pts.plane[lo, ],
               aes(x=0.008 * X1, y=0.008 * X2),
               alpha=0.3) +
  #geom_point(data=end_pts.plane,
   #          aes(x=0.008 * X1, y=0.008 * X2),
    #         size=0.5,
     #        alpha=0.1) +
  #geom_point(data=ctrs.plane,
   #          aes(x=0.008 * X1, y=0.008 * X2, size=n.post1),
    #         colour='red',
     #        shape=1,
      #       stroke=1.5) +
  #geom_point(data=ctrs.plane,
   #          aes(x=0.008 * X1, y=0.008 * X2, size=n.post2),
    #         colour='blue',
     #        shape=1,
      #       stroke=1.5) +
  geom_boxplot(data=ctrs.plane,
               aes(x=0.008 * X1,
                   y=0.008 * (min(end_pts.plane$X2) - dist_box),
                   weight=n.post1),
               width=0.008 * width_box,
               notch=TRUE,
               col='red', size=1) +
  geom_boxplot(data=ctrs.plane,
               aes(y=0.008 * X2,
                   x=0.008 * (min(end_pts.plane$X1) - dist_box),
                   weight=n.post1),
               width=0.008 * width_box,
               notch=TRUE,
               col='red', size=1) +
  geom_segment(aes(x=0.008 * (p1_x1_wm - seg_len),
                   xend=0.008 * (p1_x1_wm + seg_len),
                   y=0.008 * p1_x2_wm,
                   yend=0.008 * p1_x2_wm),
               col='red', size=3) +
  geom_segment(aes(y=0.008 * (p1_x2_wm - seg_len),
                   yend=0.008 * (p1_x2_wm + seg_len),
                   x=0.008 * p1_x1_wm,
                   xend=0.008 * p1_x1_wm),
               col='red', size=3) +
  geom_boxplot(data=ctrs.plane,
               aes(x=0.008 * X1,
                   y=0.008 * (min(end_pts.plane$X2) - 2 * dist_box),
                   weight=n.post2),
               width=0.008 * width_box,
               notch=TRUE,
               col='blue', size=1) +
  geom_boxplot(data=ctrs.plane,
               aes(y=0.008 * X2,
                   x=0.008 * (min(end_pts.plane$X1) - 2 * dist_box),
                   weight=n.post2),
               width=0.008 * width_box,
               notch=TRUE,
               col='blue', size=1) +
  geom_segment(aes(x=0.008 * (p2_x1_wm - seg_len),
                   xend=0.008 * (p2_x1_wm + seg_len),
                   y=0.008 * p2_x2_wm,
                   yend=0.008 * p2_x2_wm),
               col='blue', size=3) +
  geom_segment(aes(y=0.008 * (p2_x2_wm - seg_len),
                   yend=0.008 * (p2_x2_wm + seg_len),
                   x=0.008 * p2_x1_wm,
                   xend=0.008 * p2_x1_wm),
               col='blue', size=3) +
  geom_segment(aes(x=0.008 * p1_x1_wm,
                   xend=0.008 * p2_x1_wm,
                   y=0.008 * p1_x2_wm,
                   yend=0.008 * p2_x2_wm),
               col='yellow', size=1) +
  annotate('text',
           label=toString(paste(round(0.008 * dist_wm, 2), 'um')),
           x=0.008 * min(p1_x1_wm, p2_x1_wm),
           y=0.008 * min(p1_x2_wm, p2_x2_wm),
           size=5, col="black", hjust=1.1) +
  ggtitle(paste('Weighted median distance', post_type1, '(red) <=>', post_type2, '(blue)')) +
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.x = element_text(face="bold", color="black", size=15, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=15, angle=0),
        plot.title = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))


