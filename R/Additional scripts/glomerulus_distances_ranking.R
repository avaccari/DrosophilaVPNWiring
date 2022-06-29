# Created on Sat Jul 17 20:11:24 2021
# Name:    glomerulus_distances_ranking.R
# Purpose: Selects the top n (number of synapses) post for a given pre and
#          evaluates the median of the distances between all possible pairs
#          formed between any of the posts synapses.
#          This is similar to the average linkage distance between the two
#          synapses clusters.
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates ... :
# - ...

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
# - Evaluate distance in physical coordinates (assuming: 8x8x8 nm^3 voxels
#   Needs verification)
# - Create a 3d plot to show the split of the synapses based on the plane

# Import required libraries
library(tidyverse)
library(corrplot)
library(fields)
library(matrixStats)
library(natverse)
library(zeallot)

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
  con <- readRDS("data/hemibrain_con0.rds")
}
if (!exists('nlist')) {
  nlist <- readRDS("data/nlist1.rds")
}









###############################################################################
# Define items to analyze here
pre_type <- 'LC4'

# Top (# of synapses) of post to consider
top <- 25

# Define a plane separating lobula from glomerulus (it will be used to only
# consider the glomerulus in the analysis)
c(a, b, c, d) %<-% c(1, 0, 0, -9000)

# Plot window size
win_siz <- 500
###############################################################################








# Extract the top post in terms of synapses
top_posts <- con %>%
             filter(pre.type == pre_type, post.type != pre_type) %>%
             group_by(pre.type, post.type) %>%
             dplyr::count() %>%
             ungroup() %>%
             arrange(desc(n)) %>%
             slice(1:top) %>%
             pull(post.type)

# Evaluate all possible combinations without repetition
com <- combn(top_posts, 2)

# Extract all pre.type synapses
pre <- con %>% filter(pre.type==pre_type)

# Filter for glomerulus synapses
pre.glo <- pre %>%
           filter(pre %>%
           pull(post.coors) %>%
           xyzmatrix() %>%
           as_tibble() %>%
           mutate(GLO = a*X + b*Y + c*Z + d > 0) %>%
           select(GLO))

# Create an empty matrix to store the results
acld_median_mtrx <- matrix(0, nrow=top, ncol=top)
rownames(acld_median_mtrx) <- c(top_posts)
colnames(acld_median_mtrx) <- c(top_posts)

acld_mean_mtrx <- matrix(0, nrow=top, ncol=top)
rownames(acld_mean_mtrx) <- c(top_posts)
colnames(acld_mean_mtrx) <- c(top_posts)

median_mtrx <- matrix(0, nrow=top, ncol=top)
rownames(median_mtrx) <- c(top_posts)
colnames(median_mtrx) <- c(top_posts)

mean_mtrx <- matrix(0, nrow=top, ncol=top)
rownames(mean_mtrx) <- c(top_posts)
colnames(mean_mtrx) <- c(top_posts)

min_mtrx <- matrix(0, nrow=top, ncol=top)
rownames(min_mtrx) <- c(top_posts)
colnames(min_mtrx) <- c(top_posts)

max_mtrx <- matrix(0, nrow=top, ncol=top)
rownames(max_mtrx) <- c(top_posts)
colnames(max_mtrx) <- c(top_posts)

dunn_mtrx <- matrix(0, nrow=top, ncol=top)
rownames(dunn_mtrx) <- c(top_posts)
colnames(dunn_mtrx) <- c(top_posts)

ctrd_mtrx <- matrix(0, nrow=top, ncol=top)
rownames(ctrd_mtrx) <- c(top_posts)
colnames(ctrd_mtrx) <- c(top_posts)

# For each neuron in top_posts extract the pre synapses and evaluate dist stats
for (i in top_posts) {
  # Some info for the user
  cat(paste('Evaluating:', i, '\n'))

    # Extract pre synapses with post i
  post <- pre.glo %>% filter(post.type==i)
  post.coors <- post %>%
                pull(post.coors) %>%
                xyzmatrix()

  # Some info for the user
  cat(paste('  Number of synapses:', nrow(post.coors), '\n'))

  # Calculate distances between all synapses
  dist <- rdist(post.coors, post.coors)

  # Extract all the non-zero distances
  dist <- dist[upper.tri(dist)]

  # Calculate the median of all the distances
  med <- 0.008 * median(dist)

  # Calculate the mean of all the distances
  mean <- 0.008 * mean(dist)

  # Calculate the median of all the distances
  min <- 0.008 * min(dist)

  # Calculate the mean of all the distances
  max <- 0.008 * max(dist)

  # Some info for the user
  cat(paste('  Min distance:', min, '\n'))
  cat(paste('  Median distance:', med, '\n'))
  cat(paste('  Mean distance:', mean, '\n'))
  cat(paste('  Max distance:', max, '\n\n'))

  # Update the median matrix entry (symmetric)
  median_mtrx[i, i] <- med

  # Update the mean matrix entry (symmetric)
  mean_mtrx[i, i] <- mean

  # Update the min matrix entry (symmetric)
  min_mtrx[i, i] <- min

  # Update the max matrix entry (symmetric)
  max_mtrx[i, i] <- max
}

# For each pair of posts extract the pre synapses and evaluate dist stats
for (i in 1:ncol(com)) {
  # Extract names
  p1 <- com[1, i]
  p2 <- com[2, i]

  # Some info for the user
  cat(paste('Evaluating:', p1, 'vs.', p2, '\n'))

  # Extract pre synapses with post.type1
  post1 <- pre.glo %>% filter(post.type==p1)
  post1.coors <- post1 %>%
                 pull(post.coors) %>%
                 xyzmatrix()

  # Extract pre synapses with post.type2
  post2 <- pre.glo %>% filter(post.type==p2)
  post2.coors <- post2 %>%
                 pull(post.coors) %>%
                 xyzmatrix()

  # Some info for the user
  cat(paste('  Number of synapses: post1:', nrow(post1.coors), 'post2:', nrow(post2.coors), '\n'))



  # Average centroid linkage distance (ACLD) ----
  # Calculate the mean centroid of each post
  post1_ctrd_mean = colMeans(post1.coors)
  post2_ctrd_mean = colMeans(post2.coors)

  # Subtract the mean centroid from the opposite post
  dist1_mean <- sweep(post1.coors, 2, post2_ctrd_mean)
  dist2_mean <- sweep(post2.coors, 2, post1_ctrd_mean)

  # Calculate the distances
  dist1_mean <- sqrt(rowSums(dist1_mean * dist1_mean))
  dist2_mean <- sqrt(rowSums(dist2_mean * dist2_mean))

  # Calculate the mean of all the distances
  acld_mean <- 0.008 * mean(c(dist1_mean, dist2_mean))




  # Average centroid linkage distance (ACLD) (medians) ----
  # Calculate the median centroid of each post
  post1_ctrd_med = colMedians(post1.coors)
  post2_ctrd_med = colMedians(post2.coors)

  # Subtract the median centroid from the opposite post
  dist1_med <- sweep(post1.coors, 2, post2_ctrd_med)
  dist2_med <- sweep(post2.coors, 2, post1_ctrd_med)

  # Calculate the distances
  dist1_med <- sqrt(rowSums(dist1_med * dist1_med))
  dist2_med <- sqrt(rowSums(dist2_med * dist2_med))

  # Calculate the median of all the distances
  acld_med <- 0.008 * median(c(dist1_med, dist2_med))

  # Some info for the user
  cat(paste('  ACLD (means):', acld_mean, '\n'))
  cat(paste('  ACLD (medians):', acld_med, '\n'))

  # Distance between median centroids
  dist_ctrd <- post1_ctrd_med - post2_ctrd_med
  dist_ctrd <- 0.008 * sqrt(sum(dist_ctrd * dist_ctrd))



  # Other inter-class statistics ----
  # Calcualte the distance matrix between the two classes
  dist <- rdist(post1.coors, post2.coors)

  # Calculate the min of all the distances
  min <- 0.008 * min(dist)

  # Calculate the min of all the distances
  mean <- 0.008 * mean(dist)

  # Calculate the median of all the distances
  med <- 0.008 * median(dist)

  # Calculate the max of all the distances
  max <- 0.008 * max(dist)

  # Calculate the Dunn index
  dip1 <- min / max_mtrx[p1, p1]
  dip2 <- min / max_mtrx[p2, p2]

  # Some info for the user
  cat(paste('  Min distance:', min, '\n'))
  cat(paste('  Median distance:', med, '\n'))
  cat(paste('  Mean distance:', mean, '\n'))
  cat(paste('  Max distance:', max, '\n'))
  cat(paste('  Dunn index (p1):', dip1, '\n'))
  cat(paste('  Dunn index (p2):', dip2, '\n'))
  cat(paste('  Median centroids distance:', dist_ctrd, '\n\n'))


  # Update the ACLD median matrix entry (symmetric)
  acld_median_mtrx[p1, p2] <- acld_med
  acld_median_mtrx[p2, p1] <- acld_med

  # Update the ACLD mean matrix entry (symmetric)
  acld_mean_mtrx[p1, p2] <- acld_mean
  acld_mean_mtrx[p2, p1] <- acld_mean

  # Update the median matrix entry (symmetric)
  median_mtrx[p1, p2] <- med
  median_mtrx[p2, p1] <- med

  # Update the mean matrix entry (symmetric)
  mean_mtrx[p1, p2] <- mean
  mean_mtrx[p2, p1] <- mean

  # Update the min matrix entry (symmetric)
  min_mtrx[p1, p2] <- min
  min_mtrx[p2, p1] <- min

  # Update the max matrix entry (symmetric)
  max_mtrx[p1, p2] <- max
  max_mtrx[p2, p1] <- max

  # Update the Dunn index matrix
  dunn_mtrx[p1, p2] <- dip1
  dunn_mtrx[p2, p1] <- dip2

  # Update the median centroids matrix (symmetric)
  ctrd_mtrx[p1, p2] <- dist_ctrd
  ctrd_mtrx[p2, p1] <- dist_ctrd
}

# Show graphic results
# (There are a lot of options for this plot)
par(mfrow=c(3, 2))

corrplot(acld_mean_mtrx,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         type='upper',  # Only upper diagonal
         order='alpha',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         diag=FALSE,  # No diagonal
         tl.col='black')  # Color of the labels in black

corrplot(acld_median_mtrx,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         type='upper',  # Only upper diagonal
         order='alpha',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         diag=FALSE,  # No diagonal
         tl.col='black')  # Color of the labels in black



corrplot(mean_mtrx,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         type='upper',  # Only upper diagonal
         order='alpha',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         diag=FALSE,  # No diagonal
         tl.col='black')  # Color of the labels in black

corrplot(median_mtrx,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         type='upper',  # Only upper diagonal
         order='alpha',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         diag=FALSE,  # No diagonal
         tl.col='black')  # Color of the labels in black




corrplot(dunn_mtrx,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         type='full',  # Only upper diagonal
         order='alpha',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         diag=FALSE,  # No diagonal
         tl.col='black')  # Color of the labels in black

corrplot(ctrd_mtrx,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         type='upper',  # Only upper diagonal
         order='alpha',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         diag=FALSE,  # No diagonal
         tl.col='black')  # Color of the labels in black




