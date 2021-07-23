# Created on Sun Jul 20 10:31:02 2021
# Name:    labula_ranking.R
# Purpose: ...
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - ...
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
# Completely umbiased approach (first)
# - find the best separating line (SVM) for each post pairs with respect to pre
# - find the median line of the above
# - project on the perpendicular to the median and do box plots and evaluate
#   distance between the weighted medians
# - use these values to create a correlation matrix
# - Optimize the code. There are a lot of repeated tasks that can be optimized
#   by creating functions and storing the result in appropriate data structures

# Import required libraries
library(tidyverse)
library(natverse)
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
  con <- readRDS("data/hemibrain_con0.rds")
}
if (!exists('nlist')) {
  nlist <- readRDS("data/nlist1.rds")
}

# Source local files
source('aux_functions.R')






###############################################################################
# Define items to analyze here
pre_type <- 'LC4'

# Top (# of synapses) of post to consider
top <- 20

# Evaluate best separator using only anti-parallel?
use_anti <- TRUE
anti_threshold <- -0.3  # Max threshold allowed in gradient correlation

# Evaluate correlation matrix for all post pairs
evaluate_all <- TRUE

# Define a plane separating lobula from glomerulus (it will be used to only
# consider the lobulat in the analysis)
c(a, b, c, d) %<-% c(0.5, 0.37, -0.19, -5200)  # Original -5200; Above somas: -7000
###############################################################################







# Normalize the vector normal to the plane
w <- c(a, b, c)
w_mod <- sqrt(sum(w *w))
w_norm <- w / w_mod
w0 <- d

# Extract the top post in terms of synapses
top_posts <- con %>%
             filter(pre.type == pre_type, post.type != pre_type) %>%
             group_by(pre.type, post.type) %>%
             dplyr::count() %>%
             ungroup() %>%
             arrange(desc(n)) %>%
             slice(1:top) %>%
             pull(post.type)

# Extract all pre.type synapses
pre <- con %>% filter(pre.type==pre_type)

# Evaluate the full combinations
com_full <- combn(top_posts, 2)

# If only using anti-parallel, extract the info
if (use_anti == TRUE) {
  # Join the post in a data.frame
  for (i in 1:length(top_posts)) {
    cnt <- pre %>%
      filter(post.type == top_posts[[i]]) %>%
      group_by(pre.bodyID) %>%
      dplyr::count()
    if (i == 1) {
      posts <- cnt
    } else {
      posts <- merge(posts, cnt, by="pre.bodyID", all=TRUE)
    }
  }

  # Drop the pre.bodyID column
  posts <- posts[-1]

  # Change the names of columns to the top posts
  colnames(posts) <- c(top_posts)

  # Change NaN to zeros
  posts[is.na(posts)] <- 0

  # Evaluate the Pearson's correlation matrix
  pcorr <- cor(posts, method="pearson", use="complete.obs")

  # Create a data frame with the results
  ut <- upper.tri(pcorr)
  pcorr_df <- data.frame(row=rownames(pcorr)[row(pcorr)[ut]],
                         column=rownames(pcorr)[col(pcorr)[ut]],
                         cor=(pcorr)[ut])
  pcorr_df <- pcorr_df[order(pcorr_df$cor), ]

  # Extract the anticorrelated pairs below the threshold
  com <- t(as.matrix(pcorr_df %>%
                     filter(cor < anti_threshold) %>%
                     select(row, column)))
} else {
  # Evaluate all possible combinations without repetition
  com <- com_full
}



# Posts ----
# A matrix to store the line info
ms <- matrix(nrow=ncol(com), ncol=4)

# For each pair of posts evaluate the centroids of end points and project on
# the plane. Find the weighted centroids on the plane and the segment connecting
# these centroids. Store coordinates of the center of the segment as well and
# the coefficients of the line containing the segment.
for (c in 1:ncol(com)) {
# for (c in 1:1) {
    # Extract names
  p1 <- com[1, c]
  p2 <- com[2, c]

  # Some info for the user
  cat(paste('Evaluating:', p1, 'vs.', p2, '\n'))

  # Extract bodyIDs of pre synapses with post.type1
  post1 <- pre %>%
           filter(post.type==p1) %>%
           group_by(pre.bodyID) %>%
           dplyr::count()

  # Extract bodyIDs of pre synapses with post.type2
  post2 <- pre %>%
           filter(post.type==p2) %>%
           group_by(pre.bodyID) %>%
           dplyr::count()

  # Merge the two posts with the counts and replace NaN with zero
  post <- merge(post1,
                post2,
                by='pre.bodyID',
                all=TRUE,
                suffix=c('.post1', '.post2'))

  # Clear potential NaN
  post[is.na(post)] <- 0

  # Get the neuron list from the posts bodyIDs
  n_list <- nlist[as.character(post$pre.bodyID)]

  # Select preserved end points for each neuron
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

  # Project the preserved end points on the plane.
  end_pts.proj <- plane.proj(end_pts, w_norm, -w0)

  # Find the centroid of the end points
  center <- colMeans(end_pts.proj)

  # Calculate centroid of lobula end points
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

  # Project the centroid of the end points on the surface
  ctrs.proj <- plane.proj(ctrs, w_norm, -w0)

  # Define coordinate system in the plane
  # Use the projection of the z-axis on the plane to determine one of the new 2D
  # axes
  pointY <- c(plane.proj(t(center + c(0, 0, 2000)), w_norm, -w0))
  pointY <- pointY - center
  unitY <- - pointY / sqrt(sum(pointY * pointY))

  pointX <- cross(w_norm, pointY)
  unitX <- - pointX / sqrt(sum(pointX * pointX))

  # Centroids
  Xcm <- as.matrix(sweep(ctrs.proj, 2, center)) %*% unitX
  Ycm <- as.matrix(sweep(ctrs.proj, 2, center)) %*% unitY
  ctrs.plane <- data.frame(cbind(Xcm, Ycm))

  # Join the counts of synapses with the centroid dataset
  ctrs.plane['pre.bodyID'] <- row.names(ctrs.plane)
  ctrs.plane <- merge(ctrs.plane, post, by='pre.bodyID', all=TRUE) %>% na.omit()

  # Calculate weighted medians
  p1_x1_wm <- weightedMedian(ctrs.plane$X1, w=ctrs.plane$n.post1, interpolate=FALSE)
  p1_x2_wm <- weightedMedian(ctrs.plane$X2, w=ctrs.plane$n.post1, interpolate=FALSE)
  p2_x1_wm <- weightedMedian(ctrs.plane$X1, w=ctrs.plane$n.post2, interpolate=FALSE)
  p2_x2_wm <- weightedMedian(ctrs.plane$X2, w=ctrs.plane$n.post2, interpolate=FALSE)

  # Calculate unit vector of line connecting the two centroids
  uv <- c(p2_x1_wm, p2_x2_wm) - c(p1_x1_wm, p1_x2_wm)
  uv_mod <- sqrt(sum(uv *uv))
  uv_norm <- uv / uv_mod

  # Calculate the midpoint between the two centroid and use as origin
  origin <- 0.5 * (c(p2_x1_wm, p2_x2_wm) + c(p1_x1_wm, p1_x2_wm))

  # Store origin
  ms[c, 1:2] <- origin

  # Project the centroids onto the line
  ctrs.line <- ctrs.plane %>%
               select(-c('X1', 'X2')) %>%
               cbind(X=line.proj(ctrs.plane %>% select(X1, X2), uv_norm, origin)) %>%
               rename(X1=X.1,
                      X2=X.2)

  # If it is the first time around, plot the centers of mass
  if (c == 1) {
    plot(ctrs.plane %>% select(X1, X2), asp=1)
  }

  # Plot the points of the line
  # points(ctrs.line %>% select(X1, X2), pch=20, col='#ff0000', cex=0.5)
  reg <- lm(ctrs.line$X2 ~ ctrs.line$X1)
  if (is.na(reg$coefficients[2])) {
    abline(v=ctrs.line$X1[1], col='#ff000030')
  } else {
    abline(reg$coefficients, col='#ff000030')
  }

  # Store coefficients
  ms[c, 3:4] <- reg$coefficients
}

# Show the median line passing through the median of the center points of the
# segments connecting the centroids of each pair
# Calculate the medians of centroids and line coefficients
med <- colMedians(ms)

# Define and plot the origin
origin <- med[1:2]
points(origin[1], origin[2], pch=20, col='blue')

# Adjust the intercept to force line to pass through the median centroid
med[3] <- med[3] - ((origin[1] * med[4] + med[3]) - origin[2])

# Show the projecting line in dashed blue
abline(med[3:4], col='blue', lty='dashed')

# Evaluate the ortogonal line through the median centroid
per_coef <- -1 / med[4]
per_intr <- -origin[1] * per_coef + origin[2]

# Show the separating line in solid blue
abline(c(per_intr, per_coef), col='blue')

# Calculate unit vector of projection line
uv <- c(origin[1], origin[2]) - c(origin[1] + 100, (origin[1] + 100) * med[4] + med[3])
uv_mod <- sqrt(sum(uv *uv))
uv_norm <- uv / uv_mod

# Second round of projections ----
# Create an empty matrix to store the results
dist_mtrx <- matrix(nrow=top, ncol=top)
rownames(dist_mtrx) <- c(top_posts)
colnames(dist_mtrx) <- c(top_posts)

# Reproject the pairs on the median line evaluated above and calculate the
# distance between the weighted medians.
# Check if we are using all pairs in the evaluation
if (evaluate_all == TRUE) {
  eval_com <- com_full
} else {
  eval_com <- com
}

for (c in 1:ncol(eval_com)) {
# for (c in 1:1) {
  # Extract names
  p1 <- eval_com[1, c]
  p2 <- eval_com[2, c]

  # Some info for the user
  cat(paste('Evaluating:', p1, 'vs.', p2, '\n'))

  # Extract bodyIDs of pre synapses with post.type1
  post1 <- pre %>%
    filter(post.type==p1) %>%
    group_by(pre.bodyID) %>%
    dplyr::count()

  # Extract bodyIDs of pre synapses with post.type2
  post2 <- pre %>%
    filter(post.type==p2) %>%
    group_by(pre.bodyID) %>%
    dplyr::count()

  # Merge the two posts with the counts and replace NaN with zero
  post <- merge(post1,
                post2,
                by='pre.bodyID',
                all=TRUE,
                suffix=c('.post1', '.post2'))

  # Clear potential NaN
  post[is.na(post)] <- 0

  # Get the neuron list from the posts bodyIDs
  n_list <- nlist[as.character(post$pre.bodyID)]

  # Select preserved end points for each neuron
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

  # Project the preserved end points on the plane.
  end_pts.proj <- plane.proj(end_pts, w_norm, -w0)

  # Find the centroid of the end points
  center <- colMeans(end_pts.proj)

  # Calculate centroid of lobula end points
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

  # Project the centroid of the end points on the surface
  ctrs.proj <- plane.proj(ctrs, w_norm, -w0)

  # Define coordinate system in the plane
  # Use the projection of the z-axis on the plane to determine one of the new 2D
  # axes
  pointY <- c(plane.proj(t(center + c(0, 0, 2000)), w_norm, -w0))
  pointY <- pointY - center
  unitY <- - pointY / sqrt(sum(pointY * pointY))

  pointX <- cross(w_norm, pointY)
  unitX <- - pointX / sqrt(sum(pointX * pointX))

  # Centroids
  Xcm <- as.matrix(sweep(ctrs.proj, 2, center)) %*% unitX
  Ycm <- as.matrix(sweep(ctrs.proj, 2, center)) %*% unitY
  ctrs.plane <- data.frame(cbind(Xcm, Ycm))

  # Join the counts of synapses with the centroid dataset
  ctrs.plane['pre.bodyID'] <- row.names(ctrs.plane)
  ctrs.plane <- merge(ctrs.plane, post, by='pre.bodyID', all=TRUE) %>% na.omit()

  # Project the centroids onto the 2D line
  ctrs.line <- ctrs.plane %>%
    select(-c('X1', 'X2')) %>%
    cbind(X=line.proj(ctrs.plane %>% select(X1, X2), uv_norm, origin)) %>%
    rename(X1=X.1,
           X2=X.2)

  # Project the centroids on the 1D line
  ctrs.1d <- ctrs.line %>%
  select(-c('X1', 'X2')) %>%
  cbind(X=as.matrix(sweep(ctrs.line %>% select(X1, X2), 2, origin)) %*% uv_norm)

  # Evaluated weighted median
  p1_wm <- weightedMedian(ctrs.1d$X, w=ctrs.1d$n.post1, interpolate=FALSE)
  p2_wm <- weightedMedian(ctrs.1d$X, w=ctrs.1d$n.post2, interpolate=FALSE)

  # Distance between the weighted means
  dist <- abs(p1_wm - p2_wm)
  cat(paste('  Median distance:', dist, '\n\n'))

  # Update the distance matrix entry (symmetric)
  dist_mtrx[p1, p2] <- dist
  dist_mtrx[p2, p1] <- dist
}

# Show distance matrix
corrplot(dist_mtrx,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         type='upper',  # Only upper diagonal
         order='alpha',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         diag=FALSE,  # No diagonal
         tl.col='black')  # Color of the labels in black

