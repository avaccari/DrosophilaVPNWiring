# Created on Mon Sep 13 20:41:47 2021
# Name:    glomerulus_ranking_scan_top.R
# Purpose: 
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
#

# Import required libraries
library(tidyverse)
library(zeallot)
library(natverse)
library(e1071)
library(corrplot)

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
source("R/aux_functions.R")








###############################################################################
# Define items to analyze here
pre_type <- 'LC4'

# Top (# of synapses) of posts to consider
top_min <- 25  # minimum number of top posts to consider
top_max <- 25  # maximum number of top posts to consider

# Evaluate best separator using only anti-parallel?
use_anti <- TRUE
anti_threshold <- -0.5  # Max threshold allowed in gradient correlation

# Evaluate correlation matrix for all post pairs
evaluate_all <- TRUE

# Grab the correct plane
planes <- get.plane(pre_type, 'glomerulus')

# Plot window size
win_siz <- 1000
###############################################################################









# Identify body all posts and body IDs based on pre ----
# Extract all pre.type synapses excluding with pre
pre <- con %>% 
  filter(pre.type==pre_type, post.type != pre_type)

# Add the synapses
pre_coors <- pre %>%
  pull(post.coors) %>%
  xyzmatrix()

# Filter for glomerulus synapses using all the planes
pre.glo <- pre
for (i in 1:nrow(planes)) {
  pre.glo <- pre.glo %>%
    filter(pre.glo %>%
             pull(post.coors) %>%
             xyzmatrix() %>%
             as_tibble() %>%
             mutate(GLO = planes[i, 5] * (planes[i, 1] * X + planes[i, 2] * Y + planes[i, 3] * Z + planes[i, 4]) > 0) %>%
             select(GLO))
}

# Show filtered synapses
pre_kept <- pre.glo %>%
  pull(post.coors) %>%
  xyzmatrix()

# Find the unique posts for a given pre
all_posts <- unique(pre %>%
                      group_by(post.type) %>%
                      pull(post.type))

# Create data frame with the counts for all the posts associated with the pre
for (i in 1:length(all_posts)) {
  cnt <- pre %>%
    filter(post.type == all_posts[[i]]) %>%
    group_by(pre.bodyID) %>%
    dplyr::count()
  if (i == 1) {
    all_posts_cnt <- cnt
  } else {
    all_posts_cnt <- merge(all_posts_cnt, cnt, by="pre.bodyID", all=TRUE)
  }
}

# Zero the NA
all_posts_cnt[is.na(all_posts_cnt)] <- 0

# Change the names of columns to the top posts
colnames(all_posts_cnt) <- c('pre.bodyID', all_posts)


# Create a dataframe to store the data
rvsnp <- data.frame(top=top_min:top_max, r=0)

# Scan by varying the number of top posts to consider
for (top in top_min:top_max) {

  # Some info for the user
  cat(paste('Top neuron considered:', top))
  
  # Determine the subset of combination to analyze ----
  # Extract the top post in terms of synapses
  top_posts <- pre %>%
    group_by(post.type) %>%
    dplyr::count() %>%
    ungroup() %>%
    arrange(desc(n)) %>%
    slice(1:top) %>%
    pull(post.type)
  
  # Evaluate the combinations of the top posts
  com_full <- combn(top_posts, 2)
  
  # Extract top posts and their counts
  posts <- all_posts_cnt[, names(all_posts_cnt) %in% top_posts]
  
  # Evaluate the Spearman's correlation matrix
  pcorr <- cor(posts, method="spearman", use="complete.obs")
  
  # convert to data frame with the results and write each pair as a row entry
  ut <- upper.tri(pcorr)
  pcorr_df <- data.frame(row=rownames(pcorr)[row(pcorr)[ut]],
                         column=rownames(pcorr)[col(pcorr)[ut]],
                         cor=(pcorr)[ut])
  
  # Sort the data on the correlation values
  pcorr_df <- pcorr_df[order(pcorr_df$cor), ]
  
  # If only using anti-parallel to evaluate projection line, extract the info
  if (use_anti == TRUE) {
    
    # Extract the anticorrelated pairs below the threshold
    com <- t(as.matrix(pcorr_df %>%
                         filter(cor < anti_threshold) %>%
                         select(row, column)))
  } else {
    # Evaluate all possible combinations without repetition
    com <- com_full
  }
  
  
  
  
  # Evaluate the Posts ----
  # A matrix to store the plane info
  ms <- matrix(nrow=ncol(com), ncol=4)
  
  # For each pair of posts, evaluate the synapses with pre, the ideal plane
  # separating them and show them in 3d
  
  # Evaluate each pair
  for (c in 1:ncol(com)) {
    # for (c in 1:1) {
    # Extract names
    p1 <- com[1, c]
    p2 <- com[2, c]
    
    # Some info for the user
    # cat(paste('Evaluating:', p1, 'vs.', p2, '\n'))
    
    # Extract pre synapses with post.type1
    post1 <- pre.glo %>% filter(post.type==p1)
    post1.coors <- post1 %>%
      pull(post.coors) %>%
      xyzmatrix() %>%
      as_tibble()
    
    # Extract pre synapses with post.type2
    post2 <- pre.glo %>% filter(post.type==p2)
    post2.coors <- post2 %>%
      pull(post.coors) %>%
      xyzmatrix() %>%
      as_tibble()
    
    # Use SVM to determine best separation plane ----
    # Combine the two sinapses with different labels
    post.coors <- bind_rows(post1.coors, post2.coors, .id="label")
    
    # Change label to integer
    post.coors$label <- as.integer(as.character(post.coors$label))
  
    # Train an SVM linear model with all the data without scaling
    svm_model <- svm(label ~ .,
                     data=post.coors,
                     type='C-classification',
                     kernel='linear',
                     scale=FALSE)
  
    # Find the plane coefficients
    w <- t(svm_model$coefs) %*% svm_model$SV
    w_mod <- sqrt(sum(w *w))
    w_norm <- w/w_mod
    w0 <- svm_model$rho/w_mod
    
    # Store coefficients in the matrix
    ms[c, ] <- c(w_norm, w0)
  
    # Print the coefficient
    # cat("  Plane coefficients ax + by + cx + d = 0 (a, b, c, d):", w_norm, w0, '\n')
  }
  
  
  
  
  
  
  # Projection line construction ----
  # Image to show how the final plane is evaluated
  # Calculate the median plane and normalize
  med_plane <- colMedians(ms)
  med_plane_mod <- sqrt(sum(med_plane[1:3] * med_plane[1:3]))
  med_plane <- med_plane[1:3] / med_plane_mod
  
  
  # Find all the synapses for the evaluated posts
  # Extract pre synapses with post.type1
  post.coors <- pre.glo %>% 
    filter(post.type %in% top_posts) %>%
    pull(post.coors) %>%
    xyzmatrix() %>%
    as_tibble()
  
  center <- colMeans(post.coors)
  offset <- center %*% med_plane[1:3]

    # Show the median plane equation
  # Print the coefficient
  # cat("\nMedian plane coefficients ax + by + cx + d = 0 (a, b, c, d):", med_plane[1:3], offset, '\n')
  

  
  
  

  
  
  
  # Second round of projections ----
  # Create an empty matrix to store the results
  dist_mtrx <- matrix(0, nrow=top, ncol=top)
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
  
  
  # Project the pairs and evaluate the distance between the weighted centroids
  for (c in 1:ncol(eval_com)) {
    # Extract names
    p1 <- eval_com[1, c]
    p2 <- eval_com[2, c]
    
    # Some info for the user
    # cat(paste('Evaluating:', p1, 'vs.', p2, '\n'))
    
    # Extract pre synapses with post.type1 and project on the 1D line
    post1 <- pre.glo %>% filter(post.type==p1)
    post1.coors <- post1 %>%
      pull(post.coors) %>%
      xyzmatrix() %>%
      as_tibble()
    post1.coors.line <- as.data.frame(plane.dist(post1.coors, med_plane[1:3], center))
    
    # Extract pre synapses with post.type2 and project on 1D line
    post2 <- pre.glo %>% filter(post.type==p2)
    post2.coors <- post2 %>%
      pull(post.coors) %>%
      xyzmatrix() %>%
      as_tibble()
    post2.coors.line <- as.data.frame(plane.dist(post2.coors, med_plane[1:3], center))
    
    # Evaluated weighted median
    post1.coors.median <- median(post1.coors.line[, 1])
    post2.coors.median <- median(post2.coors.line[, 1])
    
    # Distance between the weighted means
    dist <- abs(post1.coors.median - post2.coors.median)
    # cat(paste('  Median distance:', dist, '\n\n'))
    
    # Update the distance matrix entry (symmetric)
    dist_mtrx[p1, p2] <- dist
    dist_mtrx[p2, p1] <- dist
  }
  

  
  
  
  # Check the correlation between the distance matrix and the correlation matrix
  # from the gradient ranking
  
  # Sort in alphabetical order row and columns by name for the two matrices.
  pcorr_srt <- pcorr[order(rownames(pcorr)), order(colnames(pcorr))]
  dist_mtrx_srt <- dist_mtrx[order(rownames(dist_mtrx)), order(colnames(dist_mtrx))]
  
  # Convert the matrices in dataframes
  pcorr_df <- data.frame(row=rownames(pcorr_srt)[row(pcorr_srt)[ut]],
                         column=rownames(pcorr_srt)[col(pcorr_srt)[ut]],
                         cor=(pcorr_srt)[ut])
  dist_mtrx_df <- data.frame(row=rownames(dist_mtrx_srt)[row(dist_mtrx_srt)[ut]],
                             column=rownames(dist_mtrx_srt)[col(dist_mtrx_srt)[ut]],
                             cor=(dist_mtrx_srt)[ut])
  
  # Merge the matrices by rows and columns
  merged <- merge(pcorr_df, dist_mtrx_df, by=c('row','column'))
  
  # Evaluate the correlation between the two
  corr <- cor(c(merged$cor.x), c(merged$cor.y), method='pearson')
  
  # Store in dataframe
  rvsnp$r[[top==top]] <- corr
  
  # Some info for the user
  cat(paste("\n Persons's coefficient:", corr))
}  

# Save the dataframe to disk
save(rvsnp, file=paste('output/glo_r_vs_npairs -',
                       pre_type,
                       '- top =',
                       toString(top_min),
                       '-',
                       toString(top_max),
                       '- fix thresh =',
                       toString(anti_threshold), 
                       '.Rda'))
