# Created on Mon Sep 13 20:41:47 2021
# Name:    glomerulus_ranking_scan_thresh.R
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
# - there is still a discrepancy between the results from threshold scan and 
#   actual ranking (for both glomerulus and lobula). Not extremely important, 
#   but curious why this is the case
# - And can you please also run the scan for lobula in LPLC2 with top25 posts?

# Import required libraries
library(tidyverse)
library(zeallot)
library(natverse)
library(e1071)
library(corrplot)
library(matrixStats)

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
top <- 25  # 25 for LC4 and 20 for LPLC2

# Scale SVM data
scaleSVM <- TRUE

# Evaluate best separator using only anti-parallel?
use_anti <- TRUE
# Parameters for the threshold scan
# anti_threshold_min is chosen as min in the pcorr dataframe
anti_threshold_step <- 0.01

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

# Get the post coordinates
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




# Load the data about the optimal separation planes for each pair
filename <- paste('data/', pre_type, '-glo_optimal_planes_each_pair-scaledSVM=', toString(scaleSVM), '.Rda', sep='')
if (file.exists(filename)) {
    load(filename)
  } else {
    # Evaluate the optimal plane for all the possible combinations and store in
    # a dataframe
    opt_planes <- data.frame()
    
    # Evaluate each pair
    for (c in 1:ncol(com_full)) {
      # for (c in 1:2) {
      # Extract names
      p1 <- com_full[1, c]
      p2 <- com_full[2, c]
      
      # Some info for the user
      cat(paste('Evaluating:', p1, 'vs.', p2, '\n'))
      
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
                       scale=scaleSVM)
      
      # Find the plane coefficients
      w <- t(svm_model$coefs) %*% svm_model$SV
      w_mod <- sqrt(sum(w *w))
      w_norm <- w/w_mod
      # w0 <- svm_model$rho/w_mod
      
      # Since we don't care about the offset but only the orientation of the
      # plane, force it to pass thorugh the centroids of the post coors
      center <- colMeans(post.coors %>% select('X', 'Y', 'Z'))
      w0 <- center %*% t(w_norm)

      # Store coefficients in the dataframe
      opt_planes <- rbind(opt_planes,
                          c(p1, p2, c(w_norm, w0)))
      # Some info for the user
      cat('  Optimal plane coeff:', w_norm, w0, '\n')
    }
    
    # Rename columns
    colnames(opt_planes) <- c('p1', 'p2', 'a', 'b', 'c', 'offset')
    
    # Save the data on disk
    save(opt_planes, file=filename)
}


# Sort the p1 and p2 so that p1 < p2 and then sort by p1.
# This will allow for a full merging
opt_planes[with(opt_planes, p1 > p2),1:2] <- opt_planes[with(opt_planes, p1 > p2),2:1]
opt_planes <- opt_planes[with(opt_planes, order(p1, p2)), ]









# Extract top posts and their counts
posts <- all_posts_cnt[, names(all_posts_cnt) %in% top_posts]

# Evaluate the Spearman's correlation matrix
pcorr <- cor(posts, method="spearman", use="complete.obs")

# Choose the minimum and maximum anticorrelation threshold
anti_threshold_min = ceiling(10 * log10(1/anti_threshold_step) * min(pcorr))/(10 * log10(1/anti_threshold_step))
anti_threshold_max = ceiling(10 * log10(1/anti_threshold_step) * max(pcorr))/(10 * log10(1/anti_threshold_step))

# convert to data frame with the results and write each pair as a row entry
ut <- upper.tri(pcorr)
pcorr_df <- data.frame(row=rownames(pcorr)[row(pcorr)[ut]],
                       column=rownames(pcorr)[col(pcorr)[ut]],
                       cor=(pcorr)[ut])

# Sort the data on the correlation values
pcorr_df <- pcorr_df[order(pcorr_df$cor), ]

# Change the column names
colnames(pcorr_df) <- c('p1', 'p2', 'corr')

# Sort the p1 and p2 so that p1 < p2 and then sort by p1.
# This will allow for a full merging
pcorr_df[with(pcorr_df, p1 > p2),1:2] <- pcorr_df[with(pcorr_df, p1 > p2), 2:1]
pcorr_df <- pcorr_df[with(pcorr_df, order(p1, p2)), ]


# Merge the correlation value with the optimal planes
opt_planes <- merge(opt_planes, pcorr_df, by=c('p1', 'p2'), all = TRUE)




# Create a dataframe to store the data
rvsth <- data.frame()

# Scan by varying the number of top posts to consider
for (anti_threshold in seq(anti_threshold_min, anti_threshold_max, anti_threshold_step)) {

  # Some info for the user
  cat(paste('anti_threshold considered:', anti_threshold))
  
  
  # If only using anti-parallel to evaluate projection line, extract the info
  if (use_anti == TRUE) {
    
    # Extract the anticorrelated pairs below the threshold
    com <- t(as.matrix(opt_planes %>%
                         filter(corr < anti_threshold) %>%
                         select(p1, p2)))
  } else {
    # Evaluate all possible combinations without repetition
    com <- com_full
  }
  
  # Some info for the user
  cat(paste('\n  pairs to consider:', ncol(com)))
  
  
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
   
    # Store coefficients in the matrix
    ms[c, ] <- as.numeric(opt_planes[opt_planes$p1==p1 & opt_planes$p2==p2, 3:6])
  }
  
  
  
  
  
  
  # Projection line construction ----
  # Image to show how the final plane is evaluated
  # Calculate the median plane and normalize
  med_plane <- colMedians(ms)
  med_plane_mod <- sqrt(sum(med_plane[1:3] * med_plane[1:3]))
  med_plane <- med_plane[1:3] / med_plane_mod
  
  
  # Find all the synapses for the evaluated posts
  post.coors <- pre.glo %>% 
    filter(post.type %in% top_posts) %>%
    pull(post.coors) %>%
    xyzmatrix() %>%
    as_tibble()
  
  center <- colMeans(post.coors)
  offset <- center %*% med_plane[1:3]

  # Show the median plane equation
  # Print the coefficient
  cat("\n  Median plane coefficients ax + by + cx + d = 0 (a, b, c, d):", med_plane[1:3], offset, '\n')
  

  
  
  

  
  
  
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
  dist_mtrx_srt <- dist_mtrx[order(rownames(dist_mtrx)), order(colnames(dist_mtrx))]
  
  # Convert the matrices in dataframes
  dist_mtrx_df <- data.frame(row=rownames(dist_mtrx_srt)[row(dist_mtrx_srt)[ut]],
                             column=rownames(dist_mtrx_srt)[col(dist_mtrx_srt)[ut]],
                             cor=(dist_mtrx_srt)[ut])
  
  # Change the column names
  colnames(dist_mtrx_df) <- c('p1', 'p2', 'dist')
  
  # Sort the p1 and p2 so that p1 < p2. This will allow for a full merging
  dist_mtrx_df[with(dist_mtrx_df, p1 > p2),1:2] <- dist_mtrx_df[with(dist_mtrx_df, p1 > p2),2:1]
  dist_mtrx_df <- dist_mtrx_df[with(dist_mtrx_df, order(p1, p2)), ]
  
  
  
  
  
  # Merge the matrices by rows and columns
  merged <- merge(opt_planes, dist_mtrx_df, by=c('p1','p2'))
  
  # Evaluate the correlation between the two
  corr <- cor(c(merged$corr), c(merged$dist), method='pearson')
  
  # Store in dataframe
  rvsth <- rbind(rvsth, c(anti_threshold, corr))
  
  # Some info for the user
  cat(paste("  Persons's coefficient:", corr, '\n'))
}  

# Rename columns
colnames(rvsth) <- c('anti_threshold', 'r')

# Save the dataframe to disk
save(rvsth, file=paste('output/glo_r_vs_antithresh-',
                       pre_type,
                       '-thresh=',
                       toString(anti_threshold_min),
                       '-',
                       toString(anti_threshold_max),
                       '-step=',
                       toString(anti_threshold_step),
                       '-fix_top=',
                       toString(top),
                       '-scaledSVM=',
                       toString(scaleSVM),
                       '.Rda',
                       sep=''))

# Plot the results
ggplot() +
  geom_line(data=rvsth,
             aes(x=anti_threshold,
                 y=r))
