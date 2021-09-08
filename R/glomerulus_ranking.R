# Created on Fri Jul 30 10:46:06 2021
# Name:    glomerulus_ranking.R
# Purpose: Evaluates the distance matrix of the centroid evaluated in
#          glomerulus.R evaluated according to the best separation plane
#          calculated from the top anti-parallel pairs of posts of a given pre.
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - A 3d plot identifying which synapses for a given pre are considered in the
#   analysis based on the selection planes
# - A set of 3d plots displaying the evaluated optimal plane for each post pair
# - A 3d plot showing how the median plane is evaluated
# - A set of 3d plots displaying the different pairs in relation to the median
#   plane
# - A plot showing the correlation between the distance matrix of the projected
#   centroid of each possible posts pair and the correlation between the 
#   synaptic connectivity for the same pair.
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
# - Find a way to use the scaled data with the SVM and then recover the equation
#   of the plane (especially the relationship between d and rho). This is going
#   to be more robust and determine the separating plane in more cases.
#   Since we do not care about the intercept, the method might just work without
#   too many problems. Verify that the projection on the perpendicular makes
#   sense and use that.
# - Decide if we want to drop the planes where the SVM does not converge
# - Code could be optimized by getting the various coordinates in a table and
#   then use the loop just to extract needed columns


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
pre_type <- 'LPLC2'

# Top (# of synapses) of post to consider
top <- 20  # 25 for LC4 and 20 for LPLC2

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

# Visualize all the synapses and the selection planes
open3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Add neuron skeletons (to visualize the somas)
plot3d(nlist[as.character(unique(pre %>% pull(pre.bodyID)))], soma=TRUE, col='light gray', add=TRUE)

# Add the synapses
pre_coors <- pre %>%
             pull(post.coors) %>%
             xyzmatrix()
plot3d(pre_coors, size=1, add=TRUE)

# and the planes
planes3d(a=planes[, 1], b=planes[, 2], c=planes[, 3], d=planes[, 4], alpha=0.2, add=TRUE)

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

plot3d(pre_kept, col='green', size=3, add=TRUE)

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

# Setup the 3d image
open3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))
cols <- ceiling(sqrt(ncol(com)))
rows <- ceiling(ncol(com)/cols)
mfrow3d(rows, cols, sharedMouse=TRUE)

# Evaluate each pair
for (c in 1:ncol(com)) {
  # for (c in 1:1) {
  # Extract names
  p1 <- com[1, c]
  p2 <- com[2, c]
  
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
                   scale=FALSE)

  # Find the plane coefficients
  w <- t(svm_model$coefs) %*% svm_model$SV
  w_mod <- sqrt(sum(w *w))
  w_norm <- w/w_mod
  w0 <- svm_model$rho/w_mod
  
  # Store coefficients in the matrix
  ms[c, ] <- c(w_norm, w0)

  # Print the coefficient
  cat("  Plane coefficients ax + by + cx + d = 0 (a, b, c, d):", w_norm, w0, '\n')

  
  # Plot the next pair plot
  next3d()
  plot3d(post1.coors, col='blue', add=TRUE)
  plot3d(post2.coors, col='red', add=TRUE)
  planes3d(a=w_norm[1], b=w_norm[2], c=w_norm[3], d=-w0, add=TRUE, alpha=0.2)
  legend3d
  axes3d(edges=c('x', 'y', 'z'),
         tick=FALSE,
         labels=FALSE)
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

# Plot the construction of the mean plane
# Since the offset of the plane is not relevant but we only care about the 
# orientation of the normal, the plane is forced to pass through the centroid of
# the synapses.
open3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))
plot3d(post.coors, col='green', alpha=1, add=TRUE)
center <- colMeans(post.coors)
plot3d(center[1], center[2], center[3], col='white', size=10, add=TRUE)
offset <- center %*% med_plane[1:3]
planes3d(a=ms[, 1], b=ms[, 2], c=ms[, 3], d=-ms[, 4], col='red', add=TRUE, alpha=0.2)
planes3d(a=med_plane[1], b=med_plane[2], c=med_plane[3], d=-offset, col='blue', add=TRUE)

# Show the median plane equation
# Print the coefficient
cat("\nMedian plane coefficients ax + by + cx + d = 0 (a, b, c, d):", med_plane[1:3], offset, '\n')

# Draw the line perpendicular to the mean plane an passing by the center
endpt <- rbind(center - 2000 * med_plane[1:3], center + 2000 * med_plane[1:3])
lines3d(endpt[, 1], endpt[, 2], endpt[, 3], col='black', lwd=5, add=TRUE)

# Add axes
axes3d(edges=NULL,
       tick=FALSE,
       labels=FALSE)





# Visualization of pairs with mean plane ----
# Replot with the fix plane (just for visualization. It can then be removed)
# Setup the 3d image
open3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))
cols <- ceiling(sqrt(ncol(com)))
rows <- ceiling(ncol(com)/cols)
mfrow3d(rows, cols, sharedMouse=TRUE)

# Evaluate each pair
for (c in 1:ncol(com)) {
  # for (c in 1:1) {
  # Extract names
  p1 <- com[1, c]
  p2 <- com[2, c]
  
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
  
  # Plot the next pair plot
  next3d()
  plot3d(post1.coors, col='blue', add=TRUE)
  plot3d(post2.coors, col='red', add=TRUE)
  planes3d(a=med_plane[1], b=med_plane[2], c=med_plane[3], d=-offset, add=TRUE, alpha=0.2)
  legend3d('topright', legend=paste(p1, 'vs.', p2))
  axes3d(edges=c('x', 'y', 'z'),
         tick=FALSE,
         labels=FALSE)
  
}






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
  cat(paste('Evaluating:', p1, 'vs.', p2, '\n'))
  
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
  cat(paste('  Median distance:', dist, '\n\n'))
  
  # Update the distance matrix entry (symmetric)
  dist_mtrx[p1, p2] <- dist
  dist_mtrx[p2, p1] <- dist
}

###############################################################################
# This section is needed only if you want the correlation plot to be sorted in
# the same way in all the scripts.
sort.pcorr <- order(rownames(pcorr))
pcorr <- pcorr[sort.pcorr, sort.pcorr]
pcorr.FPC <- corrMatOrder(pcorr, order='hclust')
sort.dist_mtrx <- order(rownames(dist_mtrx))
dist_mtrx <- dist_mtrx[sort.dist_mtrx, sort.dist_mtrx]
dist_mtrx_no <- dist_mtrx[pcorr.FPC, pcorr.FPC]
###############################################################################

# Show distance matrix
corrplot(0.008 * dist_mtrx_no,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         col=colorRampPalette(c("white", "darkblue"))(200),
         order='original',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         tl.col='black')  # Color of the labels in black




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


# Plot the results
ggplot() +
  xlim(-0.9, 0.9) +
  ylim (0, 7)+
  ylab("Distance between centroids in the glomerulus, um") +
  xlab("Glomerular connectivity correlation") +
  theme_classic()+
  geom_point(data=merged,
             aes(x=cor.x, y=(cor.y)*0.008),
             size=3,
             col="steelblue")+
  geom_smooth(data=merged,
              aes(x=cor.x, y=(cor.y)*0.008),
              method="lm",
              formula= 'y ~ x', 
              col="red", 
              se=TRUE) +
  annotate('text', 
           col="red", size=8,
           label=toString(paste("r =", round(corr, 2))),
           x=0.6 * max(merged$cor.x),
           y=0.8 * 0.008  * max(merged$cor.y),
           size=5)+
  theme(axis.text.x = element_text(face="bold", color="black", size=15, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=15, angle=0),
        plot.title = element_text(face="bold", color="blue", size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))




