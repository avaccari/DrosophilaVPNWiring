# Created on Sun Jul 20 10:31:02 2021
# Name:    labula_ranking.R
# Purpose: Evaluates the distance matrix of the weighted centroid evaluated in
#          lobula.R evaluated according to the best separation line calculated
#          from the top anti-parallel pairs of posts of a given pre.
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - A 3d plot identifying which end points for a given pre are considered in the
#   analysis based on the selection planes (projection plane is shown in blue)
# - A plot showing the segments connecting the weighted centroids of the top 
#   anti-parallel posts and the evaluated separationa and projection lines
# - A plot showing the midpoints of the segments above and the evaluated 
#   separation and projection lines
# - A plot showing the lines containing the segments above and the evaluated 
#   separation and projection lines
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
#

# Import required libraries
library(tidyverse)
library(natverse)
library(zeallot)
library(pracma)
library(ks)
library(matrixStats)
library(egg)
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
source('R/aux_functions.R')








###############################################################################
# Define items to analyze here
pre_type <- 'LC4'

# Top (# of synapses) of post to consider
top <- 25  # 25 for LC4 and 20 for LPLC2

# Evaluate best separator using only anti-parallel?
use_anti <- TRUE
anti_threshold <- -0.5  # Max threshold allowed in gradient correlation

# Evaluate correlation matrix for all post pairs
evaluate_all <- TRUE

# Grab the correct plane
# NOTE: the FIRST plane specified is also use as projection plane
planes <- get.plane(pre_type, 'lobula')

# Plot window size
win_siz <- 1000

###############################################################################







# Normalize the vector normal to the plane
w_norm <- planes[1, 1:3]
w0 <- planes[1, 4]



# Identify body all posts and body IDs based on pre ----
# Extract all pre.type synapses excluding with pre
pre <- con %>% 
  filter(pre.type==pre_type, post.type != pre_type)

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

# Get the list of all pre.bodyID
all_body_IDs <- all_posts_cnt[['pre.bodyID']]




# Determine coordinate system on the selected plane ----
# Get the list of all neurons defined by pre
all_nlist <- nlist[as.character(all_body_IDs)]

# Display all the neurons defined by pre
open3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Add neuron skeletons (to visualize the somas)
plot3d(all_nlist, soma=TRUE, col='light gray', add=TRUE)

# Show all endpoints
for (i in 1:length(all_nlist)) {
  neu <- all_nlist[[i]]
  ep <- neu$d[neu$EndPoints, ] %>%
    select(X, Y, Z)
  plot3d(ep, col='black', size=1, add=TRUE)
}

# and the planes
planes3d(a=planes[1, 1], b=planes[1, 2], c=planes[1, 3], d=planes[1, 4], col='blue', alpha=0.2, add=TRUE)
planes3d(a=planes[-1, 1], b=planes[-1, 2], c=planes[-1, 3], d=planes[-1, 4], alpha=0.2, add=TRUE)

# Select preserved end points for each neuron
for (i in 1:length(all_nlist)) {
  neu <- all_nlist[[i]]
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
plot3d(end_pts, col='green', size=3, add=TRUE)

# Project the preserved end points on the plane.
end_pts.proj <- plane.proj(end_pts, w_norm, -w0)

# Find the centroid of the end points
center <- colMeans(end_pts.proj)

# Define coordinate system in the plane
# Use the projection of the z-axis on the plane to determine one of the new 2D
# axes
pointY <- c(plane.proj(t(center + c(0, 0, 2000)), w_norm, -w0))
pointY <- pointY - center
unitY <- - pointY / sqrt(sum(pointY * pointY))

pointX <- cross(w_norm, pointY)
unitX <- - pointX / sqrt(sum(pointX * pointX))





# Evaluates the centroids for each body ID and project on the plane ----
# Calculate centroid of lobula end points
for (i in 1:length(all_nlist)) {
  neu <- all_nlist[[i]]
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

# SHow in red on the 3d plot
plot3d(ctrs, col='red', size=5, add=TRUE)

# Turn into a dataframe
ctrs <- data.frame(ctrs)

# Add body IDs
ctrs[['pre.bodyID']] <- lapply(all_nlist, '[[', 'bodyid')

# Merge with synapses counts and remove NA
ctrs <- merge(all_posts_cnt, ctrs, by='pre.bodyID', all=TRUE) %>%
  rename(X.3d=X,
         Y.3d=Y,
         Z.3d=Z) %>% 
  na.omit()

# Project the centroid of the end points on the plane
ctrs <- ctrs %>% 
  cbind(P=plane.proj(ctrs %>%
                       select(X.3d, Y.3d, Z.3d), w_norm, -w0)) %>%
  rename(X.proj=P.X.3d,
         Y.proj=P.Y.3d,
         Z.proj=P.Z.3d) %>% 
  na.omit()

# Calculate the coordinates on the plane
Xcm <- as.matrix(sweep(ctrs %>%
                         select(X.proj, Y.proj, Z.proj), 2, center)) %*% unitX
Ycm <- as.matrix(sweep(ctrs %>%
                          select(X.proj, Y.proj, Z.proj), 2, center)) %*% unitY
ctrs <- ctrs %>%
  cbind(X.plane=Xcm, Y.plane=Ycm) %>%
  na.omit()







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
# A matrix to store the line info
ms <- matrix(nrow=ncol(com), ncol=8)

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
  
  # Calculate weighted medians
  p1_x1_wm <- weightedMedian(ctrs$X.plane, w=ctrs[[p1]], interpolate=FALSE)
  p1_x2_wm <- weightedMedian(ctrs$Y.plane, w=ctrs[[p1]], interpolate=FALSE)
  p2_x1_wm <- weightedMedian(ctrs$X.plane, w=ctrs[[p2]], interpolate=FALSE)
  p2_x2_wm <- weightedMedian(ctrs$Y.plane, w=ctrs[[p2]], interpolate=FALSE)
  
  # Store segment points
  ms[c, 5:8] <- c(p1_x1_wm, p1_x2_wm, p2_x1_wm, p2_x2_wm)
  
  # Calculate unit vector of line connecting the two centroids
  uv <- c(p2_x1_wm, p2_x2_wm) - c(p1_x1_wm, p1_x2_wm)
  uv_mod <- sqrt(sum(uv *uv))
  uv_norm <- uv / uv_mod
  
  # Calculate the midpoint between the two centroid and use as origin
  origin <- 0.5 * (c(p2_x1_wm, p2_x2_wm) + c(p1_x1_wm, p1_x2_wm))
  
  # Store origin
  ms[c, 1:2] <- origin

  # Evaluate the projection line
  reg <- lm(c(p1_x2_wm, p2_x2_wm) ~ c(p1_x1_wm, p2_x1_wm))
  
  # Store coefficients
  ms[c, 3:4] <- reg$coefficients
}


# Evaluate the median line ----
# Calculate the medians of centroids and line coefficients
med <- colMedians(ms %>% na.omit())

# Define and plot the origin
origin <- med[1:2]

# Adjust the intercept to force line to pass through the median centroid
med[3] <- med[3] - ((origin[1] * med[4] + med[3]) - origin[2])

# Evaluate the ortogonal line through the median centroid
per_coef <- -1 / med[4]
per_intr <- -origin[1] * per_coef + origin[2]

# Calculate unit vector of projection line
uv <- c(origin[1], origin[2]) - c(origin[1] + 100, (origin[1] + 100) * med[4] + med[3])
uv_mod <- sqrt(sum(uv *uv))
uv_norm <- uv / uv_mod







# Generate the plots illustrating the process ----
# Convert the ms dataset into dataframe
ms <- data.frame(ms)

# Add columns with neuron names to ms
ms <- ms %>% cbind(p1=com[1,], p2=com[2,])

# Calculate the coordinates of the end points on the plane
Xep <- as.matrix(sweep(end_pts.proj, 2, center)) %*% unitX
Yep <- as.matrix(sweep(end_pts.proj, 2, center)) %*% unitY
end_pts.plane <- data.frame(cbind(Xep, Yep))

# Calculate and plot the convex hull for the lobula
lo <- chull(end_pts.proj)


# Segments connecting the weighted centroids
#ggplot() +             
 # coord_fixed() +
  #xlab('A-P') +
  #ylab('D-V') +
  #geom_polygon(data=end_pts.plane[lo, ],
   #            aes(x=0.008 * X1, y=0.008 * X2),
    #           alpha=0.2) +
  #geom_point(data=ctrs,
   #          aes(x=0.008 * X.plane, y=0.008 * Y.plane),
    #         shape=1) +
  #geom_segment(data=ms,
   #            aes(x=0.008 * X5, y=0.008 * X6, xend=0.008 * X7, yend=0.008 * X8),
    #           color='#ff000060') +
  #geom_text(data=ms,
   #         aes(x=0.008 * X5, y=0.008 * X6, label=p1)) +
  #geom_text(data=ms,
   #         aes(x=0.008 * X7, y=0.008 * X8, label=p2)) +
  #geom_point(data=data.frame(t(origin)),
   #          aes(x=0.008 * X1, y=0.008 * X2),
    #         color='#0000ff') +
  #ggtitle('Segments connecting weighted centroids') +
  #theme(plot.title=element_text(hjust=0.5))

# Segments connecting the weighted centroids
ggplot() +
  coord_fixed() +
  xlab('A-P axis, um') +
  ylab('D-V axis, um') +
  geom_polygon(data=end_pts.plane[lo, ],
               aes(x=0.008 * X1, y=0.008 * X2),
               alpha=0.3) +
  geom_point(data=ctrs,
             aes(x=0.008 * X.plane, y=0.008 * Y.plane),
             shape=1) +
  geom_segment(data=ms,
               aes(x=0.008 * X5, y=0.008 * X6, xend=0.008 * X7, yend=0.008 * X8),
               color='red', size=2) +
  #geom_text(data=ms,
   #         aes(x=0.008 * X5, y=0.008 * X6, label=p1)) +
  #geom_text(data=ms,
   #         aes(x=0.008 * X7, y=0.008 * X8, label=p2)) +
  geom_point(data=data.frame(t(origin)),
             aes(x=0.008 * X1, y=0.008 * X2),
             color='#0000ff') +
  geom_abline(data=data.frame(t(med)),
              aes(slope=X4, intercept=0.008 * X3),
              color='blue', size=1,
              linetype='dashed') +
  geom_abline(data=data.frame(t(c(per_intr, per_coef))),
              aes(slope=X2, intercept=0.008 * X1),
              color='blue', size=2) +
  ggtitle('Segments connectring weighted centroids (red)\n Median separation line (solid blue)\n Projection line (dashed blue)') +
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.x = element_text(face="bold", color="black", size=13, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=13, angle=0))+
  theme(axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))


# Midpoints of the segments
#ggplot() +
 # coord_fixed() +
  #xlab('A-P') +
  #ylab('D-V') +
  #geom_polygon(data=end_pts.plane[lo, ],
   #            aes(x=0.008 * X1, y=0.008 * X2),
    #           alpha=0.2) +
#  geom_point(data=ctrs,
 #            aes(x=0.008 * X.plane, y=0.008 * Y.plane),
  #           shape=1) +
  #geom_point(data=ms,
   #          aes(x=0.008 * X1, y=0.008 * X2),
    #         color='#ff0000',
     #        shape=16) +
  #geom_point(data=data.frame(t(origin)),
   #          aes(x=0.008 * X1, y=0.008 * X2),
    #         color='#0000ff') +
  #geom_abline(data=data.frame(t(med)),
   #           aes(slope=X4, intercept=0.008 * X3),
    #          color='#0000ff',
     #         linetype='dashed') +
  #geom_abline(data=data.frame(t(c(per_intr, per_coef))),
   #           aes(slope=X2, intercept=0.008 * X1),
    #          color='#0000ff') +
  #ggtitle('Midpoints of segments connectring weighted centroids') +
  #theme(plot.title=element_text(hjust=0.5))

# Lines containing the segments
ggplot() +
  coord_fixed() +
  xlab('A-P axis, um') +
  ylab('D-V axis, um') +
  geom_polygon(data=end_pts.plane[lo, ],
               aes(x=0.008 * X1, y=0.008 * X2),
               alpha=0.3) +
  geom_point(data=ctrs,
             aes(x=0.008 * X.plane, y=0.008 * Y.plane),
             shape=1) +
  geom_abline(data=ms,
              aes(slope=X4, intercept=0.008 * X3),
              color='red', size=1.5) +
  geom_point(data=data.frame(t(origin)),
             aes(x=0.008 * X1, y=0.008 * X2),
             color='#0000ff') +
  geom_abline(data=data.frame(t(med)),
              aes(slope=X4, intercept=0.008 * X3),
              color="blue", size=1.5,
              linetype='dashed') +
  geom_abline(data=data.frame(t(c(per_intr, per_coef))),
              aes(slope=X2, intercept=0.008 * X1),
              color='blue', size=2) +
  ggtitle('Lines connectring weighted centroids (red)\n median line (solid blue)\n projection line (dashed blue)') +
  theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.x = element_text(face="bold", color="black", size=13, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=13, angle=0))+
  theme(axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))






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

# Project all centroids on the median line
ctrs <- ctrs %>%
  cbind(X=line.proj(ctrs %>% 
                      select(X.plane, Y.plane), uv_norm, origin)) %>%
  rename(X.mline=X.1,
         Y.mline=X.2)

# Project the centroids on the 1D line
ctrs <- ctrs %>%
  cbind(X=as.matrix(sweep(ctrs %>%
                            select(X.mline, Y.mline), 2, origin)) %*%
          uv_norm)

# Project the pairs and evaluate the distance between the weighted centroids
for (c in 1:ncol(eval_com)) {
  # Extract names
  p1 <- eval_com[1, c]
  p2 <- eval_com[2, c]
  
  # Some info for the user
  cat(paste('Evaluating:', p1, 'vs.', p2, '\n'))
  
  # Evaluated weighted median
  p1_wm <- weightedMedian(ctrs$X, w=ctrs[[p1]], interpolate=FALSE)
  p2_wm <- weightedMedian(ctrs$X, w=ctrs[[p2]], interpolate=FALSE)
  
  # Distance between the weighted means
  dist <- abs(p1_wm - p2_wm)
  cat(paste('  Median distance:', dist, '\n\n'))
  
  # Update the distance matrix entry (symmetric)
  dist_mtrx[p1, p2] <- dist
  dist_mtrx[p2, p1] <- dist
}

# Show coefficients of the projected line used
cat('Projection line a and b:', med[3:4])


###############################################################################
# This section is needed only if you want the correlation plot to be sorted in
# the same way in all the scripts.
sort.pcorr <- order(rownames(pcorr))
pcorr <- pcorr[sort.pcorr, sort.pcorr]
pcorr.FPC <- corrMatOrder(pcorr, order='FPC')
sort.dist_mtrx <- order(rownames(dist_mtrx))
dist_mtrx <- dist_mtrx[sort.dist_mtrx, sort.dist_mtrx]
dist_mtrx_no <- dist_mtrx[pcorr.FPC, pcorr.FPC]
###############################################################################

# Show distance matrix
corrplot(0.008 * dist_mtrx_no,
         is.corr=FALSE,  # It is not a correlation matrix
         method='color',  # Color the background
         type='upper',  # Only upper diagonal
         order='original',  # Alphabetical
         addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         diag=FALSE,  # No diagonal
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
  ylim (0,50)+
  ylab("Distance between centroids in the lobula, mcm") +
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

