# Created on Tue May 11 18:04:27 2021
# Name:    glomerulus.R
# Purpose: Evaluates the spatial separation of synapses within the corresponding 
# glomerulus
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates three plots:
# - one 3D plot showing the location of the synapses for the two posts, the best
#   (or user-specified) plane separating the synapses for the two posts, and the
#   projections of the synapses on the line perpendicular to the plane passing
#   by the center of mass of the projections of the synapses on the plane.
# - one 2D plot showing the KDE, rug, and violin plot derived from the
#   projections of the synapses on the line with the origin defined by the 
#   intersection with the plane. It also shows the distance between the median 
#   of the two distributions in physical units.
# - one 2D plot showing the ECDFs derived from the projections of the synapses
#   on the line with the origin defined by theintersection with the plane and
#   the calculated Kolmogorov-Smirnov two-sample test distance.

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

# Import required libraries
library(tidyverse)
library(plotly)
library(natverse)
library(e1071)
library(zeallot)

# clean almost everything up ----
# Except the connection file if it is already loaded. 
items <- ls()
items <- items[items != 'con']  # Comment this line to reload the connection file
rm(list=items)

# Close any open plotting windows
while (dev.cur() > 1) { dev.off() }
while (rgl.cur() > 0) { rgl.close() }

# Load synapses tables ----
if (!exists('con')) {
  con <- readRDS("data/hemibrain_con0.rds")
}

# Functions ----
# Project on plane function
plane.proj <- function(points, plane.normal, plane.offset) {
  dist <- as.matrix(points) %*% plane.normal
  proj <- points - (dist - plane.offset) %*% plane.normal
}

# Project on a line function
line.proj <- function(points, line.unit.vector, origin) {
  transl <- sweep(points, 2, origin)
  proj <- as.matrix(transl) %*% line.unit.vector
  proj <- sweep(drop(proj %o% line.unit.vector), 2, origin, FUN='+')
}

# Distance from the plane function
plane.dist <- function(points, plane.normal, origin) {
  transl <- sweep(points, 2, origin)
  dist <- as.matrix(transl) %*% plane.normal
}






###############################################################################
# Define items to analyze here
pre_type <- 'LC4'
post_type1 <- 'DNp02' # Blue
post_type2 <- 'DNp05' # Red

# Define a plane separating lobula from glomerulus (it will be used to only 
# consider the glomerulus in the analysis)
c(a, b, c, d) %<-% c(1, 0, 0, -9000)

# Define the separation plane to use in the analysis. The default is the plane
# best separating DNp02 from DNp11 synapses in the LC4 glomerulus. If all values
# are zeros, the plane will be evaluated from the data: best plane separating 
# the synapses for the two items listed above.
# plane <- c(0.971986, -2.217776, 3.01801, 13742.1)  # Original found by SVM
plane <- c(-0.2462962, 0.5820498, -0.7749556, -13324.9)  # Median plane from glomerulus_ranking.R
# plane <- c(0, 0, 0, 0)

# Plot window size
win_siz <- 500
###############################################################################








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

# Extract pre synapses with post.type1
post1 <- pre.glo %>% filter(post.type==post_type1)
post1.coors <- post1 %>% 
               pull(post.coors) %>% 
               xyzmatrix() %>%
               as_tibble()
               
# Extract pre synapses with post.type2
post2 <- pre.glo %>% filter(post.type==post_type2)
post2.coors <- post2 %>% 
               pull(post.coors) %>% 
               xyzmatrix() %>%
               as_tibble()

# Check if the plane is defined
if (all(plane == 0)) {
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
} else {
  w <- plane[0:3]
  w_mod <- sqrt(sum(w *w))
  w_norm <- w/w_mod
  w0 <- plane[4]
}

# Print the coefficients
cat("Plane coefficients ax + by + cx + d = 0 (a, b, c, d):", w, w0)

# 3D analysis ----
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Plot synapses
plot3d(post1.coors, col='blue', size = 5, add=TRUE)
plot3d(post2.coors, col='red', size=5,  add=TRUE)

# Create 3d plane
planes3d(a=w_norm[1], b=w_norm[2], c=w_norm[3], d=-w0, add=TRUE)

# Project on the synapses on the plane
post1.coors.proj.plane <- plane.proj(post1.coors, w_norm, w0)
# plot3d(post1.coors.proj.plane, col='cyan', add=TRUE)
post2.coors.proj.plane <- plane.proj(post2.coors, w_norm, w0)
# plot3d(post2.coors.proj.plane, col='orange', add=TRUE)

# Find the center of mass of all projections and plot
center <- colMeans(rbind(post1.coors.proj.plane, post2.coors.proj.plane))
plot3d(center[1], center[2], center[3], col='white', size=10, add=TRUE)

# Calculate the projections onto a line perpedicular to the plane through the
# center of mass
post1.coors.proj.line <- line.proj(post1.coors, w_norm, center)
plot3d(post1.coors.proj.line, col='cyan', add=TRUE)
post2.coors.proj.line <- line.proj(post2.coors, w_norm, center)
plot3d(post2.coors.proj.line, col='orange', add=TRUE)

# Draw the line perpendicular to the plane an passing by the center
endpt <- rbind(center -2000 * w_norm, center + 2000 * w_norm)
lines3d(endpt[, 1], endpt[, 2], endpt[, 3], col='gray', add=TRUE)

# 2D analysis ----
# Project on 1D line (x-axis)
post1.coors.line <- as.data.frame(plane.dist(post1.coors, w_norm, center))
post2.coors.line <- as.data.frame(plane.dist(post2.coors, w_norm, center))

# Evaluate medians
post1.coors.median <- median(post1.coors.line[, 1])
post2.coors.median <- median(post2.coors.line[, 1])

# Show results
ggplot() + 
  geom_violin(data=post1.coors.line,
              aes(x=V1*0.008, y=matrix(0, 1, length(V1))),
              color='blue',
              trim=FALSE) +
  geom_rug(data=post1.coors.line,
           aes(x=V1*0.008, y=matrix(0, 1, length(V1))), 
           color='blue',
           sides='b') +
  geom_boxplot(data=post1.coors.line, 
               aes(x=V1*0.008, y=matrix(0, 1, length(V1))),
               width=0.3,
               notch=TRUE,
               color='blue', size=1) +
  geom_segment(aes(x=post1.coors.median*0.008, y=0, xend=post1.coors.median*0.008, yend=0.5), 
             color='black', size=1)+
  geom_violin(data=post2.coors.line,
              aes(x=V1*0.008, y=matrix(1, 1, length(V1))),
              color='red',
              trim=FALSE) +
  geom_rug(data=post2.coors.line,
           aes(x=V1*0.008, y=matrix(0, 1, length(V1))), 
           color='red',
           sides='t') +
  geom_boxplot(data=post2.coors.line, 
               aes(x=V1*0.008, y=matrix(1, 1, length(V1))),
               width=0.3,
               notch=TRUE,
               color='red', size=1) +
  theme(axis.text.x = element_text(face="bold", color="black", size=13, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=13, angle=0),
        plot.title = element_text(face="bold", color="black", size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))+
  geom_segment(aes(x=post2.coors.median*0.008, y=1, xend=post2.coors.median*0.008, yend=0.5), 
               color='black', size=1) +
  geom_segment(aes(x=post1.coors.median*0.008, y=0.5, xend=post2.coors.median*0.008, yend=0.5)) +
  #annotate('text',
   #        label=toString(paste(round(0.008*abs(post1.coors.median - post2.coors.median), 2), 'um')),
    #       x=0.5 * (post1.coors.median + post2.coors.median), 
     #      y=0.55) +
  ggtitle('Distance between centroids in the glomerulus') +
  xlab('Distance from the separation plane, um') +
  ylab('') + 
  theme(axis.text.y=element_blank(),
        
        axis.ticks.y=element_blank(),
        plot.title=element_text(hjust=0.5))

# Statistical analysis ----
# Check if th two dataset are normally distributed using Shapiro-Wilk normality
# test. Null-hypostesis: the sample distribution is normal
s1 <- shapiro.test(post1.coors.line$V1)
s2 <- shapiro.test(post2.coors.line$V1)

# If both p values are greater than 0.05, then we can consider the distribution
# normal and use a t-test otherwise, we can use the nonparametric Wilcoxon test.
# For this test, the alternative hypotesis is that one distribution is 
# stochastically greater than the other (which I think might get closer to what
# we are interested in). 
# An alternative could be the D-value from the two-sample Kolmogorov-Smirnov 
# test which measures the level of difference between the two distributions.
# In either cases, I'm not sure the p-value actually gives an accurate
# representation. Maybe the W (or D for KS) better tracks what we want to show.
if (all(c(s1$p.value, s2$p.value) >= 0.05)) {
  tt <- t.test(post1.coors.line$V1, post2.coors.line$V1)
  cat("p-value for t-test:", tt$p.value)
} else {
  wt <- wilcox.test(post1.coors.line$V1, post2.coors.line$V1)
  print(wt)
  kt <- ks.test(post1.coors.line$V1, post2.coors.line$V1)
  print(kt)
}

# Plot the ECDFs of the two posts and show the Kolmogorov-Smirnov two-sample
# distance.
# Evaluate the ECDFs
cdf1 <- ecdf(post1.coors.line[, 1])
cdf2 <- ecdf(post2.coors.line[, 1])

# Find the min and max of the post to define the x-axis range
xmin <- min(rbind(post1.coors.line, post2.coors.line))
xmax <- max(rbind(post1.coors.line, post2.coors.line))
xseq <- seq(xmin, xmax)

# Find the max distance and its location
id_max <- which.max(abs(cdf1(xseq)-cdf2(xseq)))
ks_max = xseq[id_max]
ks_dist <- abs(cdf1(ks_max) - cdf2(ks_max))


# Plot the two ECDFs and the KS2ST distance
ggplot() +
  stat_ecdf(data=post1.coors.line,
            aes(x=V1),
            color='blue', size=2) +
  geom_rug(data=post1.coors.line,
           aes(x=V1, y=matrix(0, 1, length(V1))), 
           color='blue',
           sides='b') +  
  stat_ecdf(data=post2.coors.line,
            aes(x=V1),
            color='red', size=2) +
  geom_rug(data=post2.coors.line,
           aes(x=V1, y=matrix(0, 1, length(V1))), 
           color='red',
           sides='t') +
  geom_segment(aes(x=ks_max, y=cdf1(ks_max), xend=ks_max, yend=cdf2(ks_max)), size=2) +
  annotate('text',
           label=toString(paste('D[KS] ==',round(ks_dist, 2))),
           x=0.25 * (xmax - xmin) + ks_max, 
           y=min(cdf1(ks_max) , cdf2(ks_max)),
           parse=TRUE, size=7) +
  ggtitle('Kolmogorov-Smirnov two-sample distance') +
  xlab('Samples') +
  ylab('ECDFs') +
  theme(axis.text.x = element_text(face="bold", color="black", size=13, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=13, angle=0),
        plot.title = element_text(face="bold", color="black", size=13),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))
  theme(plot.title=element_text(hjust=0.5))
  

