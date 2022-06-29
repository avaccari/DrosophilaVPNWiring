# Created on Mon Jan 24 19:33:27 2022
# Name:    labula_centroid_construction.R
# Purpose: Shows the construction for the lobula centroids
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - A 3d plot showing how the centroids are evaluated using body IDs specified
#   by the user.
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
post_type1 <- 'DNp02'  
post_type2 <- 'DNp11'  
highlight_IDs = c('1158187240', '1621357756')

# Grab the correct plane
# NOTE: the FIRST plane specified is also use as projection plane
planes <- get.plane(pre_type, 'lobula')

# Plot window size
win_siz <- 500
###############################################################################







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

# Get the desired neuron
n_highlight <- nlist[highlight_IDs]

# 3D analysis ----
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Add neuron skeletons in gray
plot3d(setdiff(n_list, n_highlight), soma=TRUE, col='light gray', add=TRUE)

# and highlighted neurons in black
plot3d(n_highlight, soma=TRUE, col='black', add=TRUE)

# and the planes
planes3d(a=planes[1, 1], b=planes[1, 2], c=planes[1, 3], d=planes[1, 4], col='blue', alpha=0.2, add=TRUE)
planes3d(a=planes[-1, 1], b=planes[-1, 2], c=planes[-1, 3], d=planes[-1, 4], alpha=0.2, add=TRUE)

# Show preserved end points for highlighted neuron
for (i in 1:length(n_highlight)) {
  neu <- n_highlight[[i]]
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
plot3d(end_pts, col='green', size=1, add=TRUE)

# Calculate centroid of highlighted neurosn end points and plot in red
for (i in 1:length(n_highlight)) {
  neu <- n_highlight[[i]]
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

# Plot the centroids in red
plot3d(ctrs, col='red', size=5, add=TRUE)

