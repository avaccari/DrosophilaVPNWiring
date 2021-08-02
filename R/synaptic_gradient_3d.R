# Created on Wed May 12 00:03:40 2021
# Name:    synaptic_gradient_3d.R
# Purpose: Colors neurons skeletons based on the synaptic gradient
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates 5 plots:
# - the first two plots are gradients for single posts based on the number of
#   synapses with pre. The color scale is 0 to the value spacified by syn_max.
#   If syn_max is zero, then the scale is 0 to the max number of synapses for
#   each post independently.
# - the following two plots are gradients for single posts based on the number
#   of synapses with pre. The color scale is nomalized between the lowest number
#   and the max number of synapses between the two posts and the pre
# - the last plot shows both posts colored according to their relative number of
#   synapses with the pre. The color scale goes 0-1 and is evaluated for each
#   post neuron using:
#   0.5 * (1 + (#post1_syn - #post2_syn)/(#post1_syn + #post2_syn))

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

# TODO:
# - include LC4 neuron without synapses (use the pre.type to identify all the
#   LC4 synapses and then use the bodyID to select the skeleton).
# - add color bar with range as legend.


# Import required libraries
library(tidyverse)
library(natverse)
library(cetcolor)

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






###############################################################################
# Define items to analyze here
pre_type <- 'LC4'
post_type1 <- 'PVLP022' # Blue
post_type2 <- 'DNp11' # Red

# Using perceptually uniform color map
# More options here:
# https://cran.r-project.org/web/packages/cetcolor/vignettes/cet_color_schemes.html
c_map_single <- 'r2'  # Map to use for individual post plot
c_size_single <- 100  # Number of colors in the map
c_map_both <- 'r2'  # Map to use for both posts together
c_size_both <- 100  # Number of colors in the map

# Define the max to use when generating single gradient post plots. The colors
# will be scaled between zero and this max. If syn_max is set to zero, then the
# max number of synapses for each post. If syn_max is set and the number of
# synapses is larger than syn_max, the color scale is capped at syn_max.
syn_max <- 40

# Plot window size
win_siz <- 500
###############################################################################





# Pre ----
# Extract all pre.type synapses
pre <- con %>% filter(pre.type==pre_type)

# Post1 ----
# Identify pre with synapses on post1 and count
post1 <- pre %>%
         filter(post.type==post_type1) %>%
         group_by(pre.bodyID) %>%
         dplyr::count()

# Evaluate max and min number of synapses
post1.syn.max <- max(post1$n)
post1.syn.min <- min(post1$n)

# Add synapses info to neuron list
for (i in 1:nrow(post1)) {
  nlist[[toString(post1[[i, 1]])]]$post1count <- post1[[i, 2]]
}

# Post2 ----
# Identify pre with synapses on post1 and count
post2 <- pre %>%
         filter(post.type==post_type2) %>%
         group_by(pre.bodyID) %>%
         dplyr::count()

# Evaluate max and min number of synapses
post2.syn.max <- max(post2$n)
post2.syn.min <- min(post2$n)

# Add synapses info to neuron list
for (i in 1:nrow(post2)) {
  nlist[[toString(post2[[i, 1]])]]$post2count <- post2[[i, 2]]
}



# Display ----
# Creating colormap
col_single <- cet_pal(c_size_single, name=c_map_single)
col_both <- cet_pal(c_size_both, name=c_map_both)

# Evaluate max and min number of synapses in the two groups and scale
# color palette accordingly. This is a full scale realtive to the two,
# not an absolute scale referenced to a max number of synapses.
post.syn.min <- min(post1.syn.min, post2.syn.min)
post.syn.max <- max(post1.syn.max, post2.syn.max)
post.syn.rng <- post.syn.max - post.syn.min


# Plots of individual gradients normalized between 0 and syn_max ----
# Gradient plot for post1 (red) - 0 to syn_max
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

if (syn_max == 0) {
  syn1.max <- post1.syn.max
} else {
  syn1.max <- syn_max
}

bodyIDs.post1 <- post1 %>% pull(pre.bodyID)

for (bodyID in bodyIDs.post1) {
  neu <- nlist[[toString(bodyID)]]

  # Get the synapses count and correct for NULL values (only connection to one
  # of the pre)
  count1 <- neu$post1count
  if (is_null(count1)) {
    count1 <- post.syn.min
  }

  # Scale the synapses count to the size of the color scale
  count1 <- min(c_size_single, round(1 + count1 * (c_size_single - 1) / syn1.max))

  # Plot the neuron
  plot3d(neu,
         soma=TRUE,
         lwd=3,
         col=col_single[count1])
}

# Gradient plot for post2 (blue) - 0 to syn_max
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

if (syn_max == 0) {
  syn2.max <- post2.syn.max
} else {
  syn2.max <- syn_max
}

bodyIDs.post2 <- post2 %>% pull(pre.bodyID)

for (bodyID in bodyIDs.post2) {
  neu <- nlist[[toString(bodyID)]]

  # Get the synapses count and correct for NULL values (only connection to one
  # of the pre)
  count2 <- neu$post2count
  if (is_null(count2)) {
    count2 <- post.syn.min
  }

  # Scale the synapses count to the size of the color scale
  count2 <- min(c_size_single, round(1 + count2 * (c_size_single - 1) / syn2.max))

  # Plot the neuron
  plot3d(neu,
         soma=TRUE,
         lwd=3,
         col=col_single[count2])
}




# Plots of individual gradient normalized between post.syn.min and post.syn.max ----
# Gradient plot for post1 (red) - post.syn.min to post.syn.max
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

bodyIDs.post1 <- post1 %>% pull(pre.bodyID)

for (bodyID in bodyIDs.post1) {
  neu <- nlist[[toString(bodyID)]]

  # Get the synapses count and correct for NULL values (only connection to one
  # of the pre)
  count1 <- neu$post1count
  if (is_null(count1)) {
    count1 <- post.syn.min
  }

  # Scale the synapses count to the size of the color scale
  count1 <- round(1 + (count1 - post.syn.min) * (c_size_single - 1) / post.syn.rng)

  # Plot the neuron
  plot3d(neu,
         soma=TRUE,
         lwd=3,
         col=col_single[count1])
}

# Gradient plot for post2 (blue) - post.syn.min to post.syn.max
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

bodyIDs.post2 <- post2 %>% pull(pre.bodyID)

for (bodyID in bodyIDs.post2) {
  neu <- nlist[[toString(bodyID)]]

  # Get the synapses count and correct for NULL values (only connection to one
  # of the pre)
  count2 <- neu$post2count
  if (is_null(count2)) {
    count2 <- post.syn.min
  }

  # Scale the synapses count to the size of the color scale
  count2 <- round(1 + (count2 - post.syn.min) * (c_size_single - 1) / post.syn.rng)

  # Plot the neuron
  plot3d(neu,
         soma=TRUE,
         lwd=3,
         col=col_single[count2])
}

# Plot joined dataset with gradient color post1 (red) -> post2 (blue)
# Join all the bodyIDs
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

bodyIDs = union(bodyIDs.post1, bodyIDs.post2)

for (bodyID in bodyIDs) {
  neu <- nlist[[toString(bodyID)]]

  # Get the synapses count and correct for NULL values (only connection to one
  # of the pre)
  count1 <- neu$post1count
  if (is_null(count1)) {
    count1 <- post.syn.min
  }

  count2 <- neu$post2count
  if (is_null(count2)) {
    count2 <- post.syn.min
  }

  # Use differences between number of synapses and scale to the size of the
  # color scale
  count <- 0.5 * (1 + (count1 - count2)/(count1 + count2))
  count <- round(1 + count * (c_size_both - 1))

  # Plot the neuron (post1: red, post2: blue)
  plot3d(neu,
         soma=TRUE,
         lwd=3,
         col=col_both[count])
}
