# Created on Mon Jun 21 19:05:27 2021
# Name:    top_synapses_highlight_3d.R
# Purpose: Highlights post neurons with highest number of synapses with pre
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates three plots:
# - one plot for each of the posts highlighting the top_neu neurons based on the
#   number of synapses. The other post neurons are plotted in gray.
# - one plot as above where the top_neu for ech post are highlighted with
#   diferent colors. The other post neurons are plotted in gray.

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
# - If neuron are within top_neu in both list they will only be displayed with
#   the first color they are plotted with (does it matter?)
# - Transparency is very slow. Change to use (or give the option to use) plotly

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
  con <- readRDS("data/hemibrain_con0.rds")
}
if (!exists('nlist')) {
  nlist <- readRDS("data/nlist1.rds")
}






###############################################################################
# Define items to analyze here
pre_type <- 'LC4'
post_type1 <- 'DNp02' # Blue
post_type2 <- 'DNp11' # Red

# Using perceptually uniform color map
# More options here:
# https://cran.r-project.org/web/packages/cetcolor/vignettes/cet_color_schemes.html
c_map <- 'd11'  # Map to use for individual post plot
c_size <- 100  # Number of colors in the map

# Define the number of top neurons from each post to highlight.
top_neu <- 1

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
         dplyr::count() %>%
         arrange(desc(n))

# Post2 ----
# Identify pre with synapses on post1 and count
post2 <- pre %>%
         filter(post.type==post_type2) %>%
         group_by(pre.bodyID) %>%
         dplyr::count() %>%
         arrange(desc(n))




# Display ----
# Creating colormap
cmap <- cet_pal(c_size, name=c_map)

# Highlight the individual post top top_neu neurons (synapses with pre) ----
# Plot top top_neu for post1 (red) with others in gray
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

for (pos in 1:nrow(post1)) {
  neu <- nlist[[toString(post1[pos, 1])]]

  # Highlight the top_neu neurons, the rest are plotted in gray
  col <- cmap[50]
  alpha <- 0.2
  if (pos <= top_neu) {
    col <- cmap[100]
    alpha <- 1.0
  }

  # Plot the neuron
  plot3d(neu,
         soma=TRUE,
         lwd=2,
         # alpha=alpha,
         col=col)
}

# Plot top top_neu for post2 (blue) with others in gray
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

for (pos in 1:nrow(post2)) {
  neu <- nlist[[toString(post2[pos, 1])]]

  # Highlight the top_neu neurons, the rest are plotted in gray
  col <- cmap[50]
  alpha <- 0.2
  if (pos <= top_neu) {
    col <- cmap[1]
    alpha <- 1.0
  }

  # Plot the neuron
  plot3d(neu,
         soma=TRUE,
         lwd=2,
         # alpha=alpha,
         col=col)
}




# Highlight the top top_neu neurons (synapses with pre) for both posts ----
# Plot joined dataset with gradient color post1 (red) -> post2 (blue)
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Plot the top top_neu from each post in red (post1) and blue (post2)
for (pos in 1:top_neu) {
  neu1 <- nlist[[toString(post1[pos, 1])]]
  neu2 <- nlist[[toString(post2[pos, 1])]]

  # Plot the neurons
  plot3d(neu1,
         soma=TRUE,
         lwd=2,
         # alpha=1.0,
         col=cmap[100])

  plot3d(neu2,
         soma=TRUE,
         lwd=2,
         # alpha=1.0,
         col=cmap[1])
}

# Plot the rest in gray
max_row = max(nrow(post1), nrow(post2))

for (pos in (top_neu + 1):max_row) {
  if (pos <= nrow(post1)) {
    neu <- nlist[[toString(post1[pos, 1])]]
    plot3d(neu,
           soma=TRUE,
           lwd=2,
           # alpha=0.2,
           col=col)
  }

  if (pos <= nrow(post2)) {
    neu <- nlist[[toString(post2[pos, 1])]]
    plot3d(neu,
           soma=TRUE,
           lwd=2,
           # alpha=0.2,
           col=col)
  }
}
