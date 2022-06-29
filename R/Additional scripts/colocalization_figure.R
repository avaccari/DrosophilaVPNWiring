# Created on Thu Nov  04 21:55:37 2021
# Name:    colocalization_figure.R
# Purpose: Create images illustrating the colocalization of selected posts
#          in the glomerulus of a given pre
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
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
library(natverse)
library(tidyverse)

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

# Load additional sekeletons
load("data/hemibrain_figure.rda")




###############################################################################
# Define items to analyze here
pre_type <- 'LC4'
post_type1 <- 'DNp02'    
post_type2 <- 'DNp11' 
post_type3 <- 'DNp04'  

pre_type_color <- 'darkgreen'
post_type1_color <- 'red'
post_type2_color <- 'purple'
post_type3_color <- 'darkgoldenrod1'

pre_type_alpha <- 1.0

# Plot window size
win_siz <- 500
###############################################################################






# Pre ----
# Extract all pre.type synapses
pre <- con %>% filter(pre.type==pre_type)

# Posts ----
# Identify pre with synapses on posts and pull bodyIDs
bodyIDs.pre <- pre %>%
  filter(post.type==post_type1 | post.type==post_type2 | post.type==post_type3) %>%
  pull(pre.bodyID) %>%
  unique()



# Plot first post
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Plot pre
plot3d(nlist[as.character(bodyIDs.pre)], soma=TRUE, alpha=pre_type_alpha, col=pre_type_color, add=TRUE)

# Plot posts
plot3d(dnp02, lwd=2.0, col=post_type1_color, add=TRUE)



# Plot second post
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Plot pre
plot3d(nlist[as.character(bodyIDs.pre)], soma=TRUE, alpha=pre_type_alpha, col=pre_type_color, add=TRUE)

# Plot posts
plot3d(dnp11, lwd=2.0, col=post_type2_color, add=TRUE)



# Plot third post
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

# Plot pre
plot3d(nlist[as.character(bodyIDs.pre)], soma=TRUE, alpha=pre_type_alpha, col=pre_type_color, add=TRUE)

# Plot posts
plot3d(dnp04, lwd=2.0, col=post_type3_color, add=TRUE)

