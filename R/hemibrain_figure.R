# Created on Thu Oct  30 18:18:11 2021
# Name:    hemibrain_figure.R
# Purpose: Create an image of the hemibrain and analyzed neurons
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
#  - A 3D plot showing the hemibrain outline with LC4, giant fiber, 
#    dnp04, and dnp11
#  - A 3D plot showing the hemibrain outline with LC4, LPLC2, and giant
#    fiber
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
r.
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
pre_type1 <- 'LC20'
pre_type2 <- 'LC28b'
post_type11 <- 'DNp02'  
post_type12 <- 'DNp11'  
post_type13 <- 'DNp04'
post_type14 <- 'Giant Fiber'
post_type21 <- 'Giant Fiber'  

# Plot window size
win_siz <- 1000
###############################################################################






# Pre ----
# Extract all pre.type synapses
pre1 <- con %>% filter(pre.type==pre_type1)
pre2 <- con %>% filter(pre.type==pre_type2)

# Posts ----
# Identify pre with synapses on posts and pull bodyIDs
bodyIDs.pre1 <- pre1 %>%
#  filter(post.type==post_type11 | post.type==post_type12 | post.type==post_type13 | post.type==post_type14) %>%
  pull(pre.bodyID) %>%
  unique()

bodyIDs.pre2 <- pre2 %>%
#  filter(post.type==post_type21) %>%
  pull(pre.bodyID) %>%
  unique()

# # Pull posts bodyIDs
# post11.ID <- pre1 %>%
#   filter(post.type==post_type11) %>%
#   pull(post.bodyID) %>%
#   unique()
# 
# post12.ID <- pre1 %>%
#   filter(post.type==post_type12) %>%
#   pull(post.bodyID) %>%
#   unique()
# 
# post13.ID <- pre1 %>%
#   filter(post.type==post_type13) %>%
#   pull(post.bodyID) %>%
#   unique()
# 
# post21.ID <- pre2 %>%
#   filter(post.type==post_type21) %>%
#   pull(post.bodyID) %>%
#   unique()


# Plot first hemibrain
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

plot3d(hemibrain.surf, alpha=0.1, col='lightblue')

# Plot pre
# Exclude a specific bodyID that we want to highlight
# bodyIDs.pre1 = bodyIDs.pre1[bodyIDs.pre1 != 1809264255]

plot3d(nlist[as.character(bodyIDs.pre1)], soma=TRUE, col='steelblue', add=TRUE)
# plot3d(nlist[as.character(bodyIDs.pre1)], soma=TRUE, col='orange', add=TRUE)

# Add the highlighted bodyID
# plot3d(nlist[['1809264255']], soma=TRUE, col='red', lwd=4.0, add=TRUE)

# Plot posts
plot3d(giantF, lwd=2.0, col='black', add=TRUE)
plot3d(dnp11, lwd=2.0, col='chartreuse3', add=TRUE)
plot3d(dnp04, lwd=2.0, col='firebrick3', add=TRUE)




# Plot second emibrain
nopen3d()
par3d('windowRect' = c(100, 100, win_siz, win_siz))

plot3d(hemibrain.surf, alpha=0.1, col='lightblue')

# Plot pres
plot3d(nlist[as.character(bodyIDs.pre1)], soma=TRUE, col='steelblue', add=TRUE)
plot3d(nlist[as.character(bodyIDs.pre2)], soma=TRUE, col='orange', add=TRUE)

# Plot posts
#plot3d(giantF, col='black', add=TRUE)



