#
# Plots neurons' skeletons colored according to the number of synapses
# Evaluates and displays the glomerulus based on the location of the
# synapses.
# Find intersections with a plane.
#

# Import required libraries
library(tidyverse)
library(nat)
library(zeallot)
library(alphashape3d)

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

# Define geometries ----
# Identify a plane separating lobula from LC4 glomerulus
c(a, b, c, d) %<-% list(1, 0, -0.60, 7500)

# Extract all LC4 synapses
lc4.all <- con %>% filter(pre.type == 'LC4')

# Restrict to synapses on the LC4 glomerulus side
lc4.glo <- lc4.all %>%
           filter(lc4.all %>%
                  pull(post.coors) %>%
                  xyzmatrix() %>%
                  as_tibble() %>%
                  mutate(GLO = a*X + b*Y + c*Z + d > 0) %>%
                  select(GLO))

# Create a concave mesh including all the synapses to identify the LC4 glomerulus
glo.mesh <- lc4.glo %>%
            pull(post.coors) %>%
            xyzmatrix() %>%
            ashape3d(., alpha = 1000) %>%
            as.mesh3d()

# Find DNP ----
# LC4 -> DNp02
# Synapses
lc4.dnp02 <- lc4.glo %>% filter(post.type == 'DNp02')
lc4.dnp02.sort <- lc4.dnp02 %>%
                  group_by(pre.bodyID, pre.type, post.type) %>%
                  count() %>%
                  arrange(desc(n))
# Skeletons
lc4.dnp02.skel <- nlist[lc4.dnp02.sort %>%
                        pull(pre.bodyID) %>%
                        as.character()]
# Color palette
lc4.dnp02.pal <- colorRampPalette(c(rgb(0, 0, 1), rgb(1, 1, 1)))(length(lc4.dnp02.skel))

# LC4 -> DNp11
# Synapses
lc4.dnp11 <- lc4.glo %>% filter(post.type == 'DNp11')
lc4.dnp11.sort <- lc4.dnp11 %>%
                  group_by(pre.bodyID, pre.type, post.type) %>%
                  count() %>%
                  arrange(desc(n))
# Skeletons
lc4.dnp11.skel <- nlist[lc4.dnp11.sort %>%
                          pull(pre.bodyID) %>%
                          as.character()]
# Color palette
lc4.dnp11.pal <- colorRampPalette(c(rgb(1, 0, 0), rgb(1, 1, 1)))(length(lc4.dnp11.skel))

# Plot results ----
# Plot the neurons
nopen3d()
par3d('windowRect' = c(100,100,1000,1000))
#plot3d(lc4.dnp02.skel, col=lc4.dnp02.pal, soma=100, add=TRUE)
#plot3d(lc4.dnp11.skel, col=lc4.dnp11.pal, soma=100, add=TRUE)
lc4.dn.skel <- union(lc4.dnp02.skel, lc4.dnp11.skel)
plot3d(lc4.dn.skel, col='lightgray', soma=100, ADD=TRUE) # add alpha=0.25 as parameter to see the transparent version

# Plot the LC4 glomerulus 3D mesh
plot3d(glo.mesh, color='lightblue', alpha=0.3, add=TRUE)

# Draw a plane separating lobula from LC4 glomerulus
planes3d(a, b, c, d=d, alpha = 0.3, add=TRUE)

# Highlight the interection between neurons and plane based on the number of synapses
# TODO: Do it using the gradent color based on the nummber of synapses
for (i in 1:length(lc4.dnp02.skel)) {
  neu <- lc4.dnp02.skel[[i]]
  int <- intersect_plane(neu, c(a, b, c, d))
  plot3d(int, col='blue', add=TRUE)
}
for (i in 1:length(lc4.dnp11.skel)) {
  neu <- lc4.dnp11.skel[[i]]
  int <- intersect_plane(neu, c(a, b, c, d))
  plot3d(int, col='red', add=TRUE)
}


# Plot the dnp02 post synapses
lc4.dnp02.coors <- lc4.dnp02 %>% pull(post.coors) %>% xyzmatrix()
plot3d(lc4.dnp02.coors, col='blue', add=TRUE)

# Plot the dnp11 post synapses
lc4.dnp11.coors <- lc4.dnp11 %>% pull(post.coors) %>% xyzmatrix()
plot3d(lc4.dnp11.coors, col='red', add=TRUE)

# Show axes and box to figure out the geometry
axes3d()
title3d(xlab='x', ylab='y', zlab='z')

