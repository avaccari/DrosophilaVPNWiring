# Created on Tue Jan 04 16:18:32 2022
# Name:    clusters_plots
# Purpose: Generate plots of the clusters obtained from the k-means.
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Based on spatial_cluster_plots.Rmd.
#
# Generates:
# - One 3D plot for each of the VPN identifying the detected clusters in
#   different colors.
#
# Copyright (c) 2022 Andrea Vaccari
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

# Clean everything up ----
# Except the connection and skeleton files if they are already loaded.
items <- ls()
# items <- items[items != 'con']  # Comment this line to reload the connection file
items <- items[items != "nlist"] # Comment this line to reload the skeleton file
rm(list = items)

# Close any open plotting windows
while (dev.cur() > 1) {
  dev.off()
}
while (rgl.cur() > 0) {
  rgl.close()
}

# Load datasets
# if (!exists('con')) {
#   con <- readRDS("data/hemibrain_con0.rds")
# }
if (!exists("nlist")) {
  nlist <- readRDS("data/nlist1.rds")
}









###############################################################################
# Define items to analyze here
# List of VPNs of interest
all_vpns <- c(
  "LC17", "LC12", "LC4", "LC6", "LPLC1", "LPLC2",
  "LC9", "LC10", "LC11", "LC13", "LC15", "LC16",
  "LC18", "LC20", "LC21", "LC22", "LC24", "LC25",
  "LC26", "LPLC4"
)

# Use fix coloring based on max number of clusters (if 0, automatical generate
# color scale)
max_clusters <- 3

# Plot window size
win_siz <- 500
###############################################################################







# Load all synapses
vpn_syn <- readRDS("output/LC_kmeans/allLC_synapses.rds")

# Generate a plot for each VPN of interest
for (vpn in all_vpns) {
  # Extract data from list
  vpn_dat <- vpn_syn[[vpn]]

  # Plot
  nopen3d()
  par3d("windowRect" = c(100, 100, win_siz, win_siz))

  # Evaluate number of clusters
  vpn_clust <- unique(vpn_dat$cluster)

  # Check if we have fix colors:
  if (max_clusters != 0) {
    # Prepare the color scale based on the number of clusters
    col_vpn <- cet_pal(max_clusters)
  } else {
    col_vpn <- cet_pal(length(vpn_clust))
  }

  # Extract and plot members of each class
  for (cl in vpn_clust) {
    # Find unique bodyids for each class
    vpn_class <- vpn_dat %>%
      filter(cluster == cl) %>%
      select("bodyid") %>%
      unique()

    # Plot the neuron in each of the class individually (we need to do this
    # beacuse we will get an error if we use nlist[vpn_class$bodyid] and
    # some of the neurons are not in nlist)
    for (bid in vpn_class$bodyid) {
      n <- nlist[[bid]]
      if (!is.null(n)) {
        plot3d(n, lw = 3, col = col_vpn[as.numeric(cl)], add = TRUE)
      }
    }
  }

  # Add boxes for LC4
  if (vpn == "LC4") {
    box3d()
  }
}
