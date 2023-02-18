# Created on Fri Apr 15 11:34:29 2022
# Name:    connectivity_bias.R
# Purpose: Evaluate the average coefficient of variation for the post synaptic
#          partners of each VPN
# Author:  Mark Dombrovski (md2ff@virginia.edu)
#
# For each of the VPN:
# - evaluates the top postsynaptic partners (number selected by user)
# - for each detected postsynaptic partner, calculates the coefficient of
#   variation of the connectivity with the VPN
# - returns the averagee of the coefficienta of variation for all the post
#   synaptic partners

#
# Copyright (c) 2022 Mark Dombrovski
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


# Import required libraries
library(tidyverse)

# Clean everything up ----
# Except the connection and skeleton files if they are already loaded.
items <- ls()
items <- items[items != "con"] # Comment this line to reload the connection file
rm(list = items)

# Close any open plotting windows
while (dev.cur() > 1) {
  dev.off()
}
while (rgl.cur() > 0) {
  rgl.close()
}

# Load datasets
if (!exists("con")) {
  con <- readRDS("data/hemibrain_con0.rds")
}




###############################################################################
# Define items to analyze here
# VPNs to consider
VPNs <- c(
  "LC4", "LC6", "LC9", "LC10", "LC11", "LC12", "LC13",
  "LC15", "LC16", "LC17", "LC18", "LC20", "LC21", "LC22",
  "LC24", "LC25", "LC26", "LPLC1", "LPLC2", "LPLC4"
)
# VPNs <- c('LC4')

# Number of top postsynaptic partner to conside
top <- 25
###############################################################################



## Functions
# Returns the top n postsynaptic partners of a given VPN
top_n <- function(vpn, top) {
  top_post <- con %>%
    filter(pre.type == vpn, post.type != vpn) %>%
    group_by(pre.type, post.type) %>%
    dplyr::count() %>%
    ungroup() %>%
    arrange(desc(n)) %>%
    slice(1:top) %>%
    pull(post.type)

  return(top_post)
}

# Return the number of connection between vpn and post
cnt <- function(post, vpn) {
  conn <- con %>%
    filter(pre.type == vpn, post.type == post) %>%
    group_by(pre.bodyID) %>%
    dplyr::count() %>%
    ungroup() %>%
    pull(n)

  return(conn)
}

## Evaluation
# Find top n postynaptic partners for all VPNs
top_posts <- as.data.frame(sapply(VPNs, top_n, top = top))

# Evaluate the connections for all the postsynaptic partners of a given VPN
for (i in colnames(top_posts)) {
  # Extract connectivity for each post
  syn <- sapply(top_posts[, i], cnt, vpn = i)
  # Evaluate standard deviation
  stdev <- sapply(syn, sd)
  # Evaluate mean
  ave <- sapply(syn, mean)
  # Evaluate coefficient of variation
  coefvar <- 100 * stdev / ave
  # Show results
  cat(i, "average CV:", mean(coefvar), "\n")
}
