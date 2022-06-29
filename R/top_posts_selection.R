# Created on Mon Aug 30 21:10:38 2021
# Name:    top_posts_selection.R
# Purpose: Evaluate how many posts are needed to provide a user selectable
#          coverage of the pre synapses
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - A scree plot with the percentage of synapses covered by each post in
#   descending order
# - A cumulative plot of the synapses covered by the posts (ordered in
#   descending order)
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

# Clean everything up ----
# Except the connection and skeleton files if they are already loaded.
items <- ls()
items <- items[items != "con"] # Comment this line to reload the connection file
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
if (!exists("con")) {
  con <- readRDS("data/hemibrain_con0.rds")
}
if (!exists("nlist")) {
  nlist <- readRDS("data/nlist1.rds")
}










###############################################################################
# Define items to analyze here
pre_type <- "LC4"

# top % of synapses required
top <- 0.82
###############################################################################







# Extract all pre.type synapses excluding with pre
pre <- con %>%
  filter(pre.type == pre_type, post.type != pre_type)

# Extract the top post in terms of synapses
top_posts <- pre %>%
  group_by(post.type) %>%
  dplyr::count() %>%
  ungroup() %>%
  arrange(desc(n))

# Normalize by the total number of synapses
top_posts$norm <- top_posts$n / sum(top_posts$n)

# Add cumulative sum
top_posts$cml <- cumsum(top_posts$norm)

# Find elements with sum below top
top_posts_select <- top_posts[top_posts$cml <= top, ]

# Create scree plot
ggplot() +
  geom_point(
    data = top_posts,
    aes(x = post.type, y = norm)
  ) +
  scale_x_discrete(limits = top_posts$post.type) +
  geom_line() +
  labs(title = paste("Scree plot: synapses percentage vs posts (", pre_type, ")")) +
  xlab("Post") +
  ylab("Percentage of synapses") +
  theme(axis.text.x = element_text(angle = 90))

# Create cumulative plot
ggplot() +
  geom_point(
    data = top_posts,
    aes(x = post.type, y = cml)
  ) +
  scale_x_discrete(limits = top_posts$post.type) +
  geom_line() +
  labs(title = paste("Cumulative plot: synapses percentage vs posts (", pre_type, ")")) +
  xlab("Post") +
  ylab("Percentage of synapses") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = top) +
  geom_point(
    data = top_posts_select,
    aes(x = post.type, y = cml),
    shape = 1,
    size = 5,
    col = "red"
  ) +
  geom_text(
    data = top_posts_select,
    aes(
      x = post.type,
      y = cml,
      label = post.type,
      hjust = "left"
    ),
    nudge_x = 5
  )
