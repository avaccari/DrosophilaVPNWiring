# Created on Thu Jul 29 16:51:37 2021
# Name:    synaptic_connectivity.R
# Purpose: Show the gradient in the neuron connectivity for the top 20 posts of
#          a given pre.
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - A plot showing the graded pre > posts synaptic connectivity for the top 
#   posts of a specific pre
#
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

# Import required libraries
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
  con <- readRDS("hemibrain_con0.rds")
}
if (!exists('nlist')) {
  nlist <- readRDS("nlist1.rds")
}









###############################################################################
# Define items to analyze here
pre_type <- 'LC4'

# Top (# of synapses) of post to consider
top <- 12
###############################################################################






# Extract all pre.type synapses excluding with pre
pre <- con %>% 
  filter(pre.type==pre_type, post.type != pre_type)

# Extract the top post in terms of synapses
top_posts <- pre %>%
  group_by(post.type) %>%
  dplyr::count() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  slice(1:top) %>%
  pull(post.type)

# Filter the data set according to the top synapses
top_cnt <- pre %>%
  filter(post.type %in% top_posts) %>%
  group_by(pre.bodyID, post.type) %>%
  dplyr::count() %>%
  arrange(desc(n))

# Add index as column
top_cnt$idx <- as.numeric(rownames(top_cnt))

# Assing a gradient color scale to the pre.bodyIDs of the post with most
# connections and use the same color for each pre.bodyID in the other posts'
# gradients representations.
top_cnt <- top_cnt %>% ungroup()

# Extract the pre.bodyIDs for the post with most connections
top_post_ids <- top_cnt[top_cnt$post.type==top_cnt[[1, 2]], ]

# Add a unique color identification to each pre.bodyID
top_post_ids <- cbind(top_post_ids, preBodyId=1:nrow(top_post_ids))

# Select only the required coloumns for the final merge
top_post_ids <- top_post_ids %>% select(pre.bodyID, preBodyId)

# Merge with top_cnt
top_cnt <- merge(top_cnt, top_post_ids)

# Plot data
ggplot() +
  xlab("Individual LC4, arranged by descending output (independent on every graph)")+
  ylab("Number of synapses")+
  theme_bw()+
  ylim(0, NA) +
  geom_point(data=top_cnt,
             aes(x=idx, y=n, col=preBodyId), 
             size=3) +
  scale_color_gradientn(colours=rainbow(71)) +
  facet_wrap(~ post.type,
             scales="free") +
  #used to have "free_x
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.text.y = element_text(face="bold", color="black", size=12, angle=0),
        plot.title = element_text(face="bold", color="blue", size=15),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        strip.text = element_text(size = 13, face="bold"),
        strip.background = element_rect(fill="#33CCFF"))
