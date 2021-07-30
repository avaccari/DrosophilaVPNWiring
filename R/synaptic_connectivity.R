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
# - 4x5 independent plots of the top 20 posts for a given pre sorted by number
#   of synapses
# con %>%
#   filter(pre.type %in% c("LC4"), post.type %in% c("SAD064","DNp11","PLP219")) %>%
#   group_by(pre.bodyID, pre.type, post.type) %>%
#   dplyr::count() %>%
#   ungroup()%>%
#   ggplot(aes(x=reorder(pre.bodyID,n), y = n))+
#   geom_point(size = 6, col="steelblue")+
#   facet_wrap(. ~ post.type, scales = "free_y")+
#   theme_bw() +
#   xlab("VPN ID, arranged by descending output")+
#   ylab("number of synapses")+
#   theme(
#     plot.title = element_text(color="black", size=14, face="bold.italic"),
#     axis.title.x = element_text(color="black", size=14, face="bold"),
#     axis.text.x = element_text(face="bold", color="black", size=10),
#     axis.title.y = element_text(color="black", size=14, face="bold"),
#     axis.text.y = element_text(face="bold", color="black", size=15),
#     strip.text = element_text(size = 20, face="bold")
#   )

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
  con <- readRDS("data/hemibrain_con0.rds")
}
if (!exists('nlist')) {
  nlist <- readRDS("data/nlist1.rds")
}









###############################################################################
# Define items to analyze here
pre_type <- 'LC4'

# Top (# of synapses) of post to consider
top <- 20
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

# Filter the dataset accroding to the top synapses
top_cnt <- pre %>%
  filter(post.type %in% top_posts) %>%
  group_by(pre.bodyID, post.type) %>%
  count() %>%
  arrange(desc(n))

# Add index as column
top_cnt$idx <- as.numeric(rownames(top_cnt))

# Plot data
ggplot() +
  xlab("VPN ID, arranged by descending output")+
  ylab("number of synapses")+
  geom_point(data=top_cnt,
             aes(x=idx, y=n)) +
  facet_wrap(~ post.type,
             scales='free_x') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

            