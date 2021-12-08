# Created on Thu Jul 08 18:29:37 2021
# Name:    gradients_correlation_pearson.R
# Purpose: Plots the number of synapses of each pre with both post against each
#          other
# Author:  Mark Dombrovski (md2ff@virginia.edu)
#
# Generates one scatter plot where the coordinates of each point represents the
# number of synapses each pre has with the two posts (one for each coordinate).
# The fit and the Pearson's correlation coefficient are shown on the plot.
#
# Copyright (c) 2021 Mark Dombrovski
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
library(tidyverse)
library(ggpubr)

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
post_type1 <- 'DNp02' # y-axis
post_type2 <- 'DNp11' # x-axis
###############################################################################





# Filter for pre and post1 and count synapses of post1 with pre
tmp0 <- con %>%
  filter(pre.type == pre_type, post.type == post_type1) %>%
  group_by(pre.bodyID, post.type) %>%
  dplyr::count()

# Filter for pre and post2 and count synapses of post2 with pre
tmp1 <- con %>%
  filter(pre.type == pre_type, post.type == post_type2) %>%
  group_by(pre.bodyID, post.type) %>%
  dplyr::count()

# Combine the datasets to plot the number of synapses against each other
merge <- merge(tmp0, tmp1, by="pre.bodyID", all=TRUE)
merge[is.na(merge)] <- 0
merge <- merge

# Calculate the Pearson's correlation coefficient
pc <- cor(merge$n.x, merge$n.y, method = "pearson", use="complete.obs")

# Normalize the counts with respect to the max between the two
count_max_x <- max(merge$n.x)
count_max_y <- max(merge$n.y)
merge$n.x.norm <- merge$n.x / count_max_x
merge$n.y.norm <- merge$n.y / count_max_y


# Generate a scatter plot. Each point is a pre and the coordinates are the
# number of synapsed with each of the posts
ggplot(merge, aes(x=n.x.norm, y=n.y.norm, col=pre.bodyID))+
  geom_point(size = 3, col="steelblue")+
  geom_smooth(method="lm",formula= 'y ~ x', col="red", se=TRUE)+
  theme_classic() +
  ylab(paste(pre_type, '>', post_type1, 'synapses'))+
  xlab(paste(pre_type, '>', post_type2, 'synapses'))+
  theme(axis.text.x = element_text(face="bold", color="black", size=13, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=15, angle=0),
        plot.title = element_text(face="bold", color="blue", size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        aspect.ratio=1)+
  annotate('text', size=5, col="red",
           label=toString(paste("Pearson's r =", round(pc, 2))),
           x=0.7 * max(merge$n.x.norm),
           y=0.2 * max(merge$n.y.norm),
           size=5)+
  xlim(0, 1)+
  ylim(0, 1)
