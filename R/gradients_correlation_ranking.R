# Created on Thu Jul 13 17:50:17 2021
# Name:    gradients_correlation_ranking.R
# Purpose: Selects the top n (number of synapses) post for a given pre and
#          evaluates the Pearson's correlation coefficient between all possible
#          pairs.
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - one correlation plot (heat map) showing Person's cross correlation
#   coefficients between between each of the to 25 (in terms of synapses)
#   posts for a specific pre.
# - one correlation plot (heat map) where the non-significant values 
#   (p < 0.05) are blanked out
# - one heat map showing the p-values
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
library(tidyverse)
library(corrplot)
library(Hmisc)

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
top <- 25 
###############################################################################





# Extract the top post in terms of synapses
top_posts <- con %>%
             filter(pre.type == pre_type, post.type != pre_type) %>%
             group_by(pre.type, post.type) %>%
             dplyr::count() %>%
             ungroup() %>%
             arrange(desc(n)) %>%
             slice(1:top) %>%
             pull(post.type)

# Select the pre
pre <- con %>%
       filter(pre.type == pre_type)

# Join the post in a data.frame
for (i in 1:length(top_posts)) {
  cnt <- pre %>%
         filter(post.type == top_posts[[i]]) %>%
         group_by(pre.bodyID) %>%
         dplyr::count()
  if (i == 1) {
    posts <- cnt
  } else {
    posts <- merge(posts, cnt, by="pre.bodyID", all=TRUE)
  }
}

# Drop the pre.bodyID column
posts <- posts[-1]

# Change the names of columns to the top posts
colnames(posts) <- c(top_posts)

# Change NaN to zeros
posts[is.na(posts)] <- 0

# Evaluate the Spearman's correlation matrix
corr <- rcorr(as.matrix(posts), type='spearman')
pcorr <- corr[[1]]
pcorr_p <- corr[[3]]
diag(pcorr_p) <- 0

# Create a data frame with the results
ut <- upper.tri(pcorr)
pcorr_df_o <- data.frame(row=rownames(pcorr)[row(pcorr)[ut]],
                         column=rownames(pcorr)[col(pcorr)[ut]],
                         cor=(pcorr)[ut])
pcorr_df <- pcorr_df_o[order(pcorr_df_o$cor), ]
pcorr_p_df_o <- data.frame(row=rownames(pcorr_p)[row(pcorr_p)[ut]],
                          column=rownames(pcorr_p)[col(pcorr_p)[ut]],
                          cor=(pcorr_p)[ut])
pcorr_p_df <- pcorr_p_df_o[order(pcorr_df_o$cor), ]

cat('Most anti-correlated:\n')
print(head(pcorr_df, 10))
cat('Most correlated:\n')
print(tail(pcorr_df, 10))

###############################################################################
# This section is needed only if you want the correlation plot to be sorted in
# the same way in all the scripts.
sort.pcorr <- order(rownames(pcorr))
pcorr_s <- pcorr[sort.pcorr, sort.pcorr]
pcorr.HC <- corrMatOrder(pcorr_s, order='hclust')
pcorr_c <- pcorr_s[pcorr.HC, pcorr.HC]
sort.pcorr_p <- order(rownames(pcorr_p))
pcorr_p_s <- pcorr_p[sort.pcorr_p, sort.pcorr_p]
pcorr_p_c <- pcorr_p_s[pcorr.HC, pcorr.HC]
###############################################################################

# Show graphic results
# (There are a lot of options for this plot)
corrplot(pcorr_c, 
         method='color',  # Color the background
         col=colorRampPalette(c("darkred", "white", "darkblue"))(200),
         order='original',  # First principal component order
         # addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         tl.col='black')  # Color of the labels in black

corrplot(pcorr_c, 
         p.mat=pcorr_p_c,  # p-values matrix
         sig.level=0.05,  # significance level
         insig='blank',  # Blank out insignificant values
         method='color',  # Color the background
         col=colorRampPalette(c("darkred", "white", "darkblue"))(200),
         order='original',  # First principal component order
         # addCoef.col='black',  # Add values in black
         number.cex=0.6,  # values size
         tl.col='black')  # Color of the labels in black

corrplot(pcorr_p_c,
         is.corr=FALSE,  # Not a correlation matrix
         method='color',  # Color the background
         col=colorRampPalette(c("darkblue", "white"))(200),
         order='original',  # First principal component order
         tl.col='black')  # Color of the labels in black
