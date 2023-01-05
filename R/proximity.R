# Created on Thu Jul 07 18:30:23 2022
# Name:    proximity.R
# Purpose: Evaluates the "overlap" between pre and post
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Based on user selection, evaluated the distance between a VPN and a post.
# It doenloads the required skeletons and saves the results.
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
library(nat)
library(neuprintr)







###############################################################################
# Define items to analyze here
pre_type <- "LC17"
post_type <- "PVLP071"

# Distance in voxel units (8nm) to use as threshold
delta <- 125

# Export to CSV
export_to_csv <- TRUE

# Define login parameters
neuprint_server <- "https://neuprint.janelia.org"
neuprint_dataset <- "hemibrain:v1.2.1"

# Use your personal token here
# neuprint_token <- <Your token here>
###############################################################################








# Connect to neuprint
conn <- neuprint_login(
    server = neuprint_server,
    token = neuprint_token,
    dataset = neuprint_dataset
)


# Grab the info for the pre and post
pre <- neuprint_search(pre_type, field = "type")
post <- neuprint_search(post_type, field = "type")

# Get the body IDs
pre_id <- pre$bodyid
post_id <- post$bodyid

# Check if available, if not download the neurons
if (file.exists(sprintf("data/neurons/%s_neu.rds", pre_type))) {
    pre_neu <- readRDS(sprintf("data/neurons/%s_neu.rds", pre_type))
} else {
    pre_neu <- neuprint_read_neurons(pre_id)
    saveRDS(pre_neu, sprintf("data/neurons/%s_neu.rds", pre_type))
}
if (file.exists(sprintf("data/neurons/%s_neu.rds", post_type))) {
    post_neu <- readRDS(sprintf("data/neurons/%s_neu.rds", post_type))
} else {
    post_neu <- neuprint_read_neurons(post_id)
    saveRDS(post_neu, sprintf("data/neurons/%s_neu.rds", post_type))
}

# Check if the score for this combination is already evaluated.
# If not, evaluate and store.
if (file.exists(
    sprintf(
        "output/overlaps/%s_%s_del%d.rds",
        pre_type,
        post_type,
        delta
    )
)) {
    os <- readRDS(
        sprintf(
            "output/overlaps/%s_%s_del%d.rds",
            pre_type,
            post_type,
            delta
        )
    )
} else {
    # Evaluate the overlap between pre and post
    os <- overlap_score(pre_neu, post_neu, delta = delta)
    saveRDS(
        os,
        sprintf(
            "output/overlaps/%s_%s_del%d.rds",
            pre_type,
            post_type,
            delta
        )
    )
}

# Check if we are exporting to CSV
if (export_to_csv == TRUE) {
    write.csv(
        os,
        file = sprintf(
            "output/overlaps/%s_%s_del%d.csv",
            pre_type,
            post_type,
            delta
        )
    )
}
