# Created on Thu Oct  7 13:14:09 2021
# Name:    gather_data.R
# Purpose: Download various datasets from the neuprint database
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
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
library(neuprintr)

# Clean everything up ----
# Except the connection and skeleton files if they are already loaded.
items <- ls()
# items <- items[items != 'con']  # Comment this line to reload the connection file
# items <- items[items != 'nlist']  # Comment this line to reload the skeleton file
rm(list = items)

# Define login parameters
neuprint_server <- "https://neuprint.janelia.org"
neuprint_dataset <- "hemibrain:v1.2.1"

# Use your personal token here
# neuprint_token <- <Your token here>

# Open connection
conn <- neuprint_login(server = neuprint_server, token = neuprint_token, dataset = neuprint_dataset)

# Get all possible bodyIDs
all_ids <- neuprint_search(".*", meta = FALSE)

# Get all synapses for all possible bodyIDs
# all_syn <- neuprint_get_synapses(all_ids)

# Get all skeletons for all possible bodyIDs
# nlist <- neuprint_read_neurons(all_ids)
# saveRDS(nlist, file='data/new_nlist.rds')

# Get few neurons of interest
dnp02 <- neuprint_read_neurons(neuprint_search("DNp02.*"), meta = FALSE)
dnp04 <- neuprint_read_neurons(neuprint_search("DNp04.*"), meta = FALSE)
dnp11 <- neuprint_read_neurons(neuprint_search("DNp11.*"), meta = FALSE)
giantF <- neuprint_read_neurons(neuprint_search("Giant.*"), meta = FALSE)

# Save them in a file
save(dnp02, dnp04, dnp11, giantF, file = "data/hemibrain_figure.rda")
