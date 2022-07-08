# Created on Wed Jun 29 18:33:41 2022
# Name:    clustering_spatial_bias.R
# Purpose: Evaluated the quality of the VPN clusters created in the cluster
#          analysis step
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - Dendritic maps showing the locations of the centroids of each of the
#   selected VPNs' neurons colored based on their cluster. The VPN can be
#   specified by the user.
# - Evaluates and displays the quality index of each of the clustering. The
#   idexes can be specified by the user
#
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
# - Add coloring based on clusters in the 3D analysis


# Import required libraries
library(clusterCrit)
library(natverse)
library(cetcolor)
library(zeallot)
library(pracma)
library(ks)
library(matrixStats)
library(egg)
library(ggpmisc)


# Clean everything up ----
# Except the connection and skeleton files if they are already loaded.
items <- ls()
items <- items[items != "con"] # Comment to reload the connection file
items <- items[items != "nlist"] # Comment to reload the skeleton file
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

# Source local files
source("R/aux_functions.R")

# Load clustering data
LC_kls <- readRDS("output/LC_kmeans/allLC_bodyID.rds")










###############################################################################
# Define items to analyze here
all_lcs <- c(
    "LC17", "LC12", "LC4", "LC6", "LPLC1", "LPLC2",
    "LC9", "LC10", "LC11", "LC13", "LC15", "LC16",
    "LC18", "LC20", "LC21", "LC22", "LC24", "LC25",
    "LC26", "LPLC4"
)

# Define clustering quality index to use
# Reference for indexes used in the package:
# https://cran.r-project.org/web/packages/clusterCrit/vignettes/clusterCrit.pdf
idx <- c("Silhouette")

# Using perceptually uniform color map
# More options here:
# https://cran.r-project.org/web/packages/cetcolor/vignettes/cet_color_schemes.html
c_map_single <- "r2" # Map to use for individual post plot
c_size_single <- 100 # Number of colors in the map
c_map_both <- "d9" # Map to use for both posts together
c_size_both <- 100 # Number of colors in the map

# Plot window size
win_siz <- 1000
###############################################################################







# Run for all the VPNs
for (pre_type in all_lcs) {

    # Grab the correct plane
    # NOTE: the FIRST plane specified is also use as projection plane
    planes <- get.plane(pre_type, "lobula")

    # Normalize the vector normal to the plane
    w_norm <- planes[1, 1:3]
    w0 <- planes[1, 4]

    # Pre ----
    # Extract all pre.type synapses
    pre <- con %>% filter(pre.type == pre_type)

    # Post ----
    # Get the neuron list from the posts bodyIDs
    n_list <- nlist[LC_kls[[pre_type]]$bodyid]

    # Drop potential NAs
    n_list <- n_list[!is.na(names(n_list))]
    # dropped <- setdiff(as.character(post$pre.bodyID), names(n_list))
    # if (length(dropped) != 0) {
    #     post <- subset(post, post$pre.bodyID != dropped)
    # }


    # 3D analysis ----
    # nopen3d()
    # par3d("windowRect" = c(100, 100, win_siz, win_siz))

    # Add neuron skeletons (to visualize the somas)
    # plot3d(n_list, soma = TRUE, col = "light gray", add = TRUE)

    # Show all endpoints
    # for (i in 1:length(n_list)) {
    #     neu <- n_list[[i]]
    #     ep <- neu$d[neu$EndPoints, ] %>%
    #         select(X, Y, Z)
    #     plot3d(ep, col = "black", size = 1, add = TRUE)
    # }

    # and the planes
    # planes3d(a = planes[1, 1], b = planes[1, 2], c = planes[1, 3], d = planes[1, 4], col = "blue", alpha = 0.2, add = TRUE)
    # planes3d(a = planes[-1, 1], b = planes[-1, 2], c = planes[-1, 3], d = planes[-1, 4], alpha = 0.2, add = TRUE)

    # Select preserved end points for each neuron
    for (i in 1:length(n_list)) {
        neu <- n_list[[i]]
        ep <- neu$d[neu$EndPoints, ]
        pts <- ep
        for (p in 1:nrow(planes)) {
            pts <- pts %>%
                mutate(LO = planes[p, 5] * (planes[p, 1] * X + planes[p, 2] * Y + planes[p, 3] * Z + planes[p, 4])) %>%
                filter(LO < 0) %>%
                select(X, Y, Z)
        }
        if (i == 1) {
            end_pts <- pts
        } else {
            end_pts <- rbind(end_pts, pts)
        }
    }

    # Show the preserved points
    # plot3d(end_pts, col = "green", size = 3, add = TRUE)

    # Project the preserved end points on the plane in yellow.
    end_pts.proj <- plane.proj(end_pts, w_norm, -w0)
    # plot3d(end_pts.proj, col = "yellow", size = 1, add = TRUE)

    # Find the centroid of the end points projections and plot in white
    center <- colMeans(end_pts.proj)
    # plot3d(center[1], center[2], center[3], col = "white", size = 10, add = TRUE)

    # Calculate centroid of lobula end points and plot in red
    for (i in 1:length(n_list)) {
        neu <- n_list[[i]]
        pts <- neu$d[neu$EndPoints, ]
        for (p in 1:nrow(planes)) {
            pts <- pts %>%
                mutate(LO = planes[p, 5] * (planes[p, 1] * X + planes[p, 2] * Y + planes[p, 3] * Z + planes[p, 4])) %>%
                filter(LO < 0) %>%
                select(X, Y, Z)
        }
        pts <- pts %>%
            colMeans() %>%
            t() %>%
            data.frame()
        pts$bodyid <- neu$bodyid
        pts$cluster <- LC_kls[[pre_type]]$cluster[LC_kls[[pre_type]]$bodyid == neu$bodyid]

        if (i == 1) {
            ctrs <- pts
        } else {
            ctrs <- rbind(ctrs, pts)
        }
    }

    # Clear rows with NaN
    ctrs <- ctrs[!is.na(ctrs[, 1]), ]

    # plot3d(ctrs, col = "red", size = 5, add = TRUE)

    # Project the centroid of the end points on the surface in blue
    ctrs.proj <- ctrs
    ctrs.proj[1:3] <- plane.proj(ctrs[1:3], w_norm, -w0)
    # plot3d(ctrs.proj, col = "blue", size = 5, add = TRUE)

    # 2D analysis ----
    # Define coordinate system in the plane
    # Use the projection of the z-axis on the plane to determine one of the new 2D
    # axes
    pointY <- c(plane.proj(t(center + c(0, 0, 2000)), w_norm, -w0))
    pointY <- pointY - center
    unitY <- -pointY / sqrt(sum(pointY * pointY))
    # arrow3d(center, center + 2000 * unitY, thikness = 1, col = "white")

    pointX <- cross(w_norm, pointY)
    unitX <- -pointX / sqrt(sum(pointX * pointX))
    # arrow3d(center, center + 2000 * unitX, thikness = 1, col = "white")

    # Project data on 2D plane ----
    # All considered endpoints
    Xep <- as.matrix(sweep(end_pts.proj, 2, center)) %*% unitX
    Yep <- as.matrix(sweep(end_pts.proj, 2, center)) %*% unitY
    end_pts.plane <- data.frame(cbind(Xep, Yep))

    # Centroids
    Xcm <- as.matrix(sweep(ctrs.proj[1:3], 2, center)) %*% unitX
    Ycm <- as.matrix(sweep(ctrs.proj[1:3], 2, center)) %*% unitY
    ctrs.plane <- cbind.data.frame(Xcm, Ycm, ctrs.proj$bodyid, ctrs.proj$cluster)
    colnames(ctrs.plane) <- c("X1", "X2", "bodyid", "cluster")

    # 2d maps of clustered centroid projections in the lobula ----
    # Create the color maps
    col_single <- cet_pal(c_size_single, name = c_map_single)

    # Evaluate parameters for the box plots
    dist_box <- 0.1 * abs(max(min(end_pts.plane$X1), min(end_pts.plane$X2)))
    width_box <- 0.75 * dist_box
    X1_rng <- abs(max(end_pts.plane$X1) - min(end_pts.plane$X1))
    X2_rng <- abs(max(end_pts.plane$X2) - min(end_pts.plane$X2))
    seg_len <- 0.05 * min(X1_rng, X2_rng)

    # Calculate the convex hull
    # For consistency use all the neuron in the pre so that it doesn't change based
    # on the chosen posts, but just on the VPN.
    hull <- pre %>%
        filter(post.type == pre_type) %>%
        group_by(pre.bodyID) %>%
        pull(pre.bodyID) %>%
        unique()

    # Get the neuron list from the posts bodyIDs
    n_list_hull <- nlist[as.character(hull)]

    # Drop potential NAs
    n_list_hull <- n_list_hull[!is.na(names(n_list_hull))]
    # dropped <- setdiff(as.character(post$pre.bodyID), names(n_list_hull))
    # post <- subset(post, post$pre.bodyID != dropped)


    # Select preserved end points for each neuron
    for (i in 1:length(n_list_hull)) {
        neu <- n_list_hull[[i]]
        ep <- neu$d[neu$EndPoints, ]
        pts <- ep
        for (p in 1:nrow(planes)) {
            pts <- pts %>%
                mutate(LO = planes[p, 5] * (planes[p, 1] * X + planes[p, 2] * Y + planes[p, 3] * Z + planes[p, 4])) %>%
                filter(LO < 0) %>%
                select(X, Y, Z)
        }
        if (i == 1) {
            end_pts_hull <- pts
        } else {
            end_pts_hull <- rbind(end_pts_hull, pts)
        }
    }

    # Project the hull end points on the 3D plane.
    end_pts_hull.proj <- plane.proj(end_pts_hull, w_norm, -w0)


    # Project the hull endpoints on the 2D plane
    Xep_hull <- as.matrix(sweep(end_pts_hull.proj, 2, center)) %*% unitX
    Yep_hull <- as.matrix(sweep(end_pts_hull.proj, 2, center)) %*% unitY
    end_pts_hull.plane <- data.frame(cbind(Xep_hull, Yep_hull))

    # Evaluate convex hull
    lo <- chull(end_pts_hull.plane)


    # Evaluate quality of clustering
    qidx <- intCriteria(
        as.matrix(ctrs.plane[1:2]),
        as.integer(ctrs.plane$cluster),
        idx
    )
    cat(toString(idx), "for", pre_type, "->", toString(qidx), "\n")

    # Plot projection of centroids
    p <- ggplot() +
        coord_fixed() +
        xlab("A-P axis, um") +
        ylab("D-V axis, um") +
        geom_polygon(
            data = end_pts_hull.plane[lo, ],
            aes(x = 0.008 * X1, y = 0.008 * X2),
            alpha = 0.3
        ) +
        geom_point(
            data = ctrs.plane,
            aes(x = 0.008 * X1, y = 0.008 * X2, col = cluster)
        ) +
        guides(colour = guide_legend(), size = guide_legend()) +
        ggtitle(paste(
            pre_type,
            "Dendritic cluster map (",
            toString(idx),
            ":",
            toString(qidx),
            ")"
        )) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
            axis.text.x = element_text(
                face = "bold",
                color = "black",
                size = 15,
                angle = 0
            ),
            axis.text.y = element_text(
                face = "bold",
                color = "black",
                size = 15,
                angle = 0
            ),
            plot.title = element_text(face = "bold", color = "black", size = 15),
            axis.title.y = element_text(size = 15),
            axis.title.x = element_text(size = 15)
        )
    print(p)
}
