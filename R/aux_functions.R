# Created on Tue Jul 20 12:40:41 2021
# Name:    aux_functions.R
# Purpose: Contains functions that are used in multiple scripts
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
# - Start with LC6

# Project on plane function
plane.proj <- function(points, plane.normal, plane.offset) {
  dist <- as.matrix(points) %*% plane.normal
  proj <- points - (dist - plane.offset) %*% drop(plane.normal)
}

# Project on a line function
line.proj <- function(points, line.unit.vector, origin) {
  transl <- sweep(points, 2, origin)
  proj <- as.matrix(transl) %*% drop(line.unit.vector)
  proj <- sweep(drop(proj %o% line.unit.vector), 2, origin, FUN='+')
}

# Distance from the plane function
plane.dist <- function(points, plane.normal, origin) {
  transl <- sweep(points, 2, origin)
  dist <- as.matrix(transl) %*% drop(plane.normal)
}

# Get plane
# Define two set o planes: one to identify the region to be considered lobula
# and one to identify the region to be considered glomerulus.
# It is possible to identify multiple planes and which side of the plane should
# be included.
# The format is a matrix where each line is a plane. The first 4 values are the
# a, b, c, and d of the plane (with a, b, and c identifying the normal to the
# plane and d the offset) and the fifth value specifies if we should preserve 
# the synapse above (1) or below (-1) the plane.
# Each plane is considered one after the other and the final result is the
# intersection of the setes of synapses preserved by each plane.
# E.g.
# planes <- matrix(c(1, 0, 0, -11000, -1,
#                    1, 0, 0, -9000, 1), ncol=5, byrow=TRUE)
get.plane <- function(pre='LC4', type='glomerulus') {
  if (type=='glomerulus' | type=='glo') {
    switch(
      pre,
      'LC4' = matrix(c(0.5, 0.3, -0.4, 1500, 1), ncol=5, byrow=TRUE),
      'LC6' = matrix(c(0.7687760,  0.5688942, -0.2921349, -9200, 1,
                       0.7687760,  0.5688942, -0.2921349, -16200, -1,
                       0.2, 0.1, 1, -23000, 1), ncol=5, byrow=TRUE),
      'LC9' = matrix(c(0.5, 0.3, -0.4, -2700, 1,
                       0.2, -0.2, 1, -16000, 1), ncol=5, byrow=TRUE),
      'LC10' = matrix(c(0.5, 0.3, -0.4, -6500, 1), ncol=5, byrow=TRUE),
      'LC11' = matrix(c(0.7687760,  0.5688942, -0.2921349, -10000, 1), ncol=5, byrow=TRUE),
      'LC12' = matrix(c(0.7687760,  0.5688942, -0.2921349, -7200, 1,
                        0.7687760,  0.5688942, -0.2921349, -11200, -1), ncol=5, byrow=TRUE),
      'LC13' = matrix(c(0.5, 0.3, -0.4, 1000, 1), ncol=5, byrow=TRUE),
      'LC15' = matrix(c(0.5, 0.3, -0.4, 1500, 1), ncol=5, byrow=TRUE),
      'LC16' = matrix(c(0.5, 0.3, -0.4, -500, 1), ncol=5, byrow=TRUE),
      'LC17' = matrix(c(0.7687760,  0.5688942, -0.2921349, -7200, 1,
                        0.7687760,  0.5688942, -0.2921349, -11200, -1), ncol=5, byrow=TRUE),
      'LC18' = matrix(c(0.5, 0.3, -0.4, 1500, 1), ncol=5, byrow=TRUE),
      'LC20' = matrix(c(0.5, 0.3, -0.4, 100, 1), ncol=5, byrow=TRUE),
      'LC21' = matrix(c(0.5, 0.3, -0.4, 500, 1), ncol=5, byrow=TRUE),
      'LC22' = matrix(c(1, 0, 0, -14500, 1), ncol=5, byrow=TRUE),
      'LC24' = matrix(c(0.5, 0.3, -0.4, 500, 1), ncol=5, byrow=TRUE),
      'LC25' = matrix(c(0.5, 0.3, -0.4, 800, 1), ncol=5, byrow=TRUE),
      'LC26' = matrix(c(0.5, 0.3, -0.4, 800, 1), ncol=5, byrow=TRUE),
      'LPLC1' = matrix(c(0.7687760,  0.5688942, -0.8921349, 6000, 1), ncol=5, byrow=TRUE),
      'LPLC2' = matrix(c(1, 0, -0.2, -3500, 1), ncol=5, byrow=TRUE),
      'LPLC4' = matrix(c(1, 0, -0.2, -3500, 1), ncol=5, byrow=TRUE)
    )
  } else if (type=='lobula' | type=='lob') {
    switch(
      pre,
      'LC4' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1), ncol=5, byrow=TRUE),
      'LC6' = matrix(c(0.80687760,  0.5688942, -0.3921349, -4000, 1), ncol=5, byrow=TRUE),
      'LC9' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1), ncol=5, byrow=TRUE),
      'LC10' = matrix(c(0.6687760,  0.5688942, -0.5921349, 1200, 1), ncol=5, byrow=TRUE),
      'LC11' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1,
                        -1.4, 0.1, 1.1, -17000, -1), ncol=5, byrow=TRUE),
      'LC12' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1), ncol=5, byrow=TRUE),
      'LC13' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5800, 1), ncol=5, byrow=TRUE),
      'LC15' = matrix(c(0.6687760,  0.4688942, -0.2921349, -3700, 1), ncol=5, byrow=TRUE),
      'LC16' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1,
                        0.1, 0.1, 1.1, -42000, 1), ncol=5, byrow=TRUE),
      'LC17' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1), ncol=5, byrow=TRUE),
      'LC18' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1,
                        0.7687760,  0.5688942, -0.0921349, -5200, -1), ncol=5, byrow=TRUE),
      'LC20' = matrix(c(0.6687760,  0.4688942, -0.3921349, -2200, 1), ncol=5, byrow=TRUE),
      'LC21' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1,
                        0.7687760,  0.7688942, 1.0921349, -29200, -1), ncol=5, byrow=TRUE),
      'LC22' = matrix(c(0.7687760,  0.5688942, -0.3921349, -3200, 1), ncol=5, byrow=TRUE),
      'LC24' = matrix(c(0.6687760,  0.4688942, -0.2521349, -5200, 1), ncol=5, byrow=TRUE),
      'LC25' = matrix(c(0.6687760,  0.4688942, -0.2521349, -5200, 1), ncol=5, byrow=TRUE),
      'LC26' = matrix(c(0.6687760,  0.4688942, -0.2521349, -5200, 1), ncol=5, byrow=TRUE),
      'LPLC1' = matrix(c(0.7687760,  0.5688942, -0.3921349, -3500, 1,
                         0.7687760,  0.5688942, -0.2921349, 500, -1,
                         0.7687760,  0.5688942, -0.8921349, 6000, 1), ncol=5, byrow=TRUE),
      'LPLC2' = matrix(c(0.7687760,  0.5688942, -0.2921349, -5200, 1,
                         0.7687760,  0.5688942, -0.2921349, -1000, -1,
                         0.7687760,  0.5688942, -0.8921349, 6300, 1), ncol=5, byrow=TRUE),
      'LPLC4' = matrix(c(0.7687760,  0.5688942, -0.2921349, -6500, 1,
                         0.7687760,  0.5688942, -0.2921349, -1000, -1,
                         0.7687760,  0.5688942, -0.9921349, 7200, 1), ncol=5, byrow=TRUE)
    )
  } else {
    print('Plane undefined!\n')
  }
}
