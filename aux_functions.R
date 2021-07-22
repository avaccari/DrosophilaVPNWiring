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
#

# Project on plane function
plane.proj <- function(points, plane.normal, plane.offset) {
  dist <- as.matrix(points) %*% plane.normal
  proj <- points - (dist - plane.offset) %*% plane.normal
}

# Project on a line function
line.proj <- function(points, line.unit.vector, origin) {
  transl <- sweep(points, 2, origin)
  proj <- as.matrix(transl) %*% line.unit.vector
  proj <- sweep(drop(proj %o% line.unit.vector), 2, origin, FUN='+')
}

# Distance from the plane function
plane.dist <- function(points, plane.normal, origin) {
  transl <- sweep(points, 2, origin)
  dist <- as.matrix(transl) %*% plane.normal
}
