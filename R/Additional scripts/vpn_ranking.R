# Created on Thu Jul 29 16:51:05 2021
# Name:    vpn_ranking.R
# Purpose: 
# Author:  Andrea Vaccari (avaccari@middlebury.edu)
#
# Generates:
# - ...
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
# - Run the cross correlation between the different VPN using one from one and 
#   the one from the other (pcorr and dist). Create a large corralation matrix
#   where the rows are distances and the columns are the pcorr. The diagonal
#   should be the best results. (VPN to use: LC4, LC6, LC9, LC10, LC11, LC12, 
#   LC13, LC15, LC16, LC17, LC18, LC20, LC21, LC22, LC24, LC25, LC26, LPLC1, 
#   LPLC2, LPLC4, LC28b) For the actual paper we only need LC4 and LPLC2. 
#   Different script.
