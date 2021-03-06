# ===============================================================================
#
# PROGRAMMERS:
#
# jean-romain.roussel.1@ulaval.ca  -  https://github.com/Jean-Romain/lidR
#
# COPYRIGHT:
#
# Copyright 2016-2018 Jean-Romain Roussel
#
# This file is part of lidR R package.
#
# lidR is free software: you can redistribute it and/or modify
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#
# ===============================================================================


#' Compute the hull of each tree.
#'
#' Compute the hull of each segmented tree. The hull can be convex, concave or a bounding box (see
#' details and references).
#'
#' The concave hull method under the hood is described in Park & Oh (2012). The function relies on
#' the \link[concaveman:concaveman]{concaveman} function.
#'
#' @template param-las
#' @param type character. Hull type. Can be 'convex', 'concave' or 'bbox'.
#' @param concavity numeric. If \code{type = "concave"}, a relative measure of concavity. 1 results
#' in a relatively detailed shape, Infinity results in a convex hull.
#' @param length_threshold numeric. If \code{type = "concave"}, when a segment length is below this
#' threshold, no further detail is added. Higher values result in simpler shapes.
#' @param attribute character. The attribute where the ID of each tree is stored. In lidR, the default is
#' "treeID".
#' @param func formula. An expression to be applied to each tree. It works like in \link{grid_metrics}
#' \link{voxel_metrics} or \link{tree_metrics} and computes, in addition to the hulls a set of metrics
#' for each tree.
#'
#' @return A \code{SpatialPolygonsDataFrame}. If a tree has less than 4 points it is not considered.
#'
#' @template LAScatalog
#' @template section-supported-option-tree_detection
#'
#' @export
#'
#' @references Park, J. S., & Oh, S. J. (2012). A new concave hull algorithm and concaveness measure
#' for n-dimensional datasets. Journal of Information science and engineering, 28(3), 587-600.
#'
#' @examples
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' poi = "-drop_z_below 0 -inside 481280 3812940 481320 3812980"
#' las = readLAS(LASfile, select = "xyz0", filter = poi)
#'
#' # NOTE: This dataset is already segmented
#'
#' #plot(las, color = "treeID", colorPalette = pastel.colors(200))
#'
#' # Only the hulls
#' convex_hulls = delineate_crowns(las)
#' plot(convex_hulls)
#'
#' # The hulls + some user-defined metrics
#' convex_hulls = delineate_crowns(las, func = ~list(Zmean = mean(Z)))
#' convex_hulls
#'
#' # The bounding box
#' bbox_hulls = delineate_crowns(las, "bbox")
#' plot(bbox_hulls)
#'
#' \dontrun{
#' # With concave hull (longer to compute)
#' concave_hulls = delineate_crowns(las, "concave")
#' plot(concave_hulls)
#'
#' spplot(convex_hulls, "ZTOP")
#' spplot(convex_hulls, "Zmean")
#' }
delineate_crowns = function(las, type = c("convex", "concave", "bbox"), concavity = 3, length_threshold = 0, func = NULL, attribute = "treeID")
{
  UseMethod("delineate_crowns", las)
}

#' @export
delineate_crowns.LAS = function(las, type = c("convex", "concave", "bbox"), concavity = 3, length_threshold = 0, func = NULL, attribute = "treeID")
{
  type <- match.arg(type)
  assert_is_a_number(concavity)
  assert_all_are_non_negative(concavity)
  assert_is_a_number(length_threshold)
  assert_all_are_non_negative(length_threshold)
  assert_is_a_string(attribute)

  # concaveman is only a 'suggested' dependency
  if (type == "concave") assert_package_is_installed("concaveman")

  # Pointer on function C style coding
  if      (type == "convex")  fhull <- stdtreehullconvex
  else if (type == "concave") fhull <- stdtreehullconcave
  else if (type == "bbox")    fhull <- stdtreehullbbox

  # Hulls computation -- aggregation by tree
  X <- Y <- Z <- NULL
  hulls <- las@data[, if (!anyNA(.BY)) fhull(X,Y,Z,.GRP, concavity, length_threshold), by = attribute]

  if (nrow(hulls) == 0)
  {
    warning("No tree found. NULL returned.", call. = FALSE)
    return(NULL)
  }

  # Convert to SpatialPolygons
  spoly <- sp::SpatialPolygons(hulls[["poly"]])
  for (i in 1:length(spoly)) spoly@polygons[[i]]@ID <- as.character(i)

  # Compute metrics
  data = hulls[,1:4]
  if (!is.null(func))
  {
    func    <- lazyeval::f_interp(func)
    call    <- lazyeval::as_call(func)
    metrics <- las@data[, if (!anyNA(.BY)) eval(call), by = attribute]
    remove  <- metrics[[attribute]] %in% hulls[[attribute]]
    metrics <- metrics[remove]
    data    <- cbind(data, metrics[,-1])
  }

  data.table::setDF(data)
  spdf <- sp::SpatialPolygonsDataFrame(spoly, data)
  sp::proj4string(spdf) <- las@proj4string
  return(spdf)
}

#' @export
delineate_crowns.LAScluster = function(las, type = c("convex", "concave", "bbox"), concavity = 3, length_threshold = 0, func = NULL, attribute = "treeID")
{
  x <- readLAS(las)
  if (is.empty(x)) return(NULL)
  metrics <- delineate_crowns(x, type, concavity, length_threshold, func, attribute)
  bbox    <- raster::extent(las)

  coords = metrics@data[,c("XTOP", "YTOP")]
  coords$ID = 1:length(metrics)
  sp::coordinates(coords) <- ~XTOP+YTOP
  coords <- raster::crop(coords, bbox)
  return(metrics[coords$ID,])
}

#' @export
delineate_crowns.LAScatalog = function(las, type = c("convex", "concave", "bbox"), concavity = 3, length_threshold = 0, func = NULL, attribute = "treeID")
{
  options <- list(need_buffer = TRUE, drop_null = TRUE, need_output_file = FALSE, automerge = TRUE)
  output  <- catalog_apply(las, delineate_crowns, type = type, concavity = concavity, length_threshold = length_threshold, func = func, attribute = attribute, .options = options)
  return(output)
}

