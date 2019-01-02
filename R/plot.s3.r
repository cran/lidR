# ===============================================================================
#
# PROGRAMMERS:
#
# jean-romain.roussel.1@ulaval.ca  -  https://github.com/Jean-Romain/lidR
#
# COPYRIGHT:
#
# Copyright 2016 Jean-Romain Roussel
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

#' Plot voxelized LiDAR data
#'
#' This function implements a 3D plot method for 'lasmetrics3d' objects
#'
#' @param x An object of the class \code{'lasmetrics3d'}
#' @param y Unused (inherited from R base)
#' @param color characters. The field used to color the points. Default is Z coordinates. Or a vector of colors.
#' @param colorPalette characters. A color palette name. Default is \code{height.colors} provided by the package lidR
#' @param bg The color for the background. Default is black.
#' @param trim numeric. Enables trimming of values when outliers break the color palette range.
#' Default is 1 meaning that the whole range of the values is used for the color palette.
#' 0.9 means that 10% of the highest values are not used to define the color palette.
#' In this case the values higher than the 90th percentile are set to the highest color. They are not removed.
#' @param \dots Supplementary parameters for \link[rgl:points3d]{points3d} if the display method is "points".
#' @examples
#' LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#' lidar = readLAS(LASfile)
#'
#' voxels = grid_metrics3d(lidar, list(Imean = mean(Intensity)))
#' plot(voxels, color = "Imean", colorPalette = heat.colors(50), trim=0.99)
#' @seealso
#' \link[lidR:grid_metrics3d]{grid_metrics3d}
#' \link[rgl:points3d]{points3d}
#' \link[lidR:height.colors]{height.colors}
#' \link[lidR:forest.colors]{forest.colors}
#' \link[grDevices:heat.colors]{heat.colors}
#' \link[grDevices:colorRamp]{colorRampPalette}
#' @export
#' @method plot lasmetrics3d
plot.lasmetrics3d = function(x, y, color = "Z", colorPalette = height.colors(50), bg = "black", trim = 1, ...)
{
  inargs <- list(...)

  inargs$col = color

  if (length(color) == 1)
  {
    if (color %in% names(x))
    {
      data = unlist(x[,color, with = FALSE])

      if (is.numeric(data))
      {
        inargs$col = set.colors(data, colorPalette, trim)
        inargs$col[is.na(inargs$col)] = "lightgray"
      }
      else if (is.character(data))
        inargs$col = data
    }
  }

  rgl::open3d()
  rgl::rgl.bg(color = bg)
  do.call(rgl::points3d, c(list(x = x$X, y = x$Y, z = x$Z), inargs))
}

#' Add a spatial object to a point cloud scene
#'
#' Add a \code{RasterLayer} object that represents a digital terrain model or a
#' \code{SpatialPointsDataFrame} that represents tree tops to a point cloud scene. To add elements
#' to a scene with a point cloud plotted with the function plot from lidR, the functions \code{add_*}
#' take as first argument the output of the plot function (see examples), because the plot function
#' does not plot the actual coordinates of the point cloud, but offsetted values. See function
#' \link[lidR:plot]{plot} and its argument \code{clear_artifacts} for more details.
#'
#' @param dtm An object of the class \code{RasterLayer}
#' @param bg The color for the background. Default is black.
#' @param \dots Supplementary parameters for \link[rgl:surface3d]{surface3d} or
#' \link[rgl:spheres3d]{spheres3d}.
#' @param x The output of the function plot used with a LAS object.
#' @param ttops A SpatialPointsDataFrame that contains tree tops coordinates.
#' @param z character. The name of the attribute that contains the height of the tree tops.
#' @param clear_artifacts logical. It is a known and documented issue that 3D visualisation with
#' \code{rgl} displays artifacts. The points and lines are inaccurately positioned in the space and thus
#' the rendering may look false or weird. This is because \code{rgl} computes with single precision \code{float}.
#' To fix this, the objects are shifted to (0,0) to reduce the number of digits needed to represent
#' their coordinates. The drawback is that the objects are not plotted at their actual coordinates.
#'
#' @name plot_3d
#' @examples
#' LASfile <- system.file("extdata", "Topography.laz", package="lidR")
#' las = readLAS(LASfile)
#'
#' dtm = grid_terrain(las, algorithm = tin())
#' ttops <- tree_detection(las, lmf(ws = 5))
#'
#' plot_dtm3d(dtm)
#'
#' x = plot(las)
#' add_dtm3d(x, dtm)
#' add_treetops3d(x, ttops)
#'
#' \dontrun{
#' library(magrittr)
#' plot(las) %>% add_dtm3d(dtm) %>% add_treetops3d(ttops)
#' }
NULL

#' @rdname plot_3d
#' @export
plot_dtm3d = function(dtm, bg = "black", clear_artifacts = TRUE, ...)
{
  rgl::open3d()
  rgl::rgl.bg(color = bg)
  shift = c(0,0)

  if (clear_artifacts)
  {
    bbox  <- raster::extent(dtm)
    shift <- c(bbox@xmin, bbox@ymin)
  }

  add_dtm3d(shift, dtm, ...)
}

#' @rdname plot_3d
#' @export
add_dtm3d = function(x, dtm, ...)
{
  args <- list(...)

  assert_is_numeric(x)
  assert_is_of_length(x, 2)

  if (!is(dtm, "RasterLayer"))
    stop("'dtm' is not RasterLayer.")

  if (is.null(args$front))
    args$front <- "lines"

  if (is.null(args$col))
    args$col <- "white"

  mx <-  t(apply(raster::as.matrix(dtm), 2, rev))
  x_ <- sort(raster::xFromCol(dtm, 1:raster::ncol(dtm))) - x[1]
  y_ <- sort(raster::yFromRow(dtm, 1:raster::nrow(dtm))) - x[2]

  args$x <- x_
  args$y <- y_
  args$z <- mx

  do.call(rgl::surface3d, args)
  return(invisible(x))
}

#' @rdname plot_3d
#' @export
add_treetops3d = function(x, ttops, z = "Z", ...)
{
  args <- list(...)

  assert_is_numeric(x)
  assert_is_of_length(x, 2)

  if (!is(ttops, "SpatialPointsDataFrame"))
    stop("'ttops' is not a SpatialPointsDataFrame")

  if (is.null(args$size))
    args$size <- 5

  if (is.null(args$col))
    args$col <- "red"

  args$add <- TRUE
  args$x   <- ttops@coords[,1] - x[1]
  args$y   <- ttops@coords[,2] - x[2]
  args$z   <- ttops@data[[z]]

  do.call(rgl::spheres3d, args)
  return(invisible(x))
}

