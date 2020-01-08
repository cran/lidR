# ===============================================================================
#
# PROGRAMMERS:
#
# jean-romain.roussel.1@ulaval.ca  -  https://github.com/Jean-Romain/lidR
#
# COPYRIGHT:
#
# Copyright 2016-2019 Jean-Romain Roussel
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

#' Plot a LAS* object
#'
#' Plot displays a 3D interactive windows-based on rgl for \link{LAS} objects\cr\cr
#' Plot displays an interactive view for \link[lidR:LAScatalog-class]{LAScatalog} objects with pan and
#' zoom capabilities based on \link[mapview:mapview-package]{mapview}. If the coordinate reference
#' system (CRS) of the \code{LAScatalog} is non empty, the plot can be displayed on top of base maps
#' (satellite data, elevation, street, and so on).\cr\cr
#' Plot displays a \link[lidR:LASheader-class]{LASheader} object exactly like it displays a LAScatalog
#' object.
#'
#'
#' @param x A \code{LAS*} object
#' @param y Unused (inherited from R base)
#' @param color characters. The attribute used to color the point cloud. Default is Z coordinates. RGB
#' is an allowed string even if it refers to three attributes simultaneously.
#' @param colorPalette characters. A vector of colors such as that generated by heat.colors,
#' topo.colors, terrain.colors or similar functions. Default is \code{"auto"} providing an automatic
#' coloring depending on the argument \code{color}
#' @param bg The color for the background. Default is black.
#' @param trim numeric. Enables trimming of values when outliers break the color palette range.
#' Every point with a value higher than \code{trim} will be plotted with the highest color.
#' @param clear_artifacts logical. It is a known and documented issue that the 3D visualisation with
#' \code{rgl} displays artifacts. The points look aligned and/or regularly spaced in some view angles.
#' This is because \code{rgl} computes with single precision \code{float}. To fix that the point
#' cloud is shifted to (0,0) to reduce the number of digits needed to represent its coordinates.
#' The drawback is that the point cloud is not plotted at its actual coordinates.
#' @param nbits integer. If \code{color = RGB} it assumes that RGB colors are coded on 16 bits as described
#' in the LAS format specification. However, this is not always respected. If the colors are stored
#' on 8 bits set this parameter to 8.
#' @param axis logical. Display axis on XYZ coordinates.
#' @param legend logical. Display a gradient color legend.
#' @param backend character. Can be \code{"rgl"} or \code{"pcv"}. If \code{"rgl"} is chosen
#' the display relies on the \code{rgl} package. If \code{"pcv"} is chosen it relies on the
#' \code{PointCloudViewer} package, which is much more efficient and can handle million of points
#' using less memory. \code{PointCloudViewer} is not available on CRAN yet and should
#' be installed from github (see. \url{https://github.com/Jean-Romain/PointCloudViewer}).
#'
#' @param mapview logical. If \code{FALSE} the catalog is displayed in a regular plot from R base.
#' @param chunk_pattern logical. Display the current chunk pattern used to process the catalog.
#'
#' @param ... Will be passed to \link[rgl:points3d]{points3d} (LAS) or \link[graphics:plot]{plot}
#' if \code{mapview = FALSE} or to \link[mapview:mapView]{mapview} if \code{mapview = TRUE} (LAScatalog).
#'
#' @examples
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- readLAS(LASfile)
#'
#' plot(las)
#' plot(las, color = "Intensity")
#'
#' # If outliers break the color range, use the trim parameter
#' plot(las, color = "Intensity", trim = 150)
#'
#' plot(las, color = "Classification")
#'
#' # This dataset is already tree segmented
#' plot(las, color = "treeID")
#'
#' # single file catalog using data provided in lidR
#' ctg = readLAScatalog(LASfile)
#' plot(ctg)
#'
#' @export
#' @method plot LAS
setGeneric("plot", function(x, y, ...)
  standardGeneric("plot"))

#' @rdname plot
setMethod("plot", signature(x = "LAS", y = "missing"), function(x, y, color = "Z", colorPalette = "auto", bg = "black", trim = Inf, backend = c("rgl", "pcv"), clear_artifacts = TRUE, nbits = 16, axis = FALSE, legend = FALSE, ...)
{
  plot.LAS(x, y, color, colorPalette, bg, trim, backend, clear_artifacts, nbits, axis, legend, ...)
})

#' @export
#' @rdname plot
setMethod("plot", signature(x = "LAScatalog", y = "missing"), function(x, y, mapview = FALSE, chunk_pattern = FALSE, ...)
{
  plot.LAScatalog(x, y, mapview, chunk_pattern, ...)
})

#' @export
#' @rdname plot
setMethod("plot", signature(x = "LASheader", y = "missing"), function(x, y, mapview = FALSE, ...)
{
  epsg <- epsg(x)
  PHB  <- x@PHB
  crs  <- tryCatch({sp::CRS(glue::glue("+init=epsg:{epsg}"))}, error = function(e) return(sp::CRS()))
  xmin <- PHB[["Min X"]]
  xmax <- PHB[["Max X"]]
  ymin <- PHB[["Min Y"]]
  ymax <- PHB[["Max Y"]]
  mtx  <- matrix(c(xmin, xmax, ymin, ymax)[c(1, 1, 2, 2, 1, 3, 4, 4, 3, 3)], ncol = 2)
  Sr   <- sp::Polygons(list(sp::Polygon(mtx)), "1")
  Sr   <- sp::SpatialPolygons(list(Sr), proj4string = crs)

  names(PHB) <- make.names(names(PHB))

  if (!is.null(PHB[["Number.of.points.by.return"]]))
  {
    PHB[["Number.of.1st.return"]] <- PHB[["Number.of.points.by.return"]][1]
    PHB[["Number.of.2nd.return"]] <- PHB[["Number.of.points.by.return"]][2]
    PHB[["Number.of.3rd.return"]] <- PHB[["Number.of.points.by.return"]][3]
    PHB[["Number.of.4th.return"]] <- PHB[["Number.of.points.by.return"]][4]
    PHB[["Number.of.5th.return"]] <- PHB[["Number.of.points.by.return"]][5]
    PHB[["Number.of.points.by.return"]] <- NULL
    PHB[["Global.Encoding"]] <- NULL
  }

  res <- new("LAScatalog")
  res@bbox <- Sr@bbox
  res@proj4string <- Sr@proj4string
  res@plotOrder <- Sr@plotOrder
  res@data <- data.table::as.data.table(PHB)
  res@polygons <- Sr@polygons

  plot.LAScatalog(res, mapview = mapview, ...)
})

plot.LAScatalog = function(x, y, mapview = FALSE, chunk_pattern = FALSE, ...)
{
  assert_is_a_bool(mapview)
  assert_is_a_bool(chunk_pattern)

  if (mapview)
  {
    if (!requireNamespace("mapview", quietly = TRUE))
    {
      message("'mapview' is required to display the LAScatalog interactively.") # nocov
      mapview <- FALSE # nocov
    }
  }

  if (mapview)
  {
    LAScatalog <- x
    mapview::mapview(LAScatalog, ...)
  }
  else if (chunk_pattern)
  {
    opt_progress(x) <- TRUE
    catalog_makecluster(x)
    return(invisible())
  }
  else
  {
    # New feature from v2.2.0 to do not process some tiles
    process <- x@data$process
    if (is.null(process)) process <- rep(TRUE, nrow(x@data))
    if (!is.logical(process)) {
      warning("The attribute 'process' of the catalog is not logical.", call. = FALSE)
      process <- rep(TRUE, nrow(x@data))
    }

    alpha   <- ifelse(process, 0.15, 0.03)
    param   <- list(...)
    xmin    <- min(x@data$Min.X)
    xmax    <- max(x@data$Max.X)
    ymin    <- min(x@data$Min.Y)
    ymax    <- max(x@data$Max.Y)
    xcenter <- (xmin + xmax)/2
    ycenter <- (ymin + ymax)/2
    col    <- grDevices::rgb(0, 0, 1, alpha = alpha)

    if (is.null(param$xlim)) param$xlim <- c(xmin, xmax)
    if (is.null(param$ylim)) param$ylim <- c(ymin, ymax)
    if (is.null(param$xlab)) param$xlab <- ""
    if (is.null(param$ylab)) param$ylab <- ""
    if (is.null(param$asp))  param$asp  <- 1
    if (!is.null(param$col)) col <- param$col

    param$col <- "white"
    param$x   <- xcenter
    param$y   <- ycenter

    op <- graphics::par(mar = c(2.5,2.5,1,1) + 0.1)

    if (is.null(param$add)) do.call(graphics::plot, param)

    graphics::rect(x@data$Min.X, x@data$Min.Y, x@data$Max.X, x@data$Max.Y, col = col)
    graphics::par(op)

    return(invisible())
  }
}

plot.LAS = function(x, y, color = "Z", colorPalette = "auto", bg = "black", trim = Inf, backend = c("rgl", "pcv"), clear_artifacts = TRUE, nbits = 16, axis = FALSE, legend = FALSE, ...)
{
  backend <- match.arg(backend)
  use_pcv <- backend == "pcv"
  use_rgl <- !use_pcv
  has_pcv <- "PointCloudViewer" %in% rownames(utils::installed.packages())
  has_col <- color %in% names(x@data)
  use_rgb <- color == "RGB"
  has_rgb <- all(c("R", "G", "B") %in% names(x@data))
  maxcol  <- 2^nbits - 1
  autocol <- all(colorPalette == "auto")
  value_index <- FALSE

  if (is.empty(x))         stop("Cannot display an empty point cloud", call. = FALSE)
  if (use_pcv & !has_pcv)  stop("'PointCloudViewer' package is needed. Please read documentation.", call. = FALSE) # nocov
  if (length(color) > 1)   stop("'color' should contain a single value.", call. = FALSE)
  if (!use_rgb & !has_col) stop("'color' should refer to an attribute of the LAS data.", call. = FALSE)
  if (use_rgb & !has_rgb)  stop("No 'RGB' attributes found.", call. = FALSE)

  if (autocol)
  {
    if (color == "Z")
      colorPalette = height.colors(50)
    else if (color == "Intensity")
      colorPalette = grDevices::heat.colors(50)
    else if (color  == "Classification")
    {
      colorPalette = lasclass.colors()
      clmin = min(x@data[["Classification"]])
      clmax = max(x@data[["Classification"]])
      trim  = min(length(colorPalette), clmax+1)
      value_index = TRUE
    }
    else if (color == "ScanAngleRank" | color == "ScanAngle")
      colorPalette = height.colors(50)
    else if (color == "ReturnNumber")
      colorPalette = rev(c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))
    else if (color == "treeID")
      colorPalette = pastel.colors(200)
    else
      colorPalette = height.colors(50)
  }

  if (use_rgb & use_pcv)
    col <- "RGB"
  else if (use_rgb & use_rgl)
    col <- grDevices::rgb(x@data[["R"]]/maxcol, x@data[["G"]]/maxcol, x@data[["B"]]/maxcol)
  else
    col <- x@data[[color]]

  args <- list(...)
  if (is.null(args$size))
    args$size <- 1.5

  if (use_rgl)
    lasplot <- .plot_with_rgl
  else
    lasplot <- .plot_with_pcv # nocov

  return(lasplot(x, bg, col, colorPalette, trim, clear_artifacts, axis, legend, args, value_index))
}

.plot_with_rgl = function(las, bg, col, pal, trim, clear_artifacts, axis, legend, args, value_index)
{
  fg   <- grDevices::col2rgb(bg)
  fg   <- grDevices::rgb(t(255 - fg)/255)
  minx <- min(las@data$X)
  miny <- min(las@data$Y)

  if (is.numeric(col))
  {
    mincol <- min(col, na.rm = TRUE)
    maxcol <- min(max(col, na.rm = TRUE), trim)
    col <- set.colors(col, pal, trim, value_index)
  }
  else if (is.character(col))
  {
    legend <- FALSE
    col <- col
  }
  else if (is.logical(col))
  {
    mincol <- 0
    maxcol <- 1
    col <- set.colors(as.numeric(col), pal)
  }

  col[is.na(col)] <- "lightgray"

  with <- c(list(x = las@data$X, y = las@data$Y, z = las@data$Z, col = col), args)

  if (clear_artifacts)
  {
    with$x <- with$x - minx
    with$y <- with$y - miny
  }

  rgl::open3d()
  rgl::rgl.bg(color = bg)
  do.call(rgl::points3d, with)

  if (axis)
  {
    rgl::axis3d("x", col = fg)
    rgl::axis3d("y", col = fg)
    rgl::axis3d("z", col = fg)
  }

  if (legend)
  {
    # nocov because this fails on some flavor on CRAN
    f <- .plot_scale_gradient(mincol, maxcol, fg, pal, bg) # nocov
    rgl::bg3d(texture = f, col = "white") # nocov
  }

  .pan3d(2)

  if (clear_artifacts)
    return(invisible(c(minx, miny)))
  else
    return(invisible(c(0,0)))
}

# nocov start
.plot_with_pcv = function(las, bg, col, pal, trim, clear_artifacts, axis, legend, args, value_index)
{
  if (is.character(col))
  {
    if (col == "RGB")
      eval(parse(text = "PointCloudViewer::plot_xyzrgb(las@data$X, las@data$Y, las@data$Z, las@data$R, las@data$G, las@data$B, args$size)"))
    else
      stop("Unexpected error.", call. = FALSE)
  }
  else
  {
    if (value_index)
      id <- col + 1
    else
    {
      if (!is.infinite(trim)) col[col > trim] <- trim
      id <- cut(col, length(pal), labels = FALSE)
    }

    eval(parse(text = "PointCloudViewer::plot_xyzcol(las@data$X, las@data$Y, las@data$Z, pal, id, args$size)"))
  }

  return(invisible(c(0,0)))
}
# nocov end

# nocov start
.plot_scale_gradient = function(min.col, max.col, text.col, scale.col, bg)
{
  f <- tempfile(fileext = ".png")
  labels <- pretty(c(min.col, max.col))
  ncol   <- length(scale.col)
  nlab   <- length(labels)
  xl <- 1 ; yb <- 1 ; xr <- 1.1 ; yt <- 2
  grDevices::png(f, 1920, 1080, bg = bg)
  graphics::layout(matrix(1:2, nrow = 1), widths = c(0.9,0.1))
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  graphics::plot(0, ann = FALSE, type = "n", axes = FALSE)
  graphics::par(mar = c(5.1, 0.5, 4.1, 0))
  graphics::plot(NA, type = "n", ann = FALSE, xlim = c(1,2), ylim = c(1,2), xaxt = "n", yaxt = "n", bty = "n")
  graphics::rect(xl, utils::head(seq(yb, yt, (yt - yb)/ncol), -1), xr, utils::tail(seq(yb, yt, (yt - yb)/ncol), -1), col = scale.col, border = NA)
  graphics::mtext(labels, side = 2, at = seq(yb, yt, length.out = nlab), las = 2, cex = 1.2, col = text.col)
  grDevices::dev.off()
  return(f)
}
# nocov end

# From rgl.setMouseCallbacks man page
# nocov start
.pan3d <- function(button, dev = rgl::rgl.cur(), subscene = rgl::currentSubscene3d(dev))
{
  start <- list()

  begin <- function(x, y)
  {
    activeSubscene <- rgl::par3d("activeSubscene", dev = dev)
    start$listeners <<- rgl::par3d("listeners", dev = dev, subscene = activeSubscene)

    for (sub in start$listeners)
    {
      init <- rgl::par3d(c("userProjection","viewport"), dev = dev, subscene = sub)
      init$pos <- c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5)
      start[[as.character(sub)]] <<- init
    }
  }

  update <- function(x, y)
  {
    for (sub in start$listeners)
    {
      init <- start[[as.character(sub)]]
      xlat <- 2*(c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5) - init$pos)
      mouseMatrix <- rgl::translationMatrix(xlat[1], xlat[2], xlat[3])
      rgl::par3d(userProjection = mouseMatrix %*% init$userProjection, dev = dev, subscene = sub )
    }
  }
  rgl::rgl.setMouseCallbacks(button, begin, update, dev = dev, subscene = subscene)
}
# nocov end
