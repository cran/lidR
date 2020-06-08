# ===============================================================================
#
# PROGRAMMERS:
#
# jean-romain.roussel.1@ulaval.ca  -  https://github.com/Jean-Romain/lidR
#
# COPYRIGHT:
#
# Copyright 2017-2018 Jean-Romain Roussel
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

#' An S4 class to represent a catalog of .las or .laz files
#'
#' A \code{LAScatalog} object is a representation of a set of las/laz files. A \code{LAScatalog} is
#' a way to manage and process an entire dataset. It allows the user to process a large area, or to
#' selectively clip data from a large area without loading all the data into computer memory.
#' A \code{LAScatalog} can be built with the function \link{readLAScatalog} and is formally an extension of
#' a \code{SpatialPolygonsDataFrame} containing extra data to allow users greater control over
#' how the dataset is processed (see details).
#'
#' A \code{LAScatalog} is formally a \code{SpatialPolygonsDataFrame} extended with new slots that
#' contain processing options. In \code{lidR}, each function that supports a \code{LAScatalog} as
#' input will respect these processing options. Internally, processing a catalog is almost always the
#' same and relies on just a few steps:
#' \enumerate{
#' \item Define chunks. A chunk is an arbitrarily-defined region of interest (ROI) of the
#' catalog. Altogether, the chunks are a wall-to-wall set of ROIs that encompass the whole dataset.
#' \item Loop over each chunk (in parallel or not).
#' \item For each chunk, load the points inside the ROI into R, run some R functions,
#' return the expected output.
#' \item Merge the outputs of the different chunks once they are all processed to build a continuous
#' (wall-to-wall) output.
#' }
#' So basically, a \code{LAScatalog} is an object that allows for batch processing but with the specificity
#' that \code{lidR} does not loop through las files, but loops seamlessly through chunks that do not
#' necessarily match with the file pattern. This way \code{lidR} can sequentially process tiny ROIs even if
#' each file may be individually too big to fit in memory. This is also why point cloud indexation
#' with \code{lax} files may significantly speed-up the processing.\cr\cr
#' It is important to note that catalogs with files that overlap each other are not natively supported
#' by \code{lidR}. When encountering such datasets the user should always filter any overlaps if
#' possible. This is possible if the overlapping points are flagged, for example in the
#' 'withheld' attribute. Otherwise \code{lidR} will not be able to process the dataset correctly.
#'
#' @slot processing_options list. A list that contains some settings describing how the catalog will be
#' processed (see dedicated section).
#'
#' @slot chunk_options list. A list that contains some settings describing how the catalog will be
#' sub-divided into chunks to be processed (see dedicated section).
#'
#' @slot output_options list. A list that contains some settings describing how the catalog will return
#' the outputs (see dedicated section).
#'
#' @slot input_options list. A list of parameters to pass to \link{readLAS} (see dedicated section).
#'
#' @section Processing options:
#' The slot \code{@processing_options} contains a \code{list} of options that determine how chunks
#' (the sub-areas that are sequentially processed) are processed.
#' \itemize{
#' \item \strong{progress}: boolean. Display a progress bar and a chart of progress. Default is TRUE.
#' Progress estimation can be enhanced by installing the package \code{progress}. See \link{opt_progress}.
#' \item \strong{stop_early}: boolean. Stop the processing if an error occurs in a chunk. If \code{FALSE}
#' the process can run until the end removing chunks that failed. Default is TRUE and the user should
#' have no reason to change this. See \link{opt_stop_early}.
#' \item \strong{wall.to.wall} logical. The catalog processing engine always guarantees to return a
#' continuous output without edge effects, assuming that the catalog is a wall-to-wall catalog. To do
#' so, some options are checked internally to guard against bad settings, such as buffer = 0 for an
#' algorithm that requires a buffer. In rare cases it might be useful to disable these controls. If
#' \code{wall.to.wall = FALSE} controls are disabled and wall-to-wall outputs cannot be guaranteed.
#' See \link{opt_wall_to_wall}
#' }
#'
#' @section Chunk options:
#' The slot \code{@clustering_options} contains a \code{list} of options that determine how chunks
#' (the sub-areas that are sequentially processed) are made.
#' \itemize{
#' \item \strong{chunk_size}: numeric. The size of the chunks that will be sequentially processed.
#' A small size allows small amounts of data to be loaded at once, saving computer memory, and vice versa.
#' The computation is usually faster but uses much more memory. If \code{chunk_size = 0} the
#' catalog is processed sequentially \emph{by file} i.e. a chunk is a file. Default is 0 i.e. by default
#' the processing engine respects the existing tiling pattern. See \link{opt_chunk_size}.
#' \item \strong{buffer}: numeric. Each chunk can be read with an extra buffer around it to ensure there is
#' no edge effect between two independent chunks and that the output is continuous. This is mandatory for
#' some algorithms. Default is 30. See \link{opt_chunk_buffer}.
#' \item \strong{alignment}: numeric. A vector of size 2 (x and y coordinates, respectively) to align the
#' chunk pattern. By default the alignment is made along (0,0), meaning that the edge of the first chunk
#' will belong on x = 0 and y = 0 and all the the other chunks will be multiples of the chunk size.
#' Not relevant if \code{chunk_size = 0}. See \link{opt_chunk_alignment}.
#' }
#'
#' @section Output options:
#' The slot \code{@output_options} contains a \code{list} of options that determine how clusters
#' (the sub-areas that are sequentially processed) are written. By "written" we mean written to files
#' or written in R memory.
#'
#' \itemize{
#' \item \strong{output_files}: string. If \code{output_files = ""} outputs are returned in R. Otherwise, if
#' \code{output_files} is a string the outputs will be written to files.
#' This is useful if the output is too big to be returned in R. A path to a templated filename
#' without extension (the engine guesses it for you) is expected. When several files are going to be
#' written a single string is provided with a template that is automatically filled. For example,
#' the following file names are possible:
#' \preformatted{
#' "/home/user/als/normalized/file_{ID}_segmented"
#' "C:/user/document/als/zone52_{XLEFT}_{YBOTTOM}_confidential"
#' "C:/user/document/als/{ORIGINALFILNAME}_normalized"
#' }
#' This option will generate as many filenames as needed with custom names for each file. The list of
#' allowed templates is described in the documentation for each function. See \link{opt_output_files}.
#' \item \strong{drivers}: list. This contains all the drivers required to seamlessly write Raster*,
#' Spatial*, LAS objects. It is recommended that only advanced users change this option. A dedicated
#' page describes the drivers in \link{lidR-LAScatalog-drivers}.
#' \item \strong{merge}: boolean. Multiple objects are merged into a single one at the end of the processing.
#' }
#'
#' @section Input options:
#' The slot \code{@input_options} contains a \code{list} of options that are passed to the function
#' \link{readLAS}. Indeed, the \code{readLAS} function is not called directly by the user but by the
#' internal processing engine. Users can propagate these options through the \code{LAScatalog} settings.
#'
#' \itemize{
#' \item \strong{select}: string. The \code{select} option. Usually this option is not respected because
#' each function knows which data must be loaded or not. This is documented in each function. See
#' \link{opt_select}.
#' \item \strong{filter}: string. The \code{filter} option. See \link{opt_filter}.
#' }
#'
#' @import data.table
#' @import methods
#' @importClassesFrom sp CRS
#' @importClassesFrom sp SpatialPolygonsDataFrame
#'
#' @export
#'
#' @useDynLib lidR, .registration = TRUE
#'
#' @examples
#' \dontrun{
#' # Build a catalog
#' ctg <- readLAScatalog("filder/to/las/files/")
#'
#' # Set some options
#' opt_filter(ctg) <- "-keep_first"
#'
#' # Summary gives a summary of how the catalog will be processed
#' summary(ctg)
#'
#' # We can seamlessly use lidR functions
#' hmean <- grid_metrics(ctg, mean(Z), 20)
#' ttops <- tree_detection(ctg, lmf(5))
#'
#' # For low memory config it is probably advisable not to load entire files
#' # and process chunks instead
#' opt_chunk_size(ctg) <- 500
#'
#' # Sometimes the output is likely to be very large
#' # e.g. large coverage and small resolution
#' dtm <- grid_terrain(ctg, 1, tin())
#'
#' # In that case it is advisable to write the output(s) to files
#' opt_output_files(ctg) <- "path/to/folder/DTM_chunk_{XLEFT}_{YBOTTOM}"
#'
#' # Raster will be written to disk. The list of written files is returned
#' # or, in this specific case, a virtual raster mosaic.
#' dtm <- grid_terrain(ctg, 1, tin())
#'
#' # When chunks are files the original names of the las files can be preserved
#' opt_chunk_size(ctg) <- 0
#' opt_output_files(ctg) <- "path/to/folder/DTM_{ORIGINALFILENAME}"
#' dtm <- grid_terrain(ctg, 1, tin())
#'
#' # For some functions, files MUST be written to disk. Indeed, it is certain that R cannot
#' # handle the entire output.
#' opt_chunk_size(ctg) <- 0
#' opt_output_files(ctg) <- "path/to/folder/{ORIGINALFILENAME}_norm"
#' opt_laz_compression(ctg) <- TRUE
#' new_ctg <- normalize_height(ctg, tin())
#'
#' # The user has access to the catalog engine through the function catalog_apply
#' output <- catalog_apply(ctg, FUN, ...)
#' }
setClass(
  Class = "LAScatalog",
  contains = "SpatialPolygonsDataFrame",
  representation(
    chunk_options = "list",
    processing_options = "list",
    output_options = "list",
    input_options = "list"
  )
)

setMethod("initialize", "LAScatalog", function(.Object)
{
  callNextMethod()

  drivers = list(
    Raster = list(
      write = raster::writeRaster,
      extension = ".tif",
      object = "x",
      path = "filename",
      param = list(format = "GTiff", NAflag = -999999)
    ),
    LAS = list(
      write = lidR::writeLAS,
      extension = ".las",
      object = "las",
      path = "file",
      param = list()
    ),
    Spatial = list(
      write = writeSpatial,
      extension = ".shp",
      object = "x",
      path = "filename",
      param = list(overwrite = FALSE)
    ),
    SimpleFeature = list(
      write = sf::st_write,
      extension = ".shp",
      object = "obj",
      path = "dsn",
      param = list(quiet = TRUE)
    ),
    DataFrame = list(
      write = data.table::fwrite,
      extension = ".txt",
      object = "x",
      path = "file",
      param = list()
    )
  )

  .Object@chunk_options <- list(
    size = 0,
    buffer = 30,
    alignment = c(0,0)
  )

  .Object@processing_options <- list(
    progress = TRUE,
    stop_early = TRUE,
    wall_to_wall = TRUE
  )

  .Object@output_options <- list(
    output_files = "",
    drivers = drivers,
    merge = TRUE
  )

  .Object@input_options <- list(
    select = "*",
    filter = "",
    alt_dir = ""
  )

  return(.Object)
})
