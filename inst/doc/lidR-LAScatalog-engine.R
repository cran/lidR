## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 2.5,
  fig.height = 2.5,
  dev.args = list(pointsize = 9)
)
knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      # record the current time before each chunk
      now <<- Sys.time()
    } else {
      # calculate the time difference after a chunk
      res <- difftime(Sys.time(), now)
      # return a character string to show the time
      #if (res > 0.1)
      #paste("<br/>========================<br/>Time for this code chunk ", options$label, " to run:", round(res,2), "<br/>========================<br/>")
    }
  }
}))
knitr::opts_chunk$set(time_it = TRUE)
#rgl::setupKnitr()
rgdal::set_thin_PROJ6_warnings(TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)

library(lidR)

## ----data, echo = FALSE-------------------------------------------------------
data = structure(list(Max.X = c(332099.99, 333600, 335099.99, 336217.52, 
332099.99, 333599.99, 335099.99, 336368.67, 332099.99, 333599.99, 
335100, 336217.52), Min.X = c(331016.91, 332100.01, 333600.01, 
335100, 331016.91, 332100, 333600, 335100, 331016.92, 332100.01, 
333600.01, 335100.01), Max.Y = c(5529993.99, 5529993.99, 5529993.99, 
5529993.99, 5528399.99, 5528399.99, 5528399.99, 5528399.99, 5526399.98, 
5526399.96, 5526399.99, 5526399.99), Min.Y = c(5528400, 5528400, 
5528400, 5528400, 5526400, 5526400, 5526400, 5526400, 5524793.5, 
5524793.5, 5524800.38, 5524793.5), Max.Z = c(53.53, 47.59, 48.66, 
49.36, 46.13, 48.16, 50.51, 50.86, 45, 74.18, 52.56, 49.33), 
    Min.Z = c(-15.95, -7.87, -3.55, -14.96, -5.94, -11.15, -5.11, 
    -4.12, -9.63, -8.27, -35.88, -20.59), filename = c("folder/file1.las", 
    "folder/file2.las", "folder/file3.las", "folder/file4.las", 
    "folder/file5.las", "folder/file6.las", "folder/file7.las", 
    "folder/file8.las", "folder/file9.las", "folder/file10.las", 
    "folder/file11.las", "folder/file12.las")), row.names = c(NA, 
-12L), class = "data.frame")

pgeom <- lapply(1:nrow(data), function(i)
{
  mtx <- matrix(c(data$Min.X[i], data$Max.X[i], data$Min.Y[i], data$Max.Y[i])[c(1, 1, 2, 2, 1, 3, 4, 4, 3, 3)], ncol = 2)
  sp::Polygons(list(sp::Polygon(mtx)),as.character(i))
})

Sr = sp::SpatialPolygons(pgeom, proj4string = sp::CRS("+init=epsg:3005"))

ctg <- new("LAScatalog")
ctg@bbox <- Sr@bbox
ctg@proj4string <- Sr@proj4string
ctg@plotOrder <- Sr@plotOrder
ctg@data <- data
ctg@polygons <- Sr@polygons

## ----setbuffer2, echo = FALSE-------------------------------------------------
opt_chunk_buffer(ctg) <- 0

## ----plotctg, fig.show='hold'-------------------------------------------------
opt_chunk_size(ctg) <- 0 # Processing by files
plot(ctg, chunk = TRUE)

opt_chunk_size(ctg) <- 1000 # Processing chunks of 1000 x 1000
plot(ctg, chunk = TRUE)

## ----setbuffer1, echo = FALSE-------------------------------------------------
opt_chunk_size(ctg) <- 0

## ----plotbuffer, fig.show='hold'----------------------------------------------
opt_chunk_buffer(ctg) <- 0 # No buffer
plot(ctg, chunk = TRUE)

opt_chunk_buffer(ctg) <- 200 # 200 m buffer
plot(ctg, chunk = TRUE)

## ----dtmnobuffer, error=TRUE--------------------------------------------------
opt_chunk_buffer(ctg) <- 0
grid_terrain(ctg, 1, tin())

## ----alignment, fig.show='hold'-----------------------------------------------
opt_chunk_size(ctg) <- 2000
opt_chunk_buffer(ctg) <- 0
plot(ctg, chunk = TRUE)

opt_chunk_size(ctg) <- 2000
opt_chunk_buffer(ctg) <- 0
opt_chunk_alignment(ctg) <- c(1000, 1000)
plot(ctg, chunk = TRUE)

## ----void, echo = FALSE, rgl=TRUE, dev='png'----------------------------------
#LASfile <- system.file("extdata", "Topography.laz", package="lidR")
#ctg = readLAScatalog(LASfile)
#opt_progress(ctg) <- FALSE
#opt_filter(ctg) <- "-keep_class 2 9"
#las = clip_circle(ctg, 273500, 5274500, 40)
#m = structure(c(0.921, -0.146, 0.362, 0, 0.386, 0.482, -0.787, 0, 
#-0.06, 0.864, 0.5, 0, 0, 0, 0, 1), .Dim = c(4L, 4L))
#plot(las)
#rgl::rgl.viewpoint(fov = 50, userMatrix = m)

## ----writeondisk, echo = FALSE, eval = FALSE----------------------------------
#  LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#  ctg2 <- readLAScatalog(LASfile)
#  opt_progress(ctg2) <- FALSE
#  opt_chunk_size(ctg2) <- 100

## ----template, eval = FALSE---------------------------------------------------
#  # Force the results to be written on disk
#  opt_output_files(ctg2) <- paste0(tempdir(), "/tree_coordinate_{XLEFT}_{YBOTTOM}")
#  trees <- find_trees(ctg2, lmf(3))
#  
#  # The output has been modified by these options and it now gives
#  # the paths to the written files (here shapefiles)
#  trees
#  #> "/tmp/RtmpJQHPNz/tree_coordinate_481200_3812900.shp" "/tmp/RtmpJQHPNz/tree_coordinate_481300_3812900.shp" "/tmp/RtmpJQHPNz/tree_coordinate_481200_3813000.shp"
#  #> [4] "/tmp/RtmpJQHPNz/tree_coordinate_481300_3813000.shp"

## ----writechm, eval = FALSE---------------------------------------------------
#  # Force the results to be written on disk
#  opt_output_files(ctg2) <- paste0(tempdir(), "/tree_coordinate_{XLEFT}_{YBOTTOM}")
#  chm <- grid_canopy(ctg2, 1, p2r())
#  
#  # Many rasters have been written on disk
#  # but a light RasterLayer has been returned anyway
#  chm
#  #> class      : RasterLayer
#  #> dimensions : 90, 90, 8100  (nrow, ncol, ncell)
#  #> resolution : 1, 1  (x, y)
#  #> extent     : 481260, 481350, 3812921, 3813011  (xmin, xmax, ymin, ymax)
#  #> crs        : +proj=utm +zone=12 +datum=NAD83 +units=m +no_defs
#  #> source     : /tmp/RtmpZVJ2hy/grid_canopy.vrt
#  #> names      : tree_coordinate_481260_3812921
#  #> values     : 0, 32.07  (min, max)

## ----clip, fig.show='hold', eval=FALSE----------------------------------------
#  opt_output_files(ctg2) <- "{tempdir()}/plot_{ID}"
#  new_ctg <- clip_circle(ctg2, x, y, 20)
#  new_ctg
#  #> class       : LAScatalog (v1.2 format 0)
#  #> extent      : 32.372, 163.136, 38.494, 198.636 (xmin, xmax, ymin, ymax)
#  #> coord. ref. : NAD83 / UTM zone 17N
#  #> area        : 3895.031 m²
#  #> points      : 44 points
#  #> density     : 8 points/m²
#  #> num. files  : 4

## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:5){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[6]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "cornflowerblue")


## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:6){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[7]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "orange")

bbox <- cl[[8]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "cornflowerblue")

## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:8){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[7]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "orange")

bbox <- cl[[9]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "red")


## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
opt_restart(ctg) <- 9
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:4){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
opt_restart(ctg) <- 1
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:length(cl)){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[7]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "orange")

bbox <- cl[[9]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "red")


## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 400
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:50){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[1]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[2]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[3]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[14]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[15]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[16]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[29]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:6) {
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

for (i in 7:11) {
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "cornflowerblue")
}

## -----------------------------------------------------------------------------
ctg$processed <- FALSE
ctg$processed[6:7] <- TRUE
plot(ctg)

## ---- echo = FALSE------------------------------------------------------------
opt_wall_to_wall(ctg) <- TRUE
opt_progress(ctg) <- FALSE

## ---- error = TRUE------------------------------------------------------------
routine <- function(chunk){ 
  las <- readLAS(chunk)
}

catalog_apply(ctg, routine)

## ----getachunk, eval=FALSE,echo=FALSE,warning=FALSE,message=FALSE,error=FALSE,results='hide',fig.keep='none'----
#  LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#  test = readLAScatalog(LASfile)
#  
#  opt_chunk_size(test) <- 150
#  opt_chunk_alignment(test) <- c(50,10)
#  opt_progress(ctg) <- FALSE
#  chunks = lidR:::catalog_makecluster(test)
#  chunk = chunks[[5]]

## ----rglbuffer, rgl = TRUE, eval = FALSE--------------------------------------
#  las <- readLAS(chunk)
#  plot(las, color = "buffer")

## ---- eval = FALSE------------------------------------------------------------
#  print(chunk)
#  #> class       : LAScluster
#  #> features    : 1
#  #> extent      : 684800, 684950, 5017810, 5017960  (xmin, xmax, ymin, ymax)
#  #> crs         : +proj=utm +zone=17 +datum=NAD83 +units=m +no_defs

## ---- warning = FALSE, eval = FALSE-------------------------------------------
#  raster::extent(chunk)
#  #> class      : Extent
#  #> xmin       : 684800
#  #> xmax       : 684950
#  #> ymin       : 5017810
#  #> ymax       : 5017960
#  sp::bbox(chunk)
#  #>         [,1]    [,2]
#  #> [1,]  684800  684950
#  #> [2,] 5017810 5017960
#  sf::st_bbox(chunk)
#  #>    xmin    ymin    xmax    ymax
#  #>  684800 5017810  684950 5017960

## ----bufferror, error = TRUE--------------------------------------------------
opt_chunk_buffer(ctg) <- 0
grid_terrain(ctg, 1, tin())

## ----routineerror, error = TRUE-----------------------------------------------
routine <- function(chunk){ 
   las <- readLAS(chunk)
   if (is.empty(las)) return(NULL)
}

options = list(need_buffer = TRUE)
catalog_apply(ctg, routine, .options = options)

## ----preparectg, echo=FALSE,warning=FALSE,message=FALSE,error=FALSE,results='hide',fig.keep='none'----
LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
ctg = readLAScatalog(LASfile)

opt_chunk_buffer(ctg) <- 10
opt_chunk_size(ctg) <- 100
opt_chunk_alignment(ctg) <- c(50,50)
opt_progress(ctg) <- FALSE

## ----applyroutine, eval = FALSE-----------------------------------------------
#  routine <- function(chunk){
#     las <- readLAS(chunk)               # read the chunk
#     if (is.empty(las)) return(NULL)     # exit if empty
#     ttop <- find_trees(las, lmf(3))     # make any computation
#     ttop <- raster::crop(ttop, extent(chunk))   # remove the buffer
#     return(ttop)
#  }
#  
#  out <- catalog_apply(ctg, routine)
#  class(out)
#  #> [1] "list"
#  print(out[[1]])
#  #> class       : SpatialPointsDataFrame
#  #> features    : 94
#  #> extent      : 481260.1, 481349.8, 3812921, 3812950  (xmin, xmax, ymin, ymax)
#  #> crs         : +proj=utm +zone=12 +datum=NAD83 +units=m +no_defs
#  #> variables   : 2
#  #> names       : treeID,     Z
#  #> min values  :      1,  2.16
#  #> max values  :    132, 32.07

## ---- eval = FALSE------------------------------------------------------------
#  out <- do.call(rbind, out)
#  print(out)
#  #> class       : SpatialPointsDataFrame
#  #> features    : 297
#  #> extent      : 481260, 481349.9, 3812921, 3813011  (xmin, xmax, ymin, ymax)
#  #> crs         : +proj=utm +zone=12 +datum=NAD83 +units=m +no_defs
#  #> variables   : 2
#  #> names       : treeID,     Z
#  #> min values  :      1,  2.16
#  #> max values  :    247, 32.07

## ----automerge, eval = FALSE--------------------------------------------------
#  options <- list(automerge = TRUE)
#  out <- catalog_apply(ctg, routine, .options = options)
#  print(out)
#  #> class       : SpatialPointsDataFrame
#  #> features    : 297
#  #> extent      : 481260, 481349.9, 3812921, 3813011  (xmin, xmax, ymin, ymax)
#  #> crs         : +proj=utm +zone=12 +datum=NAD83 +units=m +no_defs
#  #> variables   : 2
#  #> names       : treeID,     Z
#  #> min values  :      1,  2.16
#  #> max values  :    247, 32.07

