## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
rgl::setupKnitr()
library(lidR)

## ---- echo = FALSE------------------------------------------------------------
data = data.table::data.table(
  Max.X   = c(885228.88, 886993.96, 885260.93, 887025.96,
              885292.94, 887056.88, 892199.94, 893265.54, 892229.99, 893295.15,
              888759.96, 890524.95, 892259.98, 894025.98, 892289.96, 894055.93,
              888790.91, 890554.98, 888820.95, 890585.99, 892319.96, 894084.97,
              892349.89, 894114.29, 895250.23, 895094.78, 895044.96, 895053.55,
              885323.96, 887087.95, 885355.95, 887119.96, 883657.85, 885387.95,
              887150.97, 885419.98, 887182.95, 883688.44, 885442.91, 887193.9,
              888851.96, 890615.97, 888882.94, 890646.97, 892379.94, 894127.84,
              892409.97, 892676.58, 888913.92, 890676.93, 888944.86, 890707.98,
              892439.95, 894124.59, 892469.86, 894232.94, 894786.68, 888958.83,
              890713.51, 892476.43, 894239.97, 894786.07),
  Min.X   = c(885022.37,
              885204.73, 885027.52, 885229.03, 885040.86, 885261.03, 891503.09,
              892198.69, 891501.42, 892200.07, 886970.07, 888735.55, 891499.96,
              892230.05, 890521.99, 892260.01, 886994.05, 888760.09, 887026.07,
              888791.01, 890525.05, 892290.04, 890555.01, 892320.12, 894002.98,
              894026.02, 894056.02, 894085.03, 885051.45, 885293.03, 885063.29,
              885324.03, 883166.09, 885072.16, 885356.09, 883642.36, 885388.15,
              883180.23, 883658.11, 885420.02, 887057.07, 888821.02, 887088.11,
              888852.03, 890586.03, 892350.02, 890616.07, 892380.01, 887120.07,
              888883.03, 887151.11, 888914.02, 890647.06, 892410.06, 890677.07,
              892440.07, 894209.19, 887183.07, 888945.12, 890708.03, 892470.16,
              894233.07),
  Max.Y   = c(630219.48, 630214.96, 631609.95, 631604.97,
              633001.65, 632995.99, 625898.35, 625882.94, 627289.82, 627273.89,
              630174.88, 630134.94, 628681.66, 628664.99, 630094.95, 630057.95,
              631564.98, 631524.94, 632955.82, 632915.99, 631486.9, 631447.96,
              632876.93, 632838.96, 628627.89, 630019.93, 631410.97, 631740.88,
              634393.05, 634386.96, 635786.24, 635779.96, 638613.36, 637176.84,
              637169.92, 638601.99, 638560.96, 639938.36, 639926.95, 639558.31,
              634346.93, 634307.92, 635739.92, 635699.92, 634268.97, 634229.95,
              635659.89, 635622.88, 637129.84, 637089.93, 638520.94, 638481.91,
              637051.99, 637012.92, 638442.98, 638403.94, 638366.87, 639177.04,
              639133.74, 638702.56, 638702.56, 638702.56),
  Min.Y   = c(629157.18,
              629099.31, 630215.04, 630175.05, 631605.02, 631565.05, 625816.52,
              625793.6, 625883.01, 625860.81, 629036.82, 629017.72, 627274.01,
              627251.36, 628665.04, 628628.01, 630135.08, 630095.02, 631525.01,
              631487.19, 630058.02, 630020.05, 631448.08, 631411.03, 627506.32,
              628612.41, 629999.84, 631390.38, 632996.06, 632956.04, 634387.01,
              634347.01, 637939.24, 635780.07, 635740.05, 637170.11, 637130.14,
              638602.13, 638561.04, 638521.07, 632916.05, 632877.04, 634308.06,
              634269.04, 632839.06, 632801.04, 634230.04, 634223.9, 635700.07,
              635660.11, 637090.03, 637052.15, 635623.06, 635619.13, 637013.1,
              636979.71, 637259.33, 638482.01, 638443.02, 638404.08, 638367.11,
              638355.37),
  Min.Z   = c(325.12, 251.48, 244.68, 286.7, 338.86, 320.68, 118.08, 60.69, 
              -5.01, -7.58, 225.29, 252, 87.3, 41.7, 115.01, 28.77, 205.11, 
              200.85, 200.54, 169.5, 90.64, 19.72, 126.04, 28.6, 41.98, 43.15, 
              7.74, 6.7, 199.26, 190.02, 284.92, 216.16, 218.14, 318.93, 220.21, 
              218.04, 137.31, 218.13, 217.34, 190.42, 207.62, 113.84, 118.18, 
              141.52, 92, 52.5, 91.2, 77.92, 156.57, 89.53, 108.83, 93.98, 
              23.45, 1.64, 33.22, 3.29, 0.61, 108.03, 208.38, 121.18, 58.83, 
              0.95),
  Max.Z   = c(418.46, 990.54, 409.06, 1021.87, 996.42, 1005.02, 173.77, 393.97, 
              836.52, 820.98, 414.2, 936.47, 792.95, 822.51, 777.31, 837.87, 
              419.15, 741.84, 907.2, 872.27, 898.53, 822.53, 846.77, 740.65, 
              826.61, 890.21, 828.86, 680.32, 390.31, 997.2, 965.55, 969.24, 
              249.34, 849.5, 950.2, 848.64, 904.1, 880, 827.92, 888.34, 462.88, 
              906.61, 440.83, 887.34, 860.37, 747.1, 808.75, 194.76, 734.21, 
              838.34, 834.76, 758.91, 771.76, 670.1, 810.94, 761.53, 109.26, 
              303.94, 349.94, 799.8, 737.01, 593.91),
  filename = paste0("path/to/las/files/file", 1:62, ".las")
)


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

## ---- echo = FALSE------------------------------------------------------------
opt_chunk_buffer(ctg) <- 0

## ---- fig.show='hold'---------------------------------------------------------
opt_chunk_size(ctg) <- 0 # Processing by files
plot(ctg, chunk = TRUE)

opt_chunk_size(ctg) <- 900 # Processing chunks of 900 x 900
plot(ctg, chunk = TRUE)

## ---- echo = FALSE------------------------------------------------------------
opt_chunk_size(ctg) <- 0

## ---- fig.show='hold'---------------------------------------------------------
opt_chunk_buffer(ctg) <- 0 # No buffer
plot(ctg, chunk = TRUE)

opt_chunk_buffer(ctg) <- 200 # 200 m buffer
plot(ctg, chunk = TRUE)

## ---- error=TRUE--------------------------------------------------------------
opt_chunk_buffer(ctg) <- 0
grid_terrain(ctg, 1, tin())

## ---- fig.show='hold'---------------------------------------------------------
opt_chunk_size(ctg) <- 2000
opt_chunk_buffer(ctg) <- 0
plot(ctg, chunk = TRUE)

opt_chunk_size(ctg) <- 2000
opt_chunk_buffer(ctg) <- 0
opt_chunk_alignment(ctg) <- c(1000, 1000)
plot(ctg, chunk = TRUE)

## ----echo = FALSE, rgl=TRUE, dev='png'----------------------------------------
#LASfile <- system.file("extdata", "Topography.laz", package="lidR")
#ctg = readLAScatalog(LASfile)
#opt_progress(ctg) <- FALSE
#opt_filter(ctg) <- "-keep_class 2 9"
#las = clip_circle(ctg, 273500, 5274500, 40)
#m = structure(c(0.921, -0.146, 0.362, 0, 0.386, 0.482, -0.787, 0, 
#-0.06, 0.864, 0.5, 0, 0, 0, 0, 1), .Dim = c(4L, 4L))
#plot(las)
#rgl::rgl.viewpoint(fov = 50, userMatrix = m)

## ----echo = FALSE, rgl=TRUE, dev='png'----------------------------------------
LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
ctg2 <- readLAScatalog(LASfile)
opt_progress(ctg2) <- FALSE
opt_chunk_size(ctg2) <- 250

## -----------------------------------------------------------------------------
# Find the trees
trees <- find_trees(ctg2, lmf(3))

# The trees are immediately usable in subsequent analyses
trees

## -----------------------------------------------------------------------------
# Force the results to be written on disk
opt_output_files(ctg2) <- paste0(tempdir(), "/tree_coordinate_{XLEFT}_{YBOTTOM}")
trees <- find_trees(ctg2, lmf(3))

# The output has been modified by these options and it now gives
# the paths to the written files (here shapefiles)
trees

## -----------------------------------------------------------------------------
# Force the results to be written on disk
opt_output_files(ctg2) <- paste0(tempdir(), "/tree_coordinate_{XLEFT}_{YBOTTOM}")
chm <- grid_canopy(ctg2, 1, p2r())

# Many rasters have been written on disk
# but a light RasterLayer has been returned anyway
chm

## ---- echo=FALSE--------------------------------------------------------------
ctg2 = lidR:::catalog_generator(600)
opt_progress(ctg2) <- FALSE
x = c(50, 90, 150, 140)
y = c(55, 75, 80, 180)

## ---- fig.show='hold'---------------------------------------------------------
opt_output_files(ctg2) <- "{tempdir()}/plot_{ID}"
new_ctg <- clip_circle(ctg2, x, y, 20)
plot(ctg2)
points(x,y)
plot(new_ctg)

## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:18){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[19]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "cornflowerblue")


## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:22){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[23]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "orange")

bbox <- cl[[24]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "cornflowerblue")

## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:29){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[23]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "orange")

bbox <- cl[[30]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "red")


## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:length(cl)){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[23]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "orange")

bbox <- cl[[30]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "red")


## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 500
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:190){
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

bbox <- cl[[75]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[47]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[68]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[89]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")

bbox <- cl[[152]]@bbox
graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "gray")


## ---- echo=FALSE--------------------------------------------------------------
opt_chunk_size(ctg) <- 0
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- FALSE
opt_progress(ctg) <- TRUE
cl <- lidR:::catalog_makecluster(ctg)
for (i in 1:18) {
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "green3")
}

for (i in 19:22) {
  bbox <- cl[[i]]@bbox
  graphics::rect(bbox[1], bbox[2], bbox[3], bbox[4], border = "black", col = "cornflowerblue")
}

## -----------------------------------------------------------------------------
ctg$processed <- FALSE
ctg$processed[41:44] <- TRUE
plot(ctg)

## ---- echo = FALSE------------------------------------------------------------
opt_wall_to_wall(ctg) <- TRUE
opt_progress(ctg) <- FALSE

## ---- error = TRUE------------------------------------------------------------
routine <- function(chunk){ 
  las <- readLAS(chunk)
}

catalog_apply(ctg, routine)

## ---- echo=FALSE,warning=FALSE,message=FALSE,error=FALSE,results='hide',fig.keep='none'----
LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
test = readLAScatalog(LASfile)

opt_chunk_size(test) <- 150
opt_chunk_alignment(test) <- c(50,10)
opt_progress(ctg) <- FALSE
chunks = lidR:::catalog_makecluster(test)
chunk = chunks[[5]]

## ---- rgl = TRUE--------------------------------------------------------------
las <- readLAS(chunk)
plot(las, color = "buffer")

## -----------------------------------------------------------------------------
print(chunk)

## -----------------------------------------------------------------------------
extent(chunk)
bbox(chunk)

## ---- error = TRUE------------------------------------------------------------
opt_chunk_buffer(ctg) <- 0
grid_terrain(ctg, 1, tin())

## ---- error = TRUE------------------------------------------------------------
routine <- function(chunk){ 
   las <- readLAS(chunk)
   if (is.empty(las)) return(NULL)
}

options = list(need_buffer = TRUE)
catalog_apply(ctg, routine, .options = options)

## ---- echo=FALSE,warning=FALSE,message=FALSE,error=FALSE,results='hide',fig.keep='none'----
LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
ctg = readLAScatalog(LASfile)

opt_chunk_size(ctg) <- 200
opt_chunk_alignment(ctg) <- c(50,50)
opt_progress(ctg) <- FALSE

## -----------------------------------------------------------------------------
routine <- function(chunk){ 
   las <- readLAS(chunk)               # read the chunk
   if (is.empty(las)) return(NULL)     # exit if empty
   ttop <- find_trees(las, lmf(3))     # make any computation
   ttop <- crop(ttop, extent(chunk))   # remove the buffer
   return(ttop)
}

out <- catalog_apply(ctg, routine)
class(out)
print(out[[1]])

## -----------------------------------------------------------------------------
out <- do.call(rbind, out)
print(out)

## -----------------------------------------------------------------------------
options <- list(automerge = TRUE)
out <- catalog_apply(ctg, routine, .options = options)
class(out)
print(out)

## -----------------------------------------------------------------------------
find_deadtrees <- function(las, param1, param2)
{
   if (is(las, "LAScatalog"))  {
      options <- list(automerge = TRUE, need_buffer = TRUE)
      dead_trees <- catalog_apply(las, find_deadtrees, param1 = param1, param2 = param2, .options = options)
      return(dead_trees)
   }
   else if (is(las, "LAScluster")) {
      bbox <- extent(las)
      las <- readLAS(las)
      dead_trees <- find_deadtrees(las, param1, param2)
      dead_trees <- raster::crop(dead_trees, bbox)
      return(dead_trees)
   }
   else if (is(las, "LAS")) {
      # make an advanced computation
      return(dead_trees)
   }
   else {
      stop("This type is not supported.")
   }
}

