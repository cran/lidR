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
options(rmarkdown.html_vignette.check_title = FALSE)
library(lidR)

## ----echo = FALSE-------------------------------------------------------------
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

geom <- lapply(1:nrow(data), function(i)
{
  mtx <- matrix(c(data$Min.X[i], data$Max.X[i], data$Min.Y[i], data$Max.Y[i])[c(1, 1, 2, 2, 1, 3, 4, 4, 3, 3)], ncol = 2)
  sf::st_polygon(list(mtx))
})

geom <-sf::st_sfc(geom)
sf::st_crs(geom) <- 26917
data <- sf::st_set_geometry(data, geom)

ctg       <- new("LAScatalog")
ctg@data  <- data

## ----echo = FALSE-------------------------------------------------------------
bbox = as.numeric(ctg@data[42,1:4])

plot(ctg)
graphics::rect(bbox[1], bbox[3], bbox[2], bbox[4], border = "black", col = "blue")

## ----echo = FALSE-------------------------------------------------------------
bbox = as.numeric(ctg@data[42,1:4])

plot(ctg)
graphics::rect(bbox[1]+200, bbox[3]+200, bbox[2]-200, bbox[4]-200, border = "black", col = "red")
graphics::rect(bbox[2], bbox[4], bbox[1], bbox[3], border = "black", col = "blue")

## ----echo = FALSE-------------------------------------------------------------
neighbourg = ctg@data[c(19,20,23, 41, 43, 44, 45, 47),1:4]

plot(ctg)

for (i in 1:nrow(neighbourg))
{
  bbox = as.numeric(neighbourg[i,])
  graphics::rect(bbox[2], bbox[4], bbox[1], bbox[3], border = "black", col = "red")
}

## ----echo = FALSE,eval=FALSE--------------------------------------------------
# LASfile <- system.file("extdata", "Megaplot.laz", package = "lidR")
# las = readLAS(LASfile)
# f = tempfile(fileext = ".las")
# writeLAS(las, f)
# 
# 
# t1 = system.time(for (i in 1:10) readLAS(LASfile))
# t2 = system.time(for (i in 1:10) readLAS(f))
# 
# t1 = t1[3]
# t2 = t2[3]
# 
# t = c(t1,t2)
# format = c("laz", "las")
# X = data.frame(t, format)
# 
# op <- graphics::par(mar = c(4,4,1,1) + 0.1)
# barplot(X$t/min(X$t), names.arg = X$format, col = "darkred", xlab = "File format", ylab = "Relative read time", asp = 1)
# graphics::par(op)

## ----echo = FALSE-------------------------------------------------------------
X = structure(list(t = c(1.5, 0.75), format = c("laz", 
"las")), class = "data.frame", row.names = c(NA, -2L))
op <- graphics::par(mar = c(4,4,1,1) + 0.1)
barplot(X$t/min(X$t), names.arg = X$format, col = "darkred", xlab = "File format", ylab = "Relative read time", asp = 1)
graphics::par(op)

## ----echo = FALSE-------------------------------------------------------------
neighbourg = ctg@data[c(19,20,23, 41, 43, 44, 45, 47),1:4]

J = list(c(1200,1000,0,0),
         c(0,1000,0,0),
         c(0,1000,-1200,0),
         c(1200,0,0,0),
         c(1200,0,0,-1000),
         c(0,0,0,-1000),
         c(0,0,-1200,0),
         c(0,0,-1200,-1000))

plot(ctg)

for (i in 1:nrow(neighbourg))
{
  j = J[[i]]
  bbox = as.numeric(neighbourg[i,])
  graphics::rect(bbox[2] + j[1], bbox[4] + j[2], bbox[1] + j[3], bbox[3] + j[4], border = "black", col = "red")
}

## ----echo = FALSE-------------------------------------------------------------
# LASfile <- system.file("extdata", "Megaplot.laz", package = "lidR")
# las = readLAS(LASfile)
# f = tempfile(fileext = ".las")
# writeLAS(las, f)
# t1 <- t2 <- t3 <- t4 <- 0
# for(i in 1:3)
# {
#   tt1 = system.time(readLAS(f))
#   tt2 = system.time(readLAS(f, "xyz"))
#   tt3 = system.time(readLAS(LASfile))
#   tt4 = system.time(readLAS(LASfile, "xyz"))
#   
#   t1 = t1 + tt1[3]
#   t2 = t2 + tt2[3]
#   t3 = t3 + tt3[3]
#   t4 = t4 + tt4[3]
# }

t = c(0.087, 0.046, 0.182, 0.156)
format = c("las\n*", "las\nxyz", "laz\n*", "laz\nxyz")
X = data.frame(t, format)

op <- graphics::par(mar = c(4,4,1,1) + 0.1)
barplot(X$t/min(X$t), names.arg = X$format, col = "darkred", xlab = "Format and attribute selection", ylab = "Relative read time", asp = 1)
graphics::par(op)

## ----echo = FALSE-------------------------------------------------------------
Format = c("laz", "laz + lax", "las", "las + lax")
Runtime = c("40 sec", "20 sec", "10 sec", "7 sec")
Timing = data.frame(Format, Runtime)
knitr::kable(Timing)

## ----echo = FALSE-------------------------------------------------------------
Format = c("laz", "las", "las + lax")
Runtime = c("45 min", "15 min", "8 min")
Timing = data.frame(Format, Runtime)
knitr::kable(Timing)

## ----echo = FALSE-------------------------------------------------------------
Format = c("las", "las + lax")
Runtime = c("45 sec", "4 sec")
Timing = data.frame(Format, Runtime)
knitr::kable(Timing)

