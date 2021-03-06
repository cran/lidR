% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cloud_metrics.R
\name{cloud_metrics}
\alias{cloud_metrics}
\title{Compute metrics for a cloud of points}
\usage{
cloud_metrics(las, func)
}
\arguments{
\item{las}{An object of class \code{LAS}}

\item{func}{formula. An expression to be applied to the point cloud (see example)}
}
\value{
It returns a \code{list} containing the metrics
}
\description{
\code{cloud_metrics} computes a series of user-defined descriptive statistics for a LiDAR dataset.
See \link[=grid_metrics]{grid_metrics} to compute metrics on a grid. Basically there are
no predefined metrics. Users must write their own functions to create metrics (see example).
The following existing functions can serve as a guide to help users compute their own metrics:
\itemize{
\item{\link[=stdmetrics]{stdmetrics}}
\item{\link[=entropy]{entropy}}
\item{\link[=VCI]{VCI}}
\item{\link[=LAD]{LAD}}
}
}
\examples{
LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
lidar = readLAS(LASfile)

cloud_metrics(lidar, ~max(Z))
cloud_metrics(lidar, ~mean(Intensity))

# Define your own new metrics
myMetrics = function(z, i)
{
  metrics = list(
     zwimean = sum(z*i)/sum(i), # Mean elevation weighted by intensities
     zimean  = mean(z*i),       # Mean products of z by intensity
     zsqmean = sqrt(mean(z^2))  # Quadratic mean
   )

   return(metrics)
}

metrics = cloud_metrics(lidar, ~myMetrics(Z, Intensity))

# Predefined metrics
cloud_metrics(lidar, .stdmetrics)
}
\seealso{
\link[=grid_metrics]{grid_metrics}
\link[=stdmetrics]{stdmetrics}
\link[=entropy]{entropy}
\link[=VCI]{VCI}
\link[=LAD]{LAD}

Other metrics: 
\code{\link{grid_metrics}()},
\code{\link{hexbin_metrics}()},
\code{\link{point_metrics}()},
\code{\link{tree_metrics}()},
\code{\link{voxel_metrics}()}
}
\concept{metrics}
