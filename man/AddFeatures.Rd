% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_shapeprep.R
\name{AddFeatures}
\alias{AddFeatures}
\title{Add nodes or edges to a sfnetwork}
\usage{
AddFeatures(network, nodes, edges, crs, tolerance = 200, tolerance2 = 10)
}
\arguments{
\item{network}{a sfnetwork object}

\item{nodes}{a data frame with lon and lat columns, or sf object with POINT geometry}

\item{edges}{a sf object with LINESTRING geometry that includes edges to be added}

\item{crs}{a crs object, if nodes or edges has a different crs from network}

\item{tolerance}{a maximum distance for snapping a shape to the network, shapes outside this tolerance are removed}

\item{tolerance2}{a maximum distance for snapping nodes together, new nodes within this tolerance are not distinct from existing nodes}
}
\value{
A sfnetwork with the indicated nodes added
}
\description{
This function adds additional nodes or edge segments to a network. Before adding anything, it checks that the new features
are distinct from pre-existing features, and that they are relatively close to the existing network.
Note, this function has not been tested extensively.
}
