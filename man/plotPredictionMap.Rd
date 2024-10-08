% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_VASTutility.R
\name{plotPredictionMap}
\alias{plotPredictionMap}
\title{Maps, the desired data, with options for facetting}
\usage{
plotPredictionMap(
  streamvast,
  mapvar = "Density",
  facet = NA,
  FUN = "mean",
  background,
  subset,
  make.labels = F,
  xaxis.breaks = NA,
  yaxis.breaks = NA,
  palette = "turbo",
  max = Inf
)
}
\arguments{
\item{streamvast}{a streamvast object with a prediction data frame}

\item{mapvar}{the column in data to map}

\item{facet}{a column to determine the facets, or a vector of values}

\item{FUN}{a function like mean or median to summarize data}

\item{background}{a sf object to put in the background}

\item{subset}{a logical vector indicating which pred values to use}

\item{make.labels}{should the reaches be labelled, turn off for multiple years}

\item{xaxis.breaks}{Vector of x values for axis lines; useful to reduce crowding}

\item{yaxis.breaks}{Vector of y values for axis lines; useful to reduce crowding}

\item{palette}{a viridis color palette for plotting}

\item{max}{a maximum value to cap the scale; values above the max are set equal to the max}
}
\value{
A colored heatmap of the desired variable mapped over the stream network
}
\description{
Maps, the desired data, with options for facetting
}
