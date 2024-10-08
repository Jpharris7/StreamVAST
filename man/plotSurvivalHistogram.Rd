% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Survival_Functions.R
\name{plotSurvivalHistogram}
\alias{plotSurvivalHistogram}
\title{Plot a histogram of survival estimates across years and reaches}
\usage{
plotSurvivalHistogram(table, year = "all", reach = "all", title)
}
\arguments{
\item{table}{A table of survival estimates, usually from'SurvivalTable'}

\item{year}{numeric; a subset of years to include in the histogram}

\item{reach}{numeric; a subset of reach numbers to include in the histogram}

\item{title}{character; a title to pass to ggtitle()}
}
\value{
a ggplot object for the histogram
}
\description{
Plot a histogram of survival estimates across years and reaches
}
