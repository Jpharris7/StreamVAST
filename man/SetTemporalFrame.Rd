% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_VASTutility.R
\name{SetTemporalFrame}
\alias{SetTemporalFrame}
\title{This function takes a streamvast object defines a temporal sampling frame and formats everything for use with vast}
\usage{
SetTemporalFrame(
  streamvast,
  startdate = NA,
  enddate = NA,
  padzero = T,
  Time = "Year"
)
}
\arguments{
\item{streamvast}{a streamvast object, minimally with countdata, reachdata defined}

\item{startdate}{dates to truncate or expand the data in yyyy-mm-dd format}

\item{enddate}{dates to truncate or expand the data in yyyy-mm-dd format}

\item{padzero}{should survey with zero observations be added at start/end dates}

\item{Time}{character giving the timescale to use}
}
\value{
A streamvast object with formated with the appropriate temporal frame
}
\description{
This function takes a streamvast object defines a temporal sampling frame and formats everything for use with vast
}
