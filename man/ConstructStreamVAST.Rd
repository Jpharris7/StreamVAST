% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_VASTutility.R
\name{ConstructStreamVAST}
\alias{ConstructStreamVAST}
\title{ConstructStreamVAST
This function constructs objects with class "StreamVAST". Depending on the arguments supplied, it will either
initialize the object from countdata and reachdata, or it will add details to an already defined StreamVAST object.}
\usage{
ConstructStreamVAST(
  countdata,
  reachdata,
  surveydata,
  vastdata,
  vastmodel,
  countname,
  reachname,
  timetable,
  timescale,
  timerange,
  title,
  stats,
  preds,
  eval,
  spacedata,
  timedata,
  aucdata,
  unitconv,
  covariates
)
}
\arguments{
\item{countdata}{A dataframe containing count data, such as outputed from AssembleReddData}

\item{reachdata}{A sf object with data for reaches, such as from AssignReaches}

\item{surveydata}{A optional data frame or sf object with survey information, such as from MakeSurveyTracks, doesn't do much right now}

\item{vastdata}{A dataframe formatted specifically for use with VAST}

\item{vastmodel}{A fitted VAST model}

\item{countname}{A column name from countdata to use}

\item{reachname}{a column name from reachdata to use}

\item{timetable}{a dataframe converting timesteps to dates}

\item{timescale}{The temporal resolution}

\item{timerange}{A length two vector with the start and enddate for the temporal}

\item{title}{A title or name for the stream system, use for plotting}

\item{stats}{A list with summary stats, like RMSE,AIC}

\item{preds}{A dataframe of all predictions from the VAST model}

\item{eval}{A limited dataframe for evaluating model performance}

\item{spacedata}{A dataframe with results aggregated by spatial frame}

\item{timedata}{A dataframe with results aggregated by time frame}

\item{aucdata}{a dataframe with Area under the Curve calculations derived from model}

\item{covariates}{a list with elements for the covariates}

\item{StreamVAST}{An StreamVAST object to append details}

\item{vastinput}{A dataframe with spatial information for VAST, requires a specific format}

\item{vastnetwork}{A dataframe with spatial information for VAST's stream functionality}

\item{vastnetworkLL}{A dataframe with spatial information needed for VAST's stream functionality}

\item{vastsettings}{a list of settings used to run a vast model}
}
\value{
A streamvast object with appropriate formatting
}
\description{
ConstructStreamVAST
This function constructs objects with class "StreamVAST". Depending on the arguments supplied, it will either
initialize the object from countdata and reachdata, or it will add details to an already defined StreamVAST object.
}
