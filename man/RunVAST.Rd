% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_VASTutility.R
\name{RunVAST}
\alias{RunVAST}
\title{Title}
\usage{
RunVAST(streamvast, vastsettings, optimize = T, maxiter = 3, startpar)
}
\arguments{
\item{streamvast}{a streamvast object with a defined temporal frame}

\item{vastsettings}{a list of settings for a vast model}

\item{optimize}{logical, should the algorthim run multiple models with varying settings}

\item{maxiter}{a maximum number of iterations to try before giving up}

\item{startpar}{a list of starting values for the fixed parameters}
}
\value{
A streamvast object with a fitted model attached
}
\description{
Title
}
