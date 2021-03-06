% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/appflow.R
\name{writeShinyFlow}
\alias{writeShinyFlow}
\title{Write a shiny App}
\usage{
writeShinyFlow(workingfolder, app.name, run.app = TRUE)
}
\arguments{
\item{workingfolder}{the working folder}

\item{app.name}{name of the application (will be a folder name)}

\item{run.app}{should the shiny app be run?}
}
\value{
a shiny app in the app.name folder within the working folder
}
\description{
Write a shiny app that draws network connections between words.
}
