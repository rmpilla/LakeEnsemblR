% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export_config.R
\name{export_config}
\alias{export_config}
\title{Export model configuration setups}
\usage{
export_config(
  config_file,
  model = c("GOTM", "GLM", "Simstrat", "FLake"),
  folder = "."
)
}
\arguments{
\item{config_file}{name of the master LakeEnsemblR config file}

\item{model}{vector; model to export configuration file. Options include c('GOTM', 'GLM', 'Simstrat', 'FLake')}

\item{folder}{folder}
}
\description{
Create directory with file setups for each model based on a master LakeEnsemblR config file
}
\examples{


}
\author{
Tadhg Moore, Jorrit Mesman
}
\keyword{methods}
