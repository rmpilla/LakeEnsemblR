% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_stats.R
\name{calc_stats}
\alias{calc_stats}
\title{Calculate model fitness and stratification statistics}
\usage{
calc_stats(
  obs,
  model,
  depths,
  folder = ".",
  NH,
  flake_nml,
  out_time,
  par_file,
  start,
  out_hour = 0
)
}
\arguments{
\item{obs}{dataframe; observation data in the LER standardised format}

\item{model}{character; Model for which output will be analysed. Options include c('GOTM', 'GLM', 'Simstrat', 'FLake')}

\item{depths}{vector; Depths which are to be analysed in model output.}

\item{folder}{filepath; to folder which contains the model folders generated by export_config()}

\item{NH}{boolean; northern hemisphere? TRUE or FALSE. Defaults to true}

\item{flake_nml}{filepath; To FLake nml file. Only used if model = 'FLake'}

\item{out_time}{vector; of output time values to subset data by passed to `read_flake_out`.}

\item{par_file}{filepath; To Simstrat par file. Only used if model = 'Simstrat'}

\item{start}{character; of start to set as origin for Simstrat output}

\item{out_hour}{numeric; hour of output time values to subset data  to be passed to `read_flake_out`. Only used for FLake if model time step is 86400s.}
}
\description{
Calculate model fitness and stratification statistics for each model output.
}
