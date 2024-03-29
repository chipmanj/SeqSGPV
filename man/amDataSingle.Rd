% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amDataSingle.R
\name{amDataSingle}
\alias{amDataSingle}
\title{amDataSingle}
\usage{
amDataSingle(
  dataGeneration,
  dataGenArgs,
  effectGeneration,
  effectGenArgs,
  effectScale,
  allocation,
  randomize,
  modelFit,
  modelFitArgs,
  existingData
)
}
\arguments{
\item{dataGeneration}{See SeqSGPV}

\item{dataGenArgs}{See SeqSGPV}

\item{effectGeneration}{See SeqSGPV}

\item{effectGenArgs}{See SeqSGPV}

\item{effectScale}{See SeqSGPV}

\item{allocation}{See SeqSGPV}

\item{randomize}{TRUE if length(allocation) > 1}

\item{modelFit}{See SeqSGPV}

\item{modelFitArgs}{See SeqSGPV}

\item{existingData}{Previously provided data}
}
\description{
Generate a single adaptive monitoring simulation.  amDataSingle is called within SeqSGPV.
}
