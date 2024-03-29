% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.SeqSGPVeffects.R
\name{plot.SeqSGPVeffects}
\alias{plot.SeqSGPVeffects}
\title{Plot the design-based operating characteristics conditioned across a range of effects.}
\usage{
\method{plot}{SeqSGPVeffects}(
  am,
  stat,
  wait,
  steps,
  affirm,
  lag,
  N,
  xlim,
  ylim,
  log,
  addRegions = TRUE,
  addRegionLines = TRUE,
  addMain = TRUE,
  addLegend = TRUE,
  ablineH = NULL,
  ablineV = NULL,
  setMargins,
  xlab,
  ylab = "",
  ...
)
}
\arguments{
\item{am}{SeqSGPV object}

\item{stat}{Any of: rejH0, cover, bias, stopInconclusive, stopNotROPE, stopNotROME, stopInconsistent.  Use prefix lag.x when some outcomes were lagged/delayed.}

\item{wait}{Vector of possible wait times (W) before monitoring.}

\item{steps}{Vector of the number of observations (S) in between monitoring assessments.}

\item{affirm}{Vector of the number of observations required for affirming a stopping rule (A) in between raising an alert of stopping.}

\item{lag}{Vector of the number of delayed outcomes between enrolling a single observation and observing its outcome.}

\item{N}{Vector of maximum sample size. Can be set as Inf (indefinite).}

\item{xlim}{Optional limits for x-axis.}

\item{ylim}{Optional limits for y-axis.}

\item{addRegions}{Defaults to TRUE to add PRISM boundaries to plot.}

\item{addRegionLines}{Defaults to TRUE to add PRISM boundary lines to plot.}

\item{addMain}{Provides default title for figure. Defaults to TRUE.}

\item{addLegend}{Provides default legend for figure. Defaults to TRUE.}

\item{ablineH}{Add any set of horizontal reference lines.}

\item{ablineV}{Add any set of vertical reference lines.}

\item{setMargins}{Optional figure margins}

\item{ylab}{ylab}

\item{...}{additional plot parameters}
}
\description{
Condition all but two monitoring frequency parameters (wait, steps, affirm, lag, N) must have a single value.
}
