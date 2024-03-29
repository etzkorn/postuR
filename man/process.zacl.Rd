% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_zacl.R
\name{process.zacl}
\alias{process.zacl}
\title{Process Zio ZACL File}
\usage{
process.zacl(
  data,
  p = 0.95,
  k = 0.98,
  theta.star = 45,
  epoch.seconds = 60,
  nonwear.window = 3 * 94 * 60,
  nonwear.tol = 10,
  minimum.wear.bout = 94 * 60 * 24,
  cluster = "meanShift"
)
}
\arguments{
\item{data}{Data frame containing the accelerometer data. Should include variables time, x, y, z.}

\item{epoch.seconds}{Duration (seconds) of interval in which to summarise the data. Options include 10, 30, 60, 300, or "none". If epoch.seconds is "none", then data is processed at the original frequency of the data.#'}
}
\value{
A data frame with epoch level activity summaries.
For each time stamp in the resulting data frame, the summary corresponds to the following epoch.seconds
(i.e. 22/9/19 10:50:00 corresponds to the time interval from 22/9/19 10:50:00 to 22/9/19 10:54:59
if epoch.seconds=300).

\item{cluster}{Which posture cluster does an epoch belong to?}

\item{down}{Indicator for recumbent.}

\item{wear.bout}{Integer index for period of consecutive wear epochs an epoch corresponds to.}

\item{theta}{estimated inclination of the chest during the epoch.}

\item{mad}{Mean absolute deviation, a measure of activity intensity during the epoch.}

\item{x,y,z}{Median accelerations measured along each axis during an epoch.}

\item{clustering}{Clustering algorithm to create posture groups.
Either "meanShift", "ward", or "centroid".}
}
\description{
Takes the raw Zio accelerometer file and processes it into adjacent intervals of length epoch.seconds.
}
