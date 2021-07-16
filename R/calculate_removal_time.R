#' @title Calculate Removal Time of a Zio Accelerometer Device
#'
#' @description Takes a data frame containing raw triaxial accelerometer data.
#' Estimates a SINGLE removal time with accompanying inferrential metrics.
#' This function is primarily for use by process.zacl, and not directly by users.
#' @param data Data frame or matrix containing the accelerometer data.
#' Should include variables "time", "x", "y", "z".
#' We recommend the data include at least 24 hours of wear time, and data should not have any non-wear periods.
#' @param p Percentile for highest accelerations to be included in calculating the upright posture for the individual.
#' @param multiple.removals Still in development. Set to FALSE for now, for detecting only one change point.
#' @param W Window size used for detecting multiple removals. Only used when multiple.removals = T.
#' @return A data frame with statistics relevant for detecting device removal and replacement.
#'
#' * time: estimated removal time
#'
#' * r.ratio: ratio of total vs. cross.split mean resultant length,
#'
#' * r.split: mean resultant split within two splits of the data
#'
#' * r.total: marginal or overall mean resultant length (across both splits)
#'
#' * phi: angular difference between central orientation before and after split at time
#'
#' * R0, R1: Mean resultant length before and after split times n0, n1 respectively
#'
#' * r0, r1: Mean resultant length before and after split
#'
#' * mx, my, mz: x, y, z coordinates of overall central orientation (across both splits)
#'
#' * m1x, m1y, m1z: x, y, z coordinates of central orientation AFTER split
#'
#' * m0x, m0y,m0z: x, y, z coordinates of central orientation BEFORE split
#'
#' * n0, n1: number of time points before and after each split
#'
#' * p.r2: the individual's quantile of (1-r)^2 associated with percentile p.
#' (i.e. above what squared net acceleration is considered a high acceleration point).
#'
#' * start.time: start time of input accelerometer data
#'
#' * stop.time: end time of input accelerometer data
#'
#' * multiple.removals: same as input
#'
#' * W: same as input, or NA if multiple.removals = F
#'
#' @export
calculate.removal.time <- function(data, p = 0.95, multiple.removals = F, W = 500){
if(!multiple.removals){
	start.time <- min(data$time)
	stop.time <- max(data$time)
	data %>%
	na.omit %>%
	dplyr::mutate(r = sqrt(x^2 + y^2 + z^2)) %>%
	dplyr::filter(r!=0) %>%
	dplyr::mutate(x = x/r,
	       y = y/r,
	       z = z/r,
	       r2 = (1-r)^2,
	       p.r2 = quantile(r2, p)) %>%
	# only use data for which the vector length
	# is above a certain quantile
	dplyr::filter(r2 > p.r2) %>%
	dplyr::mutate(n0 = 1:n(),
	       n1 = n()-n0,
	       # overall means
	       mx = mean(x),
	       my = mean(y),
	       mz = mean(z),
	       # before-split means
	       m0x = cumsum(x)/n0,
	       m0y = cumsum(y)/n0,
	       m0z = cumsum(z)/n0,
	       # after-split means
	       m1x = (sum(x) - cumsum(x))/n1,
	       m1y = (sum(y) - cumsum(y))/n1,
	       m1z = (sum(z) - cumsum(z))/n1,
	       # mean resultant length pre/post split
	       r0 = sqrt(m0x^2 + m0y^2 + m0z^2),
	       r1 = sqrt(m1x^2 + m1y^2 + m1z^2),
	       # resultant length pre/post split
	       R0 = r0*n0,
	       R1 = r1*n1,
	       # angular difference between before/after centers
	       phi = acos((m0x*m1x + m0y*m1y + m0z*m1z)/r0/r1),
	       # net mean resultant length
	       r.total = sqrt(mean(x)^2 + mean(y)^2 + mean(z)^2),
	       r.split = (R0 + R1)/(n0+n1),
	       r.ratio = r.total/r.split) %>%
	dplyr::slice(which.max(r.split)) %>%
	dplyr::select(time, r.ratio, r.split, r.total, phi, R0, R1, r0, r1,
		  mx, my, mz, m1x, m1y, m1z,
		  m0x, m0y,m0z, n0, n1, p.r2) %>%
	dplyr::mutate(start.time = start.time,
		  stop.time = stop.time,
		  multiple.removals = multiple.removals,
		  W = NA)
}else{
	warning("Accounting for multiple removals is still in development.
	        Results from current code should not be used.")
	start.time <- min(data$time)
	stop.time <- max(data$time)
	data %>%
		na.omit %>%
		dplyr::mutate(r = sqrt(x^2 + y^2 + z^2)) %>%
		dplyr::filter(r!=0 & wear != 0) %>%
		dplyr::mutate(x = x/r,
			  y = y/r,
			  z = z/r,
			  r2 = (1-r)^2,
			  p.r2 = quantile(r2, p)) %>%
		# only use data for which the vector length
		# is above a certain quantile
		dplyr::filter(r2 > p.r2) %>%
		dplyr::mutate(# before-split running means
			  m0x = c(rep(NA, W),
			          diff(cumsum(x), lag = W)/W),
			  m0y = c(rep(NA, W),
			          diff(cumsum(y), lag = W)/W),
			  m0z = c(rep(NA, W),
			          diff(cumsum(z), lag = W)/W),
			  # after-split means
			  m1x = c(diff(cumsum(x), lag = W)/W,
			          rep(NA, W)),
			  m1y = c(diff(cumsum(y), lag = W)/W,
			          rep(NA, W)),
			  m1z = c(diff(cumsum(z), lag = W)/W,
			          rep(NA, W)),
			  # overall means
			  mx = (m0x + m1x)/2,
			  my = (m0y + m1y)/2,
			  mz = (m0z + m1z)/2,
			  # mean resultant length pre/post split
			  r0 = sqrt(m0x^2 + m0y^2 + m0z^2),
			  r1 = sqrt(m1x^2 + m1y^2 + m1z^2),
			  # resultant length pre/post split
			  R0 = r0*W,
			  R1 = r1*W,
			  # angular difference between before/after centers
			  phi = acos((m0x*m1x + m0y*m1y + m0z*m1z)/r0/r1),
			  # net mean resultant length
			  r.total = sqrt(mx^2 + my^2 + mz^2),
			  r.split = (R0 + R1)/(2*W),
			  r.ratio = r.total/r.split) %>%
		dplyr::slice(which.min(r.ratio)) %>%
		dplyr::select(time, r.ratio, r.split, r.total, phi, R0, R1, r0, r1,
			  mx, my, mz, m1x, m1y, m1z,
			  m0x, m0y, m0z, n0 = NA, n1 = NA, p.r2) %>%
		dplyr::mutate(start.time = start.time,
			  stop.time = stop.time,
			  multiple.removals = multiple.removals,
			  W = W)
}
}
