#' @title Vector Magnitude Counts (Mean Absolute Deviation)
#'
#' @description Calculate mean absolute deviation from triaxial accelerometer data for a collection of time points.
#' @param x First accelerometer channel.
#' @param y Second accelerometer channel.
#' @param z Third accelerometer channel.
#' @param roll Should mad be calculated using a sliding window?
#' @param window What number of observations should be used in the sliding window?
#' Default is 94, as Zio measures approximately 94 observations per minute.
#'
mad <- function(x, y, z, roll = FALSE, window = 94){
	if(roll == FALSE){
		m <- mean(postuR:::euclid.norm(x,y,z))
		mad.out <- mean( abs( postuR:::euclid.norm(x,y,z) - m) )
		return(mad.out)
	}else if(roll == TRUE){
		m <- runstats::RunningMean(postuR:::euclid.norm(x,y,z), W = window,circular = F )
		mad <- runstats::RunningMean( abs( postuR:::euclid.norm(x,y,z) - m), W = window,circular = F)
		return(mad)
	}
}
