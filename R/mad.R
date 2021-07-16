#' @title Vector Magnitude Counts (Mean Absolute Deviation)
#'
#' @description Calculate mean absolute deviation from triaxial accelerometer data for a collection of time points.
#' @param x First accelerometer channel.
#' @param y Second accelerometer channel.
#' @param z Third accelerometer channel.
#' @param roll Should mad be calculated using a sliding window?
#' @param window What number of observations should be used in the sliding window? Default is 90, as Zio measures approximately 90 observations per minute.
#'
mad <- function(x, y, z, roll = FALSE, window = 90){
	if(roll == FALSE){
		m <- mean(euclid.norm(x,y,z))
		mean( abs( euclid.norm(x,y,z) - m) ) %>% return()
	}else if(roll == TRUE){
		m0 <- mean(euclid.norm(x,y,z))
		m <- zoo::rollmean(euclid.norm(x,y,z), k = window, fill = m0)
		mad <- zoo::rollmean( abs( euclid.norm(x,y,z) - m), k = window, fill = 0)
		return(mad)
	}
}
