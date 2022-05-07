
#' @title Check for Non-Wear in Accelerometer Data
#'
#' @description Checks for non-wear by checking for long intervals of no change in the signal along each axis.
#' @param data Data frame containing the raw accelerometer data. Should include variables time, x, y, z.
#' @param window Number of samples for minimum contiguous period of non-wear. To be passed to accelerometry::weartime.
#' @param tol Number of allowances in any window of non-wear. To be passed to accelerometry::weartime
#' @param tol.upper Maximum allowable sequential change in any window of non-wear. To be passed to accelerometry::weartime
#' @param filter Should nonwear be removed from the data set?
#' @param minimum.wear.bout Minimum contiguous period of wear time for a period to be considered wear.
#' @importFrom rlang .data
#' @returns A data frame with variables time, x, y, z, wear, and wear.bout
#'
#' - wear = an indicator for whether a person is wearing the device.
#'
#' - wear.bout = an integer numbering which consecutive period of contiguous device wear a time point corresponds to. 0 indicates nonwear.
#'
#' @export

check.nonwear <- function(data, window = 3*94*60, tol=10, tol.upper = 0,
                          filter=F, minimum.wear.bout = 94*60*24){
    	if(filter == T){
    		data %>%
    		dplyr::mutate(change =
    		                  sqrt(
    		                      c(0, diff(.data$x))^2 +
    		                      c(0, diff(.data$y))^2 +
    		                      c(0, diff(.data$z))^2
    		                  ) > tol.upper,
    		       wear = accelerometry::weartime(counts = .data$change, window = window, tol = tol, tol_upper = 1),
    		       wear.bout = .data$wear * cumsum(.data$wear > c(0,.data$wear[-n()]))) %>%
    		group_by(wear.bout) %>%
    		mutate(bout.length = n(),
    		       wear = .data$wear &(.data$bout.length > minimum.wear.bout) ) %>%
    		ungroup %>%
    		mutate(wear.bout = .data$wear.bout*.data$wear) %>%
    		dplyr::filter(.data$wear == 1) %>%
    		dplyr::select(-.data$change, -.data$wear)
    	}else{
    		data %>%
    		dplyr::mutate(change =
    		                  sqrt(
    		                      c(0, diff(.data$x))^2 +
    		                      c(0, diff(.data$y))^2 +
    		                      c(0, diff(.data$z))^2
    		                  ) > tol.upper,
    		       wear = accelerometry::weartime(counts = .data$change, window = window, tol = tol, tol_upper = 1),
    		       wear.bout = .data$wear * cumsum(.data$wear > c(0,.data$wear[-n()]))) %>%
    		group_by(wear.bout) %>%
    		mutate(bout.length = n(),
    		       wear = .data$wear &(.data$bout.length > minimum.wear.bout) ) %>% # must be at least 30 minutes
    		ungroup %>%
    		mutate(wear.bout = .data$wear.bout*.data$wear) %>%
    		dplyr::select(-.data$change)
    	}
}

