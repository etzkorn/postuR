#' @title Inter-Daily Stability
#'
#' @description Calculate Inter-Daily Stability (IS) from one day of data using
#' a given number of epochs (usually minutes) per bin.
#' @param activity vector representing timeseries across all days of activity data.
#' Each entry corresponds to one epoch.
#' @param timestamp same length as activity vector.

iv_oneday <- function(activity_oneday, bin_epochs){
    if(length(activity_oneday)%%bin_epochs != 0){
        stop("Variability bin duration must be factor of length(activity_oneday)")
    }

    activity_oneday %*% kronecker(
        diag(length(activity_oneday)/bin_epochs),
        matrix(1/bin_epochs, nrow = bin_epochs, ncol = 1)
    )

    mean.activity = mean(activity_oneday)
    numerator = sum(diff(activity_oneday)^2)/(length(activity_oneday)-1)
    denominator = sum((activity_oneday - mean.activity)^2) / length(activity_oneday)
    return(numerator/denominator)
}
