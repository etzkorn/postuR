#' @title Inter-Daily Stability
#'
#' @description Calculate Inter-Daily Stability (IS) from one day of data using
#' a given number of epochs (usually minutes) per bin.
#' @param activity vector representing timeseries across all days of activity data.
#' Each entry corresponds to one epoch.
#' @param time same length as activity vector.
#' @param bin_epochs number of epochs per bin to collapse activity into before performing
#' IS formula.

is_alldays <- function(activity, time, bin_epochs, epoch_seconds){

    if(round(epoch_seconds * bin_epochs/60) != (epoch_seconds * bin_epochs/60)){
        stop("bin_epochs not compatible with epoch duration.")
    }

    bin_minutes = epoch_seconds*bin_epochs/60
    daily_bins = 1440/bin_minutes

    df.out <- tibble(
        time = time,
        activity = activity
    ) %>% mutate(
        date = floor_date(time, unit = "days"),
        bin = floor_date(time, unit = minutes(bin_minutes)),
        bin_group = hour(bin)*60 + minute(bin)
    ) %>% group_by(
        date, bin, bin_group
    ) %>% filter(
        sum(!is.na(activity)) == bin_epochs
    ) %>% summarise(
        bin_activity = mean(activity),
        .groups = "drop"
    ) %>% group_by(bin_group) %>% mutate(
        bin_group_count = n(),
        bin_group_mean = mean(bin_activity)
    )

    if(any(df.out$bin_group_count < 3)) return(NA)
    if(length(unique(df.out$bin_group))!= daily_bins) return(NA)

    df.out <- df.out %>% ungroup %>% mutate(
        overall_weighted_mean = sum(bin_group_mean/bin_group_count)/daily_bins
    ) %>% group_by(bin_group) %>% summarise(
        bin_group_count = bin_group_count[1],
        bin_msse = sum((bin_activity - overall_weighted_mean)^2/bin_group_count),
        bin_group_msse = (bin_group_mean[1] - overall_weighted_mean[1])^2,
        .groups = "drop"
    ) %>% ungroup %>% summarise(
        all_bin_msse = mean(bin_msse),
        all_bin_group_msse = mean(bin_group_msse),
        n_days = mean(bin_group_count),
        is = all_bin_group_msse/all_bin_msse*(n_days-1)/n_days
    ) %>% select(is) %>% unlist

    return(df.out)
}
