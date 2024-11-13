#' @title Calculate Common Activity Summaries from Summarised Data
#' @description Takes a data frame with activity measurements (MAD), and time.
#' Returns a data frame with one row with relevant global activity summaries.
#' This function was validated against `arctools::activity_stats`.
#' @param data A data frame with epoch-level activity summaries.
#' This function is designed to input the data frame returned from postuR::process.zacl directly.
#' If not, variables should include time (class dttm), mad (mean amplitude deviation),
#' wear (indicator for whether the device was actually worn during a given period),
#' and time.group (a count variable indicating the distinct wear bout).
#' @param epoch.seconds How many seconds does each row in the "data" argument correspond to?
#' @param inactive Upper bound activity intensity cutpoint for inactivity category from mean absolute deviation (milli-gravitational units).
#' @param light Upper bound activity intensity cutpoint for light activity category from mean absolute deviation (milli-gravitational units).
#' @param vlight Upper bound activity intensity cutpoint for very light activity category from mean absolute deviation (milli-gravitational units).
#' @param frag.start Hour of the day at which to start calculating fragmentation measures. Use 0 for 24 hours.
#' @param frag.stop Hour of the day at which to stop calculating fragmentation measures. Use 24 for 24 hours.
#' @returns data frame with the following variables
#'
#' - start: time of starting epoch
#'
#' - end: time of last epoch
#'
#' - epoch.seconds: duration of epoch-level summaries (same as input)
#'
#' - days: duration of recording time in days
#'
#' - wear.days: duration of wear time in days (should be less than days)
#'
#' - wear.bout.count: number of wear bouts in the case of intermittant non-wear
#'
#' - mad_0...mad_22: average mean absolute deviation of acceleration in two hour bins throughout the day.
#'
#' - mad: daily average mean absolute deviation
#'
#' - down: daily hours of time spent lying down
#'
#' - sst: daily hours of time spent sleep or sedentary (mad <= inactive)
#'
#' - lipa: daily hours of time spent in light physical activity (inactive < mad <= light)
#'
#' - mvpa: daily hours of time spent in moderate-to-vigorous physical activity (light < mad)
#'
#' @export

global.activity.summaries <- function(
        data,
        epoch.seconds = 60,
        inactive = 9.041,
        vlight = 28.187,
        light = 58.083,
        frag.start = 10,
        frag.stop = 18
){

    epoch.seconds = abs(as.numeric(difftime(data$time[1], data$time[2], unit = "secs")))

    # generate time bins
	data <-
	data %>%
	dplyr::mutate(
	    time = lubridate::as_datetime(.data$time),
		minute.bin = lubridate::hour(.data$time)*60 + lubridate::minute(.data$time),
	    hour.bin = lubridate::hour(lubridate::floor_date(.data$time, unit = lubridate::hours(2))),
	    day.bin = lubridate::floor_date(.data$time, unit = lubridate::days(1))
	)

	# pull meta data
	meta <-
		data %>%
		dplyr::summarise(
		    start = min(time),
		    end = max(time),
		    epoch.seconds = difftime(time[2], time[1], units = "secs") %>% as.numeric,
		    days = difftime(end, start, units = "days") %>% as.numeric,
		    wear.days = mean(wear)*days,
		    wear.bout.count = length(unique(wear.bout)),
		    .groups = "keep"
		)

	# generate hour-binned summaries
	hour.level <-
		data %>%
		dplyr::filter(wear == 1) %>%
	    dplyr::mutate(hour.bin = paste("mad", hour.bin, sep = "_")) %>%
		dplyr::group_by(hour.bin, .add = T) %>%
		dplyr::summarise(
		    mad = mean(mad, na.rm=T),
		    .groups = "keep") %>%
		tidyr::pivot_wider(names_from = "hour.bin", values_from = "mad") %>%
		dplyr::select(mad_0, mad_2, mad_4, mad_6, mad_8, mad_10, mad_12,
		       mad_14, mad_16, mad_18, mad_20, mad_22)

	# generate minute-binned summaries
	# average first by minute in day to adjust for differential wear periods
	global.level <- data %>%
		dplyr::filter(wear == 1) %>%
		dplyr::group_by(minute.bin, .add = T) %>%
	    #create diurnal profile
		dplyr::summarise(
		    n.day = n(),
		    down = sum(down, na.rm=T)/n.day,
		    sst = sum(mad <= inactive, na.rm = T)/n.day,
		    lipa = sum(mad <= light & mad > inactive, na.rm = T)/n.day,
		    mvpa = sum(mad > light, na.rm = T)/n.day,
		    mad = sum(mad, na.rm=T)/n.day,
		    .groups = "keep") %>%
		dplyr::ungroup(minute.bin) %>%
	    #summarize across day
		dplyr::summarise(
		    mad = sum(mad, na.rm=T)/1440*60/epoch.seconds,
		    down = sum(down, na.rm=T)/60*60/epoch.seconds,
		    sst = sum(sst, na.rm = T)/60*60/epoch.seconds,
		    lipa = sum(lipa, na.rm = T)/60*60/epoch.seconds,
		    mvpa = sum(mvpa, na.rm = T)/60*60/epoch.seconds,
		    .groups = "keep")

	# generate day-binned (M10, and L6)
	day.level <-
	    data %>%
	    dplyr::select(day.bin, minute.bin, mad, wear) %>%
	    dplyr::group_by(day.bin, .add=T) %>%
	    dplyr::filter(sum(wear, na.rm=T)>=1440*60/epoch.seconds) %>%
	    dplyr::group_modify(
	        ~full_join(tibble(minute.bin = 0:1439), ., by = "minute.bin")
	    ) %>%
	    dplyr::arrange(day.bin, minute.bin) %>%
		dplyr::mutate(
		    mad600 = runstats::RunningMean(mad, W = 600*60/epoch.seconds, circular = T),
		    mad300 = runstats::RunningMean(mad, W = 300*60/epoch.seconds, circular = T),
		    mad120 = runstats::RunningMean(mad, W = 120*60/epoch.seconds, circular = T),
		    mad60 = runstats::RunningMean(mad, W = 60*60/epoch.seconds, circular = T),
		    mad30 = runstats::RunningMean(mad, W = 30*60/epoch.seconds, circular = T),
		    mad10 = runstats::RunningMean(mad, W = 10*60/epoch.seconds, circular = T)) %>%
		dplyr::summarise(
		    valid = ifelse(sum(wear, na.rm=T)==1440*60/epoch.seconds, 1 , NA),
		    M600 = max(mad600)*valid,
		    M300 = max(mad300)*valid,
		    M120 = max(mad120)*valid,
		    M60 = max(mad60)*valid,
		    M30 = max(mad30)*valid,
		    M10 = max(mad10)*valid,
		    TM10 = which.max(mad10)*60/epoch.seconds*valid,
		    TM30 = which.max(mad30)*60/epoch.seconds*valid,
		    TM60 = which.max(mad60)*60/epoch.seconds*valid,
		    TM120 = which.max(mad120)*60/epoch.seconds*valid,
		    TM300 = which.max(mad300)*60/epoch.seconds*valid,
		    TM600 = which.max(mad600)*60/epoch.seconds*valid,
		    L600 = min(mad600)*valid,
		    L300 = min(mad300)*valid,
		    L120 = min(mad300)*valid,
		    TL600 = which.min(mad600)*60/epoch.seconds*valid,
		    TL300 = which.min(mad300)*60/epoch.seconds*valid,
		    TL120 = which.min(mad120)*60/epoch.seconds*valid,
		    iv_60 = iv_oneday(mad, 60),
		    iv_30 = iv_oneday(mad, 30),
		    iv_20 = iv_oneday(mad, 20),
		    iv_10 = iv_oneday(mad, 10),
		    iv_05 = iv_oneday(mad, 5),
		    iv_01 = iv_oneday(mad, 1),
		    .groups = "keep") %>%
		dplyr::ungroup(day.bin) %>%
		dplyr::summarise(
		    M600 = mean(M600, na.rm = T),
		    M300 = mean(M300, na.rm = T),
		    M120 = mean(M120, na.rm = T),
		    M60 = mean(M60, na.rm = T),
		    M30 = mean(M30, na.rm = T),
		    M10 = mean(M10, na.rm = T),
		    L600 = mean(L600, na.rm = T),
		    L300 = mean(L300, na.rm = T),
		    L120 = mean(L120, na.rm = T),
		    TM600 = (
		        atan2(
		            mean(sin(TM600*2*pi/1440),na.rm=T),
		            mean(cos(TM600*2*pi/1440),na.rm=T)
		            )/2/pi*1440 + 1440) %%1440, # circular mean
		    TM300 =
		        (atan2(
		            mean(sin(TM300*2*pi/1440),na.rm=T),
		            mean(cos(TM300*2*pi/1440),na.rm=T)
		            )/2/pi*1440 + 1440) %%1440,
		    TM120 =
		        (atan2(
		            mean(sin(TM120*2*pi/1440),na.rm=T),
		            mean(cos(TM120*2*pi/1440),na.rm=T)
		            )/2/pi*1440 + 1440) %%1440,
		    TM60 =
		        (atan2(
		            mean(sin(TM60*2*pi/1440),na.rm=T),
		            mean(cos(TM60*2*pi/1440),na.rm=T)
		            )/2/pi*1440 + 1440) %%1440,
		    TM30 =
		        (atan2(
		            mean(sin(TM30*2*pi/1440),na.rm=T),
		            mean(cos(TM30*2*pi/1440),na.rm=T))/2/pi*1440 + 1440) %%1440,
		    TM10 =
		        (atan2(mean(sin(TM10*2*pi/1440),na.rm=T),
		               mean(cos(TM10*2*pi/1440),na.rm=T)
		               )/2/pi*1440 + 1440) %%1440,
		    TL600 =
		        (atan2(mean(sin(TL600*2*pi/1440),na.rm=T),
		               mean(cos(TL600*2*pi/1440),na.rm=T)
		               )/2/pi*1440 + 1440) %%1440,
		    TL300 =
		        (atan2(
		            mean(sin(TL300*2*pi/1440),na.rm=T),
		            mean(cos(TL300*2*pi/1440),na.rm=T)
		            )/2/pi*1440 + 1440) %%1440,
		    TL120 =
		        (atan2(
		            mean(sin(TL120*2*pi/1440),na.rm=T),
		            mean(cos(TL120*2*pi/1440),na.rm=T)
		            )/2/pi*1440 + 1440) %%1440,
		    iv_60 = mean(iv_60),
		    iv_30 = mean(iv_30),
		    iv_20 = mean(iv_20),
		    iv_10 = mean(iv_10),
		    iv_05 = mean(iv_05),
		    iv_01 = mean(iv_01),
		    .groups = "keep")

	#Interdaily stability measures
	interdaily.stability = tibble(
	    is_60 = is_alldays(activity = data$mad, time = data$time, bin_epochs = 60, epoch_seconds = epoch.seconds),
	    is_30 = is_alldays(activity = data$mad, time = data$time, bin_epochs = 30, epoch_seconds = epoch.seconds),
	    is_20 = is_alldays(activity = data$mad, time = data$time, bin_epochs = 20, epoch_seconds = epoch.seconds),
	    is_10 = is_alldays(activity = data$mad, time = data$time, bin_epochs = 10, epoch_seconds = epoch.seconds),
	    is_05 = is_alldays(activity = data$mad, time = data$time, bin_epochs = 5, epoch_seconds = epoch.seconds),
	    is_01 = is_alldays(activity = data$mad, time = data$time, bin_epochs = 1, epoch_seconds = epoch.seconds)
	)

	#fragmentation measures
	# subset to window of time delineated by frag.start
	frag.data <- data %>%
	    dplyr::mutate(
	        wear = ifelse(
	            (minute.bin >= frag.start*60)|(minute.bin < frag.stop*60),
	            wear, 0)) %>%
	    dplyr::group_by(day.bin,.add = T) %>%
	    dplyr::group_modify(
	        ~ActFrag::fragmentation(
	            x = as.integer(.$mad >= inactive),
	            thresh = 1,
	            w = .$wear == 1,
	            metrics = "all",
	            bout.length = 1
	       ) %>%
	            dplyr::as_tibble()%>%
	            dplyr::rename(
	                mean_rest_bout = mean_r,
	                mean_active_bout = mean_a,
	                satp = SATP,
	                astp = ASTP)
	    ) %>%
	    dplyr::ungroup(day.bin) %>%
	    dplyr::summarise(
	        satp = mean(satp),
	        astp = mean(astp),
	        mean_active_bout = mean(mean_active_bout),
	        mean_rest_bout = mean(mean_rest_bout),
	        .groups = "keep")

	# merge all summaries
	if(dplyr::is.grouped_df(meta)){
	    full_join(meta, global.level, keep = F) %>%
	        full_join(hour.level, keep = F) %>%
	        full_join(day.level, keep = F) %>%
	        full_join(frag.data, keep = F)
	}else{
	    dplyr::bind_cols(meta, global.level, hour.level, day.level, interdaily.stability, frag.data)
	}
}
