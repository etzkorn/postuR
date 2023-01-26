#' @title Calculate Common Activity Summaries from Summarised Data
#' @description Takes a data frame with activity measurements (MAD), and time. Returns a data frame with one row with relevant global activity summaries.
#' @param data A data frame with epoch-level activity summaries. This function is designed to input the data frame returned from Rtbeat::process.zacl directly. If not, variables should include time (class dttm), mad (Vector magnitude counts), wear (indicator for whether the device was actually worn during a given period), and time.group (a count variable indicating the distinct wear bout).
#' @param epoch.seconds How many seconds does each row in the "data" argument correspond to?
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
#' - down: daily proportion of time spent recumbent
#'
#' @export

global.activity.summaries <- function(
        data, epoch.seconds = 60, sedentary = 0.0074, light = 0.0569){
	data <-
	data %>%

	# generate time bins
	dplyr::mutate(
	    time = lubridate::as_datetime(.data$time),
		minute.bin = lubridate::hour(.data$time)*60 + lubridate::minute(.data$time),
	    hour.bin = lubridate::hour(lubridate::floor_date(.data$time, unit = lubridate::hours(2))),
	    day.bin = lubridate::floor_date(.data$time, unit = lubridate::days(1))
	)

	# meta data
	meta <-
		data %>%
		dplyr::summarise(
		    start = min(time),
		    end = max(time),
		    epoch.seconds = difftime(time[2], time[1], units = "secs") %>% as.numeric,
		    days = difftime(end, start, units = "days") %>% as.numeric,
		    wear.days = ifelse(sum(wear) ==0, days, sum(wear)*epoch.seconds/60/1440),
		    wear.bout.count = length(unique(wear.bout)),
		    .groups = "drop"
		)

	# generate hour-binned summaries
	hour.level <-
		data %>%
		dplyr::filter(wear == 1) %>%
		dplyr::group_by(hour.bin) %>%
		dplyr::summarise(
		    mad = mean(mad, na.rm=T),
		    .groups = "drop") %>%
		dplyr::select(hour.bin, mad) %>%
		tidyr::gather(key = "var", value = "value", -hour.bin) %>%
		tidyr::unite(col = var, var, hour.bin) %>%
		dplyr::mutate(var = factor(var, ordered = T))%>%
		tidyr::spread(key = var, value = value)#%>%
		#dplyr::select(mad_0, mad_2, mad_4, mad_6, mad_8, mad_10, mad_12,
		#       mad_14, mad_16, mad_18, mad_20, mad_22)

	# generate minute-binned summaries
	# average first by minute in day to adjust for differential wear periods
	global.level <- data %>%
		dplyr::filter(wear == 1) %>%
		dplyr::group_by(minute.bin) %>%
	    #create diurnal profile
		dplyr::summarise(
		    mad = mean(mad, na.rm=T),
		    down = mean(down, na.rm=T),
		    st = mean(mad <= sedentary, na.rm = T),
		    lipa = mean(mad <= light & mad > sedentary, na.rm = T),
		    mvpa = mean(mad > light, na.rm = T),
			.groups = "drop") %>%
		dplyr::ungroup() %>%
	    #summarize across day
		dplyr::summarise(
		    mad = mean(mad, na.rm=T),
		    down = mean(down, na.rm=T)*24,
		    st = mean(st, na.rm = T)*24,
		    lipa = mean(lipa, na.rm = T)*24,
		    mvpa = mean(lipa, na.rm = T)*24)

	# generate day-binned (M10, and L6)
	# day.level <- data %>%
	#	dplyr::filter(wear == 1) %>%
	#	dplyr::group_by(day.bin) %>%
	#	dplyr::filter(n() > 10*60*60/epoch.seconds)

	#if(nrow(day.level)>=10/24*1440){
	#day.level <-
	#	dplyr::mutate(mad10 = runstats::RunningMean(mad, W = 10/24*1440, circular = T),
	#	              mad6 = runstats::RunningMean(mad, W = 6/24*1440, circular = T)) %>%
	#	dplyr::summarise(M10 = max(mad10),
	#		     TM10 = which.max(mad10),
	#		     L6 = min(mad6),
	#		     TL6 = which.min(mad6),
	#		     .groups = "drop") %>%
	#	dplyr::ungroup() %>%
	#	dplyr::summarise(M10 = mean(M10),
	#		     TM10 = atan2(mean(sin(TM10*2*pi/1440)), mean(cos(TM10*2*pi/1440))),
	#		     L6 = mean(L6),
	#		     TL6 = atan2(mean(sin(TL6*2*pi/1440)), mean(cos(TL6*2*pi/1440))),
	#		     .groups = "drop")
	#}else{
	#	day.level <- tibble(M10 = NA,TM10=NA,L6 = NA, TL6=NA)
	#}

	# fragmentation measures
	# frag <- ActFrag::fragmentation(x = as.integer(data$mad >= 0.012),
	#		  thresh = 1,
	#		  w = data$time.group != 0,
	#		  metrics = "all",
	#		  bout.length = 1) %>%
	#	dplyr::as_tibble() %>%
	#	dplyr::transmute(mean_r, mean_a, SATP, ASTP)

	# merge all summaries
	dplyr::bind_cols(meta, global.level, hour.level)
}
