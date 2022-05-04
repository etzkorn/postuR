#' @title Process Zio ZACL File
#'
#' @description Takes the raw Zio accelerometer file and processes it into adjacent intervals of length epoch.seconds.
#'
#' @param data Data frame containing the accelerometer data. Should include variables time, x, y, z.
#' @param epoch.seconds Duration (seconds) of interval in which to summarise the data. Options include 10, 30, 60, 300, or "none". If epoch.seconds is "none", then data is processed at the original frequency of the data.#'
#' @return A data frame with epoch level activity summaries.
#' For each time stamp in the resulting data frame, the summary corresponds to the following epoch.seconds
#' (i.e. 22/9/19 10:50:00 corresponds to the time interval from 22/9/19 10:50:00 to 22/9/19 10:54:59
#' if epoch.seconds=300).
#'
#' \item{cluster}{Which posture cluster does an epoch belong to?}
#'
#' \item{down}{Indicator for recumbent.}
#'
#' \item{wear.bout}{Integer index for period of consecutive wear epochs an epoch corresponds to.}
#'
#' \item{theta}{estimated inclination of the chest during the epoch.}
#'
#' \item{mad}{Mean absolute deviation, a measure of activity intensity during the epoch.}
#'
#' \item{x,y,z}{Median accelerations measured along each axis during an epoch.}
#'
#' \item{clustering}{Clustering algorithm to create posture groups.
#' Either "meanShift", "ward", or "centroid".}
#'
#' @export
process.zacl <- function(data, p = 0.95, k = 0.98, theta.star = 45, epoch.seconds = 60,
		 nonwear.window=3*94*60, nonwear.tol=10, minimum.wear.bout=94*60*24,
		 cluster = "meanShift"){
      ### Stop if epoch minutes is not a valid length
      if(!epoch.seconds %in% c(10, 30, 60, 300) & !epoch.seconds == "none"){
            stop("epoch.seconds must be 10, 30, 60, 300, or 'none'")
      }else{
      	    if(epoch.seconds == 300){
		          round.unit = lubridate::minutes(5)
	        }else round.unit = lubridate::seconds(epoch.seconds)
      }

      ### Identify non-wear (3 consecutive hours with fewer than 10 changes)
      data <- postuR::check.nonwear(data,
      							  filter = F,
      							  window = nonwear.window,
      							  tol=nonwear.tol,
      							  minimum.wear.bout = minimum.wear.bout)

      if(nrow(data) == 0){
      	stop("No wear bout is required length.")
      }

      ### FIND REMOVAL TIME POINT
      # if we deem device is device is removed within any wear bout
      # separate wear bout into two wearbout
      data <-
	      	data %>%
	      	dplyr::select(-bout.length) %>%
	      	dplyr::filter(wear.bout != 0) %>%
	      	tidyr::nest(data = c(time, x, y, z)) %>%
	      	dplyr::mutate(
	      		removal = purrr::map(data, ~postuR::calculate.removal.time(data = ., p=p)),
	      		data = purrr::map2(data, removal,
	      						   ~ dplyr::mutate(.x, wear.bout2 = (.y$r.ratio < k)* (time > .y$time)))) %>%
	      	dplyr::select(-removal) %>%
	      	tidyr::unnest(data) %>%
	      	dplyr::mutate(wear.bout = wear.bout + 0.5*wear.bout2) %>%
	      	dplyr::select(-wear.bout2) %>%
	      	dplyr::group_by(wear.bout)%>%
	      	dplyr::filter(n() > minimum.wear.bout) %>%
	      	dplyr::ungroup()

      if(nrow(data) == 0){
      	    stop("No wear bout is required length.")
      }

      ### EPOCH DATA

      min.data <-
      	data %>%
      	dplyr::mutate(time = lubridate::floor_date(time, unit = round.unit)) %>%
      	dplyr::group_by(time) %>%
      	dplyr::mutate(wear.bout = min(wear.bout)) %>%
      	dplyr::ungroup() %>%
      	dplyr::group_by(wear.bout) %>%
      	tidyr::nest() %>%
      	dplyr::mutate(
      	    # Estimate upright orientation
      	    top = purrr::map(
          		data,
          		~postuR:::find.top(.,p=p)),
      	    # Estimate Epoch-Specific Mean Absolute Deviation, Axis Medians
      	    min.data = purrr::map(
          		data,
          		~ dplyr::group_by(.,time) %>%
          			dplyr::summarise(
          				mad = mad(x,y,z),
          				x = median(x),
          				y = median(y),
          				z = median(z)) %>%
          			dplyr::mutate(
          				r = ifelse(postuR:::euclid.norm(x,y,z)==0,
          						   0.00001,
          						   postuR:::euclid.norm(x,y,z)),
          				x = x/r,
          				y = y/r,
          				z = z/r
          			) %>%
          			dplyr::arrange(time)),
            # Estimate Angle of inclination
            min.data = purrr::map2(
              	    min.data, top,
              	    ~dplyr::mutate(
              	        .x,
              	        theta = 180/pi*acos((x*.y["x"] + y*.y["y"] + z*.y["z"])/sqrt(sum(.y^2))),
              	        down0 = theta >= (theta.star))),
            # Add rotated data
            rdata = purrr::map2(
                min.data, top,
                ~postuR:::rotate.data(select(.x, x,y,z), from = .y, to = c(0,0,1)))
            )

      # derive mean shift clusters
      if(cluster=="meanShift"){
      min.data$cluster = purrr::map(min.data$min.data,
      		 ~ meanShiftR::meanShift(queryData = (.[,c("x","y","z")]) %>% as.matrix(),
      		 		bandwidth = c(0.14,0.14,0.14),
      		 		iter = 100)$assignment %>%
      		 	as.vector())
      }
      # centroid heirarchical clustering
      if(cluster=="centroid"){
      min.data$cluster = purrr::map(min.data$min.data,
      		  ~ cbind(.$x, .$y, .$z) %>%
      		  	fastcluster::hclust.vector(method = "centroid", metric = "euclidean") %>%
      		  	cutree(k = 7))
      }
      # ward heirarchical clustering
      if(cluster=="ward"){
      min.data$cluster =  purrr::map(min.data$min.data,
      		    ~ cbind(.$x, .$y, .$z) %>%
      		    	fastcluster::hclust.vector(method = "ward", metric = "euclidean") %>%
      		    	cutree(k = 5))
      }
       # adjust recumbent indicator using clusters
      min.data %>%
          dplyr::mutate(min.data = purrr::map2(
              min.data, cluster,
              ~ dplyr::mutate(.x, cluster = .y)%>%
                  dplyr::group_by(cluster) %>%
                  dplyr:: mutate(
                      p.down = mean(down0),
                      down = as.integer(p.down > .5 | down0 == 1)) %>%
                  dplyr::ungroup() %>%
                  dplyr::arrange(time))) %>%
          dplyr::select(-data, -top, -cluster)%>%
          tidyr::unnest(c(min.data, rdata, wear.bout)) %>%
          dplyr::select(-down0, -p.down)%>%
          dplyr::arrange(time)
}
