
#' @title Estimate upright position from accelerometry data
#' @description Takes a data frame with variables x,y,z and returns a vector of length 3 indicating orientation of upright position.
#' @param data a data frame with variables x,y,z

find.top <- function(data, p = 0.95){
	data %>%
		dplyr::mutate(r = sqrt(x^2 + y^2 + z^2)) %>%
		dplyr::filter(r!=0) %>%
		dplyr::mutate(x = x/r, y=y/r, z=z/r,
		       r2 = (1-r)^2) %>%
		dplyr::filter(r2 > quantile(r2, p, na.rm=T)) %>%
		dplyr::summarise(x = mean(x, na.rm=T), y=mean(y, na.rm=T), z = mean(z, na.rm=T)) %>%
		unlist
}
