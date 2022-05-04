#' @title Euclidian Norm, Radius
#'
#' @description Calculate length of a 3D vector.
#' @param x First accelerometer channel.
#' @param y Second accelerometer channel.
#' @param z Third accelerometer channel.
#' @param v optional matrix or data frame containing x,y,z measurments.

euclid.norm <- function(x,y,z, v = cbind(x,y,z)){
	apply(v, 1, FUN = function(a) sqrt(sum(a^2)))
}
