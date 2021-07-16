#' @title Rotate Triaxial Accelerometer Data
#'
#' @description Takes a data frame or data matrix with three columns, and internally generates a rotation matrix that will rotate the vector "from" to "to". The vector space orthoganol to these vectors (their cross-product) are not rotated--they become the axis of rotation.
#' @param data Data frame or matrix containing the accelerometer data. Should include ONLY variables x, y, z.
#' @param from Vector defining the starting point for the rotation.
#' @param to Vector defining the ending point for the rotation.

rotate.data <- function(data, from, to){

	# check dimensions of data and vectors
	if(ncol(data)!=3){
		cat("Data must have 3 columns.")
		return(NULL)
	}
	from <- as.vector(from)
	to <- as.vector(to)
	if(length(from)!=3 | length(to) !=3){
		cat("from and to must have length 3.")
		return(NULL)
	}

	# calculate rotation matrix
	v = pracma::cross(from,to)
	s = sum(v^2)^.5
	c = sum(from*to)
	X = c(0, -v[3], v[2],
	      v[3], 0, -v[1],
	      -v[2], v[1], 0) %>%
		matrix(nrow = 3, byrow=T)
	R = diag(3) + X + X %*% X /(1+c)

	#calculate rotation matrix and rotate
	rotated.data <- (as.matrix(data) %*% t(R)) %>% dplyr::as_tibble(.name_repair = "minimal")
	colnames(rotated.data) <- c("rx", "ry","rz")
	return(rotated.data)
}
