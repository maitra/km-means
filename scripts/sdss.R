#' SDSS
#'
#' Sloan ...
#'
#' @details
#' A dataset from the ...
#' The first column is a unique object id. The second column is the class
#' (3=galaxy, 6=star). The remaining columns are brightness (psfCounts),
#' texture, size (petroRad), shape (M_e1), shape (M_e2)
#'
#' @name sdss
#
#' @docType data
#'
#' @format A data frame with 1507 observations and 7 columns.
#
#' @keywords datasets

ihs <- function(x, theta=0)
{
	if (abs(theta) < 1e-3) {
		x
	} else {
		asinh(theta * x) / theta
	}
}

library(kmmeans)
d <- read.table("inst/extdata/sdss-all-c.data", sep=",", na.strings="?")
colnames(d) <- c('id', 'class', 'brightness', 'texture', 'size', 'shape1', 'shape2')
d$brightness <- log10(d$brightness)
d$texture <- log10(d$texture)
d$size <- ihs(d$size, theta=10)
d$shape1 <- ihs(d$shape1, theta=10)
d$shape2 <- ihs(d$shape2, theta=10)
ret <- kmmeans(d[, 3:7], 2, 1)
table(ret$partition, d[,2])
