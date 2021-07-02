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


IHS <- function(y, theta, eps = 1e-4) {
    if (theta < eps)
        y else asinh(theta * y)/theta
}

library(kmmeans)
d <- read.table("inst/extdata/sdss-all-c.data", sep=",", na.strings="?")
colnames(d) <- c('id', 'class', 'brightness', 'texture', 'size', 'Me_1', 'Me_2')

d$brightness <- log10(d$brightness)
d$texture <- log10(d$texture)
d$size <- IHS(d$size, theta=20)
d$Me_1 <- IHS(d$Me_1, theta=20)
d$Me_2 <- IHS(d$Me_2, theta=20)

ret <- kmmeans(d[, 3:7], 2, 100, 10)
print(ret)
table(ret$partition, d[,2])
