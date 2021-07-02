################################################################################
# @file kmmeans.R
#
# Purpose: R wrapper for kmmeans algorithm.
################################################################################

#' K-means clustering with missing data
#'
#' Estimate means and partition of numerical data with missing values.
#'
#' @param data		A data frame of data with possible missing values (NA).
#' @param K		Number of clusters to fit.
#' @param n.init	Number of k-means++ initializations.
#' @param kmmns.iter     Number of kmmeans iterations.
#'
#' @return A list with k-means clustering solution.
#' \itemize{
#'	\item criterion - The minimum realized objective function value.
#'	\item cluster.criteria - A vector of the minimum realized objective per
#'	cluster.
#'	\item cluster.sizes - The sizes of the clusters in the best solution.
#'	\item partition - The partition of observations into clusters in the best
#'	solution.
#'	\item modes - The estimated modes of each cluster in K x p integer matrix.
#'	\item average.criterion - The average objective function achieved across
#'	multiple initializations.
#'	\item number.initializations - The number of initializations performed.
#'	\item best.ari - If true clusters provided, then the best achieved ARI,
#'	 perhaps not from the best solution, as measured by the objective function.
#'	\item average.ari - If true clusters provided, then the average ARI achieved
#'	across all initializations.
#' }
#'
#' @references
#'
#' @examples
#'
#' @export
#' @useDynLib kmmeans

kmmeans <- function(
                    data,
                    K,
                    n.init = 1,
                    kmmns.iter=10
                    )
{
	stopifnot(is.data.frame(data))
	stopifnot(!is.na(as.integer(K) && as.integer(K)>0) & as.integer(K)>0)
	stopifnot(!is.na(as.integer(n.init)) & as.integer(n.init)>0)
        stopifnot(!is.na(as.integer(kmmns.iter)) & as.integer(kmmns.iter)>=0)

	if (!is.loaded("kmmeanspp_r", PACKAGE="kmmeans"))
		dyn.load("src/kmmeans.so")

	out <- .Call("kmmeanspp_r", as.matrix(data), as.integer(K), as.integer(n.init), as.integer(kmmns.iter))

	return(out)
}# kmmeans
