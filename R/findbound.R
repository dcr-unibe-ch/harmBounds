#' Find stopping boundary via binomial exact tests
#'
#' @param n total number of events
#' @param alpha_test nominal alpha for the binomial test
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param alternative direction of alternative, "less" or "greater"
#'
#' @return number of events in the experimental group that would lead to a stopping
#'
#' @export
#'
#' @importFrom stats binom.test
#'
#' @examples
#'	findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="greater")
#'	findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="less")
#'	
#'
findbound<-function(n, alpha_test=0.025, pH0 = 0.5, alternative="greater") {
 
	if (alternative=="greater") {
		xlim<-qbinom(p = 1-alpha_test, size = n, prob = pH0) + 1
	}
	if (alternative=="less") {
		xlim<-qbinom(p = alpha_test, size = n, prob = pH0) - 1
	}
	return(xlim)
}

