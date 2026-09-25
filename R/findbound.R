#============================
#Find single stopping boundaries
#============================

#' Find stopping boundary via binomial or beta-binomial distribution. 
#' For the latter, the overdispersion factor (or design effect)
#'	by which the variance exceeds the regular binomial variance is printed.
#'
#' @param n total number of events
#' @param alpha_test nominal alpha for the binomial test
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param alternative direction of alternative, "less" or "greater"
#' @param icc intraclass correlation if there is more than one event per patient
#'
#' @return number of events in the experimental group that would lead to a stopping
#'
#' @export
#'
#'
#' @examples
#'	findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="greater")
#'	findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="less")
#'	
#'
findbound<-function(n, alpha_test = 0.025, pH0 = 0.5, alternative = "greater", icc = NULL) {

	if (is.null(icc)) {
		res<-findboundBin(n = n, alpha_test = alpha_test, pH0 = pH0, alternative = alternative)
	} else {
		theta <- (1/icc) - 1
		
		od<-1 + (n-1)*icc
		print(paste0("Overdispersion factor: ", paste(sprintf("%0.2f",od), collapse=", ")))
		
		a <- pH0 * theta
		b <- (1 - pH0) * theta
		res<-findboundBetaBinom(n = n, alpha_test = alpha_test, a = a, b = b, pH0 = pH0, alternative = alternative)
	}
	return(res)
}

#' Find stopping boundary via binomial distribution
#'
#' @param n total number of events
#' @param alpha_test nominal alpha for the binomial test
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param alternative direction of alternative, "less" or "greater"
#'
#' @return number of events in the experimental group that would lead to a stopping
#'
#' @noRd
#'
#' @importFrom stats qbinom
#'
findboundBin<-function(n, alpha_test = 0.025, pH0 = 0.5, alternative = "greater") {
 
	if (alternative=="greater") {
		xlim<-qbinom(p = 1-alpha_test, size = n, prob = pH0) + 1
	}
	if (alternative=="less") {
		xlim<-qbinom(p = alpha_test, size = n, prob = pH0) - 1
	}
	return(xlim)
}


#' findboundBetaBinom
#' Exact single-look stopping boundary under a beta-binomial null (or
#' alternative) -- the EXACT conditional distribution of arm allocation when
#' each arm's total event count is Negative Binomial with a shared dispersion
#' parameter. 
#'
#' @param n Vector of total event counts at which a boundary is needed.
#' @param alpha_test nominal alpha for the test.
#' @param a, b beta-binomial shape parameters
#' @param pH0 proportion of events in the experimental arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param alternative "greater" or "less".
#'
#' @return vector of stopping boundaries, one per element of \code{n}
#'
#' @noRd
#'
#' @importFrom stats pbinom
#'
findboundBetaBinom <- function(n, alpha_test=0.025, a, b, pH0 = NULL, alternative="greater") {
 
	## explicit degenerate-limit handling (icc=0, theta=alpha+b=Inf)
	is.degenerate <- is.infinite(a + b)
	if (is.degenerate && is.null(pH0)) {
		stop("a+b is infinite (icc=0 / no dispersion): supply 'pH0'",
			" explicitly so the exact plain-binomial limit can be used.\n\n")
	}
 
	ri <- lapply(n, function(x) {
		possible_successes <- 0:x
		cdf <- if (is.degenerate) {
			pbinom(possible_successes, x, pH0)
		} else {
			#VGAM::pbetabinom.ab(possible_successes, x, a, b)
			pbetabinom_ab_base(possible_successes, x, a, b)
		}
 
		if (alternative == "greater") {
			en <- min(possible_successes[cdf >= (1 - alpha_test)]) + 1
		}
		if (alternative == "less") {
			en <- max(possible_successes[cdf <= alpha_test])
		}
		return(en)
	})
 
	return(unlist(ri))
}

#' pbetabinom_ab_base
#' Density for a beta-binomal distribution, same as VGAM::pbetabinom.ab(q, size, a, b)
#' @param q vector of quantiles
#' @param size 	number of trials
#' @param shape1, shape2 beta-binomial shape parameters
#'
#' @return beta binomial distribution function
#'
#' @noRd
pbetabinom_ab_base <- function(q, size, shape1, shape2) {
  sapply(q, function(qi) {
    if (qi < 0) return(0)
    if (qi >= size) return(1)
    
    k <- 0:qi
    pmf <- exp(lchoose(size, k) + lbeta(k + shape1, size - k + shape2) - lbeta(shape1, shape2))
    sum(pmf)
  })
}


#============================
#Find study schedule
#============================

#' findsched
#' 
#' Calculate testing schedule based on binomial or beta-binomial distribution. 
#' For the latter, the overdispersion factor (or design effect)
#'	by which the variance exceeds the regular binomial variance is printed.
#'
#' @param nevents Vector of event counts at which a boundary is actually needed.
#' @param pH0 proportion of events in the treatment arm under the null.
#' @param alpha_test nominal (per-look) alpha.
#' @param icc intraclass correlation if there is more than one event per patient
#'
#' @return data.frame with columns totevents, treatBound, alphaLevelBound and cutoff.
#'
#' @noRd
#'
#'
findsched<-function(nevents, pH0 = 0.5, alpha_test = 0.025, icc = NULL) {
	
	if (is.null(icc)) {
		res<-findschedBin(nevents = nevents, alpha_test = alpha_test, pH0 = pH0)
	} else {
		
		od<-1 + (nevents-1)*icc
		print(paste0("Overdispersion factor: ", paste(sprintf("%0.2f",od), collapse=", ")))
		
		theta <- (1/icc) - 1
		a <- pH0 * theta
		b <- (1 - pH0) * theta
		res<-findschedBetaBinom(nevents = nevents, alpha_test = alpha_test, a = a, b = b, pH0 = pH0)
	}
	return(res)
}


#' findschedBin
#' 
#' Calculate testing schedule based on binomial distribution 
#'
#' @param nevents Vector of event counts at which a boundary is actually needed.
#' @param pH0 proportion of events in the treatment arm under the null.
#' @param alpha_test nominal (per-look) alpha.
#'
#' @return data.frame with columns totevents, treatBound, alphaLevelBound and cutoff.
#'
#' @noRd
#'
#' @importFrom stats dbinom
#'

findschedBin<-function(nevents, pH0 = 0.5, alpha_test = 0.025) {

	bounds <- data.frame(totevents=1:max(nevents), treatBound=NA,
		alphaLevelBound=NA, cutoff=NA)

	bound <- NULL
	for (j in 1:nrow(bounds)) {

		totevents <- bounds$totevents[j]

		if (!(totevents %in% nevents)) {
			alphaVal<-0
		} else {
			alphaVal<-alpha_test
		}

		## we don't need to do the next few steps unless alphaVal is > 0
		if (alphaVal <= 0) next

		## choose the lowerBound for searching for the next cutoff value.
		if (is.null(bound)) {
			lowerBnd <- ceiling(pH0 * totevents)
		} else {
			lowerBnd <- bound
		}

		valSeq <- totevents:lowerBnd

		upperTailProbs <- cumsum(dbinom(valSeq, totevents, pH0))
		signif <- (upperTailProbs <= alphaVal)

		## if we have at least one significant value then do...
		if (isTRUE(signif[1]))	{
			## get "largest" (last) index for which signif == TRUE
			largest.index <- max(which(signif))

			## define 'bound' to be the count corresponding to
			## the 'largest.index'
			## which we have significance at per-test-level 'alphaVal'
			bound <- valSeq[largest.index]

			bounds$treatBound[j] <- bound
			bounds$alphaLevelBound[j] <- upperTailProbs[largest.index]
			bounds$cutoff[j] <- alphaVal
		}

	}
	return(bounds)
}


#' findschedBetaBinom
#'
#' The exact Beta-Binomial / Pólya urn walk: one event at a time, only searching
#' for a boundary at the requested 'nevents' checkpoints
#'
#'
#' @param nevents Vector of event counts at which a boundary is actually needed.
#' @param alpha_test nominal (per-look) alpha.
#' @param a, b beta-binomial shape parameters under the null
#' @param pH0 proportion of events in the treatment arm under the null.
#'
#' @return data.frame with columns totevents, treatBound, alphaLevelBound,
#'	cutoff -- same structure as the original.
#'
#' @noRd
findschedBetaBinom <- function(nevents, alpha_test=0.025, a, b, pH0 = NULL){
 
	Nmax <- max(nevents)
 
	## same degenerate-limit handling as pNSBetaBinom(): icc=0 (theta=Inf)
	## means a,beta are themselves Inf, so their ratio can't be recovered
	## after the fact -- pH0 must be supplied explicitly for that case.
	is.degenerate <- is.infinite(a + b)
	if (is.degenerate && is.null(pH0)) {
		stop("a+b is infinite (icc=0 / no dispersion): supply 'pH0'",
			" explicitly so the exact plain-binomial limit can be used.\n\n")
	}
	pMove <- function(i, s) {
		if (is.degenerate) return(pH0)
		(a + s) / (a + b + i - 1)
	}
 
	## running marginal distribution of s, indexed 0..(current totevents)
	p1 <- pMove(1, 0)
	dist <- c("0" = 1 - p1, "1" = p1)
 
	## create data frame to store results in
	bounds <- data.frame(totevents=1:Nmax, treatBound=NA,
		alphaLevelBound=NA, cutoff=NA)
 
	bound <- NULL
 
	for (j in 1:nrow(bounds)) {
 
		totevents <- bounds$totevents[j]
 
		## advance the walk by one event (dist already holds totevents==1's state)
		if (totevents > 1) {
			i <- totevents
			new.dist <- structure(numeric(i+1), names = as.character(0:i))
			new.dist["0"] <- (1 - pMove(i, 0)) * dist["0"]
			for (s in 1:i) {
				S <- as.character(s); S.minus.1 <- as.character(s-1)
				if (s < i) {
					new.dist[S] <- pMove(i, s-1) * dist[S.minus.1] +
									(1 - pMove(i, s)) * dist[S]
				} else {
					new.dist[S] <- pMove(i, s-1) * dist[S.minus.1]
				}
			}
			dist <- new.dist
		}
 
		if (!(totevents %in% nevents)) next
 
		## choose the lowerBound for searching for the next cutoff value
		if (is.null(bound)) {
			meanFrac <- if (is.degenerate) pH0 else a / (a + b)
			lowerBnd <- ceiling(totevents * meanFrac)
		} else {
			lowerBnd <- bound
		}
 
		valSeq <- totevents:lowerBnd
		upperTailProbs <- cumsum(dist[as.character(valSeq)])
		signif <- (upperTailProbs <= alpha_test)
 
		## if we have at least one significant value then do...
		if (isTRUE(signif[1])) {
			largest.index <- max(which(signif))
			bound <- valSeq[largest.index]
 
			bounds$treatBound[j] <- bound
			bounds$alphaLevelBound[j] <- upperTailProbs[largest.index]
			bounds$cutoff[j] <- alpha_test
		}
	}
 
	return(bounds)
}


#============================
#Beta binomial Exact test
#============================


#' betabinom_test
#'
#' Exact test for a beta-binomial outcome
#'
#' @param x number of successes
#' @param n number of observations
#' @param icc intraclass correlation
#' @param pH0 null proportion
#' @param alternative, greater or less
#'
#' @return a p-value
#'
#' @noRd
betabinom_test <- function(x, n, icc, pH0 = 0.5, alternative = "greater") {
  
  if (alternative == "greater") {
    # Sum the probabilities from the observed count up to the maximum possible events
    possible_x <- x:n
    p_value <- sum(dbetabinom_native(possible_x, n, pH0, icc))
    
  } else if (alternative == "less") {
    # Sum the probabilities from 0 up to the observed count
    possible_x <- 0:x
    p_value <- sum(dbetabinom_native(possible_x, n, pH0, icc))
    
  } else {
    stop("Alternative must be either 'greater' or 'less'")
  }
  
  return(min(p_value, 1.0))
}


#' dbetabinom_native
#'
#' Helper function
#'
#' @param x number of events
#' @param size number of evnets 
#' @param prob null proportion
#' @param pron proportion of events in the treatment arm under the null.
#' @param icc intraclass correlation
#'
#' @return vector with probabilities
#
#' @noRd
#'
#' @importFrom stats dbinom
#'
dbetabinom_native <- function(x, size, prob, icc) {
  if (icc <= 0) return(dbinom(x, size, prob))
  
  theta <- (1 / icc) - 1
  a <- prob * theta
  b <- (1 - prob) * theta
  
  log_prob <- lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b)
  return(exp(log_prob))
}
