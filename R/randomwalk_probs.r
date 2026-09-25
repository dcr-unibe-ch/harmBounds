
#' pNS
#' Helper function which creates an object containing the values of the function P(n,s)
#' defined by Breslow (1970, JASA) as, for 0 <= s <= n <= N,  the probability
#' of the binomial random walk S_n reaching S_n = s without "absorption" into
#' the rejection region (i.e. without hitting any of the stopping bounds).
#'
#'	Extended to corralated data if q is not NULL.
#'
#' @param Bounds Vector of stopping bounds after each event.
#'	For event totals where stopping is not permitted the 'Bound' should be set to NA.
#' @param pH assumed proportion of events in the treatment arm
#' @param returnPns Whether to inlcude the binomial random walk in the output
#' @param icc intraclass correlation if there is more than one event per patient
#'
#' @noRd
#'
#' @return list with stopping probabilities and boundaries
#'
#' @examples
#'
#'	pNS(Bounds=c(rep(NA,5),3), pH=0.5, returnPns=TRUE)
#'
pNS <- function(Bounds, pH = 0.5, icc = NULL, returnPns=FALSE){
    if (is.null(icc)) {
        pNSBin(Bounds = Bounds, pH = pH, returnPns = returnPns)
    } else {
		theta <- (1/icc) - 1
		a = pH * theta
		b = (1 - pH) * theta
        pNSBetaBinom(Bounds = Bounds, a = a, b = b, pH = pH, returnPns = returnPns)  
    }
}


#' pNSBin
#' Helper function which creates an object containing the values of the function P(n,s)
#' defined by Breslow (1970, JASA) as, for 0 <= s <= n <= N,  the probability
#' of the binomial random walk S_n reaching S_n = s without "absorption" into
#' the rejection region (i.e. without hitting any of the stopping bounds).
#'
#' @param Bounds Vector of stopping bounds after each event.
#'	For event totals where stopping is not permitted the 'Bound' should be set to NA.
#' @param pH assumed proportion of events in the treatment arm
#' @param returnPns Whether to inlcude the binomial random walk in the output
#'
#' @noRd
#'
#' @return list with stopping probabilities and boundaries
#'
#' @examples
#'
#'	pNSBin(Bounds=c(rep(NA,5),3), pH=0.5, returnPns=TRUE)
#'
pNSBin <- function(Bounds, pH=0.5, returnPns=FALSE){

	N <- length(Bounds)

	if (all(is.na(Bounds))) {
		stop("The vector provided for argument 'Bounds' contains only NAs.\n",
			"Exiting...\n\n")
	}

	if (any(Bounds > (1:N), na.rm=TRUE)) {
		stop("The bounds provided do not appear to correspond to the",
			" total number of events (1,2,...,N).\n", "One of more of the bounds are larger",
			"than their corresponding total.  Exiting...\n\n")
	}

	## create 'Pns' and 'Stop'. Note that Stop has initial values of 0.
	Stop <- numeric(N)
	Pns  <- vector("list", length=N)

	## 'first.bound' is the point at which stopping is first allowed
	first.bound <- which(!is.na(Bounds))[1]

	if (first.bound == 1)
	stop("Cannot stop at the first test. Why would you want to?\n\n")


	## Initialize the base level (i==1).
	Pns[[1]] <- c("0" = (1-pH), "1" = pH)

	## loop over remaining infection totals, building up the prob.s as we go
	for (i in 2:N) {

	Pns[[i]] <- structure(numeric(i+1), names= as.character(0:i))

	max.S <- ifelse(is.na(Bounds[i]), i, Bounds[i]-1)

	Pns[[i]]["0"] <- (1-pH)*Pns[[i-1]]["0"]

	for (s in 1:max.S ) {
		S         <- as.character(s)
		S.minus.1 <- as.character(s-1)

		if (s < i) {
			Pns[[i]][S] <- pH*Pns[[i-1]][S.minus.1] + (1-pH)*Pns[[i-1]][S]
		} else {
			Pns[[i]][S] <- pH*Pns[[i-1]][S.minus.1]
		}
	}
	#all probs outside the boundary are 0->trials with are stopped
	#sum(Pns[[]]) < 1 after first stop

	## 'Stop' is the prob. that we encounter a stopping bound for the first
	## time at infection count 'i'.  Only computed if stopping is possible,
	## else remains at initialized value of 0.
	if (!is.na(Bounds[i])) {

		## If this is the first infection count at which we allow stopping, then
		## the stop value must be computed allowing for people using non-standard
		## bounds (i.e. ones that that *don't* begin with bound[n] = n).
		if ((i == first.bound | is.na(Bounds[i-1])) && Bounds[i] < i ) {
			## stop prob. =
			#	having one event less in the previous round * prob of hving a further event plus
			#	no has already nee
			## that we were already at/above the number of "success" needed to stop
			## at the previous infection total (i.e. i-1) (for which we should have
			## been given a stopping bound, but weren't).
			Stop[i] <- pH*Pns[[i-1]][as.character(Bounds[i]-1)] +
						sum(Pns[[i-1]][as.character(Bounds[i]:(i-1))])
			} else {
			## if 'i' is not the the first total at which stopping is allowed, or
			## if it is, but we have a proper bound (i.e. Bounds[i] = i) then we do
			## we compute as this:
			Stop[i] <- pH*Pns[[i-1]][as.character(Bounds[i]-1)]
			}
		}
	}

	outObj <- list(
		totalStopProb = sum(Stop),
		Stop = Stop,
		Bounds  = data.frame(n=1:N, StoppingBound=Bounds),
		N = N,
		pH = pH )

	if (isTRUE(returnPns)) {
		outObj$Pns <- Pns
	}

	return(outObj)
}

#' pNSBetaBinom
#' The exact conditional distribution of arm allocation when each
#' arm's total event count is Negative Binomial with a shared dispersion
#' parameter (e.g. a Poisson-rate-per-arm model with patient-level frailty --
#' see package notes on rate vs. occurrence effects).
#'
#' Key fact used here: if X ~ NB(r1, p) and Y ~ NB(r2, p) independently (same
#' p, i.e. shared dispersion), then X | (X+Y=N) is EXACTLY Beta-Binomial(N,
#' r1, r2). Equivalently, a Beta-Binomial sequence is a PÓLYA URN: the
#' probability that event i is an experimental-arm event, given s
#' experimental-arm events among the previous i-1, is exactly
#' (a+s)/(a+b+i-1) -- a simple, state-dependent replacement for
#' the constant pH in the original pNS() recursion. No new state dimensions
#' are needed; this is the SAME O(N^2) complexity as the original pNS().
#'
#' As a, b -> infinity with a/(a+b) fixed at pH, this
#' reduces exactly to the original pNS() (see examples).
#'
#' @param Bounds Vector of stopping bounds after each event. NA where
#'	stopping is not permitted at that event count.
#' @param a,b Shape parameters of the Beta-Binomial null (or
#'	alternative) distribution.
#' @param pH assumed proportion of events in the treatment arm
#' @param returnPns Whether to include the underlying random walk in the output.
#'
#' @noRd
#'
#' @return list with stopping probabilities and boundaries
#'
#' @examples
#'
#'	## Under the null (equal rates), moderate dispersion
#'	pNSBetaBinom(Bounds=c(rep(NA,5),3), a=5, b=5)
#'
#'	## Large a,b with fixed ratio reduces to the original pNS(pH=0.5)
#'	pNSBetaBinom(Bounds=c(rep(NA,5),3), a=5000, b=5000)
#'
pNSBetaBinom <- function(Bounds, a, b, pH=NULL, returnPns=FALSE){

	N <- length(Bounds)

	if (all(is.na(Bounds))) {
		stop("The vector provided for argument 'Bounds' contains only NAs.\n",
			"Exiting...\n\n")
	}
	if (any(Bounds > (1:N), na.rm=TRUE)) {
		stop("The bounds provided do not appear to correspond to the",
			" total number of events (1,2,...,N).  Exiting...\n\n")
	}

	Stop <- numeric(N)
	Pns  <- vector("list", length=N)

	first.bound <- which(!is.na(Bounds))[1]
	if (first.bound == 1)
		stop("Cannot stop at the first test. Why would you want to?\n\n")

	## Pólya urn move probability: given s experimental-arm events among the
	## previous (i-1) draws, probability event i is experimental-arm.
	## SPECIAL CASE: as theta=a+b -> Inf (icc -> 0, no dispersion), the
	## exact limit is the CONSTANT pH = a/(a+b) -- but a,b
	## being literally Inf makes that ratio itself Inf/Inf = NaN in R, so the
	## limit must be taken explicitly (pH supplied) rather than computed from
	## already-infinite a,b.
	is.degenerate <- is.infinite(a + b)
	if (is.degenerate && is.null(pH)) {
		stop("a+b is infinite (icc=0 / no dispersion): supply 'pH'",
			" explicitly so the exact plain-binomial limit can be used.\n\n")
	}
	pMove <- function(i, s) {
		if (is.degenerate) return(pH)
		(a + s) / (a + b + i - 1)
	}

	## base level (i=1)
	p1 <- pMove(1, 0)
	Pns[[1]] <- c("0" = 1 - p1, "1" = p1)

	for (i in 2:N) {

		Pns[[i]] <- structure(numeric(i+1), names = as.character(0:i))

		max.S <- ifelse(is.na(Bounds[i]), i, Bounds[i]-1)

		Pns[[i]]["0"] <- (1 - pMove(i, 0)) * Pns[[i-1]]["0"]

		for (s in 1:max.S) {
			S         <- as.character(s)
			S.minus.1 <- as.character(s-1)

			if (s < i) {
				Pns[[i]][S] <- pMove(i, s-1) * Pns[[i-1]][S.minus.1] +
								(1 - pMove(i, s)) * Pns[[i-1]][S]
			} else {
				Pns[[i]][S] <- pMove(i, s-1) * Pns[[i-1]][S.minus.1]
			}
		}

		if (!is.na(Bounds[i])) {
			if ((i == first.bound | is.na(Bounds[i-1])) && Bounds[i] < i) {
				Stop[i] <- pMove(i, Bounds[i]-1) * Pns[[i-1]][as.character(Bounds[i]-1)] +
							sum(Pns[[i-1]][as.character(Bounds[i]:(i-1))])
			} else {
				Stop[i] <- pMove(i, Bounds[i]-1) * Pns[[i-1]][as.character(Bounds[i]-1)]
			}
		}
	}

	outObj <- list(
		totalStopProb = sum(Stop),
		Stop = Stop,
		Bounds = data.frame(n=1:N, StoppingBound=Bounds),
		N = N,
		a = a,
		b = b )

	if (isTRUE(returnPns)) {
		outObj$Pns <- Pns
	}

	return(outObj)
}