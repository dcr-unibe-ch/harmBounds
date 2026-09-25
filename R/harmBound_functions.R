# Functions for monitoring harm boundary
# --------------------------------------------------

#' Test-wise alpha necessary to control either 
#'	the family-wise type I error or the power at a specified level
#'
#'	Exactly one of alpha_total or power have to be specified.
#' The power requires the specification of an alternative via one of pH1, rrH1, orH1 or rdH1.
#'
#'	If there are several events per patient, the intraclass correlation coefficient has to be given
#'	and a beta-binomial framework is used. 
#'
#' @param nevents vector with number of events at which an interim analysis is done
#' @param alpha_total target family-wise type I error
#' @param power target power at the specified alternative
#' @param pH0 proportion of events in the treatment arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param alpha.interval Range for test-wise alpha, c(10^(-10),0.05) by default
#' @param maxevents optional maximum number of events expected for the trial (over both arms), 
#'	used to calculate the expected number of events
#' @param icc intraclass correlation if there is more than one event per patient
#' @param pH1 optional alternative, numeric vector, proportion of events in the treatment arm
#' @param rrH1 alternative specification of alternative as risk ratio (treatment / control)
#' @param orH1 alternative specification of alternative as odds ratio (treatment / control). Requires the control proportion (r0).
#' @param rdH1 alternative specification of alternative as risk difference (treatment - control). Requires the control proportion (r0).
#' @param r0 risk in the control arm. Required if the alternative is given as risk difference or odds ratio.
#' @param ... input for backward compatibility
#'
#' @return Test-wide alpha
#'
#' @export
#'
#' @importFrom stats uniroot
#' @importFrom utils capture.output
#'
#' @examples
#'	#Control overall family-wise type I error:
#'	apt<-getAlphaPerTest(nevents = c(10,50,100), alpha_total = 0.05, pH0 = 0.5)
#'	apt
#'	getHarmBound(nevents = c(10,50,100),alpha_test = apt, pH0 = 0.5)
#'
#'	#Control power assuming 80% of events in expermintal arm
#'	apt<-getAlphaPerTest(nevents = c(10,50,100), power = 0.8, pH0 = 0.5, pH1 = 0.6)
#'	apt
#'	getHarmBound(nevents = c(10,50,100),alpha_test = apt, pH0 = 0.5, pH1 = 0.6)
getAlphaPerTest <- function(nevents,
	alpha_total = NULL, power = NULL,
	pH0 = 0.5,
	alpha.interval = c(10^(-10), 1),
	maxevents = NULL,
	icc = NULL,
	pH1 = NULL, 
	rrH1 = NULL, orH1 = NULL, rdH1 = NULL,
	r0 = NULL, ...) {
	
	#totalAlpha backwards compatibility
	input<-list(...)
	if ("totalAlpha" %in% names(input) && is.null(alpha_total)) {
		warning("totalAlpha will be depreciated, please use alpha_total")	
		alpha_total<-input[["totalAlpha"]]
	}

	#either alpha_total or power
	nn <- sum(!is.null(alpha_total), !is.null(power))
	if (nn!=1) {
		stop("Either 'alpha_total' or 'power' must be specified.")
	}
	
	#check alternative if power 
	if (!is.null(power)) {
		nn<-sum(!is.null(pH1) | !is.null(rdH1) | !is.null(rrH1) | !is.null(orH1))
		if (nn!=1) {
			stop("Exactly one of pH1, rdH1, rrH1, or rrH1 must be entered with 'power'.")
		}
		if (length(pH1)>1 | length(rdH1)>1 | length(rrH1)>1 | length(orH1)>1) {
			warning("Only the first alternative is used to calibrate overall power")
		}
	}
	
    getCumAlpha <- function(alphaPerTest, nevents, pH0, ...) {
        capture.output(
			harmBounds <- getHarmBound1(
			nevents = nevents,
			alpha_test = alphaPerTest,
			pH0 = pH0,
			maxevents = maxevents,
			icc = icc,
			pH1 = pH1[1], 
			rrH1 = rrH1[1], orH1 = orH1[1], rdH1 = rdH1[1],
			r0 = r0))
			
		if (!is.null(alpha_total)) { 
			return(harmBounds$opchar[1,"cum_stop_prob"] - alpha_total)
		} else {
			return(harmBounds$opchar[2,"cum_stop_prob"] - power)
		}		
    }

	ur<-uniroot(getCumAlpha, interval = alpha.interval, tol=1e-7,
		nevents = nevents, 
		pH0 = pH0,
		maxevents = maxevents,
		icc = icc,
		pH1 = pH1, 
		rrH1 = rrH1, orH1 = orH1, rdH1 = rdH1,
		r0 = r0)
			
	return(ur$root)

}


#' Harm boundaries for safety testing
#'
#' Calculates the boundaries at each interim analysis, i.e. the number of events in the treatment arm
#' that would lead to a stopping of the trial based binomial exact tests,
#' assuming that the events should be equally distributed among both arms.
#' The indicated scenario (and all more extreme) 
#'	would lead to a rejection of H0 (equal distribution) and a stopping for safety.
#'
#'	If there are several events per patient, the intraclass correlation coefficient has to be given
#'	and a beta-binomial framework is used. The overdispersion factor (or design effect)
#'	by which the variance exceeds the regular binomial variance is printed. Note that the 
#'	effective sample size is reduced by that factor.
#'
#' The rejection region for the binomial or beta-binomial exact tests must be given for either 
#'		each test (alpha_test), overall (alpha_total, the family-wise error rate) or
#'		by the targetted power for the specified alternative.
#'
#' @param nevents vector with number of events (over both arms) at which an interim analysis is done
#' @param alpha_test the nominal alpha level to use for each test
#' @param alpha_total target overall family-wise type I error
#' @param power target power at the specified alternative
#' @param pH0 proportion of events in the treatment arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param maxevents optional maximum number of events expected for the trial (over both arms), 
#'	used to calculate the expected number of events
#' @param icc intraclass correlation if there is more than one event per patient
#' @param pH1 optional alternative, numeric vector, proportion of events in the treatment arm
#' @param rrH1 alternative specification of alternative as risk ratio (treatment / control)
#' @param orH1 alternative specification of alternative as risk ratio (treatment / control). Requires the control proportion (r0).
#' @param rdH1 alternative specification of alternative as risk difference (treatment - control). Requires the control proportion (r0).
#' @param r0 risk in the control arm. Required if the alternative is given as risk difference or odds ratio.
#' @return a list with 3 data.frames: bounds, stopprob and opchar.
#' bounds has a row for each interim analysis and columns for
#'	number of events (events),
#'	number of events in control and treatment arm that would lead to a stop
#'	(events_treatment, events_control), and the nominal alpha for each test (alpha_test).
#'  stopprob has a row for each interim analysis and columns for
#'	number of events (events),
#'	the hypothesis (pH),
#'	the stopping probability (stop_prob), and
#'	the cumulative stopping probability (cum_stop_prob)
#'	opchar has a row for each hypothesis (null plus each alternative) and columns 
#'	for the assumed proportion of events in the treatment arm (p),
#'	the cumulative stopping probabilities (cum_stop_prob) and 
#'	the expected total number of events (expected_events)
#'	for the null and each alternative.
#'
#'
#' @export
#'
#'
#' @seealso 
#' Available plot methods: \code{\link{plot.harmbound}}
#'
#' @importFrom stats dbinom qbinom
#'
#' @examples
#'
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5)
#' #adding an alternative
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, pH1=0.6)
#' 
#' #assume that a total of 150 events might occur (for the expected events)
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, pH1=0.6, maxevents=150)
#' 
#' #several alternatives
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5,
#'	pH1 = seq(0.6,0.8,by=0.05), maxevents=150)
#'
#' #using a risk ratio to specify the alternative
#' getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, rrH1=1.5, maxevents=150)
#'
#' # define the test so that an family-wise type I error of 5% is achieved
#'	getHarmBound(nevents=c(10,50,100), alpha_total=0.05, pH0=0.5)
#'
#' # define the test so that an over power of 80% is achieved
#' # needs an alternative
#'	getHarmBound(nevents=c(10,50,100), power=0.8, pH0=0.5, pH1=0.6)
getHarmBound <- function(nevents,
	alpha_test = NULL, alpha_total = NULL, power = NULL,
	pH0,
	maxevents=NULL,
	icc = NULL,
	pH1=NULL, 
	rrH1=NULL, orH1=NULL, rdH1=NULL,
	r0=NULL){
	
	#get alpha_test of not given:
	nn <- sum(!is.null(alpha_test) | !is.null(alpha_total) | !is.null(power))
	if (nn!=1) {
		stop("Exactly one of 'alpha_test', 'alpha_total' or 'power' must be specified.")
	}
	
	if (is.null(alpha_test)) {
		alpha_test<-getAlphaPerTest(nevents=nevents,
			alpha_total = alpha_total, power = power,
			pH0 = pH0,
			alpha.interval = c(10^(-10), 1),
			maxevents = maxevents,
			pH1 = pH1, 
			rrH1 = rrH1, orH1 = orH1, rdH1 = rdH1, r0 = r0)
	}
	
	res<-getHarmBound1(nevents = nevents,
		alpha_test = alpha_test,
		pH0 = pH0,
		maxevents = maxevents,
		icc = icc,
		pH1 = pH1, 
		rrH1 = rrH1, orH1 = orH1, rdH1 = rdH1, r0 = r0)
		
	return(res)
	
}

#' Harm boundaries for safety testing (main function)
#'
#' @param nevents vector with number of events (over both arms) at which an interim analysis is done
#' @param alpha_test the nominal alpha level to use for each test
#' @param pH0 proportion of events in the treatment arm under the null hypothesis,
#'	typically based on randomization ratio (e.g. 0.5 for a 1:1 randomization)
#' @param maxevents optional maximum number of events expected for the trial (over both arms), used to calculate the expected number of events
#' @param icc Intraclass correlation if there is more than one event per patient
#' @param pH1 optional alternative, numeric vector, proportion of events in the treatment arm
#' @param rrH1 alternative specification of alternative as risk ratio (treatment / control)
#' @param orH1 alternative specification of alternative as risk ratio (treatment / control). Requires the control proportion (r0).
#' @param rdH1 alternative specification of alternative as risk difference (treatment - control). Requires the control proportion (r0).
#' @param r0 risk in the control arm. Required if the alternative is given as risk difference or odds ratio.
#' @return a list with 3 data.frames: bounds, stopprob and opchar.
#' bounds has a row for each interim analysis and columns for
#'	number of events (events),
#'	number of events in control and treatment amrs that would lead to a stop
#'	(events_treatment, events_control), and the nominal alpha for each test (alpha_test).
#'  stopprob has a row for each interim analysis and columns for
#'	number of events (events),
#'	the hypothesis (pH),
#'	the stopping probability (stop_prob), and
#'	the cumulative stopping probability (cum_stop_prob)
#'	opchar has a row for each hypothesis (null plus each alternative) and columns 
#'	for the assumed proportion of events in the treatment arm (p),
#'	the cumulative stopping probabilities (cum_stop_prob) and 
#'	the expected total number of events (expected_events)
#'	for the null and each alternative.
#'
#'
#' @importFrom stats dbinom
#'
#' @noRd
#'
getHarmBound1 <- function(nevents, alpha_test, pH0,
	maxevents=NULL,
	icc = NULL,
	pH1 = NULL, 
	rrH1 = NULL, orH1 = NULL, rdH1 = NULL,
	r0 = NULL){

	#check alternative
	nn<-sum(!is.null(pH1) | !is.null(rdH1) | !is.null(rrH1) | !is.null(orH1))
	if (nn>1) {
		stop("Only one of pH1, rdH1, rrH1, or rrH1 should be entered.")
	}
	
	#check maxevents 
	if (!is.null(maxevents)) {
		if (maxevents<max(nevents)) {
			stop("maxevents has to be larger or equal than the maximum of nevents")
		}
	} else {
		maxevents<-max(nevents)
	}
	
	if (!is.null(rdH1) | !is.null(rrH1) | !is.null(orH1)) {
		
		if (!is.null(rrH1)) {
			r0<-0.5
		} else {
			if (is.null(r0)) {
				stop("r0 has to be given if the effect is given as 'orH1' or 'rdH1'.")
			}
		}
		pH1<-convertRisks(rd=rdH1, rr=rrH1, or=orH1, r0=r0, n0=(1-pH0), n1=pH0)[,"eprop"]
		pH1<-as.numeric(pH1)
	}
	
	# create boundaries
	bounds<-findsched(nevents = nevents, alpha_test =alpha_test, pH0 = pH0, icc = icc)

	#addumptions
	assumpt<-c(pH0,pH1)
	
	stopprob<-vector(length=1+length(pH1),mode="list")
	names(stopprob)<-assumpt
	cumstop<-numeric(0)
	
	i<-1
	for (i in 1:length(assumpt)) {
		
		hyp<-ifelse(i==1,"H0","H1")
					
		if (!all(is.na(bounds$treatBound))) {	
			out <- pNS(Bounds=bounds$treatBound, pH=assumpt[i], icc = icc)
			nstop<-sum(out$Bounds$n*out$Stop,(1-out$totalStopProb)*maxevents)
			outc<-data.frame(p=assumpt[i],cum_stop_prob=out$totalStopProb,
				expected_events=nstop,hyp=hyp)			
		} else {
			nr <- nrow(bounds)
			out <- list(Bounds= data.frame(n=1:nr,StoppingBound = bounds$treatBound),
				Stop = rep(0, nr),cutoff = bounds$cutoff)
			outc<-data.frame(p=pH0,cum_stop_prob=0,expected_events=maxevents,hyp=hyp)	
		}		
		
		#Boundaries
		if (i==1) {
			boundOut<-out$Bounds
			names(boundOut)[names(boundOut)=="n"] <- "events"
			names(boundOut)[names(boundOut)=="StoppingBound"] <- "events_treatment"
			boundOut$events_control <- boundOut$events - boundOut$events_treatment		
			boundOut <- cbind(boundOut,alpha_test = bounds$cutoff)
			
			stopifnot(is.na(boundOut[!boundOut$events %in% nevents,"events_treatment"]))
			boundOutna<-boundOut[boundOut$events %in% nevents, ]
			rownames(boundOutna)<-1:nrow(boundOutna)
		}
		
		#Stopping probs 
		sprobi <- data.frame(events=boundOut$events,
			pH = assumpt[i],
			hyp = hyp,
			stop_prob = out$Stop,
			cum_stop_prob=cumsum(out$Stop))
		
		stopifnot(sprobi[!sprobi$events %in% nevents,"stop_prob"]==0)
		sprobi<-sprobi[sprobi$events %in% nevents, ]
		rownames(sprobi)<-1:nrow(sprobi)
		
		stopprob[[i]]<-sprobi
		
		#summary
		if (!is.null(rrH1)) {
			outc<-cbind(outc, rr=c(1,rrH1)[i])
			outc<-outc |> dplyr::relocate(.data$rr, .after = .data$p)
		}
		if (!is.null(orH1)) {
			outc<-cbind(outc,or=c(1,orH1)[i],r0 = r0)
			outc<-outc |> 
				dplyr::relocate(.data$or, .after = .data$p) |>
				dplyr::relocate(.data$r0, .after = .data$or)
		}
		if (!is.null(rdH1)) {
			outc<-cbind(outc,rd=c(1,rdH1)[i],r0=r0)
			outc<-outc |> 
				dplyr::relocate(.data$rd, .after = .data$p) |>
				dplyr::relocate(.data$r0, .after = .data$rd)
		}
			
		cumstop<-rbind(cumstop,outc)
	
	}
	
	#combine
		
	res<-list(boundOutna,stopprob,cumstop)
	names(res)<-c("bounds","stopprob","opchar")

	class(res) <- c("harmbound", class(res))

	return(res)
}



#' Convert the proportion of events in the treatment arms to risk differences and ratios (and vice versa)
#'
#' @param eprop proportion of events in treatment arm
#' @param etotal total number of events
#' @param rd risk difference
#' @param rr risk ratio
#' @param or odds ratio
#' @param r0 risk in the control arm
#' @param n0 number of patients in the control arm
#' @param n1 number of patients in the treatment arm
#'
#' @return vector with risks in control and treatment arm (r0, r1),
#'	the risk difference (rd), risk ratio (rr) and odds ratio (or)
#'
#' @export
#'
#' @examples
#' convertRisks(eprop=0.5,etotal=100,n0=200)
#' convertRisks(eprop=0.6,etotal=100,n0=200)
#' convertRisks(rr=1.5,n0=200,r0=0.2)
#'
convertRisks<-function(eprop=NULL,etotal=NULL,
	rd=NULL,rr=NULL,or=NULL,
	r0=NULL,
	n0,n1=n0) {

	nn<-sum(!is.null(eprop) | !is.null(rd) | !is.null(rr) | !is.null(or))
	if (nn!=1) {
		stop("Only one of eprop, rd, rr or or should be given")
	}

	if (!is.null(eprop)) {

		if (is.null(etotal)) {
			stop("etotal has to be given for conversion to risk difference or ratios.")
		}

		r0<-(etotal-eprop*etotal)/n0
		r1<-eprop*etotal/n1
		rd<-r1-r0
		rr<-r1/r0
		or<-r1/(1-r1)/(r0/(1-r0))
	}

	if (!is.null(rd) | !is.null(rr) | !is.null(or)) {

		if (is.null(r0)) {
			stop("r0 has to be given for conversion to proportion of events in treatment arm.")
		}

		if (!is.null(rd)) {
			r1 <- rd + r0
		}
		if (!is.null(rr)) {
			r1 <- rr * r0
		}
		if (!is.null(or)) {
			r1 <- (or * r0/(1-r0)) / (1 + or * r0/(1-r0))
		}

		eprop <- r1*n1 / (r0*n0 + r1*n1)

		rd<-r1-r0
		rr<-r1/r0
		or<-r1/(1-r1)/(r0/(1-r0))
		etotal<-r1*n1/eprop
	}
	
	res<-cbind(eprop,etotal,n0,n1,r0,r1,rd,rr,or)
	colnames(res)<-c("eprop","etotal","n0","n1","r0","r1","rd","rr","or")
	return(res)
}

