# Harm boundaries for safety testing

Calculates the boundaries at each interim analysis, i.e. the number of
events in the intervention group that would lead to a stopping of the
trial based on a binomial exact test, assuming that the events should be
equally distributed amont both groups. The indicated scenario (and all
more extreme) would lead to a rejection of H0 (equal distribution) and a
stopping for safety.

## Usage

``` r
getHarmBound(
  nevents,
  alpha_test,
  pH0,
  maxevents = NULL,
  pH1 = NULL,
  rrH1 = NULL,
  orH1 = NULL,
  rdH1 = NULL,
  r0 = NULL,
  n = NULL
)
```

## Arguments

- nevents:

  vector with number of events (over both arms) at which an interim
  analysis is done

- alpha_test:

  the nominal alpha level to use for each test

- pH0:

  proportion of events in the intervention arm under the null
  hypothesis, typically based on randomization ratio (e.g. 0.5 for a 1:1
  randomization)

- maxevents:

  optional maximum number of events expected for the trial (over both
  arms), used to calculate the expected number of events

- pH1:

  optional alternative, numeric vector, proportion of events in the
  intervention arm

- rrH1:

  alternative specification of alternative as risk ratio (intervention /
  control)

- orH1:

  alternative specification of alternative as risk ratio (intervention /
  control). Requires the control proportion (r0).

- rdH1:

  alternative specification of alternative as risk difference
  (intervention - control). Requires the control proportion (r0) and the
  number of participants (n).

- r0:

  risk in the control group. Required if the alternative is given as
  risk difference or odds ratio.

- n:

  total number of participants. Required if the alternative is given as
  risk difference.

## Value

a list with 3 data.frames: bounds, stopprob and opchar. bounds has a row
for each interim analysis and columns for number of events (events),
number of events in control and intervention group that would lead to a
stop (events_intervention, events_control), and the nominal alpha for
each test (alpha_test). stopprob has a row for each interim analysis and
columns for number of events (events), the hypothesis (pH), the stopping
probability (stop_prob), and the cumulative stopping probability
(cum_stop_prob) opchar has a row for each hypothesis (null plus each
alternative) and columns for the assumed proportion of events in the
intervention group (p), the cumulative stopping probabilities
(cum_stop_prob) and the expected total number of events
(expected_events) for the null and each alternative.

## Examples

``` r
getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5)
#> $bounds
#>   events events_intervention events_control alpha_test
#> 1     10                   9              1      0.025
#> 2     50                  33             17      0.025
#> 3    100                  61             39      0.025
#> 
#> $stopprob
#> $stopprob$`0.5`
#>   events  pH hyp  stop_prob cum_stop_prob
#> 1     10 0.5  H0 0.01074219    0.01074219
#> 2     50 0.5  H0 0.01490030    0.02564249
#> 3    100 0.5  H0 0.01194192    0.03758441
#> 
#> 
#> $opchar
#>     p cum_stop_prob expected_events hyp
#> 1 0.5    0.03758441        98.28819  H0
#> 
#> attr(,"class")
#> [1] "harmbound" "list"     
#adding an alternative
getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, pH1=0.6)
#> $bounds
#>   events events_intervention events_control alpha_test
#> 1     10                   9              1      0.025
#> 2     50                  33             17      0.025
#> 3    100                  61             39      0.025
#> 
#> $stopprob
#> $stopprob$`0.5`
#>   events  pH hyp  stop_prob cum_stop_prob
#> 1     10 0.5  H0 0.01074219    0.01074219
#> 2     50 0.5  H0 0.01490030    0.02564249
#> 3    100 0.5  H0 0.01194192    0.03758441
#> 
#> $stopprob$`0.6`
#>   events  pH hyp stop_prob cum_stop_prob
#> 1     10 0.6  H1 0.0463574     0.0463574
#> 2     50 0.6  H1 0.2098108     0.2561682
#> 3    100 0.6  H1 0.2505214     0.5066895
#> 
#> 
#> $opchar
#>     p cum_stop_prob expected_events hyp
#> 1 0.5    0.03758441        98.28819  H0
#> 2 0.6    0.50668952        85.33730  H1
#> 
#> attr(,"class")
#> [1] "harmbound" "list"     
#assume that a total of 150 events might occur
getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, pH1=0.6, maxevents=150)
#> $bounds
#>   events events_intervention events_control alpha_test
#> 1     10                   9              1      0.025
#> 2     50                  33             17      0.025
#> 3    100                  61             39      0.025
#> 
#> $stopprob
#> $stopprob$`0.5`
#>   events  pH hyp  stop_prob cum_stop_prob
#> 1     10 0.5  H0 0.01074219    0.01074219
#> 2     50 0.5  H0 0.01490030    0.02564249
#> 3    100 0.5  H0 0.01194192    0.03758441
#> 
#> $stopprob$`0.6`
#>   events  pH hyp stop_prob cum_stop_prob
#> 1     10 0.6  H1 0.0463574     0.0463574
#> 2     50 0.6  H1 0.2098108     0.2561682
#> 3    100 0.6  H1 0.2505214     0.5066895
#> 
#> 
#> $opchar
#>     p cum_stop_prob expected_events hyp
#> 1 0.5    0.03758441        146.4090  H0
#> 2 0.6    0.50668952        110.0028  H1
#> 
#> attr(,"class")
#> [1] "harmbound" "list"     
#or several alternatives
getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5,
pH1 = seq(0.6,0.8,by=0.05), maxevents=150)
#> $bounds
#>   events events_intervention events_control alpha_test
#> 1     10                   9              1      0.025
#> 2     50                  33             17      0.025
#> 3    100                  61             39      0.025
#> 
#> $stopprob
#> $stopprob$`0.5`
#>   events  pH hyp  stop_prob cum_stop_prob
#> 1     10 0.5  H0 0.01074219    0.01074219
#> 2     50 0.5  H0 0.01490030    0.02564249
#> 3    100 0.5  H0 0.01194192    0.03758441
#> 
#> $stopprob$`0.6`
#>   events  pH hyp stop_prob cum_stop_prob
#> 1     10 0.6  H1 0.0463574     0.0463574
#> 2     50 0.6  H1 0.2098108     0.2561682
#> 3    100 0.6  H1 0.2505214     0.5066895
#> 
#> $stopprob$`0.65`
#>   events   pH hyp  stop_prob cum_stop_prob
#> 1     10 0.65  H1 0.08595444    0.08595444
#> 2     50 0.65  H1 0.43634566    0.52230009
#> 3    100 0.65  H1 0.32449164    0.84679173
#> 
#> $stopprob$`0.7`
#>   events  pH hyp stop_prob cum_stop_prob
#> 1     10 0.7  H1 0.1493083     0.1493083
#> 2     50 0.7  H1 0.6414520     0.7907603
#> 3    100 0.7  H1 0.1910207     0.9817810
#> 
#> $stopprob$`0.75`
#>   events   pH hyp  stop_prob cum_stop_prob
#> 1     10 0.75  H1 0.24402523     0.2440252
#> 2     50 0.75  H1 0.70328355     0.9473088
#> 3    100 0.75  H1 0.05210513     0.9994139
#> 
#> $stopprob$`0.8`
#>   events  pH hyp   stop_prob cum_stop_prob
#> 1     10 0.8  H1 0.375809638     0.3758096
#> 2     50 0.8  H1 0.618228306     0.9940379
#> 3    100 0.8  H1 0.005959013     0.9999970
#> 
#> 
#> $opchar
#>      p cum_stop_prob expected_events hyp
#> 1 0.50    0.03758441       146.40897  H0
#> 2 0.60    0.50668952       110.00282  H1
#> 3 0.65    0.84679173        78.10723  H1
#> 4 0.70    0.98178099        55.40060  H1
#> 5 0.75    0.99941391        42.90286  H1
#> 6 0.80    0.99999696        35.26587  H1
#> 
#> attr(,"class")
#> [1] "harmbound" "list"     
#or as risk ratio
getHarmBound(nevents=c(10,50,100), alpha_test=0.025, pH0=0.5, rrH1=1.5, maxevents=150)
#> $bounds
#>   events events_intervention events_control alpha_test
#> 1     10                   9              1      0.025
#> 2     50                  33             17      0.025
#> 3    100                  61             39      0.025
#> 
#> $stopprob
#> $stopprob$`0.5`
#>   events  pH hyp  stop_prob cum_stop_prob
#> 1     10 0.5  H0 0.01074219    0.01074219
#> 2     50 0.5  H0 0.01490030    0.02564249
#> 3    100 0.5  H0 0.01194192    0.03758441
#> 
#> $stopprob$`0.6`
#>   events  pH hyp stop_prob cum_stop_prob
#> 1     10 0.6  H1 0.0463574     0.0463574
#> 2     50 0.6  H1 0.2098108     0.2561682
#> 3    100 0.6  H1 0.2505214     0.5066895
#> 
#> 
#> $opchar
#>     p  rr cum_stop_prob expected_events hyp
#> 1 0.5 1.0    0.03758441        146.4090  H0
#> 2 0.6 1.5    0.50668952        110.0028  H1
#> 
#> attr(,"class")
#> [1] "harmbound" "list"     
```
