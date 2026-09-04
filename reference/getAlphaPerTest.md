# Test-wise alpha necessary to control either the family-wise type I error or the power at a specified level

Exactly one of totalAlpha or power have to be specified. The power
requires the specification of an alternative via one of pH1, rrH1, orH1
or rdH1.

## Usage

``` r
getAlphaPerTest(
  nevents,
  totalAlpha = NULL,
  power = NULL,
  pH0 = 0.5,
  alpha.interval = c(10^(-10), 1),
  maxevents = NULL,
  pH1 = NULL,
  rrH1 = NULL,
  orH1 = NULL,
  rdH1 = NULL,
  r0 = NULL
)
```

## Arguments

- nevents:

  vector with number of events at which an interim analysis is done

- totalAlpha:

  target overall family-wise type I error

- power:

  target power at the specified alternative

- pH0:

  proportion of events in the intervention arm under the null
  hypothesis, typically based on randomization ratio (e.g. 0.5 for a 1:1
  randomization)

- alpha.interval:

  Range for test-wise alpha, c(10^(-10),0.05) by default

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

## Value

Test-wide alpha

## Examples

``` r
#Control overall family-wise type I error:
apt<-getAlphaPerTest(nevents = c(10,50,100), totalAlpha = 0.05, pH0 = 0.5)
apt
#> [1] 0.03245429
getHarmBound(nevents = c(10,50,100),alpha_test = apt, pH0 = 0.5)
#> $bounds
#>   events events_intervention events_control alpha_test
#> 1     10                   9              1 0.03245429
#> 2     50                  33             17 0.03245429
#> 3    100                  60             40 0.03245429
#> 
#> $stopprob
#> $stopprob$`0.5`
#>   events  pH hyp  stop_prob cum_stop_prob
#> 1     10 0.5  H0 0.01074219    0.01074219
#> 2     50 0.5  H0 0.01490030    0.02564249
#> 3    100 0.5  H0 0.02085634    0.04649882
#> 
#> 
#> $opchar
#>     p cum_stop_prob expected_events hyp
#> 1 0.5    0.04649882        98.28819  H0
#> 
#> attr(,"class")
#> [1] "harmbound" "list"     

#Control power assuming 80% of events in expermintal arm
apt<-getAlphaPerTest(nevents = c(10,50,100), power = 0.8, pH0 = 0.5, pH1 = 0.6)
apt
#> [1] 0.1013193
getHarmBound(nevents = c(10,50,100),alpha_test = apt, pH0 = 0.5, pH1 = 0.6)
#> $bounds
#>   events events_intervention events_control alpha_test
#> 1     10                   8              2  0.1013193
#> 2     50                  31             19  0.1013193
#> 3    100                  57             43  0.1013193
#> 
#> $stopprob
#> $stopprob$`0.5`
#>   events  pH hyp  stop_prob cum_stop_prob
#> 1     10 0.5  H0 0.05468750     0.0546875
#> 2     50 0.5  H0 0.04648929     0.1011768
#> 3    100 0.5  H0 0.05797765     0.1591544
#> 
#> $stopprob$`0.6`
#>   events  pH hyp stop_prob cum_stop_prob
#> 1     10 0.6  H1 0.1672898     0.1672898
#> 2     50 0.6  H1 0.3260615     0.4933512
#> 3    100 0.6  H1 0.3043610     0.7977123
#> 
#> 
#> $opchar
#>     p cum_stop_prob expected_events hyp
#> 1 0.5     0.1591544        92.75366  H0
#> 2 0.6     0.7977123        68.64085  H1
#> 
#> attr(,"class")
#> [1] "harmbound" "list"     
```
