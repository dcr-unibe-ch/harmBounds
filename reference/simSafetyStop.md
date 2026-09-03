# Simulate safety stopping

The effect can be given as proportion of events in the experimental
group (pH1), the risk difference (rdH1), risk ratio (rrH1) or odds ratio
(orH1).

## Usage

``` r
simSafetyStop(
  nevents,
  pH0 = 0.5,
  alpha_test = 0.025,
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

  vector with number of events at which an interim analysis is done

- pH0:

  proportion of events in the experimental arm under the null
  hypothesis, typically based on randomization ratio (e.g. 0.5 for a 1:1
  randomization)

- alpha_test:

  nominal alpha level for binomial exact test

- pH1:

  proportion of events in the experimental arm under the alternative
  hypothesis

- rrH1:

  risk ratio (experimental / control).

- orH1:

  risk ratio (experimental / control). Requires the control proportion
  (r0).

- rdH1:

  risk difference (experimental - control). Requires the control
  proportion (r0) and the number of participants (n).

- r0:

  risk in the control group. Required if the effect is given as risk
  difference or odds ratio.

- n:

  total number of participants. Required if the effect is given as risk
  difference.

## Value

list with a dataframe with number of events in each group plus upper
limit for stopping and indicator for whether stopped, plus indicators
number of stops and time points at first stop

## Examples

``` r
set.seed(1)  
simSafetyStop(nevents=seq(10,100,by=10),pH0 = 0.5, pH1 = 0.6,alpha_test=0.025)
#> $tests
#>    nevents n1 n0 ulim   out
#> 1       10  5  5    9 FALSE
#> 2       20 10 10   15 FALSE
#> 3       30 17 13   21 FALSE
#> 4       40 23 17   27 FALSE
#> 5       50 27 23   33 FALSE
#> 6       60 35 25   39 FALSE
#> 7       70 41 29   44 FALSE
#> 8       80 46 34   50 FALSE
#> 9       90 53 37   55 FALSE
#> 10     100 57 43   61 FALSE
#> 
#> $nstop
#> [1] 0
#> 
#> $tstop
#> [1] NA
#> 

set.seed(1)  
simSafetyStop(nevents=seq(10,100,by=10),pH0 = 0.5, rrH1 = 0.6/(1-0.6), alpha_test=0.025)
#> $tests
#>    nevents n1 n0 ulim   out
#> 1       10  5  5    9 FALSE
#> 2       20 10 10   15 FALSE
#> 3       30 17 13   21 FALSE
#> 4       40 23 17   27 FALSE
#> 5       50 27 23   33 FALSE
#> 6       60 35 25   39 FALSE
#> 7       70 41 29   44 FALSE
#> 8       80 46 34   50 FALSE
#> 9       90 53 37   55 FALSE
#> 10     100 57 43   61 FALSE
#> 
#> $nstop
#> [1] 0
#> 
#> $tstop
#> [1] NA
#> 

```
