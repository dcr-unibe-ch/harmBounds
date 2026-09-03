# Test-wise alpha necessary to control the overall type I error at a specified level (0.05 by default)

Test-wise alpha necessary to control the overall type I error at a
specified level (0.05 by default)

## Usage

``` r
getAlphaPerTest(
  nevents,
  totalAlpha = 0.05,
  pH0 = 0.5,
  alpha.interval = c(10^(-10), 0.05)
)
```

## Arguments

- nevents:

  vector with number of events at which an interim analysis is done

- totalAlpha:

  Overall type I error, 0.05 by default

- pH0:

  proportion of events in the intervention arm under the null
  hypothesis, typically based on randomization ratio (e.g. 0.5 for a 1:1
  randomization)

- alpha.interval:

  Range for test-wise alpha, c(10^(-10),0.05) by default

## Value

Test-wide alpha

## Examples

``` r
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

```
