# Find stopping boundary via binomial exact tests

Find stopping boundary via binomial exact tests

## Usage

``` r
findbound(n, alpha_test = 0.025, pH0 = 0.5, alternative = "greater")
```

## Arguments

- n:

  total number of events

- alpha_test:

  nominal alpha for the binomial test

- pH0:

  proportion of events in the experimental arm under the null
  hypothesis, typically based on randomization ratio (e.g. 0.5 for a 1:1
  randomization)

- alternative:

  direction of alternative, "less" or "greater"

## Value

number of events in the experimental group that would lead to a stopping

## Examples

``` r
findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="greater")
#> [1] 15
findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="less")
#> [1] 5

```
