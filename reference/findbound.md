# Find stopping boundary via binomial or beta-binomial distribution. For the latter, the overdispersion factor (or design effect) by which the variance exceeds the regular binomial variance is printed.

Find stopping boundary via binomial or beta-binomial distribution. For
the latter, the overdispersion factor (or design effect) by which the
variance exceeds the regular binomial variance is printed.

## Usage

``` r
findbound(
  n,
  alpha_test = 0.025,
  pH0 = 0.5,
  alternative = "greater",
  icc = NULL
)
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

- icc:

  intraclass correlation if there is more than one event per patient

## Value

number of events in the experimental group that would lead to a stopping

## Examples

``` r
findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="greater")
#> [1] 15
findbound(n=20, alpha_test=0.025, pH0 = 0.5, alternative="less")
#> [1] 5

```
