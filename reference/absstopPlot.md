# Plot absolute stopping probabilities

Plot absolute stopping probabilities

## Usage

``` r
absstopPlot(harmbound)
```

## Arguments

- harmbound:

  harmbounds objects as generated using the getHarmBound function

## Value

barplot with stopping probabilities under H0 and optionally H1

## Examples

``` r
harmbound<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,pH0=0.5,pH1=c(0.6,0.7,0.8))
absstopPlot(harmbound)

```
