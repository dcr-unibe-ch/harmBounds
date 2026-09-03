# Plot cumulative stopping probabilities

Plot cumulative stopping probabilities

## Usage

``` r
cumstopPlot(harmbound)
```

## Arguments

- harmbound:

  harmbounds objects as generated using the getHarmBound function

## Value

barplot with cumulative stopping probabilities under H0 and optionally
H1

## Examples

``` r
harmbound<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,pH0=0.5,pH1=c(0.6,0.7,0.8))
cumstopPlot(harmbound)

```
