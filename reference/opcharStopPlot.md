# Plot cumulative stopping probabilities by hypothesis

Plot cumulative stopping probabilities by hypothesis

## Usage

``` r
opcharStopPlot(harmbound)
```

## Arguments

- harmbound:

  harmbounds objects as generated using the getHarmBound function

## Value

line plot with the cumulative stopping for all alternatives

## Examples

``` r
harmbound<-getHarmBound(nevents=seq(10,100,by=10),alpha_test=0.025,
pH0=0.5,pH1=seq(0.5,0.8,by=0.05),maxevents=150)
opcharStopPlot(harmbound)

```
