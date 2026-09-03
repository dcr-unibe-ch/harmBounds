# Plot method for harmbound objects produced by `getHarmBound`

Plot method for harmbound objects produced by `getHarmBound`

## Usage

``` r
# S3 method for class 'harmbound'
plot(x, which = "bounds", ...)
```

## Arguments

- x:

  harmbounds objects as generated using the getHarmBound function

- which:

  one of `"bounds"`, `"abs_stopping"`, `"cum_stopping"`, `"exp_n"`.

- ...:

  options passed to plot

## Examples

``` r
harmbound<-getHarmBound(nevents = seq(10, 100, by=10),alpha_test = 0.025,
pH0 = 0.5, pH1 = seq(0.55,0.7,by=0.05), maxevents = 150)
plot(harmbound, which = "bounds")

plot(harmbound, which = "abs_stopping")

plot(harmbound, which = "cum_stopping")

plot(harmbound, which = "opchar_stop")

plot(harmbound, which = "opchar_n")

```
