# Convert the proportion of events in the intervention groups to risk differences and ratios (and vice versa)

Convert the proportion of events in the intervention groups to risk
differences and ratios (and vice versa)

## Usage

``` r
convertRisks(
  eprop = NULL,
  etotal = NULL,
  rd = NULL,
  rr = NULL,
  or = NULL,
  r0 = NULL,
  n0,
  n1 = n0
)
```

## Arguments

- eprop:

  proportion of events in intervention group

- etotal:

  total number of events

- rd:

  risk difference

- rr:

  risk ratio

- or:

  odds ratio

- r0:

  risk in the control group

- n0:

  number of patients in the control group

- n1:

  number of patients in the intervention group

## Value

vector with risks in control and intervention group (r0, r1), the risk
difference (rd), risk ratio (rr) and odds ratio (or)

## Examples

``` r
convertRisks(eprop=0.5,etotal=100,n0=200)
#>  eprop etotal     n0     n1     r0     r1     rd     rr     or 
#>   0.50 100.00 200.00 200.00   0.25   0.25   0.00   1.00   1.00 
convertRisks(eprop=0.6,etotal=100,n0=200)
#>      eprop     etotal         n0         n1         r0         r1         rd 
#>   0.600000 100.000000 200.000000 200.000000   0.200000   0.300000   0.100000 
#>         rr         or 
#>   1.500000   1.714286 
convertRisks(rr=1.5,n0=200,r0=0.2)
#>      eprop     etotal         n0         n1         r0         r1         rd 
#>   0.600000 100.000000 200.000000 200.000000   0.200000   0.300000   0.100000 
#>         rr         or 
#>   1.500000   1.714286 
```
