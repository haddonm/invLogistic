
<!-- README.md is generated from README.Rmd. Please edit that file -->

# invLogistic

An R package to faciliate the fitting and exploration of inverse
logistic models to invertebrate tagging data.

To install this development version from [GitHub](https://github.com/)
you can use:

``` r
if (!require(devtools)){install.packages("devtools")} 

devtools::install_github("https://github.com/haddonm/invLogistic")
```

There are currently no vignettes, although each function has its own
help page.

If you were to type *methods(“plot”)* into the console you will find
both plot.bootIL and plot.IL in the list. There are help pages for
those, as well as for print.IL and summary.IL. To see the contents of
the objects output by, say, *fitIL*, you can use *str(ans)*.

Some details of the package can be found using
*packageDescription(“invLogistic”)*.

  - 22/06/2020 Some modifications to generalize the functions for use
    with invertebrates other than abalone.

Malcolm Haddon

Hobart, March 12, 2020

## 

Haddon, M., Mundy, C., and D. Tarbath (2008) Using an inverse-logistic
model to describe growth increments of blacklip abalone (*Haliotis
rubra*) in Tasmania. *Fishery Bulletin* **106**:58-71
