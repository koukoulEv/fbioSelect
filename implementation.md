Two-stage algorithm test
================

This is an R markdown document which uses a small fragment of the GDSC data and implemets the two-stage algorithm.

``` r
load("/home/koukouli/Documents/gdsc_paper/R package develop/test_gdsc_data.RData")
source("/home/koukouli/Documents/gdsc_paper/R package develop/two_stage_algorithm_functions.R")
```

Function Implementation with weights:

``` r
result = tsalgo(y = fd$Response, 
                id = fd$ID,
                t.var = fd$DosageSQ, 
                Z = model.matrix(~fd$DRUG_NAME),
                X = fd[,-c(1:5)],
                low.dim.cov.number = 1, 
                threshold.size = floor(length(unique(fd$ID))/log(length(unique(fd$ID)))),
                type.of.penalty = "grSCAD", 
                degree = 3,
                intercept = T, 
                degrFreedom = 5, 
                knots = NULL, 
                weights = fd$weights, 
                nfolds = 10,
                seed = 1000)
```

    ## Loading required package: tidyverse

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0     ✓ purrr   0.3.3
    ## ✓ tibble  3.0.0     ✓ dplyr   0.8.5
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    ## Loading required package: foreach

    ## 
    ## Attaching package: 'foreach'

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     accumulate, when

    ## Loading required package: grpreg

    ## 
    ## Attaching package: 'grpreg'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## Loading required package: doParallel

    ## Loading required package: iterators

    ## Loading required package: parallel

    ## Loading required package: splines

``` r
par(mfrow = c(2,2))
plot(seq(0,1,length = ncol(result$coefFunctions.pen.regression))^2, result$coefFunctions.pen.regression[1,], type = "l",
     xlab = "Dosage level", ylab = "Response", main = "Intercept")
plot(seq(0,1,length = ncol(result$coefFunctions.pen.regression))^2, result$coefFunctions.pen.regression[2,], type = "l",
     xlab = "Dosage level", ylab = "Response", main = "SB590885")
plot(seq(0,1,length = ncol(result$coefFunctions.pen.regression))^2, result$coefFunctions.pen.regression[3,], type = "l",
     xlab = "Dosage level", ylab = "Response", main = paste(result$active.hdcov.pen.regression[1]))
plot(seq(0,1,length = ncol(result$coefFunctions.pen.regression))^2, result$coefFunctions.pen.regression[4,], type = "l",
     xlab = "Dosage level", ylab = "Response", main = paste(result$active.hdcov.pen.regression[2]))
```

![](implementation_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
par(mfrow = c(1,1))
```
