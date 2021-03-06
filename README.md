
<!-- README.md is generated from README.Rmd. Please edit that file -->
R Package SMUT: Testing the Mediation Effect of Multiple SNPs on an Outcome Through a Mediator
==============================================================================================

<!-- badges: start -->
<!-- badges: end -->
Installation
============

You can install the released version of SMUT from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SMUT")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wjzhong/SMUT")
```

Overview
========

SMUT package has functions to test for mediation effects of multiple SNPs (G) on a continuous, binary, count or time-to-event outcome (Y) through a continuous mediator (M). Besides SNPs data, these functions can also be applied to test mediation effects in other fields. In general, SMUT package has functions to test mediation effects of multiple continuous variables (G) on a continuous, binary, count or time-to-event outcome (Y) through a continuous mediator (M).

Multi-SNP Mediation Intersection-Union Test (SMUT) : Testing mediation effects of multiple SNPs on a continuous outcome through a continuous mediator
=====================================================================================================================================================

Example of running R function SMUT
----------------------------------

The example data has a genotype matrix (G) of 100 individuals with 200 SNPs.

``` r
library(SMUT)
dim(Genotype_data)
#> [1] 100 200
Genotype_data[1:3,1:4]
#>              SNP_1 SNP_2 SNP_3 SNP_4
#> Individual_1     0     2     1     0
#> Individual_2     1     1     0     0
#> Individual_3     1     1     0     1
```

We generate one continuous mediator and one continuous outcome based on this genotype matrix, as well as two covariates.

``` r
N_individual = nrow(Genotype_data)
N_SNPs = ncol(Genotype_data)

set.seed(1)

# generate two covariates 
X1 = rnorm(N_individual, 2, 3)
X2 = sample(c(0,1), N_individual, replace = TRUE) 
X = cbind(X1, X2)

# generate coefficients: iota_M, iota_Y, beta, theta and gamma
iota_M = c(0.3,0.5)
iota_Y = c(0.2,0.6)
beta = rnorm(N_SNPs, 1, 2)
theta = 1.2
gamma = rnorm(N_SNPs, 0.5, 2)

# generate error terms
e1 = rnorm(N_individual, 0, 1)
e2 = rnorm(N_individual, 0, 1)

# generate the mediator
mediator = 1 + X %*% iota_M + Genotype_data %*% beta + e1

# generate the outcome
outcome = 2 + mediator*theta + X %*% iota_Y + Genotype_data %*% gamma + e2
```

We apply SMUT function to test the mediation effect.

``` r
result_continuous = SMUT(G = Genotype_data, mediator = mediator, outcome = outcome, covariates = X)
#> Warning: 2 SNPs with either high missing rates or no-variation are
#> excluded!
#> Warning: Genotypes of some variants are not the number of minor alleles!
#> These genotypes are flipped!
print( unlist( result_continuous ))
#>   p_value_IUT p_value_theta  p_value_beta 
#>  2.595523e-02  2.595523e-02  3.286609e-08
```

The warning messages are generated by the SKAT function in the R package SKAT. From the result, we can see that the p value of the mediation effect (*p\_value\_IUT*) is 0.02595523, which is the maximum of the p value of testing *β* (*p\_value\_beta*) and the p value of testing *θ* (*p\_value\_theta*).

References
----------

Full details of the SMUT method can be found in the manuscript:

Zhong, W., Spracklen, C.N., Mohlke, K.L., Zheng, X., Fine, J. and Li, Y., 2019. Multi-SNP mediation intersection-union test. *Bioinformatics*, 35(22), pp.4724-4729.

Generalized Multi-SNP Mediation Intersection-Union Test (GSMUT) : Testing mediation effects of multiple SNPs on a continuous, binary, count or time-to-event outcome through a continuous mediator
==================================================================================================================================================================================================

In the manuscript, GSMUT is referred to as SMUT\_GLM for a continuous, binary or count outcome and as SMUT\_PH for a time-to-event outcome.

Example of running R function GSMUT
-----------------------------------

The example data has a genotype matrix (G) of 100 individuals with 200 SNPs.

``` r
library(SMUT)
dim(Genotype_data)
#> [1] 100 200
Genotype_data[1:3,1:4]
#>              SNP_1 SNP_2 SNP_3 SNP_4
#> Individual_1     0     2     1     0
#> Individual_2     1     1     0     0
#> Individual_3     1     1     0     1
```

We generate one continuous mediator, one binary outcome and one time-to-event outcome based on this genotype matrix, as well as two covariates.

``` r
N_individual = nrow(Genotype_data)
N_SNPs = ncol(Genotype_data)

set.seed(1)
# generate two covariates 
X1 = rnorm(N_individual, 2, 3)
X2 = sample(c(0,1), N_individual, replace = TRUE) 
X = cbind(X1, X2)

# generate coefficients: iota_M, iota_Y, beta, theta and gamma
iota_M = c(0.3,0.5)
iota_Y = c(0.2,-0.6)
beta = rnorm(N_SNPs, 0, 0.5)
theta = 1
gamma = rnorm(N_SNPs, 0, 0.3)

# generate error terms
e1 = rnorm(N_individual, 0, 1)

# generate the mediator
mediator = 1 + X %*% iota_M + Genotype_data %*% beta + e1

# generate the binary outcome 
eta = 2 + mediator*theta + X %*% iota_Y + Genotype_data %*% gamma
pi = 1/(1+exp( -(eta ) ))
binary_outcome = rbinom(length(pi),size=1,prob=pi) 

# generate the time-to-event outcome based on Weibull baseline hazard
v = runif(N_individual)
lambda=0.01; rho=1; rateC=0.01
Tlat = (- log(v) / (lambda * exp( eta  )))^(1 / rho)
# censoring times
C = rexp(N_individual, rate=rateC)
# follow-up times and event indicators
time = pmin(Tlat, C)
status = as.numeric(Tlat <= C)
survival_outcome = cbind(time,status)      
colnames(survival_outcome) = c("time","status")
```

We apply GSMUT function to test the mediation effect.

``` r
result_binary = GSMUT(G = Genotype_data, mediator = mediator, outcome = binary_outcome, covariates = X, outcome_type = "binary")
#> Warning: 2 SNPs with either high missing rates or no-variation are
#> excluded!
#> Warning: Genotypes of some variants are not the number of minor alleles!
#> These genotypes are flipped!
print( unlist( result_binary ) )
#>   p_value_IUT p_value_theta     theta_hat  p_value_beta 
#>  0.0004216730  0.0002812888  0.8378383351  0.0004216730

result_survival = GSMUT(G = Genotype_data, mediator = mediator, outcome = survival_outcome, covariates = X, outcome_type = "survival")
#> Warning: 2 SNPs with either high missing rates or no-variation are
#> excluded!

#> Warning: Genotypes of some variants are not the number of minor alleles!
#> These genotypes are flipped!
print( unlist( result_survival ))
#>   p_value_IUT p_value_theta     theta_hat  p_value_beta 
#>   0.024700658   0.024700658   0.159754100   0.000421673
```

The warning messages are generated by the SKAT function in the R package SKAT. From the result, we can see that the p value of the mediation effect (*p\_value\_IUT*) for the binary outcome is 0.000421673; p value of the mediation effect for the time-to-event outcome is 0.02470066. Here the p value of the mediation effect (*p\_value\_IUT*) is the maximum of the p value of testing *β* (*p\_value\_beta*) and the p value of testing *θ* (*p\_value\_theta*). And *theta\_hat* is the point estimate of the *θ* in the outcome model.

Options for the parameter outcome\_type in the GSMUT function
-------------------------------------------------------------

| Type of outcome | outcome\_type |
|:---------------:|:-------------:|
|    continuous   |  "continuous" |
|      binary     |    "binary"   |
|      count      |    "count"    |
|  time-to-event  |   "survival"  |

References
----------

Full details of the GSMUT method can be found in the manuscript:

Zhong, W., Darville, T., Zheng, X., Fine, J. and Li, Y., 2019. Generalized Multi-SNP Mediation Intersection-Union Test. *Submitted*
