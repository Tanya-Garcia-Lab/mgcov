# mgcov: maximum likelihood estimation of Gaussian covariance matrices
The R package `mgcov` performs maximum likelihood estimation of a sparse covariance matrix for Gaussian distribution using thresholding and banding.

Relevant references can be found at:
* [A positive-definite thresholding estimator of a covariance matrix and its asymptotic efficiency] by Kim et al. (2021).
* [Estimation of a covariance matrix with zeros](https://doi.org/10.1093/biomet/asm007) by Chaudhuri et al. (2007).

## Installation

```
devtools::install_github("rakheon/mgcov", force = TRUE)
```

## Example

```
library(mgcov)
library(mvtnorm)

# data generation
p <-10; n<- 100
Sigma = diag(p); diag(Sigma[-p,-1])=0.5; diag(Sigma[-1,-p])=0.5
dat <- rmvnorm(n,mean=rep(0,ncol(Sigma)),sigma=Sigma)

# COMET (univeral thresholding) by AIC and BIC
res_uni = COmet(dat, lambda = seq(0,1,0.01))
res_uni$cov_list[[which.min(res_uni$aic)]]; res_uni$cov_list[[which.min(res_uni$bic)]]

# COMET (adaptive thresholding) by AIC and BIC
res_ada = COmet(dat, mul = 3)
res_ada$cov_list[[which.min(res_ada$aic)]]; res_ada$cov_list[[which.min(res_ada$bic)]]
```
