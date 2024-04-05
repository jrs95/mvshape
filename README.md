# mvshape <img src="man/figures/logo.png" align="right" height="139"/>
This package is used to assess the shape of the relationship between an exposure and an outcome using fractional polynomials and multivariate meta-analysis. 

## Functions
* `fracpoly`: this method fits the best-fitting fractional polynomial of degree 1 and/or 2 of the association between the outcome and exposure.  
* `mvshape`: this method performs multivariate meta-analysis of the associations of the outcome with quantiles of the exposure across studies.  
* `fracploy_plot`: this function plots the fractional polynomial model.  
* `mvshape_plot`: this function plots the mvshape model.  
* `fracpoly_mvshape_plot`: this function plots the fractional polynomial and mvshape models.  

## Installation
```
install.packages("remotes")
remotes::install_github("jrs95/mvshape")
```

## Example
```
# Library
library(mvshape)

# Data
y <- rnorm(5000)
x <- rnorm(5000, 10, 1)
c1 <- rbinom(5000, 1, 0.5)
c2 <- rnorm(5000)
study <- as.factor(c(rep("study1", 1000), rep("study2", 1000), rep("study3", 1000), rep("study4", 1000), rep("study5", 1000)))
covar <- data.frame(c1 = c1, c2 = c2, study = study)

# Analyses
fp <- fracpoly(y = y, x = x, covar = covar, family = "gaussian")
fracpoly_plot(fp)
mvma <- mvshape(y = y, x = x, covar = covar[,c("c1", "c2")], study = study, family = "gaussian")
mvshape_plot(mvma)
fracpoly_mvshape_plot(fracpoly = fp, mvshape = mvma)
```

## Citation 
Please cite this R package using the link: https://github.com/jrs95/mvshape  
