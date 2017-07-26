# sibyl

[![Build Status](https://travis-ci.org/scientific-computing-solutions/sibyl.svg?branch=master)](https://travis-ci.org/scientific-computing-solutions/sibyl)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sibyl)](https://cran.r-project.org/package=sibyl)
[![Coverage Status](https://coveralls.io/repos/github/scientific-computing-solutions/sibyl/badge.svg?branch=master)](https://coveralls.io/github/scientific-computing-solutions/sibyl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/nux777kry2ahv4ra/branch/master?svg=true)](https://ci.appveyor.com/project/bobthecat/sibyl/branch/master)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/sibyl)](https://cran.r-project.org/package=sibyl)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/sibyl)](https://cran.r-project.org/package=sibyl)


Fitting parametric semi-parametric models to time-to-event data and extrapolating 
these curves to larger survival times.

This R-package is design to simplify the statistical analyses needed to perform 
standard survival extrapolation work needed for submissions to Health Authorities 
such as NICE. The resulting reports are often tedious to produce and by standardizing 
the methods and output there is both a huge gain in manual labour and in eliminating
potential errors. The analysis is performed on one or several time-to-event variables
on either all data or subgroups. Both semi-parametric and parametric models are 
fitted to the data either as a separate model per arm or with arm as a factor; 
with or without covariates. Plots for checking model and proportional hazards 
assumptions are implemented together with survival plots. Both stratified and 
un-stratified log rank tests and restricted mean survival times (RMSTs) can be 
computed. Average survival curves can also be used to output a single curve 
representing an "average patient" based on all covariates.


## Contributors

Dalevi, Daniel (maintainer); Burkoff, Nikolas; Ouwens, Mario; Ruau, David;

## Installation

To install the development version from GitHub:
```R
install.packages("devtools")
# We spent a lot of time developing the vignettes. We recommend the read but 
# building them from source takes some time
devtools::install_github("scientific-computing-solutions/sibyl", 
                         build_vignettes = TRUE)
```
