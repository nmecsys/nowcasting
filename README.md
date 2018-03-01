
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nowcasting)](https://CRAN.R-project.org/package=nowcasting) 
[![downloads](http://cranlogs.r-pkg.org/badges/grand-total/nowcasting)](https://cran.rstudio.com/web/packages/nowcasting/index.html) 
![](http://cranlogs.r-pkg.org/badges/last-week/nowcasting?color=green)
![Build Status](https://ci.appveyor.com/api/projects/status/github/guilbran/nowcast?branch=master&svg=true)

# nowcasting
Forecasting in real time (nowcasting) the Brazilian GDP

This is a R package to facilitate the reproducibility of nowcasting in academic researches.
The backend functions of this package are based on Giannone et al. (2008) and Banbura et al. (2011) replication files.

The package is in development. Reviews, comments and pull requests are welcome.

Read more about using the package in https://goo.gl/kAL5Qs

## Installation

github
```R
devtools::install_github('nmecsys/nowcasting')
```
CRAN
```R
install.packages('nowcasting')
```
## To do

- Create functions to estimate the importance of news;

- Optimize the estimation using other Kalman Filter packages as fkf for instance;

- Class output to use plot() instead of nowcast.plot().



