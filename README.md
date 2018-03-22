
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nowcasting)](https://CRAN.R-project.org/package=nowcasting) 
[![downloads](http://cranlogs.r-pkg.org/badges/grand-total/nowcasting)](https://cran.rstudio.com/web/packages/nowcasting/index.html) 
![](http://cranlogs.r-pkg.org/badges/last-week/nowcasting?color=green)
![Build Status](https://ci.appveyor.com/api/projects/status/github/guilbran/nowcast?branch=master&svg=true)

# nowcasting
An R Package for Forecasting Models with Real-Time Data.

The **nowcasting** package contains useful tools for developing real-time forecasting models so that researchers can find in this package a simple and practical way to reproduce several of these models. In this version of the package we present three methods, based on seminal articles in this literature: *Giannone et al. 2008* and *Bańbura et al. 2011*. The package also provides vintages (information observed at the time of publication) of Brazilian economic data.


**The package is in development. Reviews, comments and pull requests are welcome.**


## Installation

github
```R
devtools::install_github('nmecsys/nowcasting')
```
CRAN
```R
install.packages('nowcasting')
```

## How to use nowcasting package

Two examples of how the nowcasting package can be used are discussed below. In the first example we use the dataset in the classic paper by *Giannone et al. (2008)*, which uses a panel of time series for the US economy, and in the second we show how to obtain vintages for the main series for the Brazilian economy.

### nowcasting US GDP

The dataset used in *Giannone et al. (2008)* can be added to the environment with the following command.

```{r warning=FALSE}
library(nowcasting)
data(USGDP)
```

USDGP is a list with two data frames:

* `USGDP$base`: this is an unbalanced panel with the model variables;
* `USGDP$legend`: this contains the legend with the relevant information for each variable in the `USGDP$base` dataframe.

The first step is to define the variables $y_{t}$ (quarterly) and $x_{t}$ (monthly). To identify the variable that represents the GDP in the dataset, view the legend in `USGDP$legend`. Then the GDP is selected in the dataset and transformed into a quarterly variable by `month2qtr` function. The explanatory variables are all the other variables in the dataset.

```{r warning=FALSE}
gdp <- month2qtr(x = USGDP$base[,"RGDPGR"])
```

The second step covers the treatment of the explanatory variables, that will be done by `Bpanel` function. This function creates a balanced panel using an unbalanced panel as input. The missing observations and outliers are substituted using the outlier correction methodology in *Giannone et al. (2008)* and there are some transformation options that make the variables stationary. The object `USGDP$legend` contains all the transformations used in *Giannone et al. (2008)*. The `Bpanel` function also allows monthly variables to be aggregated to represent quarterly quantities.

```{r warning=FALSE}
gdp_position <- which(colnames(USGDP$base) == "RGDPGR")
base <- Bpanel(base = USGDP$base[,-gdp_position], 
               trans = USGDP$legend$Transformation[-gdp_position], aggregate = TRUE)
```
Once these variables have been treated, the `nowcast` function can be used to estimate the model parameters according to the estimation method selected, the number *r* of dynamic factors, the lag order of the factors *p* and the number *q* of shocks in the factors. In this example we use the *Two-Stage - quarterly factors* method. The arguments *r*, *p* and *q* were defined according to *Giannone et al. (2008)*.

```{r warning=FALSE}
nowcastUSGDP <- nowcast(y = gdp, x = base, r = 2, p = 2, q = 2, method = '2sq')
```

To see the **results** of the forecasts in the object `nowcastUSGDP`, use the following command.

```{r warning=FALSE}
# y forecasts
tail(nowcastUSGDP$yfcst,8)

# the regression between y and its factors can be accessed using `$reg`.
summary(nowcastUSGDP$reg)

# the results related to the estimation of factors 
tail(nowcastUSGDP$factors$dynamic_factors) # factors
head(nowcastUSGDP$factors$Lambda) # Lambda matrix
nowcastUSGDP$factors$A # A matrix
nowcastUSGDP$factors$BB # BB': u's variance covariance matrix (factor equation)
diag(nowcastUSGDP$factors$Psi) # Psi: epsilon's variance covariance matrix (x equation)

# the forecasts of the explanatory variables are in `$xfcst`.
tail(nowcastUSGDP$xfcst[,1:5]) # x forecasts (first 5 variables)
```

The **graphs** available with the `nowcast.plot` function allow the forecast of the variable of interest and the estimated factors to be easily viewed.

```{r warning=FALSE}
 # y fcst
nowcast.plot(nowcastUSGDP, type = "fcst")

# factors
nowcast.plot(nowcastUSGDP, type = "factors") 

 # how much of the variability in the dataset is explained by each factor 
nowcast.plot(nowcastUSGDP, type = "eigenvalues")

# importance of each variable in the first factor
nowcast.plot(nowcastUSGDP, type = "eigenvectors") 
```

### vintages

Vintage in this context refers to the dataset observed on a specific date. It is useful to assess the performance of models with information published at different times and so reproduce the real past.

The **nowcasting** package contains a dataset of very important Brazilian economic time series, which can be accessed with the `RTDB` function. If no argument is specified, the codes of all the available series are returned. The codes are the same as those used on the Central Bank of Brazil platform to download time series.

```{r warning=FALSE}
# the available series
head(RTDB())

# serie 1: availables vintages
head(RTDB(series_code = 1)) 

# vintage 2017-04-04: availables series
 head(RTDB(vintage = "2017-04-04")) 

# serie 1, vintage 2017-04-04: data
tail(RTDB(series_code = 1, vintage = "2017-04-04")) 
```

The `PRTDB` function is intended to simulate the vintages of any dataset. The function only excludes observations from the time series based on lag information provided by the user and simulates what would be observed on the reference date.

```{r warning=FALSE}
# BRGDP data (last six observations)
tail(BRGDP)

# BRGDP data observed on the reference date (2017-10-01)
tail(PRTDB(mts = BRGDP, delay = c(1,30,60,90,20,10,30,60), vintage = "2017-10-01"))
```

