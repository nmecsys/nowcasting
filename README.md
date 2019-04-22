
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nowcasting)](https://CRAN.R-project.org/package=nowcasting) 
[![downloads](http://cranlogs.r-pkg.org/badges/grand-total/nowcasting)](https://cran.rstudio.com/web/packages/nowcasting/index.html) 
![](http://cranlogs.r-pkg.org/badges/last-week/nowcasting?color=green)
![Build Status](https://ci.appveyor.com/api/projects/status/github/guilbran/nowcast?branch=master&svg=true)

# nowcasting
An R Package for Forecasting Models with Real-Time Data.

The **nowcasting** package contains useful tools for using dynamic factor models. In this version of the package we present three methods, based on the articles of *Giannone et al. 2008* and *Bańbura et al. 2011*. Furthermore, the package offers auxiliary functions to treat variables, constuct vintages, visualize results, etc.


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

## How to use the nowcasting package

Two examples of how the nowcasting package can be used are discussed below. In the first example we use the dataset from *Giannone et al. (2008)*, which uses a time series panel for the US economy. In the second example we show how to make forecasts based on pseudo real time vintages using Brazilian data and how to use information criteria for determining the number of factors and shocks to use in the model.

### nowcasting US GDP 1/2

The dataset used in *Giannone et al. (2008)* can be added to the environment with the following command.

```{r warning=FALSE}
library(nowcasting)
data(USGDP)
```

USDGP is a list with two data frames:

* `USGDP$base`: this is an unbalanced panel with the model´s variables;
* `USGDP$legend`: this contains the legend with the relevant information for each variable in the `USGDP$base` dataframe.

In order to find the GDP time series in the dataset one can use the legend in `USGDP$legend`. We transform GDP into a quarterly variable by using the `month2qtr` function. The explanatory variables are all the remaining variables from the dataset.

```{r warning=FALSE}
gdp <- month2qtr(x = USGDP$base[,"RGDPGR"], reference_month = 3)
```

The second step consists in treating the explanatory variables. This can be done by using the `Bpanel` function. This function creates a balanced panel using an unbalanced panel as input. The missing observations and outliers are substituted using the outlier correction methodology from *Giannone et al. (2008)*. The function includes most usual transformations to obtain stationary variables. The object `USGDP$legend` contains all the transformations used in *Giannone et al. (2008)*. The `Bpanel` function also allows monthly variables to be aggregated to represent quarterly quantities.

```{r warning=FALSE}
gdp_position <- which(colnames(USGDP$base) == "RGDPGR")
base <- Bpanel(base = USGDP$base[,-gdp_position], 
               trans = USGDP$legend$Transformation[-gdp_position], aggregate = TRUE)
```
Once these variables have been treated, the `nowcast` function can be used to estimate the model parameters according to the estimation method selected, the number *r* of dynamic factors, the lag order of the factors *p* and the number *q* of shocks to the factors. For this example we use the *Two-Stage - quarterly factors* method. The arguments *r*, *p* and *q* were defined according to *Giannone et al. (2008)*.

```{r warning=FALSE}
nowcastUSGDP <- nowcast(y = gdp, x = base, r = 2, p = 2, q = 2, method = '2s')
```
The in sample evaluation from *Giannone et al. (2008)* could be reproduced by looking at the ACF of the residuals of the model specified above.

```{r warning=FALSE}
res <- ts(nowcastUSGDP$reg$residuals, start = start(gdp), frequency = 4)
acf(window(res, start = c(1985,1), end = c(2004,4)))
```
The **results** can be accessed from the object `nowcastUSGDP`.

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

The **graphs** available with the `nowcast.plot` function allow the to visualize some results of interest.

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
### nowcasting US GDP 2/2

In this example we work with the data the Federal Reserve of New York uses in its weekly nowcasting report. The explanatory variables are mixed frequencies including both monthly and quarterly series. 

```{r warning=FALSE}
library(nowcasting)
data(NYFED)
```
Similarly to the previous working example, the object *NYFED* contains all the necessary information to run the `nowcast` function. The block structure, the transformations to make the variables stationary and the frequencies of the variables can be loaded as illustrated below. 

```{r warning=FALSE}
base <- NYFED$base
blocks <- NYFED$blocks$blocks
trans <- NYFED$legend$Transformation
frequency <- NYFED$legend$Frequency
```
The same setting as the NY FED is used. We therefore limit the number of factors, r, per block to one and define the factor process as a VAR(1). The algorithm displays the convergence of the loglikelihood function every 5 iterations.

```{r warning=FALSE}
nowEM <- nowcast(y = "GDPC1", x = x, r = 1, p = 1, method = "EM", blocks = blocks, frequency = frequency)
```
The forecasts can be visualized using the function `nowcast.plot` as illustrated below.

```{r warning=FALSE}
nowcast.plot(nowEM)
```

### Nowcasting Brazilian GDP using vintages

A vintage is a dataset observed on a specific date. The latter is useful to evaluate the out of sample performance of our model. The **nowcasting** package contains a dataset of Brazilian economic time series.

```{r warning=FALSE}
library(nowcasting)
data(BRGDP)
```

The `PRTDB` function is intended to construct pseudo real time vintages of any dataset. The function excludes observations from time series based on the lag information provided by the user and simulates what would be observed on the reference date. In this case we have construct a 10 year database ending in 2015-06-01.

```{r warning=FALSE}
vintage <- PRTDB(mts = BRGDP$base, delay = BRGDP$delay, vintage = "2015-06-01")
base <- window(vintage, start = c(2005,06), frequency = 12)
x <- Bpanel(base = base, trans = BRGDP$trans)
```

The variable to be forecast is then made stationary. We also use the `month2qtr` function to cast GDP as a quaterly variable.

```{r warning=FALSE}
GDP <- base[,which(colnames(base) == "PIB")]
GDP_qtr <- month2qtr(x = GDP,reference_month = 3)
y <- diff(diff(GDP_qtr), 4)
```
Information criteria can be used in order to help determine the number of factors *r* and shocks to the factors *q* that the model should have. 

```{r warning=FALSE}
ICR1 <- ICfactors(x = x, type = 1)
ICR2 <- ICfactors(x = x, type = 2)
ICQ1 <- ICshocks(x = x, r = 2, p = 2)
```
The user is now ready to forecast the variable of interest. The summary of the regression can be accessed as illustrated below.

```{r warning=FALSE}
now <- nowcast(y = y, x = x, r = 2, q = 2 , p = 2)
summary(now$reg)
```
Finally the in- and out of sample forecasts for this particular vintage can be visualized using the `nowcast.plot` function.

```{r warning=FALSE}
nowcast.plot(now, type = "fcst")
```



