# fnets

Contains methods for network estimation and forecasting for high-dimensional time series under a factor-adjusted VAR model. See 

> _fnets_: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling

by Dom Owens, Haeran Cho and Matteo Barigozzi [arXiv:2301.11675](https://arxiv.org/abs/2301.11675) accompanying the R package, and

> _FNETS_: Factor-adjusted network estimation and forecasting for high-dimensional time series

by Matteo Barigozzi, Haeran Cho and Dom Owens [arXiv:2201.06110](https://arxiv.org/abs/2201.06110) for details of the methodology.


## Installation

To install `fnets` from CRAN:

```
install.packages("fnets")
```


To install from GitHub:

```
devtools::install_github("https://github.com/Dom-Owens-UoB/fnets")
```

## Usage

We can generate an example dataset used in the above paper for simulation studies, by separately generating the factor-driven common component and the idiosyncratic VAR process as
```
set.seed(123)
n <- 500
p <- 50
common <- sim.unrestricted(n, p)
idio <- sim.var(n, p)
x <- common$data + idio$data
```

Fit a factor-adjusted VAR model with `q = 2` factors and `lasso` for VAR transition matrix estimation
```
out <- fnets(x, q = 2, var.order = 1, var.method = "lasso", do.lrpc = FALSE)
```

Plot the Granger network induced by the estimated VAR transition matrices:
```
plot(out, type = "granger", display = "network")
```

Estimate and plot the partial-correlation and long-run partial correlation-based networks:
```
plrpc <- par.lrpc(out)
out$lrpc <- plrpc
out$lrpc.method <- 'par'
plot(out, type = "lrpc", display = "heatmap")
```

Estimate the (long-run) partial correlation-based networks directly using `fnets`:
```
out <- fnets(x, q = 2, var.order = 1, var.method = "lasso", do.lrpc = TRUE)
```

Forecast `n.ahead` steps:
```
pr <- predict(out, n.ahead = 1, common.method = "restricted")
pr$forecast
```






