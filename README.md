# fnets

Factor-adjusted Network Analysis and Forecasting for High-dimensional Time Series

## Installation

To install `fnets` from Github:

```
devtools::install_github("https://github.com/Dom-Owens-UoB/fnets")
```

## Usage

We can generate example data with a common and idiosyncratic component via:
```
set.seed(222)
n <- 200
p<- 100
chi <- sim.factor.M1(n,p)
xi <- sim.idio(n,p)
sample.data <- chi$data + xi$data
```

Fit the model by:
```
model <- fnets(sample.data, q=2, idio.method = "lasso")
```

Predict by:
```
pr <- predict(model,sample.data, common.method = "static")
cpre <- common.predict(model,sample.data, common.method = "static")
ip <- idio.predict(model,sample.data, cpre)
```

Plot the model:
```
plot(model, type = "heatmap")
```
![Granger](figures/model.pdf)

Estimate and plot the long-run partial correlation network:
```
net <- param.lrpc(model, sample.data)
plot(net, type="heatmap")
```
![Omega](figures/omega.pdf)
![Delta](figures/delta.pdf)





