% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fnets.R
\name{predict.fnets}
\alias{predict.fnets}
\title{Forecasting by fnets}
\usage{
\method{predict}{fnets}(
  object,
  newdata = NULL,
  n.ahead = 1,
  fc.restricted = TRUE,
  r = c("ic", "er"),
  ...
)
}
\arguments{
\item{object}{\code{fnets} object}

\item{newdata}{input time series matrix; by default, uses input to \code{object}.
Valid only for the case where \code{newdata} is modelled as a VAR process without any factors}

\item{n.ahead}{forecasting horizon}

\item{fc.restricted}{whether to forecast using a restricted or unrestricted, blockwise VAR representation of the common component}

\item{r}{number of static factors, or a string specifying the factor number selection method when \code{fc.restricted = TRUE};
 possible values are:
\describe{
   \item{\code{"ic"}}{ information criteria of Alessi, Barigozzi & Capasso (2010)}
   \item{\code{"er"}}{ eigenvalue ratio of Ahn & Horenstein (2013)}
}}

\item{...}{not used}
}
\value{
a list containing
\item{forecast}{ forecasts for the given forecasting horizon}
\item{common.pred}{ a list containing forecasting results for the common component}
\item{idio.pred}{ a list containing forecasting results for the idiosyncratic component}
\item{mean.x}{ \code{mean.x} argument from \code{object}}
}
\description{
Produces forecasts of the data for a given forecasting horizon by
separately estimating the best linear predictors of common and idiosyncratic components
}
\examples{
out <- fnets(data.restricted, q = 2, do.lrpc = FALSE, var.args = list(n.cores = 2))
pre.unr <- predict(out, fc.restricted = FALSE)
pre.res <- predict(out, fc.restricted = TRUE)
}
\seealso{
\link[fnets]{fnets}
}
