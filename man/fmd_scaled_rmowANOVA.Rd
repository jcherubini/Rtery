% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmd_scaled_rmowANOVA.R
\name{fmd_scaled_rmowANOVA}
\alias{fmd_scaled_rmowANOVA}
\title{Allometric scaling for flow-mediated dilation: Repeated measures one-way ANOVA}
\usage{
fmd_scaled_rmowANOVA(dat)
}
\arguments{
\item{dat}{data frame object; does not have to be named 'dat', but the columns that correspond to peak artery diameter, baseline artery diameter, and group (or condition) must be labelled "dpeak", "dbase", and "group", respectively.}
}
\value{
This function returns the following:
\item{model.coef}{A dataframe continaing coefficients from the linear model}
\item{main.effects}{A dataframe containing main effects and statistical contrasts}
\item{transformed.emmeans}{Backtransformed estimated marginal means and model standard error}
\item{plot}{Graphical representation using \code{ggplot2} syntax; contains means and 95\% confidence intervals}
}
\description{
This function calculates allometrically-scaled flow-mediated dilation responses for repeated measures one-way ANOVA study designs, and returns  comparisons between groups. The function will also return back-transformed means, standard errors, and 95\% confidence intervals. Users must label the data set columns as "dpeak", "dbase", and "group" respectively. This \code{Rtery} function also requires that data be arranged in long format.
}
\examples{
Create a sample data frame (simulated data; as such, values may not correspond with physiological norms) with peak and baseline artery diameters from sample participants:

dat <- data.frame(pid = as.factor(rep(1:10, 3)), dpeak = rnorm(30, 4.27, 1.12), dbase = rnorm(30, 4.16, 1.21)*0.8, group = as.factor(c(rep("NS", 10), rep("SD", 10), rep("SR", 10))))
dat <- dat |>
dplyr::arrange(pid)
res <- fmd_scaled_rmowANOVA(dat)

}
\keyword{FMD;}
\keyword{allometric}
\keyword{scaled;}
