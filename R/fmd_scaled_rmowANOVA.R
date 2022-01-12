#' Allometric scaling for flow-mediated dilation: Repeated measures one-way ANOVA
#'
#' This function calculates allometrically-scaled flow-mediated dilation responses for repeated measures one-way ANOVA study designs, and returns  comparisons between groups. The function will also return back-transformed means, standard errors, and 95% confidence intervals. Users must label the data set columns as "dpeak", "dbase", and "group" respectively. This `Rtery` function also requires that data be arranged in long format.
#' @param dat data frame object; does not have to be named 'dat', but the columns that correspond to peak artery diameter, baseline artery diameter, and group (or condition) must be labelled "dpeak", "dbase", and "group", respectively.
#' @return This function returns the following:
#'       \item{model.coef}{A dataframe continaing coefficients from the linear model}
#'       \item{main.effects}{A dataframe containing main effects and statistical contrasts}
#'       \item{transformed.emmeans}{Backtransformed estimated marginal means and model standard error}
#'       \item{plot}{Graphical representation using `ggplot2` syntax; contains means and 95% confidence intervals}
#' @keywords FMD; scaled; allometric
#' @export
#' @examples
#' Create a sample data frame (simulated data; as such, values may not correspond with physiological norms) with peak and baseline artery diameters from sample participants:
#' library(tidyverse)
#' dat <- tibble::tibble(pid = as.factor(rep(1:10, 3)), dpeak = rnorm(30, 4.27, 1.12), dbase = rnorm(30, 4.16, 1.21)*0.8, group = as.factor(c(rep("NS", 10), rep("SD", 10), rep("SR", 10))))
#' dat <- dat %>%
#' arrange(pid)
#' fmd_scaled_rmowANOVA(dat)

# ----- WRITING THE FUNCTION ----- #

fmd_scaled_rmowANOVA <- function(dat){

  log_dbase <- log(dat[["dbase"]])
  log_dpeak <- log(dat[["dpeak"]])
  difflogs <- log_dpeak - log_dbase
  df <- tibble::tibble(log_dbase = log_dbase, log_dpeak = log_dpeak, difflogs = difflogs)
  newdat <- cbind(dat, df)

  fit <- lm(log_dpeak ~ log_dbase, data = df)

  if(fit$coefficients[2] != 1 & confint(fit, level = 0.95)[2,2] < 1){
    # Assign to mixed effects model
    fit2 <- lmerTest::lmer(difflogs ~ log_dbase + group + (1|pid), data = newdat)
    lmcoef.aov <- anova(fit2, type = "III", ddf = "Kenward-Roger")
    em.means <- emmeans::emmeans(fit2, "group")
    em.means.cont.p <- pairs(em.means, adjust = "tukey")
    em.means.cont.ci <- summary(em.means.cont.p, infer = c(TRUE, FALSE))

    # Now create a table for transformed emmeans...
    tf.em.means <- data.frame("Group" = (as.data.frame(em.means)[1]), "Estimated Marginal Means" = ((exp(as.data.frame(em.means)[2])-1)*100), "Std. Error" = ((exp(as.data.frame(em.means)[3])-1)*100), "LL" = ((exp(as.data.frame(em.means)[5])-1)*100), "UL" = ((exp(as.data.frame(em.means)[6])-1)*100))
    tf.em.means <- tf.em.means %>%
      rename("Group" = group, "Est.MarginalMeans" = emmean, "Std.Error" = SE, "LL" = lower.CL, "UL" = upper.CL)

    # Create a table for transformed emmeans stats...
    tf.stats <- data.frame("contrast" = (as.data.frame(em.means.cont.p)[,1]), "std.error" = ((exp(as.data.frame(em.means.cont.p)[,3])-1)*100), "df" = (as.data.frame(em.means.cont.p)[,4]), "t.stat" = (as.data.frame(em.means.cont.p)[,5]), "p-value" = (as.data.frame(em.means.cont.p)[,6]), "lower.95CI" = ((exp(as.data.frame(em.means.cont.ci)[,5])-1)*100), "upper.95CI" = ((exp(as.data.frame(em.means.cont.ci)[,6])-1)*100))

    plot <- ggplot2::ggplot(tf.em.means, ggplot2::aes(x = Group, y = Est.MarginalMeans, group = 1)) +
      ggplot2::geom_errorbar(aes(ymin = LL, ymax = UL), width = 0, size =1.5, colour = "dodgerblue1", alpha = 0.5) +
      ggplot2::geom_point(size = 2) +
      ggplot2::ylab("Scaled FMD") +
      ggplot2::xlab("Group") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(face = "bold", size = 12),
                     axis.title.x = element_text(face = "bold", size = 16),
                     axis.title.y = element_text(face = "bold", size = 16),
                     axis.text.y = element_text(face = "bold", size = 12))

    print("-------------------------------------------------------------------")
    print("Linear Model Coefficients")
    print(lmcoef.aov)
    print("-------------------------------------------------------------------")
    print("Backtransformed estimated Marginal Means")
    print(tf.em.means)
    print("-------------------------------------------------------------------")

    # Create a list of output objects
    value <- list(
      model.coef = fit2,
      main.effects = lmcoef.aov,
      transformed.emmeans = tf.em.means,
      plot = plot,
      contrasts = tf.stats
    )
    attr(value, "class") <- "fmd_scaled_rmowANOVA"
    value

  }
  else {
    print("Allometric scaling not required because: (i) the unstandardized regression coefficient (beta) does not deviate from 1; (ii) the 95% CI was calculated to have an upper limit greater than 1")
  }
}
