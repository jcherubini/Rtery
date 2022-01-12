#' Allometric scaling for flow-mediated dilation: Independent sample t-tests
#'
#' This function calculates allometrically-scaled flow-mediated dilation responses for independent sample t-tests and returns the corresponding comparisons between groups. The function will also return back-transformed means, standard errors, and 95% confidence intervals.
#' @param dat data frame object; does not have to be named 'dat', but the columns that correspond to peak artery diameter, baseline artery diameter, and group (or condition) must be labeled "dpeak", "dbase", and "group", respectively.
#' @keywords FMD; scaled; allometric
#' @export
#' @examples
#' Simulated data with peak and baseline artery diameters from 12 participants:
#'
#' EXAMPLE 1: Allometric scaling not required
#' dat <- data.frame(participant = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12"), dpeak = c(4.1, 3.2, 6.5, 5.9, 4.3, 2.1, 4.0, 6.3, 3.3, 4.9, 5.2, 7.1), dbase = c(4.0, 2.9, 6.0, 5.7, 4.1, 2.0, 3.7, 6.2, 3.2, 4.3, 5.1, 6.9), group = as.factor(c("NS", "NS", "NS", "NS", "SD", "SD", "SD", "SD", "SR", "SR", "SR", "SR")))
#' res <- fmd_scaled_ttest(dat)
#'
#' EXAMPLE 2: Allometric scaling required and calculated:
#' dat <- data.frame(dpeak = rnorm(30, 4.27, 1.12), dbase = rnorm(30, 4.16, 1.21)*0.8, group = as.factor(c(rep("NS", 10), rep("SD", 10), rep("SR", 10))))
#' res <- fmd_scaled_ttest(dat)
#'
#' @import dplyr
#' @import emmeans
#' @import car
#' @import ggplot2

fmd_scaled_ttest <- function(dat){
# Take log of data:
  log_dbase <- log(dat[["dbase"]])
  log_dpeak <- log(dat[["dpeak"]])
  difflogs <- log_dpeak - log_dbase
  df <- data.frame(log_dbase = log_dbase, log_dpeak = log_dpeak, difflogs = difflogs)
  newdat <- cbind(dat, df)

  fit <- lm(log_dpeak ~ log_dbase, data = df)

  if(fit$coefficients[2] != 1 & confint(fit, level = 0.95)[2,2] < 1){
    # Assign to univariate general linear model
    fit2 <- lm(difflogs ~ group + log_dbase, data = newdat)
    lmcoef.aov <- car::Anova(fit2, type = "III", ddf = "Kenward-Roger")
    fit3 <- emmeans(fit2, "group")
    emmeans.cont.p <- pairs(fit3, adjust = "tukey") #pval
    emmeans.cont.ci <- summary(emmeans.cont.p, infer = c(TRUE, FALSE)) #ci

    # Table for backtransformed estimated marginal means and measure of error:
    tf.em.means <- as.data.frame(c("group" = (as.data.frame(fit3)[1]), (exp((summary(fit3)[2]))-1)*100, (exp((summary(fit3)[3]))-1)*100, (exp((summary(fit3)[5]))-1)*100, (exp((summary(fit3)[6]))-1)*100))

    tf.em.means <- tf.em.means %>%
      rename("Group" = group.group,
             "EstMarginalMeans" = emmean,
             "SE" = SE,
             "LL" = lower.CL,
             "UL" = upper.CL) #rename columns

    # Table for contrasts
    tf.stats <- data.frame("contrast" = (as.data.frame(emmeans.cont.p)[,1]),
                           "std.error" = ((exp(as.data.frame(emmeans.cont.p)[,3])-1)*100),
                           "df" = (as.data.frame(emmeans.cont.p)[,4]),
                           "t.stat" = (as.data.frame(emmeans.cont.p)[,5]),
                           "p-value" = (as.data.frame(emmeans.cont.p)[,6]),
                           "lower.95CI" = ((exp(as.data.frame(emmeans.cont.ci)[,5])-1)*100),
                           "upper.95CI" = ((exp(as.data.frame(emmeans.cont.ci)[,6])-1)*100))

    plot <- ggplot(tf.em.means, aes(x = Group, y=EstMarginalMeans, group = 1)) +
      geom_errorbar(aes(ymin = LL, ymax = UL), width = 0, size =1.5, colour = "dodgerblue1", alpha = 0.5) +
      geom_point(size = 2) +
      ylab("Scaled FMD") +
      xlab("Group") +
      theme_bw() +
      theme(axis.text.x = element_text(face = "bold", size = 12),
                     axis.title.x = element_text(face = "bold", size = 16),
                     axis.title.y = element_text(face = "bold", size = 16),
                     axis.text.y = element_text(face = "bold", size = 12))

    print("------------------------------------------------------")
    print("Linear Model Coefficients")
    print(lmcoef.aov)
    print("------------------------------------------------------")
    print("Backtransformed estimated Marginal Means")
    print(tf.em.means)
    print("------------------------------------------------------")

    # Create a list of output objects
    value <- list(
      model.coef = lmcoef.aov,
      transformed.emmeans = tf.em.means,
      plot = plot,
      contrasts = tf.stats
    )
    attr(value, "class") <- "fmd_scaled_ttest"
    value

  }
  else {
    print("Allometric scaling not required because: (i) the unstandardized regression coefficient (beta) does not deviate from 1; (ii) the 95% CI was calculated to have an upper limit greater than 1")
    }
}

