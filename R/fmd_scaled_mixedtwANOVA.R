#' Allometric scaling for flow-mediated dilation: Mixed two-way ANOVA
#'
#' This function calculates allometrically-scaled flow-mediated dilation responses for a mixed two-way ANOVA study design and returns pairwise comparisons between groups with a conservative Tukey post-hoc correction. The function will also return transformed means, standard errors, and 95% confidence intervals.
#'
#' @param dat Data frame object; does not have to be named 'dat', but the columns that correspond to peak artery diameter, baseline artery diameter, group (or condition), and time must be labeled "dpeak", "dbase", "group", and "time", respectively.
#'
#' @return This function returns the following:
#'       \item{modelcoef}{A dataframe that contains the model coeficients from a linear mixed effects model}
#'       \item{tf.emmeans}{Backtransformed estimated marginal means}
#'       \item{maingroup}{Frequentist comparisons accross groups or conditions}
#'       \item{maintime}{Frequentist comparisons accross times or conditions.}
#'       \item{interactiongroup}{Linear and quadratic contrasts of group for each time. From the `emmeans` package.}
#'       \item{interactiontime}{Linear and quadratic contrasts of time for each group. From the `emmeans` package.}
#'       \item{all.comparisons}{Pairwise comparisons accross main effects and interactions. Treated with Tukey post-hot correction method and Kenward-Roger degrees of freedom method. Output from the `emmeans` package.}
#'       \item{plot}{Graphical representation using `ggplot` syntax}
#'
#' @keywords FMD; scaled; allometric
#' @export
#' @examples
#' Simulated exercise intervention data consisting of an independent grouping (CON, MICT, SIT) x repeated measures at different times (0, 6 weeks):
#'
#' dat <- data.frame(dbase = rnorm(300, 4.1, 1.2), dpeak = rnorm(300, 4.3, 1.1), group = as.factor(sort(rep(c("CON", "SIT", "MICT"), 100))), time = as.factor(rep(sort(rep(c("0wk", "6wk"), 50)), 3)), pid = as.factor(c(rep(1:50, 2), rep(51:100, 2), rep(101:150, 2))))
#' dat <- dat |>
#' dplyr::arrange(pid)
#'
#' res <- fmd_scaled_mixedtwANOVA(dat)
#'
#' @import dplyr
#' @import lmerTest
#' @import emmeans
#' @import ggplot2

fmd_scaled_mixedtwANOVA <- function(dat){
  log_dbase <- log(dat[["dbase"]])
  log_dpeak <- log(dat[["dpeak"]])
  difflogs <- log_dpeak - log_dbase
  df <- data.frame(log_dbase = log_dbase,
                   log_dpeak = log_dpeak,
                   difflogs = difflogs)
  newdat <- cbind(dat, df)

  fit <- lm(log_dpeak ~ log_dbase, data = df)

  if(fit$coefficients[2] != 1 & confint(fit, level = 0.95)[2,2] < 1){

    # Assign to mixed effects model
    fit2 <- lmer(difflogs ~ log_dbase + group*time + (1|pid), data = newdat)
    lmcoef.aov <- anova(fit2, type = "III", ddf = "Kenward-Roger")
    em.means <- emmeans(fit2, "group", "time")
    main.group.p <- emmeans(fit2, pairwise ~ group | time)
    main.group.ci <- summary(main.group.p, infer = c(TRUE, FALSE))
    main.time.p <- emmeans(fit2, pairwise ~ time | group)
    main.time.ci <- summary(main.time.p, infer = c(TRUE, FALSE))
    em.means.cont.p <- pairs(em.means, adjust = "tukey")
    em.means.cont.ci <- summary(em.means.cont.p, infer = c(TRUE, FALSE))
    interaction.group <- contrast(main.group.p[[1]], interaction = c("poly", "consec"), by = NULL)
    interaction.time <- contrast(main.time.p[[1]], interaction = c("poly", "consec"), by = NULL)
    allcomps <- emmeans(fit2, pairwise ~ time * group)

    # Now create a table for transformed emmeans...
    tf.em.means <- data.frame("Group" = (as.data.frame(em.means)[1]),
                              "Time" = (as.data.frame(em.means)[2]),
                              "Estimated Marginal Means" = ((exp(as.data.frame(em.means)[3])-1)*100),
                              "Std. Error" = ((exp(as.data.frame(em.means)[4])-1)*100),
                              "LowerCI" = ((exp(as.data.frame(em.means)[6])-1)*100),
                              "UpperCI" = ((exp(as.data.frame(em.means)[7])-1)*100))

    tf.em.means <- tf.em.means %>%
      rename("Group" = group,
             "Time" = time,
             "Est.MarginalMeans" = emmean,
             "Std.Error" = SE,
             "LowerCI" = lower.CL,
             "UpperCI" = upper.CL)

    # Create a table for transformed emmeans stats
    tf.stats <- data.frame("contrast" = (as.data.frame(em.means.cont.p)[,1]),
                           "time" = (as.data.frame(em.means.cont.p)[,2]),
                           "std.error" = ((exp(as.data.frame(em.means.cont.p)[,4])-1)*100),
                           "df" = (as.data.frame(em.means.cont.p)[,5]),
                           "t.ratio" = (as.data.frame(em.means.cont.p)[,6]),
                           "p-value" = (as.data.frame(em.means.cont.p)[,7]),
                           "lower.95CI" = ((exp(as.data.frame(em.means.cont.ci)[,6])-1)*100),
                           "upper.95CI" = ((exp(as.data.frame(em.means.cont.ci)[,7])-1)*100))

    # Create a table for main group contrasts
    maingroupp <- as.data.frame(main.group.p$contrasts)
    maingroupci <- as.data.frame(main.group.ci$contrasts)
    main.group.contrasts.table <- data.frame("contrast" = as.data.frame(maingroupp[1]),
                                             "time" = as.data.frame(maingroupp[2]),
                                             "SE" = (exp(as.data.frame(maingroupp[4]))-1)*100,
                                             "df" = as.data.frame(maingroupp[5]),
                                             "t.ratio" = as.data.frame(maingroupp[6]),
                                             "p.value" = (as.data.frame(maingroupp[7])),
                                             "lower.95CI" = (exp(as.data.frame(maingroupci[6]))-1)*100,
                                             "upper.95CI" = (exp(as.data.frame(maingroupci[7]))-1)*100)

    # Create a table for main time contrasts
    maintimep <- as.data.frame(main.time.p$contrasts)
    maintimeci <- as.data.frame(main.time.ci$contrasts)
    main.time.contrasts.table <- data.frame("contrast" = as.data.frame(maintimep[1]),
                                            "time" = as.data.frame(maintimep[2]),
                                            "SE" = (exp(as.data.frame(maintimep[4]))-1)*100,
                                            "df" = as.data.frame(maintimep[5]),
                                            "t.ratio" = as.data.frame(maintimep[6]),
                                            "p.value" = (as.data.frame(maintimep[7])),
                                            "lower.95CI" = (exp(as.data.frame(maintimeci[6]))-1)*100,
                                            "upper.95CI" = (exp(as.data.frame(maintimeci[7]))-1)*100)

    plot <-
      ggplot(data = tf.em.means, aes(x = Group, y=Est.MarginalMeans, group = 1)) +
      geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0, size =1.5, colour = "dodgerblue1", alpha = 0.6) +
      geom_point(size = 2) +
      ylab("Scaled FMD") +
      xlab("Group") +
      facet_wrap(~ Time) +
      theme_bw()

    print("-------------------------------------------------------------------")
    print("Linear Model Coefficients")
    print(lmcoef.aov)
    print("-------------------------------------------------------------------")
    print("Backtransformed estimated Marginal Means")
    print(tf.em.means)
    print("-------------------------------------------------------------------")

    # Create a list of output objects:
    value <- list(
      model.coef = lmcoef.aov,
      contrast.time = main.time.contrasts.table,
      contrast.group = main.group.contrasts.table,
      interaction.group = interaction.group,
      interaction.time = interaction.time,
      transformed.emmeans = tf.em.means,
      all.comparisons = allcomps$contrasts,
      plot = plot
    )
    attr(value, "class") <- "fmd_scaled_mixedtwANOVA"
    value

  }
  else {
    print("Allometric scaling not required because: (i) the unstandardized regression coefficient (beta) does not deviate from 1; (ii) the 95% CI was calculated to have an upper limit greater than 1")
  }
}
