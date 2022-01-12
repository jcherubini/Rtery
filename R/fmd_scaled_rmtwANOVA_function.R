#' Allometric scaling for flow-mediated dilation: Two-way repeated measure ANOVA
#'
#' This function calculates allometrically-scaled flow-mediated dilation responses for two-way repeated measure ANOVA study designs and returns pairwise comparisons between groups with a Tukey post-hoc correction. The function will also return backtransformed means, standard errors, and 95% confidence intervals.
#'
#' @param dat data frame object; does not have to be named 'dat', but the columns that correspond to peak artery diameter, baseline artery diameter, and group (or condition) must be labeled "dpeak", "dbase", and "group", and "time", respectively.
#' @keywords FMD; scaled; allometric
#' @export
#' @examples
#' Simulated data comparing FMD before and after normal sleep (NS), sleep deprivation (SD), and sleep restriction (SR), among 50 participants:
#' dat <- data.frame(dbase = rnorm(300, 4.1, 1.2), dpeak = rnorm(300, 4.3, 1.1), group = as.factor(rep(c("NS", "SR", "SD"), 100)), time = as.factor(rep(sort(rep(c("pre", "post"), 50)), 3)), pid = c(rep(1:50, 6)))
#'
#' dat <- dat |>
#' dplyr::arrange(pid)
#'
#' res <- fmd_scaled_rmtwANOVA(dat)
#'
#' @import dplyr
#' @import geepack
#' @import emmeans
#' @import ggplot2

fmd_scaled_rmtwANOVA <-function(dat){
  log_dbase <- log(dat[["dbase"]])
  log_dpeak <- log(dat[["dpeak"]])
  difflogs <- log_dpeak - log_dbase
  df <- data.frame(log_dbase = log_dbase, log_dpeak = log_dpeak, difflogs = difflogs)
  newdat <- cbind(dat, df)

  fit <- lm(log_dpeak ~ log_dbase, data = df)

  if(fit$coefficients[2] != 1 & confint(fit, level = 0.95)[2,2] < 1){

    # Assign to gee model
    library(geepack)
    fit_gee <- geepack::geeglm(difflogs ~ group + time + group*time + log_dbase,
                               data = newdat,
                               id = pid,
                               corstr = "exchangeable")
    geecoef <- summary(fit_gee)
    maineffects <- anova(fit_gee)
    em.means <- emmeans(fit_gee, "group", "time")
    main.group.p <- emmeans(fit_gee, pairwise ~ group | time)
    main.group.ci <- summary(main.group.p, infer = c(TRUE, FALSE))
    main.time.p <- emmeans(fit_gee, pairwise ~ time | group)
    main.time.ci <- summary(main.time.p, infer = c(TRUE, FALSE))
    em.means.cont.p <- pairs(em.means, adjust = "tukey")
    em.means.cont.ci <- summary(em.means.cont.p, infer = c(TRUE, FALSE))
    interaction.group <- contrast(main.group.p[[1]], interaction = c("poly", "consec"), by = NULL)
    interaction.time <- contrast(main.time.p[[1]], interaction = c("poly", "consec"), by = NULL)
    allcomps <- emmeans(fit_gee, pairwise ~ time * group)

    # Create a table for transformed emmeans
    tf.em.means <- data.frame("Group" = (as.data.frame(em.means)[1]),
                              "Time" = (as.data.frame(em.means)[2]),
                              "Estimated Marginal Means" = ((exp(as.data.frame(em.means)[3])-1)*100),
                              "Std. Error" = ((exp(as.data.frame(em.means)[4])-1)*100),
                              "LL" = ((exp(as.data.frame(em.means)[6])-1)*100),
                              "UL" = ((exp(as.data.frame(em.means)[7])-1)*100))

    tf.em.means <- tf.em.means %>%
      rename("Group" = group,
             "Est.MarginalMeans" = emmean,
             "Std.Error" = SE,
             "LL" = asymp.LCL,
             "UL" = asymp.UCL)

    # Create a table for transformed emmeans stats...
    tf.stats <- data.frame("contrast" = (as.data.frame(em.means.cont.p)[,1]),
                           "time" = (as.data.frame(em.means.cont.p)[,2]),
                           "std.error" = ((exp(as.data.frame(em.means.cont.p)[,4])-1)*100),
                           "df" = (as.data.frame(em.means.cont.p)[,5]),
                           "z.ratio" = (as.data.frame(em.means.cont.p)[,6]),
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

    print("-------------------------------------------------------------------")
    print("GEE effects")
    print(maineffects)
    print("-------------------------------------------------------------------")
    print("Backtransformed estimated marginal means")
    print(tf.em.means)
    print("-------------------------------------------------------------------")


    plot <- ggplot(data = tf.em.means, aes(x = Group, y = Est.MarginalMeans, group = 1)) +
      geom_errorbar(aes(ymin = LL, ymax = UL), width = 0, size =1.5, colour = "dodgerblue1", alpha = 0.6) +
      geom_point(size = 2) +
      ylab("Scaled FMD") +
      xlab("Group") +
      facet_wrap(~ time) +
      theme_bw()

    # Create a list of output objects
    value <- list(
      model.coef = geecoef,
      main.effects = maineffects,
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

