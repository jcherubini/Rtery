# Rtery: An R package for endothelial function analysis
<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

`Rtery` is an R package that hosts several functions to support statistical analyses related to arterial structure and function. `Rtery` was created with a user-friendly objective and specifically designed for users whose native programming language may not be similar to R syntax. However, even those familiar with R syntax may find convenience in some of the functions offered by `Rtery`!

<p align="center">
  <img src="https://github.com/jcherubini/Rtery/blob/main/Figures/RteryLogo.png" width="175" height="175">
</p>
<p align="center">
  <em>Figure 1: The Rtery emblem.</em>
</p>

You can install `Rtery` into your local R environment using the following code:

    devtools::install_github("jcherubini/Rtery")

## **Functionality and capabilities**

`Rtery` supports a variety of functions related to the analysis of endothelial function and arterial caliber. These functions can be summarized into the following domains:
- Absolute and relative flow-mediated dilation
- Allometric scaling of flow-mediated dilation
- Calculating flow-mediated dilation normalized to shear stimuli (coming soon!)
- Hemodynamic variables (coming soon!)

#### Flow-mediated dilation (FMD)
Endothelial function can be estimated by calculating the flow-mediated dilation (FMD) of an artery in response to arterial shear (Thijssen et al. 2019). FMD can be calculated using absolute values, and expressed as the difference between peak arterial diameter (`dpeak`) and baseline arterial diameter (`dbase`): `FMD = dpeak - dbase`. FMD can also be calculated as a percent-change (%FMD) relative to baseline diameter. The %FMD can be calculated as the difference between peak diameter and baseline diameter, relative to the baseline diameter of an artery: `%FMD = ((dpeak - dbase)/dbase)*100`.

The `Rtery::fmd()` function can calculate both absolute and relative FMD for given 'dpeak' and 'dbase' values, and concatenate the values as a column into a data frame for subsequent analysis. Users can specify the type of FMD that is to be calculated using the `type` argument: `Rtery::fmd(data, type = c("absolute"))`. Several calculations of FMD can thereby be obtained using the `Rtery` package. 

#### FMD normalized to shear stimulus: 
Baseline arterial diameters can differ between individuals. Thus, arteries can experience different magnitudes of shear stimuli depending on arterial diameter. Therefore, it is reccomended to normalize %FMD to the magnitude of shear stimulus in an artery (Pyke and Tschakovsky, 2005). The `Rtery` package will soon conveniently offer the function `Rtery::fmd_shear()` to normalize %FMD to given shear rates, and concatenate the normalized FMD values as a column into a data frame for analysis. 

#### Allometric scaling of FMD:

The %FMD is a simple and canonical method to estimate endothelial function. However,as Atkinson and Batterham (2013) note, statisticians caution against using change scores to inform physiological parameters. Covarying by baseline arterial diameter is an insufficient treatment to surmount the baseline diameter dependency of %FMD scores (Atkinson and Batterham, 2013). As a consequence, allometric scaling of %FMD is sometimes necessary and scaled FMD responses may indeed be a preferable alternative to the traditional %FMD calculation. 

Allometric scaling may be technically challenging to perform and can vary across study designs and statistical software, and it can be difficult to ascertain when allometric treatments should be applied. `Rtery` offers a quick and convenient method to first, test for the necessity of allometric scaling and, if necessary, calculate allometrically-scaled FMD responses for a variety of different study designs using functions from the `lmerTest` (Kuznetsova et al. 2020) and `emmeans` (Lenth, 2022) packages (Figure 2).

<p align="center">
  <img src="https://github.com/jcherubini/Rtery/blob/main/Figures/RteryProcess.png" width="350" height="175" alt="Alternative text caption here">
  </p>
<p align="center">
  <em>Figure 2: Rtery uses a three-step computational process for allometric treatment of FMD data.</em>
</p>

Consider the following experimental design with data simulated for three independent groups: Normal sleep (NS), sleep deprivation (SD), and sleep restriction (SR). A one-way ANOVA study design may be used to compare the %FMD scores across three groups, and thus, the `Rtery::fmd_scaled_owANOVA()` (owANOVA = **o**ne-**w**ay **an**alysis **o**f **va**riance) function would be used to examine statistical differences across the groups using frequentist methods.

    set.seed(2021)
    dat <- data.frame(pid = as.factor(rep(1:30)), dpeak = rnorm(30, 4.27, 1.12), dbase = rnorm(30, 4.16, 1.21)*0.8, group = as.factor(c(rep("NS", 10), rep("SR", 10), rep("SD", 10))))  
    res <- Rtery::fmd_scaled_owANOVA(dat)
    
R will compute the differences in scaled FMD for all three groups and return the corresponding model output, as well as the estimated marginal means, p-value, model standard error, and 95% confidence interval for each condition.

Using functions provided by the `lmerTest` package, `Rtery` can also handle data that violate assumptions of independence and include instances of repeated measures. For instance, consider an exercise training intervention that compares %FMD across three difference exercise interventions, sprint interval training (SIT), moderate intensity continuous training (MICT), and a control group (CON). A mixed two-way ANOVA (group x time) may be used to treat this data. The corresponding `Rtery` call would therefore be `fmd_scaled_mixedtwANOVA()`. Again, simulating data:

    dat <- data.frame(dbase = rnorm(300, 4.1, 1.2), dpeak = rnorm(300, 4.3, 1.1), group = as.factor(sort(rep(c("CON", "SIT", "MICT"), 100))), time = as.factor(rep(sort(rep(c("0wk", "6wk"), 50)), 3)), pid = as.factor(c(rep(1:50, 2), rep(51:100, 2), rep(101:150, 2))))
    dat <- dat |>
      dplyr::arrange(pid)
    res <- Rtery::fmd_scaled_mixedtwANOVA(dat)
    
All necessary statistics and contrasts will be generated and reported in a table, again, courtesy of the `emmeans` and `lmerTest` packages. Additional contrast outputs can be accessed using the `$` operator after calling the function, `res <- fmd_scaled_mixedtwANOVA`, into a results object.

## Additional hemodynamic features soon to be included in `Rtery`

The `Rtery` package will include functions that calculate multiple hemodynamic parameters provided that the user inputs the constituent variables. For instance, `Rtery` may calculate arterial distensibility `dist()`, pulse wave velocity `pwv()`, and pulse pressure `pp()`. These functions may be helpful for individuals who wish to generate a list of parameters and concatenate into a data frame for analysis.  

## Future directions for the `Rtery` package
`Rtery` users can expect a maintained and reflexive analysis experience. The `Rtery` development team is currently working on several extensions to the package such as advanced interaction analyses through the `emmeans` package, and a Shiny web application hosted by the R ecosystem. As a consequence, `Rtery` will retain its experimental status until the development team is sufficiently satisfied to pursue a stable status in the future. 

Interested in joining the `Rtery` team, contributing to the package, or sharing an idea? All are encouraged to submit a pull request or contact jcherubini29@gmail.com to inquire further.

## References

1. Thijssen DH, Bruno RM, van Mil AC, Holder SM, Faita F, Greyling A, Zock PL, Taddei S, Deanfield JE, Luscher T, Green DJ. **Expert consensus and evidence-based recommendations for the assessment of flow-mediated dilation in humans.** European heart journal. 2019 Aug 7;40(30):2534-47.
2. Pyke KE, Tschakovsky ME. **The relationship between shear stress and flow‐mediated dilatation: implications for the assessment of endothelial function.** The Journal of physiology. 2005 Oct 1;568(2):357-69.
3. Atkinson G, Batterham AM. **The percentage flow-mediated dilation index: a large-sample investigation of its appropriateness, potential for bias and causal nexus in vascular medicine.** Vascular medicine. 2013 Dec;18(6):354-65.
4. Lenth RV. **emmeans: Estimated Marginal Means, aka Least-Squares Means.** R
  package version 1.7.2. 2022: https://CRAN.R-project.org/package=emmeans
5. Kuznetsova A, Brockhoff PB, Christensen RHB. **lmerTest Package: Tests in Linear
Mixed Effects Models.** Journal of Statistical Software, 2017 82(13), 1-26. doi:
10.18637/jss.v082.i13 (URL: https://doi.org/10.18637/jss.v082.i13). 

--------------------------------------------------------------------------------------------------------------------------------------------------------

*Disclaimer: All `Rtery` functions are currently in beta-testing and may not yet provide optimal computation and structure. Please contact the package creator, Joshua Cherubini (jcherubini29@gmail.com), for any questions related to `Rtery` use or development. Images created using Biorender.*
