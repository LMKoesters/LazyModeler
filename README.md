# LazyModeler

`LazyModeler` is a statistical package for the programming language R that allows users to easily perform regression modeling. It includes removal of autocorrelated variables, choice between several types of (non)linear regression models, standard stepwise model simplification, various model quality checks, plotting of coefficient estimates and relationships, and output generation.

# Overview and major functions

`LazyModeler` automatizes all necessary steps needed for use of (non)linear regression models. It comprises three major functions that are included within the main function `optimize_model`.

The first major function `remove_autocorrelations` checks for any autocorrelations (\|r\| \> 0.7) (Dormann et al. 2013) given a list of variables sorted by relevance. Automatic removal of these autocorrelations is possible through the use of a function parameter. Removal will follow the order of the list of variables, ensuring that the user's expertise on the importance of features is respected. A named list is returned with a) a vector containing all removed predictors, and b) a dataframe listing autocorrelations and information on deleted variables.

The main function provides the model formula to the second major function `simplify_model`. If autocorrelations were detected, the formula is updated accordingly. The regression model is then calculated. Options for the models are: `lm`, `glm`, `lmer`, `glmer`, `gam`, or `nlmer`, with all possible distributions of the response variable being allowed. Stepwise backward simplification or forward model selection takes place using an iterative process where each time the metric(s) specified by the user are applied on the model to check whether further simplification/selection is needed. Main variables are kept when they are involved in interactions. Options for the metrics are: `aov`, `aic`, `aicc`, or `bic`. The final model is returned to the main function alongside its metadata as well as simplification history if requested by the user.

Using the third major function `fancy_plotting`, the final model then undergoes multiple visualization steps. Plots to assess model quality are created using the standard plot function available through base R, or model check included in the `performance` R package (Lüdecke et al. 2021). Furthermore, the script produces regression, box, or violin plots for each numerical or categorical coefficient as well as plots depicting effects sizes and estimates. All generated plots are returned to the user within a named list. The main function additionally returns the output of both the model simplification/selection and autocorrelation functions as well as the summary of the final model.

`LazyModeler` makes use of the R package `corrplot` (Wei and Simko 2021) to calculate correlations between variables, `lme4` (Bates et al. 2024) and `lmerTest` (Kuznetsova, Brockhoff, and Christensen 2017) for regression modeling, `tidyverse` (Wickham
et al. 2019) for data handling, and `MuMIn` (Bartoń 2024) for calculation of AICc scores. For generation of plots visualizing regression, effect size, and estimates, the script further leverages `tidyverse` and color palettes included in the `colorspace` (Zeileis et al. 2020) and `viridis` (Garnier et al. 2024) R packages.

# Example

``` r

# import example data
data(plants)

# check data structure
str(plants)
summary(plants)

# testing dataset (subset) based on Karbstein et al. 2021
#(https://onlinelibrary.wiley.com/doi/10.1111/mec.15919)

results_example <- optimize_model(plants, quote(sexual_seed_prop ~
altitude + latitude_gps_n + longitude_gps_e + (solar_radiation +
annual_mean_temperature + isothermality)^2 + I(isothermality^2) +
habitat + ploidy),  autocorrelation_cols = c("solar_radiation",
"annual_mean_temperature", "isothermality", "altitude",
"latitude_gps_n", "longitude_gps_e"), automatic_removal=TRUE,
autocorrelation_threshold = 0.8, correlation_method="spearman",
model_type = "glm", model_family = "quasibinomial",
assessment_methods=c("anova"), simplification_direction="backward",
omit.na="overall", scale_predictor=TRUE,
plot_quality_assessment="performance", round_p=3,
cor_use="complete.obs", plot_relationships=TRUE, jitter_plots=TRUE,
plot_type="violinplot",  stat_test="wilcox",
backward_simplify_model=TRUE, trace=TRUE)
```

![Navigating through the output. For example, (a) simply click on dataframe button highlighted with a red arrow to (b) illustrate the final model output.](paper/assets/figure1.png)
Navigating through the output. For example, (a) simply click on dataframe button highlighted with a red arrow to (b) illustrate the final model output.

![(a) Model quality check and (b,c) exemplary output plots of significant relationships.](paper/assets/figure2.png)
(a) Model quality check and (b,c) exemplary output plots of significant relationships.


# References
Bartoń, Kamil. 2024. *MuMIn: Multi-Model Inference*.
<https://CRAN.R-project.org/package=MuMIn>.

Bates, D., M. Maechler, B. Bolker, and S. Walker. 2024.
“<span class="nocase">lme4 - Linear mixed-effects models using ’Eigen’
and S4</span>.” <https://github.com/lme4/lme4/>.

Dormann, Carsten F, Jane Elith, Sven Bacher, Carsten Buchmann, Gudrun
Carl, Gabriel Carré, Jaime R. García Marquéz, et al. 2013.
“<span class="nocase">Collinearity: a review of methods to deal with it
and a simulation study evaluating their performance</span>.” *Ecography*
36 (1): 27–46. <https://doi.org/10.1111/j.1600-0587.2012.07348.x>.

Garnier, Simon, Ross, Noam, Rudis, Robert, Camargo, et al. 2024.
*<span class="nocase">viridis(Lite)</span> - Colorblind-Friendly Color
Maps for r*. <https://doi.org/10.5281/zenodo.4679423>.

Lüdecke, Daniel, Mattan S. Ben-Shachar, Indrajeet Patil, Philip
Waggoner, and Dominique Makowski. 2021. “Performance: An r Package for
Assessment, Comparison and Testing of Statistical Models.” *Journal of
Open Source Software* 6 (60): 3139.
<https://doi.org/10.21105/joss.03139>.

Kuznetsova, Alexandra, Per B. Brockhoff, and Rune H. B. Christensen.
2017. “<span class="nocase">lmerTest</span> Package: Tests in Linear
Mixed Effects Models.” *Journal of Statistical Software* 82 (13): 1–26.
<https://doi.org/10.18637/jss.v082.i13>.

Wei, Taiyun, and Viliam Simko. 2021. *R Package ’Corrplot’:
Visualization of a Correlation Matrix*.
<https://github.com/taiyun/corrplot>.

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019.
“Welcome to the <span class="nocase">tidyverse</span>.” *Journal of Open
Source Software* 4 (43): 1686. <https://doi.org/10.21105/joss.01686>.

Zeileis, Achim, Jason C. Fisher, Kurt Hornik, Ross Ihaka, Claire D.
McWhite, Paul Murrell, Reto Stauffer, and Claus O. Wilke. 2020.
“<span class="nocase">colorspace</span>: A Toolbox for Manipulating and
Assessing Colors and Palettes.” *Journal of Statistical Software* 96
(1): 1–49. <https://doi.org/10.18637/jss.v096.i01>.

