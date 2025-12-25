<div align="center" style="display: flex; justify-content: space-between; align-items: center;">
  <a href="https://github.com/LMKoesters/LazyModeler">
    <img src="https://github.com/LMKoesters/LazyModeler/blob/main/logo_lazymodeler.png?raw=true" 
         alt="LazyModeler Logo" 
         width="400" 
         style="border-radius: 10px;">
  </a>
</div>

# LazyModeler

`LazyModeler` is a statistical package for the programming language R that enables users to easily perform regression modeling. It includes removal of autocorrelated variables, choice between several types of (non)linear regression models, standard stepwise model simplification, various model quality checks, plotting of coefficient estimates and relationships, and output generation.

# Statement of need

Statistical modeling describes the process of finding a mathematical function with specific statistical assumptions that best fits the observed data (Crawley, 2007, 2015; Henley et al., 2020).
This process attempts, in practice, to find a (causal) relationship between a dependent response variable `y` and an independent predictor variable `x` for any postulated hypothesis. For statistical inference and graphics in science, the programming environment R (R Core Team, 2024) has become highly popular.

Our R package `LazyModeler` enables users to automatically remove autocorrelated variables, choose between several types of (non)linear regression models (e.g., LM, GLM, LMER, GLMER, GAM, or NLMER), perform stepwise model simplification, check model quality, plot coefficient estimates and relationships, and generate the output of the final model.

# Overview and major functions

`LazyModeler` automatizes all necessary steps needed for use of (non)linear regression models. It comprises three major functions that are included within the main function `optimize_model`.

The first major function `remove_autocorrelations` checks for any autocorrelations (\|r\| \> 0.7) (Dormann et al. 2013) given a list of variables sorted by relevance. Automatic removal of these autocorrelations is possible through the use of a function parameter. Removal will follow the order of the list of variables, ensuring that the user's expertise on the importance of features is respected. A named list is returned with a) a vector containing all removed predictors, and b) a dataframe listing autocorrelations and information on deleted variables.

The main function provides the model formula to the second major function `simplify_model`. If autocorrelations were detected, the formula is updated accordingly. The regression model is then calculated. Options for the models are: `lm`, `glm`, `lmer`, `glmer`, `gam`, or `nlmer`, with all possible distributions of the response variable being allowed. Stepwise backward simplification or forward model selection takes place using an iterative process where each time the metric(s) specified by the user are applied on the model to check whether further simplification/selection is needed. Main variables are kept when they are involved in interactions. Options for the metrics are: `aov`, `aic`, `aicc`, or `bic`. The final model is returned to the main function alongside its metadata as well as simplification history if requested by the user.

Using the third major function `fancy_plotting`, the final model then undergoes multiple visualization steps. Plots to assess model quality are created using the standard plot function available through base R, or model check included in the `performance` R package (Lüdecke et al. 2021). Furthermore, the script produces regression, box, or violin plots for each numerical or categorical coefficient as well as plots depicting effects sizes and estimates. All generated plots are returned to the user within a named list. The main function additionally returns the output of both the model simplification/selection and autocorrelation functions as well as the summary of the final model.

`LazyModeler` makes use of the R package `corrplot` (Wei and Simko 2021) to calculate correlations between variables, `lme4` (Bates et al. 2024) and `lmerTest` (Kuznetsova, Brockhoff, and Christensen 2017) for regression modeling, `tidyverse` (Wickham
et al. 2019) for data handling, and `MuMIn` (Bartoń 2024) for calculation of AICc scores. For generation of plots visualizing regression, effect size, and estimates, the script further leverages `tidyverse` and color palettes included in the `colorspace` (Zeileis et al. 2020) and `viridis` (Garnier et al. 2024) R packages.

# How to install
You have two options to install `LazyModeler`. You can either install through GitHub using the `remotes` package.
``` r
remotes::install_github("LMKoesters/LazyModeler")
```

Alternatively, you need to download the tarball from GitHub and then install using `install.packages`.
``` r
install.packages("PATH_TO_TARBALL/LazyModeler-v.0.2.0.tar.gz", repos = NULL, type="source")
```

# Example

``` r

# import example data
data(plants)

# check data structure
str(plants)
summary(plants)

# testing dataset (subset) based on Karbstein et al. 2021
# (https://onlinelibrary.wiley.com/doi/10.1111/mec.15919)

results_example <- optimize_model(plants,
    quote(sexual_seed_prop ~
    altitude + latitude_gps_n + longitude_gps_e +
    (solar_radiation + annual_mean_temperature +
    isothermality)^2 + I(isothermality^2) +
    habitat + ploidy),
    autocorrelation_cols = c("solar_radiation",
    "annual_mean_temperature", "isothermality", "altitude",
    "latitude_gps_n", "longitude_gps_e"),
    automatic_removal=TRUE,
    autocorrelation_threshold = 0.8,
    correlation_method="spearman",
    model_type = "glm",
    model_family = "quasibinomial",
    assessment_methods=c("anova"),
    simplification_direction="backward",
    omit.na="overall",
    scale_predictor=TRUE,
    plot_quality_assessment="performance",
    round_p=3,
    cor_use="complete.obs",
    plot_relationships=TRUE,
    jitter_plots=TRUE,
    plot_type="violinplot",
    stat_test="wilcox",
    backward_simplify_model=TRUE,
    trace=TRUE)
```

![Navigating through the output. For example, (a) simply click on dataframe button highlighted with a red arrow to (b) illustrate the final model output.](paper/assets/figure1.png)
Navigating through the output. For example, (a) simply click on dataframe button highlighted with a red arrow to (b) illustrate the final model output.

![(a) Model quality check and (b,c) exemplary output plots of significant relationships.](paper/assets/figure2.png)
(a) Model quality check and (b,c) exemplary output plots of significant relationships.

# Community guidelines
We welcome contributions, feedback, and suggestions to improve this project.
- Encountered a bug or unexpected behavior? Open an issue on GitHub. Just make sure that your issue hasn't been reported yet by checking existing issues before opening a new one.
- Contributions (e.g., code improvements, new features, documentation) are welcome via pull requests. When contributing, describe your changes clearly and provide sufficient context to help us understand your work.

# Important note
The model selection procedures implemented in LazyModeler are provided for convenience and exploratory analysis, and reflect practices recommended in widely used applied statistics literature [@Crawley2007; @Crawley2015]. Users should be aware, however, that statistical inference reported from a model chosen in a data-driven way may be anti-conservative (e.g., p-values may appear smaller than they truly are, confidence intervals narrower). This issue is known as post-selection inference (PSI). Specialized methods have been developed to address it, for instance [@Lee2016], but they are not yet broadly applicable across the full range of model classes supported by LazyModeler. We have implemented PSI for (generalized) linear regression models based on the 'selcorr' R package [@Cattaneo2021], but users are free to use the retained model from LazyModeler for more sophisticated PSI analyses.


# References
Bartoń, K. (2024). *MuMIn: Multi-Model Inference*.
<https://CRAN.R-project.org/package=MuMIn>.

Bates, D., Maechler, M., Bolker, B., and Walker,S. (2024).
“<span class="nocase">lme4 - Linear mixed-effects models using ’Eigen’
and S4</span>.” <https://github.com/lme4/lme4/>.

Bauer, M., & Albrecht, H. (2020). Vegetation monitoring in a 100-year-old calcareous grassland
reserve in Germany. Basic and Applied Ecology, 42, 15–26. <https://doi.org/10.1016/j.baae.2019.11.003>.

Cai, L., Kreft, H., Taylor, A., Denelle, P., Schrader, J., Essl, F., Kleunen, M. van, Pergl, J.,
Pyšek, P., Stein, A., Winter, M., Barcelona, J. F., Fuentes, N., Inderjit, Karger, D. N.,
Kartesz, J., Kuprijanov, A., Nishino, M., Nickrent, D., … Weigelt, P. (2023). Global models
and predictions of plant diversity based on advanced machine learning techniques. New
Phytologist, 237(4), 1432–1445. <https://doi.org/10.1111/nph.1853>.

Cattaneo, M. (2021). selcorr: Post-Selection Inference for Generalized Linear Models. R package version 1.0. 
https://CRAN.R-project.org/package=selcorr

Crawley, M. J. (2007). The R Book (p. 942). John Wiley & Sons, Ltd. <https://doi.org/10.1251002/9780470515075126>

Crawley, M. J. (2015). Statistics: an introduction using R (sec. ed., p. 339). John Wiley & Sons. ISBN: 1118448960

Dormann, C. F., Elith, J., Bacher, S., Buchmann, C.,
Carl, G., Carré, G., García Marquéz, J. R, et al. (2013).
“<span class="nocase">Collinearity: a review of methods to deal with it
and a simulation study evaluating their performance</span>.” *Ecography*
36 (1): 27–46. <https://doi.org/10.1111/j.1600-0587.2012.07348.x>.

Forstmeier, W., & Schielzeth, H. (2011). Cryptic multiple hypotheses testing in linear models:
overestimated effect sizes and the winner’s curse. Behavioral Ecology and Sociobiology,
65(1), 47–55. <https://doi.org/10.1007/s00265-010-1038-5>.

Garnier, S., Ross, N., Rudis, B., Filipovic-Pierucci, A., et al. (2024).
*<span class="nocase">viridis(Lite)</span> - Colorblind-Friendly Color
Maps for r*. <https://doi.org/10.5281/zenodo.4679423>.

Hastie, T. (2023). gam: Generalized Additive Models. <https://cran.r-project.org/web/>.

Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B.
(2017). “<span class="nocase">lmerTest</span> Package: Tests in Linear
Mixed Effects Models.” *Journal of Statistical Software* 82 (13): 1–26.
<https://doi.org/10.18637/jss.v082.i13>.

Lee, J. D., Sun, D. L., Sun, Y., & Taylor, J. E. (2016). Exact post-selection inference, 
with application to the lasso. *Annals of Statistics* 44: 907–27. 
<https://doi.org/10.1214/15-AOS1371>.

Li, J. (2023). Overview of high dimensional linear regression models. Theoretical and Natural
Science, 5(1), 656–661. <https://doi.org/10.54254/2753-8818/5/20230427>.

Lüdecke, D., Ben-Shachar, M. S., Patil, I.,
Waggoner, P., & Makowski, D. (2021). “Performance: An r Package for
Assessment, Comparison and Testing of Statistical Models.” *Journal of
Open Source Software* 6 (60): 3139.
<https://doi.org/10.21105/joss.03139>.

R Core Team. (2024). R: a language and environment for statistical computing. 
R Foundation for Statistical Computing. <http://www.r-project.org/>

Römermann, C., Bucher, S. F., Hahn, M., & Bernhardt-Römermann, M. (2016). Plant
functional traits – fixed facts or variable depending on the season? Folia Geobotanica,
51(2), 143–159. <https://doi.org/10.1007/s12224-016-9250-3>.

Schielzeth, H. (2010). Simple means to improve the interpretability of regression coefficients.
Methods in Ecology and Evolution, 1(2), 103–113. <https://doi.org/10.1111/j.2041-210X.2010.00012.x>

Schielzeth, H., Dingemanse, N. J., Nakagawa, S., Westneat, D. F., Allegue, H., Teplitsky, C.,
Réale, D., Dochtermann, N. A., Garamszegi, L. Z., & Araya‐Ajoy, Y. G. (2020). Robustness
of linear mixed‐effects models to violations of distributional assumptions. Methods in
Ecology and Evolution, 11(9), 1141–1152. <https://doi.org/10.1111/2041-210X.13434>.

Wei, T., & Viliam, S. (2021). *R Package ’Corrplot’:
Visualization of a Correlation Matrix*.
<https://github.com/taiyun/corrplot>.

Wicke, S., Müller, K. F., DePamphilis, C. W., Quandt, D., Bellot, S., & Schneeweiss, G. M. (2016).
Mechanistic model of evolutionary rate variation en route to a nonphotosynthetic lifestyle in plants.
Proceedings of the National Academy of Sciences of the United States of America, 113(32), 9045-9050.
<https://doi.org/10.1073/pnas.1607576113>.

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019.
“Welcome to the <span class="nocase">tidyverse</span>.” *Journal of Open
Source Software* 4 (43): 1686. <https://doi.org/10.21105/joss.01686>.

Zeileis, A., Fisher, J. C., Hornik, K., Ihaka,R.,
McWhite, C. D., Murrell, P., Stauffer, R., & Wilke, C. O. 2020.
“<span class="nocase">colorspace</span>: A Toolbox for Manipulating and
Assessing Colors and Palettes.” *Journal of Statistical Software* 96
(1): 1–49. <https://doi.org/10.18637/jss.v096.i01>.

