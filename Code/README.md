# multidim-wp1 Code

## Recommended code workflow and script overview

1.  [*PrepareData.R*](PrepareData.R) - from Living Planet Database data (available at: <https://www.livingplanetindex.org/data_portal>), clean, cross reference and filter time series prior to model fitting.

2.  [*GlobalModel.R*](GlobalModel.R) - using the cleaned data, fit generalised linear mixed effects model between log abundance and threat combinations across all appropriate time series.

3.  [*RealmModels.R*](RealmModels.R) - fit generalised linear mixed effects model between log abundance and threat combinations across time series subset by realm (freshwater, marine or terrestrial).

4.  [*TaxonModels.R*](TaxonModels.R) - fit generalised linear mixed effects model between log abundance and threat combinations across time series subset by taxa (amphibian, bird, fish, mammal or reptile).

5.  [*SupplementaryModels.R*](SupplementaryModels.R) - fit generalised linear mixed effects model between log abundance and threat combinations across time series subset by time series length (10 or 20 time points).

6.  [*ModelDiagnostics.R*](ModelDiagnostics.R) - assess model fits using the normality of residuals, posterior predictive checks and residual autocorrelation.

7.  [*Figure1.R*](Figure1.R) - generate Figure 1: a map of time series locations and data coverage.

8.  [*Figure2.R*](Figure2.R) - generate Figure 2: coefficient estimates from the [*global model*](GlobalModel.R).

9.  [*Figure3.R*](Figure3.R) - generate Figure 3: estimate proportion of interactive vs additive threat influences. A diagrammatic example of the process is available in the [*supplementary information*](https://github.com/duncanobrien/multiple-threats/tree/main/Supplementary/supplementary_materials.pdf).

10. [*Figure4.R*](Figure4.R) - generate Figure 4: boxplot of counterfactual effect of removing each threat.

11. [*Counterfactual.R*](Counterfactual.R) - duplication of [*Figure4.R*](Figure4.R).

12. [*ReportedStatistics.R*](ReportedStatistics.R) - script to calculate the statistics reported in the main text.

13. [*SupplementaryFigures.R*](SupplementaryFigures.R) - supplementary figures.

14. [*SupplementaryTables.R*](SupplementaryTables.R) - supplementary tables.

## Supporting functions

The [utilities folder](https://github.com/duncanobrien/multiple-threats/tree/master/Code/utils) contains supporting functions used in the primary workflow scripts.

The [archive folder](https://github.com/duncanobrien/multiple-threats/tree/master/Code/archive) contains previous script versions in case historic code is useful.
