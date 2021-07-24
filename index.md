# Ensemble Random Forests
## R package for Ensemble Random Forests

This R package implements Ensemble Random Forests from [Siders et al. (2020)](https://www.int-res.com/abstracts/esr/v43/p183-197/). Ensemble Random Forests is exactly as it sounds, an ensemble model built from multiple Random Forests. The theory is that Random Forests perform exceptionally well as machine learning tools for categorical data when the data is *balanced*, i.e. that the proportion of data in each category is close to equal. There many instances where a *balanced* dataset occurs and there have been various methods to help the Random Forests model acheive good performance. The motivating example for developing the Ensemble Random Forests model comes from protected species interactions in fisheries, which by design, are often rare or hyper-rare. In this case, many of the costs of the balancing adjustments to Random Forests to tend to outweigh the benefits. Enter Ensemble Random Forests. By using an ensemble of Random Forests, the effects of poorly learning the majority class can be ameliorated to some extent and can facilitate running a Random Forests approach in hyper-rare instances. 

### Tutorial

Running ERFs on a given dataset is easy. The function `erf()` will take a given dataset in R `data.frame` format, amend it for modeling using `erf_data_prep()` and `erf_formula_prep()`, run each RF in the ensemble using `rf_ens_fn()`, and return a fitted ERF object. This object can then be passed to various output functions: `erf_plotter()` and ... to visualize and summarize.

##### To generate a new dataset:
Just run `data_sim()` and a new dataset with the default parameters will be generated. The examples in `data_sim()` provide a few code snippets to quickly visualize the simulated dataset. 

##### Or use the provided simulated dataset:
The provided dataset is a `list` object that contains a `data.frame` of the sampled locations, the beta coefficients of the logistic model used to predict the probability of occurrence, and a `raster` `brick` object containing the gridded covariates, log-odds of occurrence, and probabilities of occurrence. 

###### The simulated covariates

![](reference/figures/simulated_covariates.png)

###### The log-odds and the probability of occurrence with the observed presences (white dots)

![](reference/figures/simulated_real.png)

##### Check model performance:
We can check model performance using the `rocr_ens` function. This calculates a battery of performance metrics. This function works on any set of predictions (ranging from (0,1)) and any set of observations (as a `factor`). We can test this on our simulated data. 

![](reference/figures/simulated_roc.png)

