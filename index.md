# Ensemble Random Forests
## R package for Ensemble Random Forests

This R package implements Ensemble Random Forests from [Siders et al. (2020)](https://www.int-res.com/abstracts/esr/v43/p183-197/). Ensemble Random Forests is exactly as it sounds, an ensemble model built from multiple Random Forests. The theory is that Random Forests perform exceptionally well as machine learning tools for categorical data when the data is *balanced*, i.e. that the proportion of data in each category is close to equal. There many instances where a *balanced* dataset occurs and there have been various methods to help the Random Forests model acheive good performance. The motivating example for developing the Ensemble Random Forests model comes from protected species interactions in fisheries, which by design, are often rare or hyper-rare. In this case, many of the costs of the balancing adjustments to Random Forests to tend outweigh the benefits. Enter Ensemble Random Forests. By using an ensemble of Random Forests, the effects of poorly learning the majority class can be ameliorated to some extent and can facilitate running a Random Forests approach in hyper-rare instances. 

## Support
The initial work for this package was funded by the [Western Pacific Regional Fishery Management Council](https://www.wpcouncil.org), the [National Oceanic and Atmospheric Adminstration Pacific Islands Fisheries Science Center](https://www.fisheries.noaa.gov/about/pacific-islands-fisheries-science-center), and the [National Oceanic and Atmospheric Adminstration Pacific Islands Regional Office](https://www.fisheries.noaa.gov/about/pacific-islands-regional-office)under Award No: NA15NMF4410008 and NA15NMF4410066. Current development and support is under Award No: NA20NMF4410015. 

Additional simulation code to replicate [Siders et al. (2020)](https://www.int-res.com/abstracts/esr/v43/p183-197/) can be found on [Open Science Framework](https://osf.io/q9wfn/).


