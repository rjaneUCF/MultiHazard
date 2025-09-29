
# MultiHazard <img src="https://user-images.githubusercontent.com/15319503/195926656-9d3d37b5-86ab-4d4b-9e6d-3c70d5399c73.png" align="right" height="200"/>

## Overview

The `MultiHazard` package provides tools for stationary multivariate
statistical modeling, for example, to estimate joint occurrence
probabilities of MULTIple co-occurring HAZARDs. 

<br>

## Installation

Install the latest version of this package by entering the following in
R:

``` r
install.packages("remotes")
remotes::install_github("rjaneUCF/MultiHazard")
```

<br>

## Features

The package contains
functions for pre-processing data, including imputing missing values,
detrending and declustering time series as well as analyzing pairwise
correlations over a range of lags. Functionality is built in to
implement the conditional sampling - copula theory approach described in
[Jane et al. (2020)](https://doi.org/10.5194/nhess-20-2681-2020)
including the automated threshold selection method from [Solari et
al. (2017)](https://doi.org/10.1002/2016WR019426). There is a function
that calculates joint probability contours using the method of
overlaying conditional contours given in [Bender et
al. (2016)](https://doi.org/10.1080/02626667.2015.1052816) and extracts
design events such as the “most likely” event or an ensemble of possible
design events. The package also includes methods from [Murphy-Barltrop
et al. (2023)](https://doi.org/10.1002/env.2797) and [Murphy-Barltrop et
al. (2024)](https://doi.org/10.1007/s10687-024-00490-4) for deriving
isolines using the [Heffernan and Tawn
(2004)](https://doi.org/10.1111/j.1467-9868.2004.02050.x) \[HT04\]
and [Wadsworth and Tawn (2013)](https://doi.org/10.3150/12-BEJ471)
\[WT13\] models, together with a novel bootstrap procedure for
quantifying sampling uncertainty in the isolines. Three higher
dimensional approaches — standard (elliptic/Archimedean) copulas, Pair
Copula Constructions (PCCs) and a conditional threshold exceedance
approach (HT04) — are coded. Finally, the package can be implemented to
derive temporally coherent extreme events comprising a hyetograph and
water level curve for simulated peak rainfall and peak sea level events,
as outlined in (Report).

## Learn More

For detailed tutorials and examples, see the package [vignette](https://github.com/rjaneUCF/MultiHazard/tree/master/vignettes)

## Citation

If you use this package, please cite:

> Jane, R., Wahl, T., Peña, F., Obeysekera, J., Murphy-Barltrop, C.,
> Ali, J., Maduwantha, P., Li, H., and Malagón Santos, V. (under review)
> MultiHazard: Copula-based Joint Probability Analysis in R. Journal of
> Open Source Software. \[under revision\]

<br>

## Community guidelines

Contributions to the `MultiHazard` package are welcome! Please feel free
to submit issues or pull requests on GitHub.

<br>

## Applications of package

This package has been used in several peer-reviewed publications:

> Li, H., Jane, R. A., Eilander, D., Enríquez, A. R., Haer, T., and
> Ward, P. J. (2025). Assessing the spatial correlation of potential
> compound flooding in the United States, EGUsphere \[preprint\],
> <https://doi.org/10.5194/egusphere-2025-2993>.

> Amorim, R., Villarini, G., Kim, H., Jane, R., and Wahl, T. (2025). A
> Practitioner’s approach to process-driven modeling of compound
> rainfall and storm surge extremes for coastal Texas, J. Hydrol. Eng.,
> 30(5), 04025025, <https://doi.org/10.1061/JHYEFF.HEENG-648>.

> Maduwantha, P., Wahl, T., Santamaria-Aguilar, S., Jane, R., Booth, J.
> F., Kim, H., and Villarini, G. (2024). A multivariate statistical
> framework for mixed storm types in compound flood analysis, Nat.
> Hazards Earth Syst. Sci., 24, 4091–4107,
> <https://doi.org/10.5194/nhess-24-4091-2024>.

> Nasr, A. A., Wahl, T., Rashid, M. M., Jane, R., Camus, P. and Haigh,
> I. D. (2023). Temporal changes in dependence between compound coastal
> and inland flooding drivers around the contiguous United States
> coastline, Weather Clim. Extrem., 41, 100594,
> <https://doi.org/10.1016/j.wace.2023.100594>.

> Kim, H., Villarini, G., Jane, R., Wahl, T., Misra, S., and Michalek,
> A. (2023). On the generation of high‐resolution probabilistic design
> events capturing the joint occurrence of rainfall and storm surge in
> coastal basins, Int. J. Climatol, 43(2), 761-771,
> <https://doi.org/10.1002/joc.7825>.

> Kim, T., Villarini, G., Kim, H., Jane, R., and Wahl, T. (2023). On the
> compounding of nitrate loads and discharge, J. Environ. Qual., 52,
> 706–717. <https://doi.org/10.1002/jeq2.20458>.

> Peña, F., Obeysekera, J., Jane R., Nardi, F., Maran, C., Cadogan, A.,
> de Groen, F., and Melesse, A. (2023). Investigating compound flooding
> in a low elevation coastal karst environment using multivariate
> statistical and 2D hydrodynamic modeling, Weather Clim. Extrem., 39,
> 100534. <https://doi.org/10.1016/j.wace.2022.100534>.

> Jane, R., Cadavid, L., Obeysekera, J., and Wahl, T. (2020).
> Multivariate statistical modelling of the drivers of compound flood
> events in South Florida, Nat. Hazards Earth Syst. Sci., 20, 2681–2699,
> <https://doi.org/10.5194/nhess-20-2681-2020>.



