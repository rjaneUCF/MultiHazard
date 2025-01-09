---
title: '`MultiHazard`: Copula-based Joint Probability Analysis in R'
tags:
- R
- copula
- isolines
- flood risk analysis
- multivariate statistics
date: "1 December 2024"
output:
  pdf_document: default
  html_document:
    df_print: paged
authors:
- name: Robert A. Jane
  orcid: "0000-0003-4194-852X"
  corresponding: true
  equal-contrib: true
  affiliation: 1
- name: Thomas Wahl
  orcid: "0000-0003-3643-5463"
  equal-contrib: true
  affiliation: 1
- name: Francisco Pena
  orcid: "0000-0002-1779-9926"
  equal-contrib: true
  affiliation: 2
- name: Jayantha Obeysekera
  orcid: "0000-0002-6203-1190"
  equal-contrib: true
  affiliation: 3
- given-names: Callum Murphy-Barltrop
  orchid: "0000-0002-3479-2902"
  affiliation: 4
- name: Javed Ali
  orcid: "0000-0002-6203-1190"
  equal-contrib: true
  affiliation: 1
- name: Pravin Maduwantha
  orcid: "0000-0002-9819-5170"
  equal-contrib: true
  affiliation: 1
- name: Huazhi Li
  orchid: "0000-0001-9589-2918"
  affiliation: 5
affiliations:
- name: University of Central Florida, USA
  index: 1
  ror: 00hx57361
- name: South Florida Water Management District, USA
  index: 2
- name: Florida International University, USA
  index: 3
- name: TU Dresden, DE
  index: 4
- name: Vrije Universiteit Amsterdam, NL
  index: 5
bibliography: paper.bib
aas-doi: "10.3847/xxxxx <- update this with the DOI from AAS once you know it."
aas-journal: "Journal of Open Source Software"
---

# Summary

Compound events occur when combinations of drivers and/or hazards contribute to a societal/environmental impact [@Zscheischler:2020]. Even if none of the individual drivers or hazards are extreme, their combination can produce extreme impacts. Assessing the potential for compound extreme events is therefore critical for effective risk management and mitigation planning. To determine the probability of compound events, statistical models are applied to time series data of the drivers or hazards, typically as the first step in the risk-analysis modeling chain. 

The `MultiHazard` R package is designed to enable practitioners to estimate the likelihood of compound events. Although the methods in the package are well-established in the scientific literature, they are not widely adopted by the engineering community despite guidelines increasingly mandating the estimation of compound event likelihoods. Functions within the package are designed to allow practitioners to apply their best judgement in making subjective choices. Inputs are time series representing the drivers/hazards; these may be historical observations or numerically generated with models (for the past or future). Outputs are: 

-	Estimates of the joint return periods for specific combinations of drivers/hazards.
-	Isolines with uncertainty bounds containing drivers/hazards with a specified joint return period along with the “most-likely” or an ensemble of events sampled along an isoline.
-	Synthetic sets of events where the peak magnitude of at least one driver/hazard is extreme.


# Statement of need

The `MultiHazard` package was developed in collaberation with the South Florida Water Management District (SFWMD) to improve level of service assessments for coastal infrastructure affected by both inland and coastal drivers [@Jane:2020]. Initially, the package was created to implement the conditional sampling - copula theory approach. In two-sided conditional sampling a driver is conditioned to be extreme and paired with the maximum value of the other driver within a specified lag-time. The process is repeated with the drivers reversed yielding two conditional samples. The best fitting of 40 copulas are tested to model the dependence between drivers in each conditional sample. The isoline corresponding to a user specified return period is given by the outer envelope when overlaying the (conditional) contours from the copula model fit to each sample @Bender:2016, see \autoref{fig:isolines}a,b. To obtain a single design event, the probability of events on the isolines is calculated using an empirical density estimate, selecting the event with the highest likelihood. Some experts recommend sampling an ensemble of events from the isoline to account for uncertainty in design event selection which is also possible in the package. The conditional sampling – copula theory approach has continued to gain traction [@Kim:2023; @Maduwantha:2024] and in a review of the best-available, actionable science was highlighted as an approach that Federal agencies in the United States may wish to develop detailed technical guidance on how to use it meet their needs [@FFRMS:2023]. 

The package possesses the functionaity to compute the desired isoline by interpolating joint return periods calculated over a user-defined grid which can result in smoother curves than overlaying the partial isolines. Isoline generation for time series comprising events from multiple (independent) populations e.g. different types of storms is faciliated by multiplying independent annual-non exceedance probabilities associated with each poopulation as outlined in @Maduwantha:2024. Also coded in the package are the methods proposed in @Barltrop:2023 to derive isolines using the HT04 [@Heffernan:2004] and WT13 [@Wadsworth:2013] models, along with the novel bootstrap procedure for calculating sample uncertainties, see \autoref{fig:isolines}c,d. The HT04 and WT13 models allow for the analysis of compound events without assuming any copula form. 

![The 100-year isoline at case study site derived using two methods. Two-sided consitoina sample - copula theory method: (a) Partial isolines from the samples conditioned on rainfall (red line) and water level (blue line) whose outer envelopen when overlaid gives (b) the full isoline. Heffernan and Tawn model method: (c) Sample uncertainties are calculated along 100 rays eminating from the point comprising the mimum obseved values of both drivers. (d) The Isoline (red line) and associated 90% confidence interval (dashed black lines) obtained via a block bootstrap procedure . \label{fig:isolines}](Figure_1.png)

Compound events are of increasing concern for entities responsible for managing flood risk. `MultiHazard` is designed as a comprehensive user-friendly tool for copula-based joint probability analysis in R. As such, the package provides functions for pre-processing data including imputing missing values, detrending and declustering time series, as well as exploratory data analysis, e.g., analyzing pairwise correlations over a range of time-lags between the two drivers/hazards. The package also contains an automated threshold selection method for the Generalized Pareto distribution [@Solari:2017] and approaches for robustly capturing the dependence structure when there are more than two relevant drivers/hazards,namely, standard (elliptic/Archimedean) copulas, Pair Copula Constructions (PCCs) and the HT04 model. For the analysis undertaken for the SFWMD, the higher dimensional approaches enabled groundwater level to be included in the analysis. More recently, an approach to generate time varying synthetic events, i.e., hyetographs and hydrographs, was added. Time varying conditions are a prerequisite for non-steady state hydrodynamic modeling. 

# Related packages

Several packages employ copula-based approaches to derive isolines. The MATLAB toolbox `MvCAT` (Multivariate Copula Analysis Toolbox) [@Sadegh:2017] utilizing up to 26 copula families to model the dependence structure between a pair of random variables. `MvCAT` includes fewer copula families than MultiHazard and does not consider the two-sided conditional sampling - copula theory methodology but rather uses a one-sided sampling approach. On the other hand, `MvCAT` adopts multiple criteria to select among the candidate copula families and a Bayesian framework to account for the uncertainty range for the copula parameters. `ReturnCurves` [@Andre:2024] is a recently released R package that generates bivariate isolines based on the angular dependence function (@Wadsworth:2013 model) [@Barltrop:2023]. The bootstrap procedure in the `ReturnCurves`package used to assess the sampling uncertainty of the estimated isolines is implemented in `MultiHazard`. 

# Acknowledgements

The Authors thank the South Florida Water Management District for funding the development of this package through several projects over the past few years. Robert Jane and Thomas Wahl acknowledge financial support from the USACE Climate Preparedness and Resilience Community of Practice.

# References
