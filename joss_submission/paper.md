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
- name: "Robert A. Jane"
  orcid: "0000-0003-4194-852X"
  corresponding: true
  equal-contrib: true
  affiliation: 1
- name: Thomas Wahl
  orcid: "0000-0003-3643-5463"
  equal-contrib: true
  affiliation: 1
- name: Luis Cadavid
  equal-contrib: true
  affiliation: 2
- name: Jayantha Obeysekera
  orcid: "0000-0002-6203-1190"
  equal-contrib: true
  affiliation: 3
- name: Javed Ali
  orcid: "0000-0002-6203-1190"
  equal-contrib: true
  affiliation: 1
- name: Author with no affiliation
  affiliation: 2
- given-names: Ludwig
  dropping-particle: van
  surname: Beethoven
  affiliation: 3
bibliography: paper.bib
aas-doi: "10.3847/xxxxx <- update this with the DOI from AAS once you know it."
aas-journal: "Journal of Open Source Software"
affiliations:
- name: University of Central Florida, USA
  index: 1
  ror: 00hx57361
- name: South Florida Water Management District, USA
  index: 2
- name: Florida International University, USA
  index: 3
---

# Summary

Compound events occur when combinations of drivers and/or hazards contribute to a societal/environmental impact [@Zscheischler:2020]. Even if none of the individual drivers or hazards are extreme, their combination can produce extreme effects. Robust estimates of the potential for compound events are critical for effective risk management and mitigation planning. To determine the probability of compound events, statistical models are applied to time series data of the drivers or hazards, typically as the first step in the risk-analysis modeling chain. 

The `MultiHazard` R package is designed to enable practitioners to estimate the likelihood of compound events. Although the methods in the package are well-established in scientific literature, they are not widely adopted by the engineering community despite guidelines increasingly mandating the estimation of compound event likelihoods. Functions within the package are designed to allow practitioners to apply their best judgement in making subjective choices. Inputs are time series representing the drivers/hazards, they may be historic or numerically generated. Outputs are: 

-	Estimates of the joint return periods for specific combinations of drivers/hazards.
-	Isolines containing drivers/hazards with a specified joint return period along with the “most-likely” or an ensemble of events sampled along an isoline.
-	Synthetic sets of events events where the peak magnitude of at least one driver/hazard is extreme.


# Statement of need

Compound events are increasingly concern for entities responsible for managing flood risk. The MultiHazard package was developed for the South Florida Water Management District to improve level of service assessments for the coastal infrastructure affected by both inland and coastal drivers [@Jane:2020]. Initially, the package was created to implement the conditional sampling-copula theory approach. In two-sided conditional sampling a driver is conditioned to be extreme and paired with the maximum value of the other driver within a specified lag-time. Process is represented with the drivers reversed yielding two conditional samples. The best fitting of 40 copulas are tested to model the dependence between drivers in each conditional sample. The isoline corresponding to a user specified return period is given by the outer envelope when overlaying the (conditional) contours from the copula model fit each sample @Bender:2016, see Figure 1. To obtain a single design event, the probability of events on the isolines is calculated using an empirical density estimate, selecting the event with the highest probability. Some experts recommend sampling an ensemble of events from the isoline to account for uncertainty in design event selection which is also possible in the package. The conditional sampling – copula theory approach has continued to gain traction [@Kim:2023; @Maduwantha:2024] and in a review of the best-available, actionable science as gaining increasing acceptance in the scientific community [@FFRMS:2023] was highlighted as an approach that Federal agencies in the United States may wish to develop detailed technical guidance on how to use it meet their needs. 

`MultiHazard` is a designed as a comprehensive user-friendly tool for copula-based joint probability analysis in R. As such the package provides functions for pre-processing data including imputing missing values, detrending and declustering time series as well as exploratory data analysis e.g., analyzing pairwise correlations over a range of lags. The package also contains an automated threshold selection method for the Generalized Pareto distribution [@Solari:2017] and approaches for robustly capturing the dependence structure when there are more than two relevant drivers/hazards,namely, standard (elliptic/Archimedean) copulas, Pair Copula Constructions (PCCs) and the conditional threshold exceedance approach of @Heffernan:2004. For the analysis undertaken for the SFWMD, the higher dimensional approaches enabled groundwater level to be included in the analysis.  More recently, an approach to generate time varying synthetic events i.e., hyetographs and hydrographs was added. Time varying condition are a prerequisite for non-steady state hydrodynamic modeling. 

# Related packages

Several packages employ copula-based approaches to derive isolines. The MATLAB toolbox `MvCAT` (Multivariate Copula Analysis Toolbox) [@Sadegh:2017] utilizing up to 26 copula families to model the dependence structure between a pair of random variables. `MvCAT` has fewer copula families than MultiHazard and does not explicitly implement the two-sided conditional sampling -copula theory methodology. On the other hand, `MvCAT` adopts multiple criteria to select among the candidate copula families and a Bayesian framework to account for the uncertainty range for the copula parameters. `ReturnCurves` is a recently released R package that generates bivariate isolines based on the angular dependence function. The bootstrap procedure in the `ReturnCurve`package used to assess the sampling uncertainty of the estimated isolines is implemented in `MultiHazard`. 

# Acknowledgements

The Authors thank the South Florida Water Management District for funding the development of this package through several projects over the past few years. Robert Jane and Thomas Wahl acknowledge financial support from the USACE Climate Preparedness and Resilience Community of Practice

# References
