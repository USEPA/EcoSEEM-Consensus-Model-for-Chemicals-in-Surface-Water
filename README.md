# EcoSEEM-Consensus-Model-for-Chemicals-in-Surface-Water

U.S. EPA’s Center for Computational Toxicology and Exposure provides tools to rapidly generate quantitative toxicity, human exposure, and internal dose estimates. The SEEM (Systematic Empirical Evaluation of Models) meta-modeling approach uses forward-backward inference to identify when complex models add more signal than noise. SEEM has been used to successfully develop exposure predictions built from different data streams and model results that correspond to human biomonitoring data. The current project, EcoSEEM, extends this approach to surface water concentrations. We assess the correlation between measured values from surface water monitoring sites and results from openly-available mechanistic models that provide nationwide surface water concentration estimates for many chemicals. Seasonal concentrations for hundreds of chemicals from 1998 to 2018 inform estimates of likely loading values and the predictivity of each model. The three high-throughput models evaluated cover both far- and near-field chemical releases, and predict chemical fate based on various physical chemical properties. The overall association between model results and observations is used to create a consensus regression model for screening-level concentration estimates (and estimate uncertainties) of human and ecological contact with the thousands of chemicals which may be present in surface water, with or without monitoring data, which can contribute to the prioritization of safety evaluations for a broad range of chemicals. This abstract does not necessarily represent the views or policies of the U.S. Environmental Protection Agency. 
## Characterizing surface water concentrations of hundreds of organic chemicals in United States for environmental risk prioritization

### Supplemental Information

* observation_data/load_water_data.py
* observation_data/sayre_water_data.R
* observation_data/smaller_df_21july2020.zip

### Abstract

#### Background

Thousands of chemicals are observed in freshwater, typically at trace levels. Measurements are collected for different purposes, so sample characteristics vary. Due to inconsistent data availability for exposure and hazard, it is complex to prioritize which chemicals may pose risks. 

#### Objective

We evaluated the influence of data curation and statistical practices aggregating surface water measurements of organic chemicals into exposure distributions intended for prioritizing based on nation-scale potential risk.

#### Methods

The Water Quality Portal includes millions of observations describing over 1700 chemicals in 93% of hydrologic subbasins across the United States. After filtering to maintain quality and applicability while including all possible samples, we compared concentrations across sample types. We evaluated statistical methods to estimate per-chemical distributions for chosen samples. Overlaps between resulting exposure ranges and distributions representing no-effect concentrations for multiple freshwater species were used to rank estimated chemical risks for further assessment.

#### Results 

We make screening-level estimates of surface water chemical concentration for 186 organic chemicals. From the original set, the chemical space decreased primarily due to censored values. In the final set of 1.5 million samples, one chemical exceeded the 5th percentile of no-effect concentrations for the most delicate freshwater species with the median environmental concentration (the highest priority risk condition identified here), and a further 29 chemicals were identified for possible further evaluation of environmental occurrence. We further identify 169 chemicals where all measurements were censored but with the range comparison method, risk could still be prioritized.

### Authors

* Risa R. Sayre [sayre.risa@epa.gov]
* R. Woodrow Setzer [setzer.woodrow@epa.gov]
* Marc L. Serre [marc_serre@unc.edu]
* John F. Wambaugh [wambaugh.john@epa.gov] 
