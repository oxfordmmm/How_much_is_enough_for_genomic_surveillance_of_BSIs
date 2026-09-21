# How much is enough for genomic surveillance of *Escherichia* and *Klebsiella* BSIs in England?

This repository contains code to estiamte the same size of bloostream infection (BSI)-associated *E. coli* and *Klebsiella* isolates required for genomic sequencing to achieve a specified sampling coverage using Bayesian analysis combined with power calculation. These analyses use data from the NEKSUS (National *E. coli* and *Klebsiella* bloodstream infection and CPE UK Surveillance) Study. Data for this study are underembargo until April 2027. The results of this analysis are presented as a poster at ESCMID 2026. 

# Repository contents
## RShiny application
The `RShiny` folder contains the code to easily apply this Bayesian sample size estimation method to your own data. Some example datasets are included. Please follow [think link](https://connect.posit.cloud/dorottyanagy96/content/01a0c3dd-90d2-e597-58a5-51391a468387) for the deloyed version of the sample size estimator tool. The README.md within the RShiny folder contains further information on how to use the tool.

## R scripts
The `R_scripts` folder containd R analysis scripts for generating fastBAPS clusters, performing comparison of ecological estimators to Bayesian analysis, and performing Bayesian estimation of sample size for standard and hierarchical models.

## Citation

If you use this application or its methodology, please cite:

**Nagy et al., 2026. *How much is enough? Optimising sampling frames for genomic surveillance of Escherichia coli and Klebsiella spp. bloodstream infections – a retrospective study.***



For data files associated with the publication, please see: [https://10.6084/m9.figshare.32326584](https://doi.org/10.6084/m9.figshare.32326584).
