# Genomic surveillance sampling-frame Shiny app

This Shiny application implements a Bayesian approach to sample-size estimation for genomic surveillance of bacterial pathogens. It is designed to help estimate the sampling effort required to capture a specified proportion of genomic diversity at a given level of certainty, given input sampling data.

The application accompanies the methods described in:

**Nagy et al., 2026. *How much is enough? Optimising sampling frames for genomic surveillance of Escherichia coli and Klebsiella spp. bloodstream infections – a retrospective study.***

## Usage
Please go to [this link](https://connect.posit.cloud/dorottyanagy96/content/01a0c3dd-90d2-e597-58a5-51391a468387) to use the deployed app in Posit Cloud Connect.

The application provides two analysis frameworks.

### Isolate-level features

The **Isolate-level features** tab is intended for mutually exclusive characteristics, where each isolate belongs to a single category, such as multilocus sequence type (MLST), or fastBAPS cluster.

Input data should contain:

* `mlst_profile` – the feature or category name (can be modified in 'feature column' on the first tab)
* `count` – the number of isolates belonging to that category

A template CSV is available from within the application.

### Sub-isolate-level features

The **Sub-isolate level features** tab is intended for characteristics that are not mutually exclusive, where an isolate may contain zero, one, or multiple features. Examples include antimicrobial resistance (AMR) genes and plasmids.

Input data should be supplied as a presence/absence matrix, with one isolate per row and one feature per column. Feature names are automatically obtained from the CSV column headers.

A template CSV is available from within the application.

## Getting started

You can use the application by:

* entering data directly into the provided fields;
* using one of the example datasets; or
* uploading your own CSV file following the appropriate input template.

The application estimates the posterior distribution of feature frequencies, including allowance for features that have not yet been observed, and uses these estimates to explore the sampling effort required to capture genomic diversity.

## R packages

The application requires the following R packages:

* `shiny`
* `dplyr`
* `tidyr`
* `ggplot2`
* `plotly`
* `readr`
* `htmltools`
* `matrixStats`

To install the required packages locally:

```r
install.packages(c(
  "shiny",
  "dplyr",
  "tidyr",
  "ggplot2",
  "plotly",
  "readr",
  "htmltools",
  "matrixStats"
))
```

## Running the application locally

From the `RShiny` directory, run:

```r
shiny::runApp(".")
```

Alternatively, from the root of this repository:

```r
shiny::runApp("RShiny")
```

## Citation

If you use this application or its methodology, please cite:

**Nagy et al., 2026. *How much is enough? Optimising sampling frames for genomic surveillance of Escherichia coli and Klebsiella spp. bloodstream infections – a retrospective study.***


## Hosting

This app is deployed using Posit Connect Cloud. The app itself does not contain or depend on a personal Shiny Server URL.
