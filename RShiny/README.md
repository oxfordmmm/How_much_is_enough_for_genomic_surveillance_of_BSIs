# Genomic surveillance sampling-frame Shiny app

This Shiny application implements the Bayesian analysis-based sample size estiamtion for genomic surveillance of bacterail pathogens. 

## Usage

The `Isolate-level features` tab lets you analyse features that are mutually exclusive in isolates (e.g: MLST), while the second `Sub-isolate level features` tab is for features that are not mutually exclusive (e.g.: carriage of antimicrobial resistance [AMR] genes).

You can either enter you won data in the provided fields, use one of the example datasets, or upload your own csv data matching the available template. 


## R packages

The app requires:

- shiny
- dplyr
- tidyr
- ggplot2
- plotly
- readr
- htmltools
- matrixStats

Install once with:

```r
install.packages(c("shiny", "dplyr", "tidyr", "ggplot2", "plotly", "readr", "htmltools", "matrixStats"))
```

Run locally with:

```r
shiny::runApp(".")
```

## Hosting

This app is deployed using Posit Connect Cloud. The app itself does not contain or depend on a personal Shiny Server URL.
