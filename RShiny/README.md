# Genomic surveillance sampling-frame Shiny app

This Shiny application implements the Bayesian sampling-frame analysis described by the supplied manuscript code and the additional prior-equating rules supplied with the app.

## Main changes

- Isolate-level default feature column is `mlst_profile`.
- Isolate-level observed features all receive the same `alpha_prior`.
- The CRP estimates the next-feature novelty probability; `set_alpha_novel_sum_from_crp()` converts this to novel prior mass. Following the supplied manuscript commands, the CRP-derived mass is rounded upward with `ceiling()`, then `alpha_novel_num = alpha_novel_sum / alpha_prior`, so each novel feature has exactly the same prior as an observed feature. With `alpha_prior = 1`, the number of novel features equals the rounded novel prior mass.
- Sub-isolate-level input has no feature-column selector. `sample` is used as the isolate ID when present; otherwise the first non-numeric column is used as the ID and every remaining CSV header is treated as a gene/feature.
- Sub-isolate observed features all receive the same `beta_prior`.
- The app derives `u_hat` from singleton feature mass (`number of features present in exactly one isolate / total feature occurrences`) and then applies the supplied equation `beta_novel_num = ceiling((u_hat * K) / (1 - u_hat))`; `beta_novel_sum = beta_novel_num * beta_prior`. This keeps each novel feature's beta prior equal to the observed-feature prior.
- The sub-isolate posterior uses the supplied `subisolate_bayes_beta()` formulation.
- The exact unrounded 99% detection grid is used for both tabs: `n_grid <- c(seq(0, 10000, by=1), seq(10005, 50000, by=5), seq(50010, 100000, by=10))` and `f_grid_99 <- 1 - ((1 - 0.99)^(1/n_grid))`.
- Sensitivity-grid functions have been removed from the app.

## Important reproducibility note

The supplied materials define the CRP-derived novel mass calculation explicitly. They do not define `u_hat` in the sub-isolate snippet, although they supply the downstream equation using `u_hat`. The app therefore makes the singleton-mass definition explicit so the calculation is inspectable rather than hidden. If your manuscript's earlier code used a different `u_hat`, that helper should be replaced with that exact definition before treating the app as numerically identical to the publication.

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

Use a GitHub repository as the source of truth and deploy the app to Posit Connect Cloud. The app itself does not contain or depend on a personal Shiny Server URL.
