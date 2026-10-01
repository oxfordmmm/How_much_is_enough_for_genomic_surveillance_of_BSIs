# EpiSENTRY genomic surveillance sampling-frame Shiny app

EpiSENTRY is an R Shiny application for estimating the sample size required to capture a specified proportion of population genomic diversity in pathogen genomic surveillance. It uses Bayesian analysis of user-provided survey or sample data and supports both mutually exclusive isolate-level features and non-mutually-exclusive sub-isolate-level features.

The application accompanies:

[**Nagy et al., 2026. _How much is enough? Optimising sampling frames for genomic surveillance of Escherichia coli and Klebsiella spp. bloodstream infections – a retrospective study._**](https://www.medrxiv.org/content/10.64898/2026.09.28.26364154v1)

## Why use this tool?

A conventional sample-size or power calculation can be used to estimate the number of samples needed to detect a feature at a specified minimum frequency. However, this requires the user to choose a minimum frequency threshold and generally considers one feature at a time.

EpiSENTRY is designed for genomic surveillance data in which many genomic features may be relevant simultaneously. The method:

- considers multiple genomic features at once;
- accounts for feature categories that may be present in the underlying population but were not observed in the initial sample;
- incorporates uncertainty in estimated feature frequencies or prevalences; and
- relates sample size to **coverage**, defined here as the proportion of the population represented by features detected at least once in the sample.

For comparison, a conventional sample-size calculator is available from [ClinCalc](https://clincalc.com/stats/samplesize.aspx).

## Usage

### 1. Select the appropriate feature tab

Use **Isolate-level features** when each isolate or sampling unit can have only one value for the feature of interest. Examples include MLST, cgMLST, fastBAPS clusters, or another mutually exclusive lineage/category label.

Use **Sub-isolate-level features** when each isolate or sampling unit can contain zero, one, or multiple features. Examples include antimicrobial resistance genes, virulence genes, plasmids, or other presence/absence genomic features.

### 2. Provide input data from an initial genomic survey

Use data from the population of interest or, where appropriate, from a sufficiently similar population.

For either tab you can:

1. select **Upload CSV** and upload a `.csv` file;
2. download the relevant template CSV from the app, complete it, and upload it; or
3. select **Paste CSV text** and paste comma-separated data directly into the text box.

If both an uploaded file and pasted text are present, the uploaded file is used.

Example data for _E. coli_ and _Klebsiella_ are available in the associated GitHub repository.

### 3. Select the column names and analysis parameters

Each input field in the app has a short description. The default target coverage and target feature richness are both 0.80.

### 4. Press **Run analysis**

The analysis can take a few minutes, particularly for large datasets or analyses using many posterior draws. A progress indicator is displayed while the analysis is running.

### 5. Download results if required

The numerical results can be downloaded as a CSV and the plots can be downloaded as a PDF.

## Input data requirements

### Isolate-level features

The isolate-level input represents a frequency table in which each row contains a genomic feature category and the number of isolates in that category.

By default, the app expects these column names:

- `mlst_profile` — feature/category label;
- `count` — number of isolates observed with that feature/category.

The names can be changed in the **Feature column** and **Count column** fields in the app.

Example:

```csv
mlst_profile,count
ST131,120
ST73,55
ST95,31
```

Requirements:

- the feature column may contain character or numeric labels;
- blank or missing feature labels are removed;
- the count column must contain numeric, finite, non-negative values;
- repeated feature labels are permitted and are automatically combined by summing their counts;
- additional columns may be present but are not used by the isolate-level analysis;
- the total number of observed isolates is calculated as the sum of the count column.

### Sub-isolate-level features

The sub-isolate input is a sample-by-feature matrix. Each row represents one isolate or sampling unit and each genomic feature is represented by a separate column.

Example:

```csv
sample,blaCTX-M-15,tetA,sul1
ISO001,1,1,0
ISO002,0,1,0
ISO003,1,0,1
```

#### Isolate ID column

The app identifies the isolate/sample ID column automatically:

1. if a column called `sample` is present, it is used as the ID column;
2. otherwise, the first non-numeric column is used;
3. if all columns are numeric, the first column is treated as the ID column.

The ID column itself is not included as a genomic feature.

#### Feature columns

Every remaining column is treated as a genomic feature and its column header becomes the feature name.

Feature values must be numeric and finite. The app converts values to presence/absence as follows:

- values `> 0` are treated as **present** (`1`);
- values `<= 0` are treated as **absent** (`0`).

For standard presence/absence data, using explicit `0` and `1` values is recommended.

Do not include metadata columns such as species, hospital, date, location, or patient group among the feature columns unless they are intentionally being analysed as numeric presence/absence features. If such metadata are needed in the source dataset, prepare a separate input CSV containing only the isolate ID and feature columns before uploading it to the app.

## Outputs and interpretation

The app reports posterior feature-frequency/prevalence estimates, estimated unseen feature categories, and sample-size estimates for the selected coverage and feature-richness targets.

### Coverage

**Coverage** is the proportion of the population represented by features detected at least once in the sample. The coverage plots show:

- how estimated population coverage changes with sample size;
- the posterior population mass contributed by features above each frequency threshold; and
- uncertainty in the sample size required to achieve the selected target coverage.

### Feature richness

**Feature richness** is the proportion of modelled unique feature categories expected to be detected. In the current implementation, the modelled feature total includes both features observed in the input dataset and the predicted novel feature categories.

The richness plots show:

- how the proportion of modelled feature categories expected to be detected changes with sample size;
- the proportion of modelled feature categories above each frequency threshold; and
- uncertainty in the sample size required to achieve the selected target richness.

### Observed versus posterior prevalence

The observed-versus-posterior plot compares the prevalence/frequency calculated directly from the input data with the posterior median and 95% credible interval for each observed feature. Points close to the diagonal indicate close agreement between the observed and posterior estimates.

## Statistical implementation

For isolate-level features, the app uses a Chinese restaurant process to estimate the probability of a previously unseen feature category and a Dirichlet posterior for feature frequencies. Observed and predicted novel features are assigned the same per-feature `alpha` prior.

For sub-isolate-level features, the app estimates unseen-feature probability using the Good-Turing estimator used in the application and applies a feature-wise beta posterior. Observed and predicted novel features are assigned the same per-feature `beta` prior.

The sample-size grid corresponds to 99% probability of detecting a feature and explores sample sizes from 0 to 100,000 using the grid implemented in `app.R`.

## R packages

The app requires:

- `shiny`
- `dplyr`
- `tidyr`
- `ggplot2`
- `plotly`
- `readr`
- `htmltools`
- `matrixStats`

Install the required packages with:

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

## Running locally

From the `RShiny` directory:

```r
shiny::runApp(".")
```

Or from the root of this repository:

```r
shiny::runApp("RShiny")
```

## Citation

If you use the application or its methodology, please cite:

**Nagy et al., 2026. _How much is enough? Optimising sampling frames for genomic surveillance of Escherichia coli and Klebsiella spp. bloodstream infections – a retrospective study._**
