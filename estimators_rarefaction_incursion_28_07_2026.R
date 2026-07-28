#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ANALYSIS OF ECOLOGICAL ESTIMATORS, RAREFACTION, LINEAGE INCURSION AND COVERAGE ESTIMATORS ####
# for NEKSUS study overall 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages
install.packages("tidyverse")
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
install.packages("ggplot2")
install.packages("purrr")
install.packages("bayesboot")
install.packages("ggh4x")
#install.packages("grafify")
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(purrr)
library(grafify)
library(tibble)
library(stringr)
library(scales)
library(forcats)
library(patchwork)
library(gtsummary)
library(gt)
library(ggh4x)

# set workign directory
setwd("~")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load and prepare data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ecoli_bsi_samples_metadata <- read.csv("model_results/ecoli_bsi_samples_metadata.csv")
kleb_bsi_samples_metadata <- read.csv("model_results/kleb_bsi_samples_metadata.csv")

# load amrfinder and plasmid data 
## this df already had plasmid "community_subcommunity" (pling) and "rep_types_whole_plasmid" (mob-suite)
ecoli_bsi_amrfinder_metadata <- read.csv("model_results/ecoli_bsi_amrfinder_metadata.csv")
kleb_bsi_amrfinder_metadata <- read.csv("model_results/kleb_bsi_amrfinder_metadata.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. ECOLOGICAL ESTIMATORS iNEXT AND PRESEQR ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 1.a preseqR ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.ai E. coli preseqR ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ecoli BSI upsampling with preseq
#install.packages("PreseqR") # this doesn't work as removed from CRAN :(
install.packages("fs")
find.package("cli")
packageVersion("cli")
remove.packages("cli")
install.packages("cli")
install.packages("devtools")
devtools::install_github("smithlabcode/PreseqR")
library(preseqR)


# if the above not working try:
# 1. Make sure you have a clean R session
#.rs.restartR()   # RStudio; or restart R manually

# 2. Install the missing package(s)
#install.packages("polynom")

# 3. Reinstall preseqR from GitHub (ask to install dependencies too)
# Use remotes (lighter) or devtools; remotes is recommended for single installs
#if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
#remotes::install_github("smithlabcode/PreseqR", dependencies = TRUE)
#library(preseqR)


# Count profile frequencies - use MLSTs from kleborate
profile_counts <- ecoli_bsi_samples_metadata |>
  mutate(mlst_profile = as.character(escherichia__mlst_achtman__ST)) |>
  group_by(mlst_profile) |>
  summarise(n = n()) |>
  group_by(n) |>
  summarise(no_mlsts = n_distinct(mlst_profile))
#View(profile_counts)
profile_counts <- as.matrix(profile_counts)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Species richness (i.e. species accumulation curve (SAC) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Question: how many unique MLSTs occur at least r times? (as a function of sample size)

# define current and estiamted upsampling sample sizes for E.coli BSIs:
current_n <- nrow(ecoli_bsi_samples_metadata) #1471
estimated_total_n <- current_n / 0.65 #2263 = capture proporiton from SGSS
estimated_national_n <- 42224  # E. coli BSIs occuring in 23/24 FY

#unique mlsts observed currently
actual_unique_mlsts <- length(unique(ecoli_bsi_samples_metadata$escherichia__mlst_achtman__ST)) # 263

# define sample sizes
relative_sample_sizes <- c(seq(0,30, 0.1), (estimated_total_n / current_n), (estimated_national_n / current_n))
relative_sample_sizes <- sort(relative_sample_sizes)
#length(relative_sample_sizes)
actual_sample_sizes <- relative_sample_sizes * current_n
actual_sample_sizes <- sort(actual_sample_sizes)
# define parameters, including range of r = minimum number of times it is important to detect an MLST).
r_values <- c(seq(1)) 
bootstrap_times <- 1000
conf_level <- 0.95

# use lapply to apply preseqR.rSAc.bootstrap function over all values of r_values
results_list <- lapply(r_values, function(rval) {
  message(sprintf("Running preseqR for r = %s ...", rval))
  sac <- preseqR.rSAC.bootstrap(n = profile_counts, r = rval, mt = 20,
                                size=SIZE.INIT, mu=MU.INIT, times = bootstrap_times, conf = conf_level)
  
  tibble(r = rval, 
         relative_sample_size = relative_sample_sizes,
         actual_sample_size = actual_sample_sizes,
         est = sac$f(relative_sample_sizes),
         se = sac$se(relative_sample_sizes),
         lower = sac$lb(relative_sample_sizes),
         upper = sac$ub(relative_sample_sizes)
  )
})

all_results_df <- bind_rows(results_list)
# factor r for plotting
all_results_df <- all_results_df |>
  mutate(r = factor(r, levels = sort(unique(r))))
#View(all_results_df)
# save
write.csv(all_results_df, "rarefaction/ecoli_bsi_mlst_preseqR.csv", row.names = FALSE)
#~~~~~~~~~~~#
# Entry point for plotting
#all_results_df <- read.csv("rarefaction/ecoli_bsi_mlst_preseqR_by.csv")
#~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Coverage (aka Sample coverage) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Question: what is the proportion of the bacterial population belong to/have a feature that has already been observed in the sample
# aka: what is the probability of observing an MLST at least r times? (as a function of sample size)

# set r = 1 as threshold of detection
r_values <- c(1)

results_list_sample_cov  <- lapply(r_values, function(rval) {
  message(sprintf("Running preseqR sample coverage for r = %s ...", rval))
  sac <- preseqR.sample.cov.bootstrap(n = profile_counts, r = rval, mt = 20,
                                      times = bootstrap_times, conf = conf_level)
  
  tibble(r = rval, 
         relative_sample_size = relative_sample_sizes ,
         actual_sample_size = actual_sample_sizes ,
         est = sac$f(relative_sample_sizes ),
         se = sac$se(relative_sample_sizes ),
         lower = sac$lb(relative_sample_sizes ),
         upper = sac$ub(relative_sample_sizes )
  )
})

# bind results and factor r
all_results_sample_cov_df <- bind_rows(results_list_sample_cov)
all_results_sample_cov_df <- all_results_sample_cov_df |>
  mutate(r = factor(r, levels = sort(unique(r))))
#View(all_results_sample_cov_df)

# save
write.csv(all_results_sample_cov_df, "rarefaction/ecoli_bsi_mlst_preseqR_coverage.csv", row.names = FALSE)
#~~~~~~~~~~~#
# Entry point for plotting
#all_results_sample_cov_df <- read.csv("rarefaction/ecoli_bsi_mlst_preseqR_coverage.csv")
#~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.ai Klebsiella preseqR ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# count mlst profile frequencies
profile_counts <- kleb_bsi_samples_metadata |>
  mutate(klebsiella_mlst_ST  = as.character(klebsiella_mlst_ST )) |> # use kleborate MLST assignation
  group_by(klebsiella_mlst_ST ) |>
  summarise(n = n()) |>
  group_by(n) |>
  dplyr::summarise(no_mlsts = n_distinct(klebsiella_mlst_ST ))
#View(profile_counts)
profile_counts <- as.matrix(profile_counts)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Species richness (i.e. species accumulation curve (SAC) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define current and estiamted upsampling sample sizes for Kleb BSIs:
current_n <- nrow(kleb_bsi_samples_metadata) #468
estimated_total_n <- current_n / 0.65 #720 - capture fraction from SGSS
estimated_national_n <- 13078 # from 23/24 FY

#unique mlsts observed currently
actual_unique_mlsts <- length(unique(kleb_bsi_samples_metadata$klebsiella_mlst_ST)) # 295

# define sample sizes
relative_sample_sizes <- c(seq(0,30, 0.1), (estimated_total_n/current_n), (estimated_national_n / current_n))
relative_sample_sizes <- sort(relative_sample_sizes)
length(relative_sample_sizes)
actual_sample_sizes <- relative_sample_sizes * current_n
actual_sample_sizes <- sort(actual_sample_sizes)
# define parameters, including r = minimum number of times it is important to detect an MLST.
r_values <- c(1)
bootstrap_times <- 1000
conf_level <- 0.95

# use lapply to apply preseqR.rSAc.bootstrap function over all values of r_values
results_list <- lapply(r_values, function(rval) {
  message(sprintf("Running preseqR for r = %s ...", rval))
  sac <- preseqR.rSAC.bootstrap(n = profile_counts, r = rval, mt = 20,
                                size=SIZE.INIT, mu=MU.INIT, times = bootstrap_times, conf = conf_level)
  
  tibble(r = rval, 
         relative_sample_size = relative_sample_sizes,
         actual_sample_size = actual_sample_sizes,
         est = sac$f(relative_sample_sizes),
         se = sac$se(relative_sample_sizes),
         lower = sac$lb(relative_sample_sizes),
         upper = sac$ub(relative_sample_sizes)
  )
})

all_results_df <- bind_rows(results_list)
# factor r for plotting
all_results_df <- all_results_df |> mutate(r = factor(r, levels = sort(unique(r))))
#View(all_results_df)
# save
write.csv(all_results_df, "rarefaction/kleb_bsi_mlst_preseqR.csv", row.names = FALSE)
#~~~~~~~~~~~#
# Entry point for plotting
#all_results_df <- read.csv("rarefaction/kleb_bsi_mlst_preseqR.csv")
#~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Coverage (aka Sample coverage) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
r_values <- c(1)

results_list_sample_cov  <- lapply(r_values, function(rval) {
  message(sprintf("Running preseqR sample coverage for r = %s ...", rval))
  sac <- preseqR.sample.cov.bootstrap(n = profile_counts, r = rval, mt = 20,
                                      times = bootstrap_times, conf = conf_level)
  
  tibble(r = rval, 
         relative_sample_size = relative_sample_sizes ,
         actual_sample_size = actual_sample_sizes ,
         est = sac$f(relative_sample_sizes ),
         se = sac$se(relative_sample_sizes ),
         lower = sac$lb(relative_sample_sizes ),
         upper = sac$ub(relative_sample_sizes )
  )
})

all_results_sample_cov_df <- bind_rows(results_list_sample_cov)
all_results_sample_cov_df <- all_results_sample_cov_df |> mutate(r = factor(r, levels = sort(unique(r))))
#View(all_results_sample_cov_df)
#save
write.csv(all_results_sample_cov_df, "rarefaction/kleb_bsi_mlst_preseqR_sample_cov.csv", row.names = FALSE)
#~~~~~~~~~~~#
# Entry point for plotting
#all_results_sample_cov_df <- read.csv("rarefaction/kleb_bsi_mlst_preseqR_sample_cov.csv")
#~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 1.b iNEXT ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.bi E. coli iNEXT ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
install.packages("iNEXT")
library(iNEXT)

# Input: Abundance data (when counting how may of MLSTs present) 
length(unique(ecoli_bsi_samples_metadata$escherichia__mlst_achtman__ST)) # 263 unique mlsts
# prepare input data as matrix with rownames as mlst_profiles and column names as sites, with first column being overall. 
ecoli_bsi_mlst_abundance_df <- ecoli_bsi_samples_metadata |>
  group_by(escherichia__mlst_achtman__ST, region) |>
  summarise(n = n()) |>
  pivot_wider(names_from = region, values_from = n, values_fill = 0) |>
  mutate(Total = sum(Midlands, `North East A`, `North East B`, `North West`, `South West`, 
                     `East`, `London`,`South East A`, `South East B`, `South East C`)) |>
  dplyr::select(escherichia__mlst_achtman__ST, Total, everything())
ecoli_bsi_mlst_abundance_df <- data.frame(ecoli_bsi_mlst_abundance_df)
rownames(ecoli_bsi_mlst_abundance_df) <- ecoli_bsi_mlst_abundance_df[,1]
ecoli_bsi_mlst_abundance_df <- ecoli_bsi_mlst_abundance_df[,-1]
#View(ecoli_bsi_mlst_abundance_df)
sum(ecoli_bsi_mlst_abundance_df$Total) # check

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Species Richness (Species accumulation curve) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# get vector of total N sample sizes for which to estimate diversity
# define current and estiamted upsampling sample sizes for E.coli BSIs:
current_n <- nrow(ecoli_bsi_samples_metadata) #1471
estimated_total_n <- current_n / 0.65 #2263
estimated_national_n <- 42224

#unique mlsts observed currently
actual_unique_mlsts <- length(unique(ecoli_bsi_samples_metadata$escherichia__mlst_achtman__ST)) # 287

# define sample sizes
relative_sample_sizes <- c(seq(0,30, 0.1), (estimated_total_n/current_n), (estimated_national_n / current_n))
relative_sample_sizes <- sort(relative_sample_sizes)
#length(relative_sample_sizes)
actual_sample_sizes <- relative_sample_sizes * current_n
actual_sample_sizes <- sort(actual_sample_sizes)
str(actual_sample_sizes)
# make integer
actual_sample_sizes <- as.integer(round(actual_sample_sizes))

# for different diveristy measures; q = 0 -> richness; q = 1 -> Shannon diversity; q = 2 -> simpron diversity
ecoli_bsi_iNEXT <- iNEXT(ecoli_bsi_mlst_abundance_df, q=c(0,1,2), datatype="abundance", size = actual_sample_sizes)
# inspect
#class(ecoli_bsi_iNEXT)
#str(ecoli_bsi_iNEXT)
ecoli_bsi_iNEXT$DataInfo
ecoli_bsi_iNEXT$iNextEst
ecoli_bsi_iNEXT$iNextEst$size_based # CIs obtained for fixed sample size
ecoli_bsi_iNEXT$iNextEst$coverage_based # CIs obtained for coverage
ecoli_bsi_iNEXT$AsyEst

# save
saveRDS(ecoli_bsi_iNEXT, "rarefaction/ecoli_bsi_mlst_iNEXT.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#ecoli_bsi_iNEXT <- readRDS("rarefaction/ecoli_bsi_mlst_iNEXT.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# covert p to factor for plotting
ecoli_bsi_iNEXT_size_based <- as.data.frame(ecoli_bsi_iNEXT$iNextEst$size_based)
ecoli_bsi_iNEXT_size_based$Order.q <- factor(ecoli_bsi_iNEXT_size_based$Order.q, levels = sort(unique(ecoli_bsi_iNEXT_size_based$Order.q)))
ecoli_bsi_iNEXT_coverage_based <- as.data.frame(ecoli_bsi_iNEXT$iNextEst$coverage_based)
ecoli_bsi_iNEXT_coverage_based$Order.q <- factor(ecoli_bsi_iNEXT_coverage_based$Order.q)
ecoli_bsi_asy <- as.data.frame(ecoli_bsi_iNEXT$AsyEst)
ecoli_bsi_asy$Diversity <- factor(ecoli_bsi_asy$Diversity)
#View(ecoli_bsi_iNEXT_size_based)


# split off total data
ecoli_bsi_iNEXT_size_based_total <- ecoli_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total" & m != 0 & m != 1)

ecoli_bsi_iNEXT_coverage_based_total <- ecoli_bsi_iNEXT_coverage_based |>
  dplyr::filter(Assemblage == "Total" & m != 0 & m != 1 )

ecoli_bsi_asy <- ecoli_bsi_asy |>
  mutate(Order.q = case_when(Diversity == "Species richness" ~ 0,
                             Diversity == "Shannon diversity" ~ 1,
                             Diversity == "Simpson diversity" ~ 2))
ecoli_bsi_asy$Order.q <- factor(ecoli_bsi_asy$Order.q, levels = sort(unique(ecoli_bsi_asy$Order.q)))
ecoli_bsi_asy_total <- ecoli_bsi_asy |>
  dplyr::filter(Assemblage == "Total")

# save
write.csv(ecoli_bsi_iNEXT_size_based_total, "rarefaction/ecoli_bsi_iNEXT_size_based_total.csv",  row.names = FALSE)
write.csv(ecoli_bsi_iNEXT_coverage_based_total, "rarefaction/ecoli_bsi_iNEXT_coverage_based_total.csv",  row.names = FALSE)
write.csv(ecoli_bsi_asy_total, "rarefaction/ecoli_bsi_asy_total.csv",  row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.bii Klebsiella iNEXT ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Abundance data (when counting how may of MLSTs present) 
#  sample size based sampling curve 
length(unique(kleb_bsi_samples_metadata$klebsiella_mlst_ST)) # 295 unique mlsts
# prepare input data as matrix with rownames as mlst_profiles and column names as sites, with first column being overall. 
kleb_bsi_mlst_abundance_df <- kleb_bsi_samples_metadata |>
  group_by(klebsiella_mlst_ST, region) |>
  summarise(n = n()) |>
  pivot_wider(names_from = region, values_from = n, values_fill = 0) |>
  mutate(Total = sum(Midlands, `North East A`, `North East B`, `North West`, `South West`, 
                     `East`, `London`,`South East A`, `South East B`, `South East C`)) |>
  dplyr::select(klebsiella_mlst_ST, Total, everything())
kleb_bsi_mlst_abundance_df <- data.frame(kleb_bsi_mlst_abundance_df)
rownames(kleb_bsi_mlst_abundance_df) <- kleb_bsi_mlst_abundance_df[,1]
kleb_bsi_mlst_abundance_df <- kleb_bsi_mlst_abundance_df[,-1]
#View(kleb_bsi_mlst_abundance_df)
sum(kleb_bsi_mlst_abundance_df$Total) # 468 check

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Species Richness (Species accummulation curve) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# get vector of sample sizes for which to estimate diversity
# define current and estiamted upsampling sample sizes for Kleb BSIs:
current_n <- nrow(kleb_bsi_samples_metadata) #468
estimated_total_n <- current_n / 0.65 #720
estimated_national_n <- 13078

#unique mlsts observed currently
actual_unique_mlsts <- length(unique(kleb_bsi_samples_metadata$klebsiella_mlst_ST)) # 297

# define sample sizes
relative_sample_sizes <- c(seq(0,30, 0.1), (estimated_total_n/current_n), (estimated_national_n / current_n))
relative_sample_sizes <- sort(relative_sample_sizes)
#length(relative_sample_sizes)
actual_sample_sizes <- relative_sample_sizes * current_n
actual_sample_sizes <- sort(actual_sample_sizes)
str(actual_sample_sizes)
# make integer
actual_sample_sizes <- as.integer(round(actual_sample_sizes))

kleb_bsi_iNEXT <- iNEXT(kleb_bsi_mlst_abundance_df, q=c(0,1,2), datatype="abundance", size = actual_sample_sizes)
# inspect
#class(kleb_bsi_iNEXT)
#str(kleb_bsi_iNEXT)
kleb_bsi_iNEXT$DataInfo
kleb_bsi_iNEXT$iNextEst
kleb_bsi_iNEXT$iNextEst$size_based # CIs obtained for fixed sample size
kleb_bsi_iNEXT$iNextEst$coverage_based # CIs obtained for coverage
kleb_bsi_iNEXT$AsyEst

# save
saveRDS(kleb_bsi_iNEXT, "rarefaction/kleb_bsi_mlst_iNEXT.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#kleb_bsi_iNEXT <- readRDS("rarefaction/kleb_bsi_iNEXT.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# covert p to factor for plotting
kleb_bsi_iNEXT_size_based <- as.data.frame(kleb_bsi_iNEXT$iNextEst$size_based)
kleb_bsi_iNEXT_size_based$Order.q <- factor(kleb_bsi_iNEXT_size_based$Order.q, levels = sort(unique(kleb_bsi_iNEXT_size_based$Order.q)))
kleb_bsi_iNEXT_coverage_based <- as.data.frame(kleb_bsi_iNEXT$iNextEst$coverage_based)
kleb_bsi_iNEXT_coverage_based$Order.q <- factor(kleb_bsi_iNEXT_coverage_based$Order.q)
kleb_bsi_asy <- as.data.frame(kleb_bsi_iNEXT$AsyEst)
kleb_bsi_asy$Diversity <- factor(kleb_bsi_asy$Diversity)
#View(kleb_bsi_iNEXT_size_based)

# split into regional and total data
kleb_bsi_iNEXT_size_based_total <- kleb_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total" & m != 0 & m != 1)

kleb_bsi_iNEXT_coverage_based_total <- kleb_bsi_iNEXT_coverage_based |>
  dplyr::filter(Assemblage == "Total" & m != 0 & m != 1 )

kleb_bsi_asy <- kleb_bsi_asy |>
  mutate(Order.q = case_when(Diversity == "Species richness" ~ 0,
                             Diversity == "Shannon diversity" ~ 1,
                             Diversity == "Simpson diversity" ~ 2))
kleb_bsi_asy$Order.q <- factor(kleb_bsi_asy$Order.q, levels = sort(unique(kleb_bsi_asy$Order.q)))
kleb_bsi_asy_total <- kleb_bsi_asy |>
  dplyr::filter(Assemblage == "Total")

# save
write.csv(kleb_bsi_iNEXT_size_based_total, "rarefaction/kleb_bsi_iNEXT_size_based_total.csv", row.names = FALSE)
write.csv(kleb_bsi_iNEXT_coverage_based_total, "rarefaction/kleb_bsi_iNEXT_coverage_based_total.csv",  row.names = FALSE)
write.csv(kleb_bsi_asy_total, "rarefaction/kleb_bsi_asy_total.csv",  row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Compare Preseq and iNEXT Sample coverage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load ecoli preseq data
ecoli_bsi_mlst_preseq <- read.csv("rarefaction/ecoli_bsi_mlst_preseqR.csv")
# load ecoli iNEXT data
ecoli_bsi_iNEXT_size_based_total <- read.csv("rarefaction/ecoli_bsi_iNEXT_size_based_total.csv")
# load kleb preseq data
kleb_bsi_mlst_preseq <- read.csv("rarefaction/kleb_bsi_mlst_preseqR.csv")
# load kleb iNEXT data
kleb_bsi_iNEXT_size_based_total <- read.csv("rarefaction/kleb_bsi_iNEXT_size_based_total.csv")

# prep data
ecoli_bsi_mlst_preseq <- ecoli_bsi_mlst_preseq |>
  dplyr::filter(r == 1) |>
  mutate(Method = "preseqR",
         Genus = "Escherichia") |>
  dplyr::select(Genus, Method, r, actual_sample_size, est, lower, upper)

ecoli_bsi_iNEXT_size_based_total <- ecoli_bsi_iNEXT_size_based_total |>
  dplyr::filter(Assemblage == "Total", Order.q == 0) |>
  mutate(Method = "iNEXT", 
         r = "iNEXT",
         Genus = "Escherichia") |>
  rename(actual_sample_size = m ,
         est = qD,
         lower = qD.LCL, 
         upper = qD.UCL) |>
  dplyr::select(c(Genus, Method, r, actual_sample_size, est, lower, upper))


kleb_bsi_mlst_preseq <- kleb_bsi_mlst_preseq |>
  dplyr::filter(r == 1) |>
  mutate(Method = "preseqR",
         Genus = "Klebsiella") |>
  dplyr::select(Genus, Method, r, actual_sample_size, est, lower, upper)

kleb_bsi_iNEXT_size_based_total <- kleb_bsi_iNEXT_size_based_total |>
  dplyr::filter(Assemblage == "Total", Order.q == 0) |>
  mutate(Method = "iNEXT", 
         r = "iNEXT",
         Genus = "Klebsiella") |>
  rename(actual_sample_size = m ,
         est = qD,
         lower = qD.LCL, 
         upper = qD.UCL) |>
  dplyr::select(c(Genus, Method, r, actual_sample_size, est, lower, upper))

# rbind
combined_ecological_estimators <- rbind(ecoli_bsi_mlst_preseq, ecoli_bsi_iNEXT_size_based_total,
                                        kleb_bsi_mlst_preseq, kleb_bsi_iNEXT_size_based_total)
#View(combined_ecological_estimators)
#colnames(combined_ecological_estimators)

# plot
make_panel <- function(dat, genus_name, neksus_x, total_bsi_x) {
  panel_max <- max(c(dat$est, dat$upper), na.rm = TRUE)
  ggplot(dat, aes(x = actual_sample_size, y = est, colour = Method, fill = Method)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
    geom_vline(xintercept = neksus_x, linetype = "longdash", colour = "grey20", linewidth = 0.8) +
    geom_vline(xintercept = total_bsi_x, linetype = "longdash", colour = "grey20", linewidth = 0.8) +
    annotate(
      "text", x = neksus_x, y = panel_max * 1.05,
      label = "NEKSUS sample", angle = 90,
      vjust = -0.5, hjust = 0.5, size = 3.5
    ) +
    annotate(
      "text", x = total_bsi_x, y = panel_max * 1.05,
      label = "Total BSIs in 23/24 FY", angle = 90,
      vjust = -0.5, hjust = 0.5, size = 3.5
    ) +
    labs(
      title = genus_name,
      x = "Sample Size",
      y = "Species Richness"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      
    )
}

p_ecoli <- make_panel(
  dat = combined_ecological_estimators |> filter(Genus == "Escherichia"),
  genus_name = "Escherichia",
  neksus_x = 1471,
  total_bsi_x = 42224
)

p_kleb <- make_panel(
  dat = combined_ecological_estimators |> filter(Genus == "Klebsiella"),
  genus_name = "Klebsiella",
  neksus_x = 468,
  total_bsi_x = 13078
)

p <- p_ecoli + p_kleb + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p
#save
ggsave(filename = "rarefaction/combined_ecoli_kleb_bsi_mlst_iNEXT_vs_preseq_count.png", plot = p, width = 14, height = 7, units = "in", dpi = 300 )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. RAREFACTION OF NEKSUS DATA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 2.a MLST rarefaction ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make alternative function that calculates species richness, shannon and simpson diversity
rarefy_species <- function(df,
                           feature = "mlst",
                           sample_col = NULL,
                           sample_sizes = NULL,
                           max_samples = NULL,
                           reps = 100,
                           replacement = TRUE,
                           n_cores = NULL,
                           seed = NULL,
                           verbose = TRUE) {
  
  # df: data.frame with one row per individual (isolate)
  # feature: column name (string) containing the species/MLST feature
  # sample_col: optional column name for isolate/sample id (not used for sampling here)
  # sample_sizes: vector of sample sizes to evaluate (overrides max_samples)
  # max_samples: if sample_sizes NULL, use 1:max_samples (default = nrow(df))
  # reps: number of replicates per sample size
  # replacement: sample with replacement? (TRUE/FALSE)
  # n_cores: number of cores to use; default = detectCores() - 1 (at least 1)
  # seed: optional seed for reproducibility
  # returns: data.frame with one row per sample size and columns:
  #   richness_est, richness_lcl, richness_ucl,
  #   shannon_est, shannon_lcl, shannon_ucl,
  #   simpson_est, simpson_lcl, simpson_ucl,
  #   coverage_est, coverage_lcl, coverage_ucl
  
  # Input checks
  if (!is.data.frame(df)) stop("df must be a data.frame")
  if (!(feature %in% names(df))) stop("feature must be a column name in df")
  N_total <- nrow(df)
  if (is.null(max_samples)) max_samples <- N_total
  if (is.null(sample_sizes)) sample_sizes <- seq_len(max_samples)
  if (any(sample_sizes < 1)) stop("sample_sizes must be >= 1")
  if (!replacement && any(sample_sizes > N_total)) {
    stop("Some requested sample_sizes exceed number of rows in df while replacement = FALSE")
  }
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  } else {
    n_cores <- max(1, as.integer(n_cores))
  }
  if (!is.null(seed)) set.seed(seed)
  
  # helper: compute metrics for a sample of features
  compute_metrics_from_features <- function(features_vec) {
    features_vec <- features_vec[!is.na(features_vec)]
    N <- length(features_vec)
    if (N == 0) {
      return(list(richness = 0, shannon = NA_real_, simpson = NA_real_, coverage = NA_real_))
    }
    counts <- as.integer(table(features_vec))
    richness <- length(counts)
    p <- counts / sum(counts)
    shannon <- - sum(p * log(p))
    if (N <= 1) {
      simpson <- NA_real_
    } else {
      simpson <- 1 - sum(counts * (counts - 1)) / (N * (N - 1))
    }
    f1 <- sum(counts == 1)
    coverage <- 1 - (f1 / N)
    list(richness = richness, shannon = shannon, simpson = simpson, coverage = coverage)
  }
  
  # worker for a single sample size s
  run_reps_for_s <- function(s) {
    reps_results <- vector("list", reps)
    for (i in seq_len(reps)) {
      idx <- if (replacement) {
        sample.int(N_total, size = s, replace = TRUE)
      } else {
        sample.int(N_total, size = s, replace = FALSE)
      }
      feats <- df[[feature]][idx]
      reps_results[[i]] <- compute_metrics_from_features(feats)
    }
    richness_v <- vapply(reps_results, function(x) x$richness, numeric(1))
    shannon_v <- vapply(reps_results, function(x) x$shannon, numeric(1))
    simpson_v <- vapply(reps_results, function(x) x$simpson, numeric(1))
    coverage_v <- vapply(reps_results, function(x) x$coverage, numeric(1))
    
    summarise_metric <- function(v) {
      v_non_na <- v[!is.na(v)]
      if (length(v_non_na) == 0) {
        return(c(est = NA_real_, lcl = NA_real_, ucl = NA_real_))
      }
      est <- mean(v_non_na)
      lcl <- as.numeric(quantile(v_non_na, probs = 0.025, na.rm = TRUE, type = 6))
      ucl <- as.numeric(quantile(v_non_na, probs = 0.975, na.rm = TRUE, type = 6))
      c(est = est, lcl = lcl, ucl = ucl)
    }
    
    r_rich <- summarise_metric(richness_v)
    r_shan <- summarise_metric(shannon_v)
    r_simp <- summarise_metric(simpson_v)
    r_cov  <- summarise_metric(coverage_v)
    
    data.frame(
      sample_size = s,
      richness_est = r_rich["est"], richness_lcl = r_rich["lcl"], richness_ucl = r_rich["ucl"],
      shannon_est  = r_shan["est"],  shannon_lcl  = r_shan["lcl"],  shannon_ucl  = r_shan["ucl"],
      simpson_est  = r_simp["est"],  simpson_lcl  = r_simp["lcl"],  simpson_ucl  = r_simp["ucl"],
      coverage_est = r_cov["est"],   coverage_lcl = r_cov["lcl"],   coverage_ucl = r_cov["ucl"],
      stringsAsFactors = FALSE
    )
  }
  
  # run in parallel
  if (.Platform$OS.type != "windows") {
    if (verbose) message("Using mclapply with ", n_cores, " cores")
    results_list <- parallel::mclapply(sample_sizes, FUN = run_reps_for_s, mc.cores = n_cores)
    results_df <- do.call(rbind, results_list)
  } else {
    if (verbose) message("Using parLapply with ", n_cores, " cores (Windows fallback)")
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, varlist = c("df", "feature", "N_total", "reps", "replacement",
                                            "compute_metrics_from_features", "run_reps_for_s"),
                            envir = environment())
    results_list <- parallel::parLapply(cl, sample_sizes, fun = run_reps_for_s)
    parallel::stopCluster(cl)
    results_df <- do.call(rbind, results_list)
  }
  
  # Pivot to long/tidy format
  # Desired Metrics and column mapping
  long_df <- results_df |>
    tidyr::pivot_longer(
      cols = -sample_size,
      names_to = c("Metric", ".value"),
      names_pattern = "(.*)_(est|lcl|ucl)"
    ) |>
    mutate(
      Metric = dplyr::case_when(
        Metric == "richness" ~ "Species richness",
        Metric == "shannon"  ~ "Shannon diversity",
        Metric == "simpson"  ~ "Simpson diversity",
        Metric == "coverage" ~ "Sample coverage",
        TRUE ~ Metric
      )
    ) |>
    dplyr::select(sample_size, Metric, est, lcl, ucl) |>
    dplyr::arrange(sample_size, Metric)
  
  rownames(long_df) <- NULL
  return(long_df)
}



# Run rarefaction
max_samples <- length(unique(ecoli_bsi_samples_metadata$isolateid)) # 1471
n_reps <- 1000
# with replacement
# df must have column "mlst" (or change feature param)
ecoli_bsi_mlst_rarefaction_replacement <- rarefy_species(df = ecoli_bsi_samples_metadata,
                                                         feature = "escherichia__mlst_achtman__ST", max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = 12, 
                                                         #seed = 42
)

# run rarefactionwithout replacement
ecoli_bsi_mlst_rarefaction_no_replacement <- rarefy_species(df = ecoli_bsi_samples_metadata,
                                                            feature = "escherichia__mlst_achtman__ST", max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = 12, 
                                                            #seed = 42
)

# Clean, rescale and combine tables
# rbind 2 dfs after adding column to identify replacement vs no replacement
ecoli_bsi_mlst_rarefaction_replacement <- ecoli_bsi_mlst_rarefaction_replacement |>
  mutate(Replacement = "With replacement")
ecoli_bsi_mlst_rarefaction_no_replacement <- ecoli_bsi_mlst_rarefaction_no_replacement |>
  mutate(Replacement = "No replacement")
ecoli_bsi_mlst_rarefaction_combined <- rbind(ecoli_bsi_mlst_rarefaction_replacement,
                                             ecoli_bsi_mlst_rarefaction_no_replacement)

# add proportion of diversity captured to output df
total_mlst_diversity <- length(unique(ecoli_bsi_samples_metadata$escherichia__mlst_achtman__ST)) #263
max_shannon <- round(max(ecoli_bsi_mlst_rarefaction_replacement[ecoli_bsi_mlst_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

ecoli_bsi_mlst_rarefaction_combined <- ecoli_bsi_mlst_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_mlst_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor)|>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_mlst_rarefaction_combined)

#reorder metircs:
ecoli_bsi_mlst_rarefaction_combined$Metric <- factor(ecoli_bsi_mlst_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
ecoli_bsi_mlst_rarefaction_combined$Replacement <- factor(ecoli_bsi_mlst_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# save output table
write.csv(ecoli_bsi_mlst_rarefaction_combined, file = "rarefaction/ecoli_bsi_mlst_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#ecoli_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_mlst_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# repeat for klebsiella
# Run rarefaction
max_samples <- length(unique(kleb_bsi_samples_metadata$isolateid)) # 468
n_reps <- 1000
# with replacement
# df must have column "mlst" (or change feature param)
kleb_bsi_mlst_rarefaction_replacement <- rarefy_species(df = kleb_bsi_samples_metadata,
                                                         feature = "klebsiella_mlst_ST ", max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = 12, 
                                                         #seed = 42
)

# run rarefactionwithout replacement
kleb_bsi_mlst_rarefaction_no_replacement <- rarefy_species(df = kleb_bsi_samples_metadata,
                                                            feature = "klebsiella_mlst_ST ", max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = 12, 
                                                            #seed = 42
)

# Clean, rescale and combine tables
# rbind 2 dfs after adding column to identify replacement vs no replacement
kleb_bsi_mlst_rarefaction_replacement <- kleb_bsi_mlst_rarefaction_replacement |>
  mutate(Replacement = "With replacement")
kleb_bsi_mlst_rarefaction_no_replacement <- kleb_bsi_mlst_rarefaction_no_replacement |>
  mutate(Replacement = "No replacement")
kleb_bsi_mlst_rarefaction_combined <- rbind(kleb_bsi_mlst_rarefaction_replacement,
                                             kleb_bsi_mlst_rarefaction_no_replacement)

# add proportion of diversity captured to output df
total_mlst_diversity <- length(unique(kleb_bsi_samples_metadata$escherichia__mlst_achtman__ST)) #263
max_shannon <- round(max(kleb_bsi_mlst_rarefaction_replacement[kleb_bsi_mlst_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

kleb_bsi_mlst_rarefaction_combined <- kleb_bsi_mlst_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_mlst_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor)|>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_mlst_rarefaction_combined)

#reorder metircs:
kleb_bsi_mlst_rarefaction_combined$Metric <- factor(kleb_bsi_mlst_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
kleb_bsi_mlst_rarefaction_combined$Replacement <- factor(kleb_bsi_mlst_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# save output table
write.csv(kleb_bsi_mlst_rarefaction_combined, file = "rarefaction/kleb_bsi_mlst_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#kleb_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_mlst_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 2.b Key Gene Rarefaction ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rarefy_genes <- function(df,
                         gene_col = "gene",
                         sample_col = "sample",
                         sample_sizes = NULL,
                         max_samples = NULL,
                         reps = 100,
                         replacement = TRUE,
                         n_cores = NULL,
                         seed = NULL,
                         verbose = TRUE) {
  
  if (!is.data.frame(df)) stop("df must be a data.frame")
  if (!(gene_col %in% names(df))) stop("gene_col must be a column name in df")
  if (!(sample_col %in% names(df))) stop("sample_col must be a column name in df")
  
  # Unique sampling units
  samples <- unique(df[[sample_col]])
  samples <- samples[!is.na(samples)]
  N_units <- length(samples)
  if (N_units == 0) stop("No valid sample IDs found in sample_col")
  
  if (is.null(max_samples)) max_samples <- N_units
  if (is.null(sample_sizes)) sample_sizes <- seq_len(max_samples)
  
  if (any(sample_sizes < 1)) stop("sample_sizes must be >= 1")
  if (!replacement && any(sample_sizes > N_units)) {
    stop("Some requested sample_sizes exceed number of unique samples while replacement = FALSE")
  }
  
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  } else {
    n_cores <- max(1, as.integer(n_cores))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # Precompute: for each sample, the unique genes present (excluding NA)
  # This makes sampling fast.
  # gene_list[[sample_id]] = character vector of unique genes in that sample
  df2 <- df[, c(sample_col, gene_col)]
  df2 <- df2[!is.na(df2[[sample_col]]), , drop = FALSE]
  
  # Remove NA genes for the gene sets; samples with only NA remain as empty vectors
  df_non_na <- df2[!is.na(df2[[gene_col]]), , drop = FALSE]
  
  if (nrow(df_non_na) == 0) {
    # all gene entries NA -> no genes anywhere
    gene_list <- setNames(vector("list", length(samples)), samples)
    for (s in samples) gene_list[[as.character(s)]] <- character(0)
  } else {
    # unique gene presence per sample
    df_non_na <- unique(df_non_na)
    gene_list <- split(df_non_na[[gene_col]], df_non_na[[sample_col]])
    # ensure all samples exist (including those with zero genes)
    missing_samples <- setdiff(as.character(samples), names(gene_list))
    if (length(missing_samples) > 0) {
      for (s in missing_samples) gene_list[[s]] <- character(0)
    }
  }
  
  # Metric computation from a sampled set of sample IDs
  compute_metrics_from_sample_ids <- function(sample_ids) {
    # combine gene presence across sampled samples
    genes_all <- unlist(gene_list[as.character(sample_ids)], use.names = FALSE)
    
    # gene occurrence counts across sampled samples (presence/absence per sample already)
    # If sampling WITH replacement, the same sample may appear multiple times; we treat that as repeated draws.
    # So genes from repeated sample draws count repeatedly, consistent with replacement sampling.
    if (length(genes_all) == 0) {
      # no genes observed in sampled set
      return(list(richness = 0, shannon = NA_real_, simpson = NA_real_, coverage = NA_real_))
    }
    
    counts <- as.integer(table(genes_all))
    richness <- length(counts)
    
    N <- sum(counts)
    p <- counts / N
    shannon <- -sum(p * log(p))
    
    if (N <= 1) {
      simpson <- NA_real_
    } else {
      simpson <- 1 - sum(counts * (counts - 1)) / (N * (N - 1))
    }
    
    # Good's coverage on gene-occurrence sample:
    f1 <- sum(counts == 1)
    coverage <- 1 - (f1 / N)
    
    list(richness = richness, shannon = shannon, simpson = simpson, coverage = coverage)
  }
  
  summarise_metric <- function(v) {
    v_non_na <- v[!is.na(v)]
    if (length(v_non_na) == 0) return(c(est = NA_real_, lcl = NA_real_, ucl = NA_real_))
    est <- mean(v_non_na)
    lcl <- as.numeric(quantile(v_non_na, probs = 0.025, na.rm = TRUE, type = 6))
    ucl <- as.numeric(quantile(v_non_na, probs = 0.975, na.rm = TRUE, type = 6))
    c(est = est, lcl = lcl, ucl = ucl)
  }
  
  # Worker: run reps replicates for a single sample size
  run_reps_for_s <- function(s) {
    reps_results <- vector("list", reps)
    for (i in seq_len(reps)) {
      chosen <- if (replacement) {
        sample(samples, size = s, replace = TRUE)
      } else {
        sample(samples, size = s, replace = FALSE)
      }
      reps_results[[i]] <- compute_metrics_from_sample_ids(chosen)
    }
    
    richness_v <- vapply(reps_results, function(x) x$richness, numeric(1))
    shannon_v  <- vapply(reps_results, function(x) x$shannon, numeric(1))
    simpson_v  <- vapply(reps_results, function(x) x$simpson, numeric(1))
    coverage_v <- vapply(reps_results, function(x) x$coverage, numeric(1))
    
    r_rich <- summarise_metric(richness_v)
    r_shan <- summarise_metric(shannon_v)
    r_simp <- summarise_metric(simpson_v)
    r_cov  <- summarise_metric(coverage_v)
    
    data.frame(
      sample_size = s,
      richness_est = r_rich["est"], richness_lcl = r_rich["lcl"], richness_ucl = r_rich["ucl"],
      shannon_est  = r_shan["est"],  shannon_lcl  = r_shan["lcl"],  shannon_ucl  = r_shan["ucl"],
      simpson_est  = r_simp["est"],  simpson_lcl  = r_simp["lcl"],  simpson_ucl  = r_simp["ucl"],
      coverage_est = r_cov["est"],   coverage_lcl = r_cov["lcl"],   coverage_ucl = r_cov["ucl"],
      stringsAsFactors = FALSE
    )
  }
  
  # Parallel over sample sizes
  if (.Platform$OS.type != "windows") {
    if (verbose) message("Using mclapply with ", n_cores, " cores")
    results_list <- parallel::mclapply(sample_sizes, FUN = run_reps_for_s, mc.cores = n_cores)
    wide_df <- do.call(rbind, results_list)
  } else {
    if (verbose) message("Using parLapply with ", n_cores, " cores (Windows fallback)")
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(
      cl,
      varlist = c("samples", "gene_list", "reps", "replacement",
                  "compute_metrics_from_sample_ids", "summarise_metric",
                  "run_reps_for_s", "sample_col", "gene_col"),
      envir = environment()
    )
    results_list <- parallel::parLapply(cl, sample_sizes, fun = run_reps_for_s)
    parallel::stopCluster(cl)
    wide_df <- do.call(rbind, results_list)
  }
  
  # Long-form output
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(tidyr)
  library(dplyr)
  
  long_df <- wide_df |>
    tidyr::pivot_longer(
      cols = -sample_size,
      names_to = c("Metric", ".value"),
      names_pattern = "(.*)_(est|lcl|ucl)"
    ) |>
    dplyr::mutate(
      Metric = dplyr::case_when(
        Metric == "richness" ~ "Species richness",
        Metric == "shannon"  ~ "Shannon diversity",
        Metric == "simpson"  ~ "Simpson diversity",
        Metric == "coverage" ~ "Sample coverage",
        TRUE ~ Metric
      )
    ) |>
    dplyr::select(sample_size, Metric, est, lcl, ucl) |>
    dplyr::arrange(sample_size, Metric)
  
  rownames(long_df) <- NULL
  long_df
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# filter amrfinder df for AMR genes only
ecoli_bsi_arg <- ecoli_bsi_amrfinder_metadata |>
  filter(Type == "AMR"| is.na(Type))
#table(ecoli_bsi_arg$Type, useNA = "ifany")
#table(ecoli_bsi_amrfinder_metadata$Type, useNA = "ifany")
# check all isolates present in bith df
length(unique(ecoli_bsi_amrfinder_metadata$sample))
length(unique(ecoli_bsi_arg$sample))

# define run params
max_samples <- length(unique(ecoli_bsi_amrfinder_metadata$sample)) # 1471
n_reps <- 1000
n_cores <- 12

# Run
ecoli_bsi_arg_rarefaction_replacement <- rarefy_genes(ecoli_bsi_arg, sample_col = "sample", gene_col = "Element.symbol",
                                                      max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = n_cores #, seed = 1
)
ecoli_bsi_arg_rarefaction_no_replacement <- rarefy_genes(ecoli_bsi_arg, sample_col = "sample", gene_col = "Element.symbol",
                                                         max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = n_cores #, seed = 1
)

# Clean, rescale and combine tables
# rbind 2 dfs after adding column to identify replacement vs no replacement
ecoli_bsi_arg_rarefaction_replacement <- ecoli_bsi_arg_rarefaction_replacement |>
  mutate(Replacement = "With replacement")
ecoli_bsi_arg_rarefaction_no_replacement <- ecoli_bsi_arg_rarefaction_no_replacement |>
  mutate(Replacement = "No replacement")
ecoli_bsi_arg_rarefaction_combined <- rbind(ecoli_bsi_arg_rarefaction_replacement,
                                            ecoli_bsi_arg_rarefaction_no_replacement)


# combine dfs
# add proportion of diversity captured to output df
total_ARG_diversity <- length(unique(ecoli_bsi_amrfinder_metadata[!is.na(ecoli_bsi_amrfinder_metadata$Element.symbol), ]$Element.symbol)) #331
max_shannon <- round(max(ecoli_bsi_arg_rarefaction_replacement[ecoli_bsi_arg_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

ecoli_bsi_arg_rarefaction_combined <- ecoli_bsi_arg_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_ARG_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor)|>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_arg_rarefaction_combined)
#reorder metircs:
ecoli_bsi_arg_rarefaction_combined$Metric <- factor(ecoli_bsi_arg_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
ecoli_bsi_arg_rarefaction_combined$Replacement <- factor(ecoli_bsi_arg_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# save output table
write.csv(ecoli_bsi_arg_rarefaction_combined, file = "rarefaction/ecoli_bsi_arg_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#ecoli_bsi_arg_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_arg_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# filter amrfinder df for AMR genes only
kleb_bsi_arg <- kleb_bsi_amrfinder_metadata |>
  filter(Type == "AMR"| is.na(Type))
#table(kleb_bsi_arg$Type, useNA = "ifany")
#table(kleb_bsi_amrfinder_metadata$Type, useNA = "ifany")
# check all isolates present in bith df
length(unique(kleb_bsi_amrfinder_metadata$sample))
length(unique(kleb_bsi_arg$sample))

# define run params
max_samples <- length(unique(kleb_bsi_amrfinder_metadata$sample)) # 468
n_reps <- 1000
n_cores <- 12

# Run
kleb_bsi_arg_rarefaction_replacement <- rarefy_genes(kleb_bsi_arg, sample_col = "sample", gene_col = "Element.symbol",
                                                      max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = n_cores #, seed = 1
)
kleb_bsi_arg_rarefaction_no_replacement <- rarefy_genes(kleb_bsi_arg, sample_col = "sample", gene_col = "Element.symbol",
                                                         max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = n_cores #, seed = 1
)

# Clean, rescale and combine tables
# rbind 2 dfs after adding column to identify replacement vs no replacement
kleb_bsi_arg_rarefaction_replacement <- kleb_bsi_arg_rarefaction_replacement |>
  mutate(Replacement = "With replacement")
kleb_bsi_arg_rarefaction_no_replacement <- kleb_bsi_arg_rarefaction_no_replacement |>
  mutate(Replacement = "No replacement")
kleb_bsi_arg_rarefaction_combined <- rbind(kleb_bsi_arg_rarefaction_replacement,
                                            kleb_bsi_arg_rarefaction_no_replacement)


# combine dfs
# add proportion of diversity captured to output df
total_ARG_diversity <- length(unique(kleb_bsi_amrfinder_metadata[!is.na(kleb_bsi_amrfinder_metadata$Element.symbol), ]$Element.symbol)) #331
max_shannon <- round(max(kleb_bsi_arg_rarefaction_replacement[kleb_bsi_arg_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

kleb_bsi_arg_rarefaction_combined <- kleb_bsi_arg_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_ARG_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor)|>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_arg_rarefaction_combined)
#reorder metircs:
kleb_bsi_arg_rarefaction_combined$Metric <- factor(kleb_bsi_arg_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
kleb_bsi_arg_rarefaction_combined$Replacement <- factor(kleb_bsi_arg_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# save output table
write.csv(kleb_bsi_arg_rarefaction_combined, file = "rarefaction/kleb_bsi_arg_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#kleb_bsi_arg_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_arg_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 2.c Plasmid subcommunities Rarefaction ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# check column
ecoli_bsi_amrfinder_metadata$community_subcommunity
# define run params
max_samples <- length(unique(ecoli_bsi_amrfinder_metadata$sample)) # 1471
n_reps <- 1000
n_cores <- 12

# Run
ecoli_bsi_pling_rarefaction_replacement <- rarefy_genes(ecoli_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "community_subcommunity",
                                                        max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = n_cores #, seed = 1
)
ecoli_bsi_pling_rarefaction_no_replacement <- rarefy_genes(ecoli_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "community_subcommunity",
                                                           max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = n_cores #, seed = 1
)

# Clean, rescale and combine tables
# rbind 2 dfs after adding column to identify replacement vs no replacement
ecoli_bsi_pling_rarefaction_replacement <- ecoli_bsi_pling_rarefaction_replacement |>
  mutate(Replacement = "With replacement")
ecoli_bsi_pling_rarefaction_no_replacement <- ecoli_bsi_pling_rarefaction_no_replacement |>
  mutate(Replacement = "No replacement")
ecoli_bsi_pling_rarefaction_combined <- rbind(ecoli_bsi_pling_rarefaction_replacement,
                                              ecoli_bsi_pling_rarefaction_no_replacement)


# combine dfs
# add proportion of diversity captured to output df
total_pling_diversity <- length(unique(ecoli_bsi_amrfinder_metadata[!is.na(ecoli_bsi_amrfinder_metadata$community_subcommunity), ]$community_subcommunity)) #331
max_shannon <- round(max(ecoli_bsi_pling_rarefaction_replacement[ecoli_bsi_pling_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

ecoli_bsi_pling_rarefaction_combined <- ecoli_bsi_pling_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_pling_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor)|>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_pling_rarefaction_combined)
#reorder metircs:
ecoli_bsi_pling_rarefaction_combined$Metric <- factor(ecoli_bsi_pling_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
ecoli_bsi_pling_rarefaction_combined$Replacement <- factor(ecoli_bsi_pling_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# save output table
write.csv(ecoli_bsi_pling_rarefaction_combined, file = "rarefaction/ecoli_bsi_pling_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#ecoli_bsi_pling_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_pling_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#

# Klebsiella plasmids
# check column
kleb_bsi_amrfinder_metadata$community_subcommunity
# define run params
max_samples <- length(unique(kleb_bsi_amrfinder_metadata$sample)) # 1471
n_reps <- 1000
n_cores <- 12

# Run
kleb_bsi_pling_rarefaction_replacement <- rarefy_genes(kleb_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "community_subcommunity",
                                                        max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = n_cores #, seed = 1
)
kleb_bsi_pling_rarefaction_no_replacement <- rarefy_genes(kleb_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "community_subcommunity",
                                                           max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = n_cores #, seed = 1
)

# Clean, rescale and combine tables
# rbind 2 dfs after adding column to identify replacement vs no replacement
kleb_bsi_pling_rarefaction_replacement <- kleb_bsi_pling_rarefaction_replacement |>
  mutate(Replacement = "With replacement")
kleb_bsi_pling_rarefaction_no_replacement <- kleb_bsi_pling_rarefaction_no_replacement |>
  mutate(Replacement = "No replacement")
kleb_bsi_pling_rarefaction_combined <- rbind(kleb_bsi_pling_rarefaction_replacement,
                                              kleb_bsi_pling_rarefaction_no_replacement)


# combine dfs
# add proportion of diversity captured to output df
total_pling_diversity <- length(unique(kleb_bsi_amrfinder_metadata[!is.na(kleb_bsi_amrfinder_metadata$community_subcommunity), ]$community_subcommunity)) #331
max_shannon <- round(max(kleb_bsi_pling_rarefaction_replacement[kleb_bsi_pling_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

kleb_bsi_pling_rarefaction_combined <- kleb_bsi_pling_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_pling_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor)|>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_pling_rarefaction_combined)
#reorder metircs:
kleb_bsi_pling_rarefaction_combined$Metric <- factor(kleb_bsi_pling_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
kleb_bsi_pling_rarefaction_combined$Replacement <- factor(kleb_bsi_pling_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# save output table
write.csv(kleb_bsi_pling_rarefaction_combined, file = "rarefaction/kleb_bsi_pling_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#kleb_bsi_pling_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_pling_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Combined plot with E.coli and Klebsiella ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# for MLST, ARGs and plasmid PLING subcommunities
# load data
#ecoli_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_mlst_rarefaction_combined.csv")
#kleb_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_mlst_rarefaction_combined.csv")
#ecoli_bsi_arg_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_arg_rarefaction_combined.csv")
#kleb_bsi_arg_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_arg_rarefaction_combined.csv")
#ecoli_bsi_pling_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_pling_rarefaction_combined.csv")
#kleb_bsi_pling_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_pling_rarefaction_combined.csv")

# combine per genomic feature and add feature column
combined_bsi_mlst_rarefaction <- rbind(ecoli_bsi_mlst_rarefaction_combined, kleb_bsi_mlst_rarefaction_combined)
combined_bsi_mlst_rarefaction <- combined_bsi_mlst_rarefaction |>
  mutate(Feature = "MLST")

combined_bsi_arg_rarefaction <- rbind(ecoli_bsi_arg_rarefaction_combined, kleb_bsi_arg_rarefaction_combined)
combined_bsi_arg_rarefaction <- combined_bsi_arg_rarefaction |>
  mutate(Feature = "AMR genes")

combined_bsi_pling_rarefaction <- rbind(ecoli_bsi_pling_rarefaction_combined, kleb_bsi_pling_rarefaction_combined)
combined_bsi_pling_rarefaction <- combined_bsi_pling_rarefaction |>
  mutate(Feature = "PLING plasmid subcommunity")

# combine all 3 dfs
combined_bsi_rarefaction <- rbind(combined_bsi_mlst_rarefaction, combined_bsi_arg_rarefaction)
combined_bsi_rarefaction <- rbind(combined_bsi_rarefaction, combined_bsi_pling_rarefaction)


# set order of factors
combined_bsi_rarefaction <- combined_bsi_rarefaction |>
  mutate(Feature = case_when(Feature == "AMR genes" ~ "ARG",
                            TRUE ~ Feature))
combined_bsi_rarefaction$Feature <- factor(combined_bsi_rarefaction$Feature, levels = c("MLST", "ARG", "PLING plasmid subcommunity"))

combined_bsi_rarefaction$Replacement <- factor(combined_bsi_rarefaction$Replacement, levels = c("With replacement", "No replacement"))
combined_bsi_rarefaction <- combined_bsi_rarefaction |>
  mutate(Metric = case_when(Metric == "Sample coverage" ~ "Coverage",
                            TRUE ~ Metric))
combined_bsi_rarefaction$Metric <- factor(combined_bsi_rarefaction$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Coverage"))

combined_bsi_rarefaction <- combined_bsi_rarefaction |>
  mutate(ReplacementGenus = interaction(Replacement, Genus, sep = "_"))

combined_bsi_rarefaction <- combined_bsi_rarefaction |>
  filter(!is.na(est)) |>
  filter(!(Metric == "Simpson diversity" & est <=0.85 ))
nrow(combined_bsi_rarefaction)

pal <- c(
  "No replacement_Escherichia"   = "seagreen3",
  "No replacement_Klebsiella"   = "darkorange",
  "With replacement_Escherichia"= "grey90",
  "With replacement_Klebsiella" = "grey90"
)

# plot combined data
#NOTE plot not in the best form due to axes 
p <- ggplot(data = combined_bsi_rarefaction, aes(x = sample_size, y = est, fill = ReplacementGenus, colour = ReplacementGenus)) +
  geom_line( linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3, colour = NA) +
  scale_colour_manual(
    name = "Sampling scheme",
    values = pal,
    labels = c(
      "No replacement_Escherichia"    = "No replacement – Escherichia",
      "No replacement_Klebsiella"    = "No replacement – Klebsiella",
      "With replacement_Escherichia" = "With replacement – Escherichia",
      "With replacement_Klebsiella"  = "With replacement – Klebsiella"
    )) +
  scale_fill_manual(
    name = "Sampling scheme",
    values = pal,
    labels = c(
      "No replacement_Escherichia"    = "No replacement – Escherichia",
      "No replacement_Klebsiella"    = "No replacement – Klebsiella",
      "With replacement_Escherichia" = "With replacement – Escherichia",
      "With replacement_Klebsiella"  = "With replacement – Klebsiella"
    )) +
  facet_wrap(Feature ~ Metric, scales = "free_y", ncol = 4) +
  
  #facet_grid( rows = vars(Metric),
  # cols = vars(Feature),
  #scales = "free",   # <<< free y (per row)
  #switch = "y" ) +
  
  ggh4x::facetted_pos_scales(
    y = list(
      Metric == "Simpson" ~ scale_y_continuous(limits = c(0.5, 1.0))
    )
  ) +
  theme_minimal(base_size = 16) +
  labs(title = "E. coli and Klebsiella BSI Rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14),
    ## Top strip (Metric)
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 16,
                                #face = "bold"
    ),
    ## Left strip (Feature) — vertical like y-axis
    strip.text.y.left = element_text(angle = 90, size = 16, 
                                     #face = "bold"
    ),
    ## Reduce visual clutter
    panel.spacing = unit(1, "lines")
  )
print(p)
# save
ggsave("rarefaction/combined_bsi_rarefaction.png", plot = p, width = 15, height = 10.5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Rarefaction summary table ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define rarefaction result tables
rarefaction_tables <- list(
  MLST = ecoli_bsi_mlst_rarefaction,
  Gene_AMR = gene_rarefaction_AMR_only,
  Plasmid_subcommunities = ecoli_bsi_plasmid_subcommunity_rarefaction
)

# Define thresholds of interest
thresholds <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)

# Function to get sample size required to reach each threshold
get_sample_sizes <- function(df, thresholds) {
  map_dbl(thresholds, function(thresh) {
    df |>
      filter(prop_diversity >= thresh) |>
      slice_head(n = 1) |>
      pull(sample_size) |>
      (\(x) if (length(x) == 0) NA else x)()
  })
}


# Custom names for the metrics
metric_labels <- c(
  MLST = "MLST",
  Gene_AMR = "ARG",
  Plasmid_subcommunities = "Plasmid Subcomm."
)

# Build wide-format summary table with metrics as rows and thresholds as columns
rarefaction_summary_wide <- map_dfr(
  names(rarefaction_tables),
  function(metric) {
    df <- rarefaction_tables[[metric]]
    sizes <- get_sample_sizes(df, thresholds)
    max_div <- round(max(df$prop_diversity, na.rm = TRUE), 3)
    tibble(
      Metric = metric_labels[[metric]],  # Apply custom label
      `Max. Prop. Diversity` = max_div,
      !!!set_names(sizes, as.character(thresholds))
    )
  }
)

View(rarefaction_summary_wide)

rarefaction_summary_wide$`Max. Prop. Diversity` <- round(rarefaction_summary_wide$`Max. Prop. Diversity`, 3)

threshold_cols <- intersect(names(rarefaction_summary_wide), as.character(thresholds))
rarefaction_summary_wide <- rarefaction_summary_wide |>
  rename_with(
    .fn = ~ paste0(as.numeric(.) * 100, "%"),
    .cols = all_of(threshold_cols)
  )


# Clean table with gt
summary_table <- rarefaction_summary_wide |>
  gt(rowname_col = "Metric") |>
  tab_header(
    title = "E. coli BSI Rarefaction Summary Table",
    subtitle = "Sample size needed to reach diversity thresholds"
  ) |>
  fmt_missing(everything(), missing_text = "--")

print(summary_table)
# Save to Word or HTML
gtsave(summary_table, "rarefaction/ecoli_bsi_rarefaction_summary.docx")
gtsave(summary_table, "rarefaction/ecoli_bsi_rarefaction_summary.html")
write.csv(rarefaction_summary_wide, file = "rarefaction/ecoli_bsi_rarefaction_summary.csv")


# kleb
# Define your rarefaction result tables
rarefaction_tables <- list(
  MLST = kleb_mlst_rarefaction,
  Gene_AMR = gene_rarefaction_AMR_only,
  Plasmid_subcommunities = kleb_bsi_plasmid_subcommunity_rarefaction
)

# Define thresholds of interest
thresholds <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

# Function to get sample size required to reach each threshold
get_sample_sizes <- function(df, thresholds) {
  map_dbl(thresholds, function(thresh) {
    df |>
      filter(prop_diversity >= thresh) |>
      slice_head(n = 1) |>
      pull(sample_size) |>
      (\(x) if (length(x) == 0) NA else x)()
  })
}

# Custom names for the metrics
metric_labels <- c(
  MLST = "MLST",
  Gene_AMR = "ARG",
  Plasmid_subcommunities = "Plasmid PLING Subcommunities"
)

# Build wide-format summary table with metrics as rows and thresholds as columns
rarefaction_summary_wide <- map_dfr(
  names(rarefaction_tables),
  function(metric) {
    df <- rarefaction_tables[[metric]]
    sizes <- get_sample_sizes(df, thresholds)
    max_div <- round(max(df$prop_diversity, na.rm = TRUE), 3)
    tibble(
      Metric = metric_labels[[metric]],  # Apply custom label
      `Max. Prop. Diversity` = max_div,
      !!!set_names(sizes, as.character(thresholds))
    )
  }
)
#View(rarefaction_summary_wide)
rarefaction_summary_wide$`Max. Prop. Diversity` <- round(rarefaction_summary_wide$`Max. Prop. Diversity`, 3)
threshold_cols <- intersect(names(rarefaction_summary_wide), as.character(thresholds))
rarefaction_summary_wide <- rarefaction_summary_wide |>
  rename_with(
    .fn = ~ paste0(as.numeric(.) * 100, "%"),
    .cols = all_of(threshold_cols)
  )


# Clean table with gt
summary_table <- rarefaction_summary_wide |>
  gt(rowname_col = "Metric") |>
  tab_header(
    title = "Klebsiella BSI Rarefaction Summary Table",
    subtitle = "Sample size needed to reach diversity thresholds"
  ) |>
  fmt_missing(everything(), missing_text = "--")

print(summary_table)
# Save to Word or HTML
gtsave(summary_table, "rarefaction/kleb_bsi_rarefaction_summary.docx")
gtsave(summary_table, "rarefaction/kleb_bsi_rarefaction_summary.html")
write.csv(rarefaction_summary_wide, file = "rarefaction/kleb_bsi_rarefaction_summary.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. NOVEL LINEAGE INCURSION ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make function tp simulate species incursion 
# assuming that new species will replace existing species completely at random, which is unlikely to be accurate at higher frequencies, but may hold for lower frequencies.
simulate_species_incursion <- function(df,
                                       feature = "mlst",
                                       sample_col = "sample",
                                       sample_sizes = NULL,
                                       max_samples = NULL,
                                       k_interval = 5,
                                       reps = 10,
                                       new_species_label = "NEW_SPECIES",
                                       n_cores = NULL,
                                       seed = NULL,
                                       verbose = TRUE) {
  # Input checks
  if (!is.data.frame(df)) stop("df must be a data.frame")
  if (!(feature %in% names(df))) stop("feature must be a column name in df")
  if (!(sample_col %in% names(df))) stop("sample_col must be a column name in df")
  # unique samples and rows
  df <- df[!is.na(df[[sample_col]]), , drop = FALSE]
  unique_samples <- unique(df[[sample_col]])
  N_total_units <- length(unique_samples)
  if (N_total_units == 0) stop("No samples found in sample_col")
  
  if (is.null(max_samples)) max_samples <- N_total_units
  if (is.null(sample_sizes)) {
    # default: seq(100, max_samples, by=100) but ensure last value = max_samples
    seqs <- if (max_samples <= 100) unique(c(max_samples)) else seq(100, max_samples, by = 100)
    if (tail(seqs, 1) != max_samples) seqs <- unique(c(seqs, max_samples))
    sample_sizes <- seqs
  }
  if (any(sample_sizes < 1)) stop("sample_sizes must be >= 1")
  if (any(sample_sizes > N_total_units)) stop("sample_sizes cannot exceed number of unique samples")
  reps <- as.integer(reps)
  if (reps < 1) stop("reps must be >= 1")
  
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  } else {
    n_cores <- max(1, as.integer(n_cores))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # For speed: create a named vector feature_by_sample where each sample has its feature
  # If a sample appears multiple times in df (shouldn't), pick the first or ensure df is one-row-per-sample.
  # We'll force one row per sample by taking the first occurrence per sample.
  df_one <- df[!duplicated(df[[sample_col]]), , drop = FALSE]
  # feature value may be NA (sample with no feature)
  feature_by_sample <- setNames(as.character(df_one[[feature]]), df_one[[sample_col]])
  
  # metrics computation for a vector of features (one per sampled unit)
  compute_metrics_from_feature_vec <- function(feat_vec) {
    # feat_vec: character vector length = sample_size, may contain NA if no gene in sample
    feat_vec <- feat_vec[!is.na(feat_vec)]  # drop NA features (samples with no gene)
    N <- length(feat_vec)
    if (N == 0) {
      return(list(richness = 0, shannon = NA_real_, simpson = NA_real_, coverage = NA_real_, new_freq = 0))
    }
    counts <- as.integer(table(feat_vec))
    richness <- length(counts)
    p <- counts / sum(counts)
    shannon <- -sum(p * log(p))
    if (sum(counts) <= 1) {
      simpson <- NA_real_
    } else {
      simpson <- 1 - sum(counts * (counts - 1)) / (sum(counts) * (sum(counts) - 1))
    }
    f1 <- sum(counts == 1)
    coverage <- 1 - (f1 / sum(counts))
    new_freq <- if (new_species_label %in% names(counts)) counts[[as.character(new_species_label)]] / sum(counts) else 0
    list(richness = richness, shannon = shannon, simpson = simpson, coverage = coverage, new_freq = new_freq)
  }
  
  # worker for a single sample_size s
  run_for_sample_size <- function(s) {
    if (verbose) message(sprintf("Processing sample_size = %d ...", s))
    # select a random subset of samples of size s (no replacement)
    sampled_units <- sample(unique_samples, size = s, replace = FALSE)
    # vector of feature values for these sampled units
    base_feat_vec <- as.character(feature_by_sample[as.character(sampled_units)])
    names(base_feat_vec) <- sampled_units
    
    # For each k = 0..s (number to replace), run 'reps' replicates
    ks <- c(seq(0,s, by = k_interval))
    
    # We'll build results as a list of data.frames per k
    out_list_k <- vector("list", length(ks))
    names(out_list_k) <- as.character(ks)
    
    for (k in ks) {
      # for k==0, we can compute once or repetitively; still do reps for consistent variance
      replicate_results <- vector("list", reps)
      for (rep_i in seq_len(reps)) {
        if (k == 0) {
          # no replacement: use base vector
          feat_vec <- base_feat_vec
        } else {
          # choose k distinct sample indices to replace
          replaced <- sample(seq_along(sampled_units), size = k, replace = FALSE)
          feat_vec <- base_feat_vec
          feat_vec[replaced] <- new_species_label
        }
        metrics <- compute_metrics_from_feature_vec(feat_vec)
        replicate_results[[rep_i]] <- c(metrics, list(sample_size = s, k = k, replicate = rep_i))
      }
      # convert replicate_results to data.frame
      df_k <- do.call(rbind, lapply(replicate_results, function(x) {
        data.frame(
          sample_size = x$sample_size,
          k = x$k,
          replicate = x$replicate,
          Species_richness = x$richness,
          Shannon_diversity = x$shannon,
          Simpson_diversity = x$simpson,
          Sample_coverage = x$coverage,
          New_species_frequency = x$new_freq,
          stringsAsFactors = FALSE
        )
      }))
      out_list_k[[as.character(k)]] <- df_k
    } # end ks
    
    # bind all k rows
    df_out_s <- do.call(rbind, out_list_k)
    # return replicate-level wide df for this sample_size
    df_out_s
  } # end run_for_sample_size
  
  # Parallel over sample_sizes
  if (.Platform$OS.type != "windows") {
    if (verbose) message("Using mclapply with ", n_cores, " cores")
    results_list <- parallel::mclapply(sample_sizes, FUN = run_for_sample_size, mc.cores = n_cores)
    replicates_df <- do.call(rbind, results_list)
  } else {
    if (verbose) message("Using parLapply with ", n_cores, " cores (Windows fallback)")
    cl <- parallel::makeCluster(n_cores)
    # export needed objects and functions
    parallel::clusterExport(cl, varlist = c("unique_samples", "feature_by_sample", "reps",
                                            "new_species_label", "compute_metrics_from_feature_vec", "run_for_sample_size"),
                            envir = environment())
    results_list <- parallel::parLapply(cl, sample_sizes, fun = run_for_sample_size)
    parallel::stopCluster(cl)
    replicates_df <- do.call(rbind, results_list)
  }
  
  # Make sure types are correct
  replicates_df$sample_size <- as.integer(replicates_df$sample_size)
  replicates_df$k <- as.integer(replicates_df$k)
  replicates_df$replicate <- as.integer(replicates_df$replicate)
  
  # long form of replicates
  library(tidyr)
  replicates_long <- replicates_df |>
    tidyr::pivot_longer(
      cols = c(Species_richness, Shannon_diversity, Simpson_diversity, Sample_coverage),
      names_to = "Metric",
      values_to = "value"
    ) |>
    dplyr::arrange(sample_size, k, replicate, Metric)
  
  # summary across replicates: mean and percentile CI per (sample_size, k, Metric)
  summary_df <- replicates_long |>
    dplyr::group_by(sample_size, k, Metric) |>
    dplyr::summarise(
      est = mean(value, na.rm = TRUE),
      lcl = as.numeric(quantile(value, probs = 0.025, na.rm = TRUE, type = 6)),
      ucl = as.numeric(quantile(value, probs = 0.975, na.rm = TRUE, type = 6)),
      .groups = "drop"
    ) |>
    dplyr::mutate(new_species_frequency = k / sample_size) |>
    dplyr::select(sample_size, k, new_species_frequency, everything()) |>
    dplyr::arrange(sample_size, k, Metric)
  
  return(summary_df)
}

#set parameters
max_samples <- length(unique(ecoli_bsi_samples_metadata$isolateid))
n_reps <- 1000
n_cores <- 12

# run for e. coli
ecoli_bsi_mlst_incursion <- simulate_species_incursion(ecoli_bsi_samples_metadata,
                                                       feature = "escherichia__mlst_achtman__ST", sample_col = "isolateid",
                                                       sample_sizes = NULL, 
                                                       max_samples = max_samples,
                                                       k_interval = 10,
                                                       reps = n_reps, new_species_label = "NEW_SPECIES",
                                                       n_cores = n_cores, seed = NULL, verbose = TRUE)


# run for Kleb
max_samples <- length(unique(kleb_bsi_samples_metadata$isolateid))
kleb_bsi_mlst_incursion <- simulate_species_incursion(kleb_bsi_samples_metadata,
                                                      feature = "klebsiella_mlst_ST", sample_col = "isolateid",
                                                      sample_sizes = NULL, 
                                                      max_samples = max_samples,
                                                      k_interval = 10,
                                                      reps = n_reps, new_species_label = "NEW_SPECIES",
                                                      n_cores = n_cores, seed = NULL, verbose = TRUE)


ecoli_bsi_mlst_incursion <- ecoli_bsi_mlst_incursion |> mutate(Genus = "Escherichia")
kleb_bsi_mlst_incursion <- kleb_bsi_mlst_incursion |> mutate(Genus = "Klebsiella")

combined_bsi_mlst_incursion <- rbind(ecoli_bsi_mlst_incursion, kleb_bsi_mlst_incursion)
#set stratifying features as factors:
# replace sample coverage with coverage
combined_bsi_mlst_incursion <- combined_bsi_mlst_incursion |>
  mutate(Metric = case_when(Metric == "Sample_coverage" ~ "Coverage",
                            TRUE ~ Metric))
combined_bsi_mlst_incursion$Metric <- factor(combined_bsi_mlst_incursion$Metric, levels = c("Simpson_diversity", "Shannon_diversity", "Species_richness", "Coverage"))
combined_bsi_mlst_incursion$sample_size <- as.factor(combined_bsi_mlst_incursion$sample_size)

# save
write.csv(combined_bsi_mlst_incursion, "rarefaction/combined_bsi_mlst_incursion.csv", row.names = FALSE)

# REPROCESS: create ordered sample_size factor and per-genus colour palettes
# assume df is your data and genus_pal_list is defined as before
df <- combined_bsi_mlst_incursion
df <- df |> mutate(sample_size = as.integer(as.character(sample_size)))
sample_levels <- sort(unique(df$sample_size))

# base palettes (adjust as you like)
genus_pal_list <- list(
  "Escherichia" = c("#c7f0c7", "#8fdc8f", "#57c257", "#2f962f", "#1b5e1b"),
  "Klebsiella"  = c("#ffe6cc", "#ffc38f", "#ff9a3d", "#e86d00", "#b04a00")
)

make_pal <- function(base_cols, n) {
  colorRampPalette(base_cols)(n)
}

# build lookup grid
lookup <- expand.grid(
  Genus = sort(unique(as.character(df$Genus))),
  sample_size = sample_levels,
  stringsAsFactors = FALSE) |> 
  arrange(Genus, sample_size)

# create colors per genus
lookup2 <- lookup |>
  group_by(Genus) |>
  summarise(n_needed = n(), .groups = "drop") |>
  rowwise() |>
  mutate(
    base_cols = {
      # safe grab of palette; fallback to greys if not defined
      g <- as.character(Genus)
      pal <- genus_pal_list[[g]]
      if (is.null(pal)) pal <- c("#bdbdbd", "#252525")
      list(pal)
    },
    pal_vec = list(make_pal(base_cols[[1]], n_needed))
  ) |>
  select(Genus, pal_vec) |>
  unnest_longer(pal_vec) |>      # expands one row per palette color
  group_by(Genus) |>
  mutate(sample_size = sample_levels) |>  # ensure sample_size ordering inside group
  ungroup() |>
  rename(col = pal_vec) |>
  select(Genus, sample_size, col)

# if you prefer fill slightly lighter, create fill by alpha
#library(scales)
lookup2 <- lookup2 |>
  mutate(fill = scales::alpha(col, 0.5))

# sanity check: lookup2 should have one row per (Genus, sample_size)
stopifnot(nrow(lookup2) == length(unique(df$Genus)) * length(sample_levels))

#make lookup2 sample sizes as factors
lookup2$sample_size <- factor(lookup2$sample_size, levels = sample_levels)

# merge back to main df
df2 <- df |>
  mutate(sample_size = factor(sample_size, levels = sample_levels)) |>
  left_join(lookup2, by = c("Genus", "sample_size"))

# create key for mapping in scale_manual
df2 <- df2 |> mutate(key = paste(Genus, as.character(sample_size), sep = "_"))

cols_vec <- setNames(df2 |> distinct(key, col) |> pull(col),
                     df2 |> distinct(key, col) |> pull(key))
fills_vec <- setNames(df2 |> distinct(key, fill) |> pull(fill),
                      df2 |> distinct(key, fill) |> pull(key))
# Optional: nicer labels for legend (just show numeric sample_size)
#legend_labels <- setNames(as.character(lookup$sample_size), lookup$key)

# PLOT 
p <- ggplot(data = df2, 
            aes(x = new_species_frequency, y = est, 
                colour = key, fill = key, group = key)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.35, colour = NA) +
  geom_line(linewidth = 1) +
  
  # Use our manual scales (keys map to colours)
  scale_colour_manual(name = "Sample size",
                      values = cols_vec,
                      # labels = legend_labels,
                      breaks = names(cols_vec)) +
  scale_fill_manual(name = "Sample size",
                    values = fills_vec,
                    #labels = legend_labels,
                    breaks = names(fills_vec),
                    guide = "none") + # hide fill legend to avoid duplicate entries
  scale_x_log10() +
  
  facet_wrap(Genus ~ Metric, scales = "free", ncol = 4) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11),
    legend.position = "right",
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(x = "New species frequency", y = "Estimate")

print(p)
# save
ggsave("rarefaction/combined_bsi_mlst_incursion.png", plot = p, width = 15, height = 10.5, units = "in", dpi = 300)

#str(cobmined_bsi_mlst_incursion$sample_size)
##unique(cobmined_bsi_mlst_incursion$sample_size)
#sample_levels
#cobmined_bsi_mlst_incursion$sample_size <- factor(cobmined_bsi_mlst_incursion$sample_size, levels  = rev(sample_levels))

# plot
p <- ggplot(data = combined_bsi_mlst_incursion, aes(x = new_species_frequency, y = est, colour = sample_size, fill = sample_size)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3, colour = NA) +
  scale_fill_viridis_d(option = "mako") +
  scale_colour_viridis_d(option = "mako") +
  scale_x_log10() +
  facet_wrap(Genus ~ Metric, scales= "free", ncol = 4) +
  labs(x = "Frequency of novel species" ,
       y = "Diversity metric value",
        ) +
  theme_minimal(base_size = 16)

print(p)
# save
ggsave("rarefaction/combined_bsi_mlst_incursion_logscale.png", plot = p, width = 15, height = 10, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Coverage Estimators ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## load E. coli Bayesian estimation data
ecoli_mass_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
ecoli_mass_df_long <- ecoli_mass_df |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(Est = median_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, Est, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "Bayesian 90%",
                               Estimator == "min_sample_95" ~ "Bayesian 95%",
                               Estimator == "min_sample_99" ~ "Bayesian 99%", 
                               TRUE ~ Estimator))

# load rarefaction data 
ecoli_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_mlst_rarefaction_combined.csv")
ecoli_bsi_mlst_rarefaction_sc <- ecoli_bsi_mlst_rarefaction_combined |>
  filter(Metric == "Sample coverage" & Replacement == "No replacement") |>
  rename(Est = est, q2.5 = lcl, q97.5 = ucl) |>
  mutate(Estimator = "rarefaction") |>
  dplyr::select(Estimator, sample_size, Est, q2.5, q97.5)

# load E. coli preseqR coverage data
ecoli_bsi_preseq_sam_cov <- read.csv("rarefaction/ecoli_bsi_mlst_preseqR_coverage.csv")
ecoli_bsi_preseq_sam_cov_filt <- ecoli_bsi_preseq_sam_cov |>  
  filter(r == 1) |>
  rename(sample_size = actual_sample_size, Est = est, q2.5 = lower, q97.5 = upper) |>
  mutate(Estimator = "preseqR") |>
  dplyr::select(Estimator, sample_size, Est, q2.5, q97.5)
#View(ecoli_bsi_preseq_sam_cov_filt)

#iNEXT data
ecoli_bsi_iNEXT_size_based_total <- read.csv("rarefaction/ecoli_bsi_iNEXT_size_based_total.csv")
ecoli_bsi_iNEXT_size_based_total <- ecoli_bsi_iNEXT_size_based_total |>
  dplyr::filter(Assemblage == "Total" & Order.q == 0 ) |>
  rename( sample_size = m, Est = SC, q2.5 = SC.LCL, q97.5 = SC.UCL) |>
  mutate(Estimator = "iNEXT") |>
  dplyr::select(Estimator, sample_size, Est, q2.5, q97.5)
#View(ecoli_bsi_iNEXT_size_based_total)

# bind all rows together
ecoli_bsi_combined_sample_cov_estimators <- rbind(ecoli_mass_df_long,
                                                  ecoli_bsi_mlst_rarefaction_sc,
                                                  ecoli_bsi_preseq_sam_cov_filt,
                                                  ecoli_bsi_iNEXT_size_based_total)

# set order of Estimators:
ecoli_bsi_combined_sample_cov_estimators$Estimator <- factor(ecoli_bsi_combined_sample_cov_estimators$Estimator, 
                                                             levels = c("Bayesian 99%",
                                                                        "Bayesian 95%",
                                                                        "Bayesian 90%",
                                                                        "rarefaction",
                                                                        "iNEXT",
                                                                        "preseqR"))
#add colour key
sam_cov_colours <- c(
  `Bayesian 99%` = "seagreen4",
  `Bayesian 95%` = "seagreen3",
  `Bayesian 90%` = "seagreen1",
  rarefaction = "black", 
  iNEXT = "grey40",
  preseqR = "grey70"
)

# define max sample size
max_actual_sample_coverage <- ecoli_bsi_combined_sample_cov_estimators |>
  filter(Estimator == "rarefaction") |>
  filter(sample_size == max(sample_size)) |>
  pull(Est)

# plot
sample_coverage_comparison_plot <- ggplot(ecoli_bsi_combined_sample_cov_estimators, aes(y = Est, x = sample_size, colour = Estimator, fill = Estimator)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
  # add line for 80%
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  # add solid line for actual collected sample
  #geom_hline(yintercept = max_actual_sample_coverage, linetype = "solid") +
  #annotate("label", y = 0.72, x = 10000, label = "80% sample coverage", angle = 0, vjust = -0.1, hjust = 0.1, size = 4, colour = "black", fill = "white") +
  #annotate("label", y = 0.92, x = 10000, label = "NEKSUS sample collection", angle = 0, vjust = -0.1, hjust = 0.5,  size = 4, colour = "black", fill = "white") +
  scale_fill_manual(name = "Estimator", values = sam_cov_colours) +
  scale_colour_manual(name = "Estimator", values = sam_cov_colours) +
  scale_y_continuous() +
  scale_x_log10(limits = c(10, 10000),
                breaks = c(10, 100, 1000, 10000)) +
  labs(y = "Coverage",
       x = "Sample size",
       #title = "Comparison of sample coverage estimators for E. coli BSI MLSTs",
  ) +
  theme_minimal() +
  guides(
    colour = guide_legend(ncol = 2),
    fill = guide_legend(ncol = 2)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    plot.margin = margin(5.5, 20, 5.5, 5.5)  # space for vertical labels
  )
sample_coverage_comparison_plot
# save
ggsave("rarefaction/ecoli_bsi_mlst_sample_coverage_comparison_plot.png", plot = sample_coverage_comparison_plot, width = 6, height = 6, units = "in", dpi = 300)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## load Klebsiella Bayesian estimation data
kleb_mass_df <- read.csv("model_results/kleb_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
kleb_mass_df_long <- kleb_mass_df |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(Est = median_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, Est, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "Bayesian 90%",
                               Estimator == "min_sample_95" ~ "Bayesian 95%",
                               Estimator == "min_sample_99" ~ "Bayesian 99%", 
                               TRUE ~ Estimator))

# load rarefaction data 
kleb_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_mlst_rarefaction_combined.csv")
kleb_bsi_mlst_rarefaction_sc <- kleb_bsi_mlst_rarefaction_combined |>
  filter(Metric == "Sample coverage" & Replacement == "No replacement") |>
  rename(Est = est, q2.5 = lcl, q97.5 = ucl) |>
  mutate(Estimator = "rarefaction") |>
  dplyr::select(Estimator, sample_size, Est, q2.5, q97.5)

# load Klebsiella preseqR coverage data
kleb_bsi_preseq_sam_cov <- read.csv("rarefaction/kleb_bsi_mlst_preseqR_sample_cov.csv")
kleb_bsi_preseq_sam_cov_filt <- kleb_bsi_preseq_sam_cov |>  
  filter(r == 1) |>
  rename(sample_size = actual_sample_size, Est = est, q2.5 = lower, q97.5 = upper) |>
  mutate(Estimator = "preseqR") |>
  dplyr::select(Estimator, sample_size, Est, q2.5, q97.5)
#View(kleb_bsi_preseq_sam_cov_filt)

#iNEXT data
kleb_bsi_iNEXT_size_based_total <- read.csv("rarefaction/kleb_bsi_iNEXT_size_based_total.csv")
kleb_bsi_iNEXT_size_based_total <- kleb_bsi_iNEXT_size_based_total |>
  dplyr::filter(Assemblage == "Total" & Order.q == 0 ) |>
  rename( sample_size = m, Est = SC, q2.5 = SC.LCL, q97.5 = SC.UCL) |>
  mutate(Estimator = "iNEXT") |>
  dplyr::select(Estimator, sample_size, Est, q2.5, q97.5)
#View(kleb_bsi_iNEXT_size_based_total)

# bind all rows together
kleb_bsi_combined_sample_cov_estimators <- rbind(kleb_mass_df_long,
                                                  kleb_bsi_mlst_rarefaction_sc,
                                                  kleb_bsi_preseq_sam_cov_filt,
                                                  kleb_bsi_iNEXT_size_based_total)

# set order of Estimators:
kleb_bsi_combined_sample_cov_estimators$Estimator <- factor(kleb_bsi_combined_sample_cov_estimators$Estimator, 
                                                             levels = c("Bayesian 99%",
                                                                        "Bayesian 95%",
                                                                        "Bayesian 90%",
                                                                        "rarefaction",
                                                                        "iNEXT",
                                                                        "preseqR"))
sam_cov_colours <- c(
  `Bayesian 99%` = "darkorange4",
  `Bayesian 95%` = "darkorange3",
  `Bayesian 90%` = "darkorange1",
  rarefaction = "black", 
  iNEXT = "grey40",
  preseqR = "grey70"
)

# define max sample size
max_actual_sample_coverage <- kleb_bsi_combined_sample_cov_estimators |>
  filter(Estimator == "rarefaction") |>
  filter(sample_size == max(sample_size)) |>
  pull(Est)

# plot
sample_coverage_comparison_plot <- ggplot(kleb_bsi_combined_sample_cov_estimators, aes(y = Est, x = sample_size, colour = Estimator, fill = Estimator)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
  # add line for 80%
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  # add solid line for actual collected sample
  #geom_hline(yintercept = max_actual_sample_coverage, linetype = "solid") +
  #annotate("label", y = 0.72, x = 10000, label = "80% sample coverage", angle = 0, vjust = -0.1, hjust = 0.1, size = 4, colour = "black", fill = "white") +
  #annotate("label", y = 0.92, x = 10000, label = "NEKSUS sample collection", angle = 0, vjust = -0.1, hjust = 0.5,  size = 4, colour = "black", fill = "white") +
  scale_fill_manual(name = "Estimator", values = sam_cov_colours) +
  scale_colour_manual(name = "Estimator", values = sam_cov_colours) +
  scale_y_continuous() +
  scale_x_log10(limits = c(10, 10000),
                breaks = c(10, 100, 1000, 10000)) +
  labs(y = "Coverage",
       x = "Sample size",
       #title = "Comparison of sample coverage estimators for E. coli BSI MLSTs",
  ) +
  theme_minimal() +
  guides(
    colour = guide_legend(ncol = 2),
    fill = guide_legend(ncol = 2)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    plot.margin = margin(5.5, 20, 5.5, 5.5)  # space for vertical labels
  )
sample_coverage_comparison_plot
# save
ggsave("rarefaction/kleb_bsi_mlst_sample_coverage_comparison_plot.png", plot = sample_coverage_comparison_plot, width = 6, height = 6, units = "in", dpi = 300)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# get summary table of median, and 95% CIS fo sample sizes needed to detect samples at minimum frequency x, and how this translates to population covered.
#View(ecoli_bsi_combined_sample_cov_estimators)
ecoli_bsi_sample_coverage_summary_table <- ecoli_bsi_combined_sample_cov_estimators |>
  group_by(Estimator)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# helper: for a single estimator df and a threshold, find mins for a given column
min_sample_at_or_above <- function(df, colname, thr) {
  # return NA if no rows meet the condition
  res <- df |>
    filter(!is.na(.data[[colname]])) |>
    filter(.data[[colname]] >= thr) |>
    summarize(min_ss = if (n() == 0) NA_real_ else min(sample_size, na.rm = TRUE)) |>
    pull(min_ss)
  if (length(res) == 0) NA_real_ else res
}

# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_sample_coverage_summary_table <- ecoli_bsi_combined_sample_cov_estimators |>
  group_by(Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      est_ss  <- min_sample_at_or_above(df, "Est", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(est_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(est_ss)) "NA" else formatC(est_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             est_ss = est_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Estimator, threshold_label, cell, est_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, est_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_sample_coverage_summary_table_cells_only <- ecoli_bsi_sample_coverage_summary_table |>
  select(Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_sample_coverage_summary_table_cells_only)
#View(ecoli_bsi_sample_coverage_summary_table_cells_only)

write.csv(ecoli_bsi_sample_coverage_summary_table_cells_only, "rarefaction/ecoli_bsi_sample_coverage_summary_table_cells_only.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# get summary table of median, and 95% CIS fo sample sizes needed to detect samples at minimum frequency x, and how this translates to population covered.
#View(kleb_bsi_combined_sample_cov_estimators)
kleb_bsi_sample_coverage_summary_table <- kleb_bsi_combined_sample_cov_estimators |>
  group_by(Estimator)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# helper: for a single estimator df and a threshold, find mins for a given column
min_sample_at_or_above <- function(df, colname, thr) {
  # return NA if no rows meet the condition
  res <- df |>
    filter(!is.na(.data[[colname]])) |>
    filter(.data[[colname]] >= thr) |>
    summarize(min_ss = if (n() == 0) NA_real_ else min(sample_size, na.rm = TRUE)) |>
    pull(min_ss)
  if (length(res) == 0) NA_real_ else res
}

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_sample_coverage_summary_table <- kleb_bsi_combined_sample_cov_estimators |>
  group_by(Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      est_ss  <- min_sample_at_or_above(df, "Est", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(est_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(est_ss)) "NA" else formatC(est_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             est_ss = est_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Estimator, threshold_label, cell, est_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, est_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_sample_coverage_summary_table_cells_only <- kleb_bsi_sample_coverage_summary_table |>
  select(Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_sample_coverage_summary_table_cells_only)
#View(kleb_bsi_sample_coverage_summary_table_cells_only)

write.csv(kleb_bsi_sample_coverage_summary_table_cells_only, "rarefaction/kleb_bsi_sample_coverage_summary_table_cells_only.csv", row.names = FALSE)