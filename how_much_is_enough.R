#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# BAYESIAN BOOTSTRAPPING AND POWER CALCULATION ####
# FOR ESTIMATING SAMPLE SIZE REQUIRED FOR GENOMIC SURVEILLANCE 
# INFORMED BY DATA FROM THE NEKSUS STUDY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. SETUP ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Install and load packages
install.packages("tidyverse")
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
install.packages("ggplot2")
install.packages("purrr")
install.packages("bayesboot")
install.packages("grafify")
install.packages("gt")
install.packages("gtsummary")
install.packages("viridis")
install.packages("forcats")
install.packages("scales")
install.packages("gtools")
install.packages("posterior")
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
install.packages("loo")


library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(purrr)
library(grafify)
library(tibble)
library(stringr)
library(patchwork)
library(gtsummary)
library(gt)
library(bayesboot)
library(viridis)
library(forcats)
library(scales)
library(gtools)
library(posterior)
library(cmdstanr)
library(loo)


# set workign directory
setwd("~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/")

# read-in cleaned/deduplicated data (under embargo until April 2027)

# load E.coli and Kleb MLST and fastbaps cluster data
ecoli_bsi_samples_metadata <- read.csv("rarefaction/ecoli_bsi_samples_metadata.csv")
kleb_bsi_samples_metadata <- read.csv("rarefaction/kleb_bsi_samples_metadata.csv")
#colnames(ecoli_bsi_samples_metadata)
#colnames(kleb_bsi_samples_metadata)


# load amrfinder and plasmid data 
## this df already had plasmid "community_subcommunity" (pling) and "rep_types_whole_plasmid" (mob-suite)
ecoli_bsi_amrfinder_metadata <- read.csv("rarefaction/ecoli_bsi_amrfinder_metadata.csv")
kleb_bsi_amrfinder_metadata <- read.csv("rarefaction/kleb_bsi_amrfinder_metadata.csv")
#colnames(ecoli_bsi_amrfinder_metadata)
#colnames(kleb_bsi_amrfinder_metadata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Explore data: ####
# Frequency distribution histograms of MLSTs, fastbaps clusters, plasmids and ARGs 
ecoli_samples_df <- ecoli_bsi_samples_metadata
ecoli_amr_df     <- ecoli_bsi_amrfinder_metadata
klebsiella_samples_df <- kleb_bsi_samples_metadata   
klebsiella_amr_df     <- kleb_bsi_amrfinder_metadata
colnames(ecoli_bsi_samples_metadata)
# Column names 
mlst_col_ecoli <- "escherichia__mlst_achtman__ST"
mlst_col_kleb  <- "klebsiella_mlst_ST"   
fb_cols <- c("Level.1", "Level.2", "Level.3")
arg_symbol_col <- "Element.symbol"
arg_type_col   <- "Type"   # filter Type == "AMR"
subcommunity_col <- "community_subcommunity"
contig_id_col    <- "Contig.id"
sample_id_cols   <- c("sample","run","Contig.id")  # used to deduplicate contig entries per sample

# Colors
genus_cols <- c(Escherichia = "seagreen3", Klebsiella = "darkorange")

n_bins <- 30

# helper functions 
# given a samples df and the column name for MLST, compute counts per MLST
get_mlst_counts <- function(samples_df, mlst_col, genus_name) {
  samples_df %>%
    filter(!is.na(.data[[mlst_col]])) %>%
    group_by(!!sym(mlst_col)) %>%
    summarise(count = n(), .groups = "drop") %>%
    transmute(metric = "MLSTs", genus = genus_name, value = count)
}

# given samples df and fastbaps level column name compute counts
get_fb_counts <- function(samples_df, level_col, genus_name, level_label) {
  samples_df %>%
    filter(!is.na(.data[[level_col]])) %>%
    group_by(!!sym(level_col)) %>%
    summarise(count = n(), .groups = "drop") %>%
    transmute(metric = paste0("FastBAPS clusters ", level_label), genus = genus_name, value = count)
}

# ARG counts (count occurrences of each gene across samples)
get_arg_counts <- function(amr_df, genus_name) {
  amr_df %>%
    filter(!is.na(.data[[arg_symbol_col]]), .data[[arg_type_col]] == "AMR") %>%
    group_by(!!sym(arg_symbol_col)) %>%
    summarise(count = n(), .groups = "drop") %>%
    transmute(metric = "AMR genes", genus = genus_name, value = count)
}

# plasmid subcommunity counts (deduplicate contigs per sample then count number of contigs per subcommunity)
get_plasmid_subcommunity_counts <- function(amr_df, genus_name) {
  amr_df %>%
    filter(!is.na(.data[[subcommunity_col]])) %>%
    # dedupe by sample/run/contig so contig not double counted
    group_by(!!!syms(sample_id_cols)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(!!sym(subcommunity_col)) %>%
    summarise(count = n(), .groups = "drop") %>%
    transmute(metric = "Plasmid subcommunities", genus = genus_name, value = count)
}

# call functions
# Escherichia
ecoli_mlst   <- get_mlst_counts(ecoli_samples_df, mlst_col_ecoli, "Escherichia")
ecoli_fb_l1  <- get_fb_counts(ecoli_samples_df, fb_cols[1], "Escherichia", "L1")
ecoli_fb_l2  <- get_fb_counts(ecoli_samples_df, fb_cols[2], "Escherichia", "L2")
ecoli_fb_l3  <- get_fb_counts(ecoli_samples_df, fb_cols[3], "Escherichia", "L3")
ecoli_args   <- get_arg_counts(ecoli_amr_df, "Escherichia")
ecoli_plasm  <- get_plasmid_subcommunity_counts(ecoli_amr_df, "Escherichia")

# Klebsiella
kleb_mlst   <- get_mlst_counts(klebsiella_samples_df, mlst_col_kleb, "Klebsiella")
kleb_fb_l1  <- get_fb_counts(klebsiella_samples_df, fb_cols[1], "Klebsiella", "L1")
kleb_fb_l2  <- get_fb_counts(klebsiella_samples_df, fb_cols[2], "Klebsiella", "L2")
kleb_fb_l3  <- get_fb_counts(klebsiella_samples_df, fb_cols[3], "Klebsiella", "L3")
kleb_args   <- get_arg_counts(klebsiella_amr_df, "Klebsiella")
kleb_plasm  <- get_plasmid_subcommunity_counts(klebsiella_amr_df, "Klebsiella")

# bind everything
all_counts <- bind_rows(
  ecoli_mlst, ecoli_fb_l1, ecoli_fb_l2, ecoli_fb_l3, ecoli_args, ecoli_plasm,
  kleb_mlst,  kleb_fb_l1,  kleb_fb_l2,  kleb_fb_l3,  kleb_args,  kleb_plasm
)
#View(all_counts)


# label factor order for nice faceting
metric_levels <- c("MLSTs", "FastBAPS clusters L1", "FastBAPS clusters L2", "FastBAPS clusters L3",  "Plasmid subcommunities", "AMR genes")
all_counts$metric <- factor(all_counts$metric, levels = metric_levels)

# plot
# 1) Counts histogram (counts on x, number of clusters/STs on y) - only fastbaps L3
counts_plot <- ggplot(all_counts[!all_counts$metric %in% c("FastBAPS clusters L1", "FastBAPS clusters L2"), ], aes(x = value, fill = genus, colour = NULL)) +
  geom_histogram(bins = n_bins, position = "identity", alpha = 0.6, closed = "right") +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  scale_fill_manual(values = genus_cols, name = "Genus") +
  labs(x = "Number of isolates per unit", y = "Number of unique units") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")
counts_plot
ggsave("panel_counts_histograms_four_panel.png", counts_plot, width = 6, height = 4.5, dpi = 300)

# 2) Frequency histogram (counts on x, number of clusters/STs on y)
freq_plot_by_genus <- ggplot(all_counts[!all_counts$metric %in% c("FastBAPS clusters L1", "FastBAPS clusters L2"), ], aes(x = value, fill = genus, group = genus)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins = n_bins,
                 position = "identity", alpha = 0.6, colour = NA) +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(Escherichia = "seagreen3", Klebsiella = "darkorange")) +
  labs(x = "Size (number of isolates per unit)", y = "Proportion of units (within genus)", fill = "Genus") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))
freq_plot_by_genus
ggsave("panel_freq_histograms_four_panel.png", freq_plot_by_genus, width = 6, height = 4.5, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. BAYESIAN BOOTSTRAPPING FOR MLST AND FASTBAPS CLUSTERS - OVERALL NATIONAL COUNTS ####
# estimate posterior frequency distributions of MLSTs and fastbaps clusters in the E. coli and Klebsiella BSI population in England, using a Bayesian bootstrap approach with a Dirichlet prior. This will allow us to estimate the total frequency mass of rare types (e.g. those with frequency <0.1%) and thus inform sample size calculations for genomic surveillance.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 1.1 MLSTs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Write bayesboot function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTION: bayesboot_mlst
# Inputs:
#   df: data.frame with columns 'mlst_profile' and 'count'
#   alpha_named: named numeric vector of alpha pseudo-counts for each observed MLST
#                (names must match df$mlst_profile). May be 0 for some types.
#   alpha_novel: numeric pseudo-count for the grouped "novel" category (can be non-integer)
#   B: number of posterior draws
#   use_exact_dirichlet: if TRUE, use direct Dirichlet sampling (supports non-integer alphas)
#                         if FALSE, attempt to use bayesboot by adding integer pseudo-observations
# Returns:
#   A list with:
#     draws: matrix B x (K+1) of posterior frequency draws (columns named by mlst + "NOVEL")
#     summary_df: data.frame summarizing mean, sd, quantiles per category
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# function #
bayesboot_mlst <- function(df,
                           feature_col = "mlst_profile",
                           alpha_named = NULL,
                           alpha_novel = 1,
                           B = 1000,
                           use_exact_dirichlet = TRUE,
                           seed = 2026) {
  # dependencies: dplyr, tidyr, bayesboot
  stopifnot(requireNamespace("dplyr", quietly = TRUE),
            requireNamespace("tidyr", quietly = TRUE))
  set.seed(seed)
  
  # basic checks
  stopifnot(is.data.frame(df))
  stopifnot(is.character(feature_col) && length(feature_col) == 1)
  stopifnot(all(c(feature_col, "count") %in% colnames(df)))
  
  # collapse to unique feature rows (group by suppied column name)
  df <- df |> 
    dplyr::group_by(.data[[feature_col]]) |> 
    dplyr::summarise(count = sum(.data[["count"]]), .groups = "drop") 
  
  # K names and counts
  K <- nrow(df)
  names_counts <- as.character(df[[feature_col]])
  n_k <- df$count
  N <- sum(n_k)
  
  # default alphas if not provided = 1 for each observed
  if(is.null(alpha_named)) {
    alpha_named <- rep(1, K)
    names(alpha_named) <- names_counts
  } else {
    if(!all(names_counts %in% names(alpha_named))) {
      stop("alpha_named must have names for every observed MLST in df$mlst_profile")
    }
    # re-order to match df
    alpha_named <- alpha_named[names_counts]
  }
  
  # If using exact Dirichlet sampling (recommended), draw directly with rgamma:
  if(use_exact_dirichlet) {
    Kplus <- K + 1
    colnames_out <- c(names_counts, "NOVEL")
    draws <- matrix(NA_real_, nrow = B, ncol = Kplus)
    colnames(draws) <- colnames_out
    
    for(b in seq_len(B)) {
      # sample gamma variables with shape = n_k + alpha_k
      shapes_obs <- n_k + as.numeric(alpha_named) # add one to everyhitng?? why
      g_obs <- rgamma(K, shape = shapes_obs, rate = 1)
      # novel mass
      g_novel <- rgamma(1, shape = alpha_novel, rate = 1)
      G <- c(g_obs, g_novel)
      p <- G / sum(G)
      draws[b, ] <- p
    }
    
    summary_df <- as.data.frame(t(apply(draws, 2, function(x) {
      c(mean = mean(x), sd = sd(x), q2.5 = quantile(x, 0.025), q97.5 = quantile(x, 0.975))
    })))
    summary_df[[feature_col]] <- rownames(summary_df)
    rownames(summary_df) <- NULL
    
    return(list(draws = draws, summary_df = summary_df))
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Otherwise: use bayesboot by creating pseudo-observations to represent alpha.
  # NOTE: bayesboot uses Dirichlet(1,...,1) over observations; replicating a
  # MLST alpha times approximates adding alpha pseudo-counts BUT this requires
  # integer alphas. We'll round down then use remainder as fractional (approx).
  # This path is provided for users who *must* use bayesboot object output.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Check integrality of alpha_named and alpha_novel
  if(any(alpha_named < 0) || alpha_novel < 0) stop("alpha values must be non-negative")
  if(any(alpha_named != floor(alpha_named)) || alpha_novel != floor(alpha_novel)) {
    stop("When use_exact_dirichlet = FALSE, alpha_named and alpha_novel must be integer-valued (or you can set use_exact_dirichlet = TRUE).")
  }
  
  # expand observed isolates into a vector of length N
  obs_vector <- rep(df[[feature_col]], times = df$count)
  
  # construct pseudo-observations for alpha
  pseudo_vector <- unlist(mapply(function(name, a) {
    if(a>0) rep(name, a) else character(0)
  }, name = names_counts, a = alpha_named, SIMPLIFY = FALSE), use.names = FALSE)
  
  if(alpha_novel > 0) pseudo_vector <- c(pseudo_vector, rep("NOVEL", alpha_novel))
  
  combined_vector <- c(obs_vector, pseudo_vector)
  combined_df <- data.frame(feature = combined_vector, stringsAsFactors = FALSE)
  
  # Define statistic function that takes data and weights (use.weights = TRUE)
  # and returns vector of weighted proportions for each category (observed + NOVEL)
  stat_fn <- function(dat, weights) {
    # dat is the data frame passed to bayesboot; we expect a column 'mlst'
    wdf <- data.frame(feature = dat$feature, w = weights)
    res <- wdf |> 
      dplyr::group_by(feature) |> 
      dplyr::summarise(m = sum(w), .groups = "drop") 
    # ensure we return columns in consistent order: observed names then NOVEL (if present)
    out_names <- c(names_counts, "NOVEL")
    res_wide <- res |> tidyr::pivot_wider(names_from = feature, values_from = m, values_fill = 0)
    # fill missing columns
    missing_cols <- setdiff(out_names, colnames(res_wide))
    for (mc in missing_cols) res_wide[[mc]] <- 0
    # return vector in order
    as.numeric(res_wide[out_names])
  }
  
  # run bayesboot
  if (!requireNamespace("bayesboot", quietly = TRUE)) {
    stop("Package 'bayesboot' is required for use_exact_dirichlet = FALSE path. Please install it or set use_exact_dirichlet = TRUE.")
  }
  bb <- bayesboot::bayesboot(combined_df, statistic = stat_fn, use.weights = TRUE, R = B)
  
  # bb is a data.frame with B rows and K+? columns (as returned)
  # Ensure columns match and convert to matrix
  draws <- as.matrix(bb)
  colnames(draws) <- colnames(bb)
  
  # summarise
  summary_mat <- as.data.frame(t(apply(draws, 2, function(x) {
    c(mean = mean(x), sd = sd(x), q2.5 = stats::quantile(x, 0.025), q97.5 = stats::quantile(x, 0.975))
  })))
  
  summary_df <- as.data.frame(summary_mat)
  summary_df[[feature_col]] <- rownames(summary_df)
  rownames(summary_df) <- NULL
  
  return(list(draws = draws, summary_df = summary_df, bayesboot_obj = bb))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTION: compute_mass_curve
# For each posterior draw, compute for a grid of thresholds f the total mass
# of categories with p_k >= f. Return a summary df with mean and 95% CI.
# Inputs:
#   draws: matrix B x (K+1) with columns for each MLST and "NOVEL"
#   f_grid: numeric vector of thresholds between 0 and 1 (default log-spaced)
# Returns:
#   data.frame with columns: f, mean_mass, sd_mass, q2.5, q97.5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
compute_mass_curve <- function(draws, f_grid = NULL) {
  B <- nrow(draws)
  if(is.null(f_grid)) {
    # denser near small f values (log-spaced)
    f_grid <- unique(c(0, round(10^seq(-5, 0, length.out = 200), 8)))
    f_grid <- sort(f_grid)
  }
  out <- map_dfr(f_grid, function(f) {
    masses <- rowSums(draws * (draws >= f))
    tibble(f = f,
           mean_mass = mean(masses),
           sd = sd(masses),
           q2.5 = quantile(masses, 0.025),
           q97.5 = quantile(masses, 0.975))
  })
  out
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.1a E. coli MLSTs #### 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare data - E. coli
ecoli_bsi_mlst_df <- ecoli_bsi_samples_metadata |>
  dplyr::rename(kleborate_mlst = escherichia__mlst_achtman__ST) |>
  dplyr::group_by(kleborate_mlst) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(kleborate_mlst, count) |>
  dplyr::arrange(count, kleborate_mlst)
#View(ecoli_bsi_mlst_df)
#table(ecoli_bsi_mlst_df$count)
#nrow(ecoli_bsi_mlst_df) # 263


# estimates for alpha novel
# prior for alpha_novel
# df has columns mlst_profile and count
df <- ecoli_bsi_mlst_df
N <- sum(df$count)
f1 <- sum(df$count == 1)  # 161       # number of singletons
q_hat <- f1 / N  # Good-Turing first-order - proportion of singletons

K <- length(df$count) # 263
alpha_named <- rep(1, K)                   # set uninformative priors
A_obs <- sum(alpha_named)         # sum of per-MLST prior pseudo-counts you plan to use # 278 # 1 per MLST profile
alpha_novel <- (q_hat / (1 - q_hat)) * (N + A_obs) # 213
alpha_novel_rounded <- round(alpha_novel)

#alternative alpha-novels
alpha_novel_null <-  0 # (no unseen mass)
alpha_0.5 <- 0.5
alpha_1 <- 1
alpha_novel_GT_anchored = (f1 / N ) * sum(alpha_named) # → Good–Turing–anchored prior  # 34
alpha_novel_f1 <- f1 #(Good–Turing) # 161
alpha_novel_2f1 = 2 * f1  #(more conservative) # 364


# Preferred / exact approach (supports non-integer alphas):
ecoli_bsi_mlst_bayesboot_exact <- bayesboot_mlst(ecoli_bsi_mlst_df,
                                                 feature_col = "kleborate_mlst",
                                                 alpha_named = NULL, alpha_novel = alpha_1, 
                                                 B = 10000, use_exact_dirichlet = TRUE)
#ecoli_bsi_mlst_bayesboot_non_exact <- bayesboot_mlst(ecoli_bsi_mlst_df, feature_col = "kleborate_mlst", alpha_named = NULL, alpha_novel = alpha_novel_rounded, B = 10000, use_exact_dirichlet = FALSE)

# quick look at per-category posterior means:
#View(ecoli_bsi_mlst_bayesboot_exact$summary_df)
#print(ecoli_bsi_mlst_bayesboot_non_exact$summary_df)
#print(res$draws)

# save
saveRDS(ecoli_bsi_mlst_bayesboot_exact, "rarefaction/ecoli_bsi_mlst_bayesboot_exact.rds")
#saveRDS(ecoli_bsi_mlst_bayesboot_non_exact, "rarefaction/ecoli_bsi_mlst_bayesboot_non_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#ecoli_bsi_mlst_bayesboot_exact <- readRDS("rarefaction/ecoli_bsi_mlst_bayesboot_exact.rds")
#ecoli_bsi_mlst_bayesboot_non_exact <- readRDS("rarefaction/ecoli_bsi_mlst_bayesboot_non_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#

# compute mass >= f curve
#f_grid <- c(0, seq(0.0001, 0.001, by = 0.0001), seq(0.002, 0.01, by = 0.001), seq(0.02, 1, by = 0.01), seq(0.02, 1, by = 0.01))
f_grid <- c(0, seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
ecoli_mass_df <- compute_mass_curve(ecoli_bsi_mlst_bayesboot_exact$draws, f_grid = f_grid)
ecoli_mass_df_non_exact <- compute_mass_curve(ecoli_bsi_mlst_bayesboot_non_exact$draws, f_grid = f_grid)
#View(ecoli_mass_df)
#View(ecoli_mass_df_non_exact)

# add sample sizes
ecoli_mass_df <- ecoli_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
         min_sample_95 = log(1-0.95)/log(1-f),
         min_sample_99 = log(1-0.99)/log(1-f),
         Genus = "Escherichia",
         Exact = "Exact") |>
  filter(f != 0)
# add sample sizes to non-exact
ecoli_mass_df_non_exact <- ecoli_mass_df_non_exact |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
         min_sample_95 = log(1-0.95)/log(1-f),
         min_sample_99 = log(1-0.99)/log(1-f),
         Genus = "Escherichia",
         Exact = "Non-exact") |>
  filter(f != 0)

# save  mass df
write.csv(ecoli_mass_df, "ecoli_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
write.csv(ecoli_mass_df_non_exact, "ecoli_bsi_kelborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#ecoli_mass_df <- read.csv("ecoli_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv")
#ecoli_mass_df_non_exact <- read.csv("ecoli_bsi_kelborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.1b Klebsiella MLSTs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare data - Kleb
kleb_bsi_mlst_df <- kleb_bsi_samples_metadata |>
  dplyr::rename(kleborate_mlst = klebsiella_mlst_ST) |>
  dplyr::group_by(kleborate_mlst) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(kleborate_mlst, count) |>
  dplyr::arrange(count, kleborate_mlst)
#View(kleb_bsi_mlst_df)
# check MLST count distribution
table(kleb_bsi_mlst_df$count)
nrow(kleb_bsi_mlst_df) # 297
# estimates for alpha novel
# prior for alpha_novel
# df has columns mlst_profile and count
df <- kleb_bsi_mlst_df
N <- sum(df$count) # 468
f1 <- sum(df$count == 1)  # 227        # number of singletons
q_hat <- f1 / N  #0.485 =  Good-Turing first-order - proportion of singletons

K <- length(df$count) #297 # number of unique MLSTs
alpha_named <- rep(1, K)                   # set uninformative priors
A_obs <- sum(alpha_named)    #297     # sum of per-MLST prior pseudo-counts you plan to use # 278 # 1 per MLST profile
alpha_novel <- (q_hat / (1 - q_hat)) * (N + A_obs) # 720.56
alpha_novel_rounded <- round(alpha_novel)

#alternative alpha-novels
alpha_novel_null <-  0 # (no unseen mass)
alpha_0.5 <- 0.5
alpha_1 <- 1
alpha_novel_GT_anchored = (f1 / N ) * sum(alpha_named) # 144.0577 → Good–Turing–anchored prior  
alpha_novel_f1 <- f1 #(Good–Turing) # 227 = count of singletons - approximation if small prior alpha
alpha_novel_2f1 = 2 * f1  #(more conservative) # 454


# Preferred / exact approach (supports non-integer alphas):
kleb_bsi_mlst_bayesboot_exact <- bayesboot_mlst(kleb_bsi_mlst_df,
                                                feature_col = "kleborate_mlst",
                                                alpha_named = NULL, alpha_novel = alpha_1, 
                                                B = 10000, use_exact_dirichlet = TRUE)
kleb_bsi_mlst_bayesboot_non_exact <- bayesboot_mlst(kleb_bsi_mlst_df,
                                                    feature_col = "kleborate_mlst",
                                                    alpha_named = NULL, alpha_novel = alpha_1, 
                                                    B = 10000, use_exact_dirichlet = FALSE)

print(kleb_bsi_mlst_bayesboot_exact$summary_df)
print(kleb_bsi_mlst_bayesboot_non_exact$summary_df)

# save
saveRDS(kleb_bsi_mlst_bayesboot_exact, "rarefaction/kleb_bsi_kleborate_mlst_bayesboot_exact.rds")
saveRDS(kleb_bsi_mlst_bayesboot_non_exact, "rarefaction/kleb_bsi_kleborate_mlst_bayesboot_non_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#kleb_bsi_mlst_bayesboot_exact <- readRDS("rarefaction/kleb_bsi_kleborate_mlst_bayesboot_exact.rds")
#kleb_bsi_mlst_bayesboot_non_exact <- readRDS("rarefaction/kleb_bsi_kleborate_mlst_bayesboot_non_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#

# compute mass >= f curve
#f_grid <- c(0, seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
kleb_mass_df <- compute_mass_curve(kleb_bsi_mlst_bayesboot_exact$draws, f_grid = f_grid)
kleb_mass_df_non_exact <- compute_mass_curve(kleb_bsi_mlst_bayesboot_non_exact$draws, f_grid = f_grid)

# add sample sizes
kleb_mass_df <- kleb_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Klebsiella",
                Exact = "Exact") |>
  filter(f != 0)

# add sample sizes to non-exact
kleb_mass_df_non_exact <- kleb_mass_df_non_exact |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Klebsiella",
                Exact = "Non-exact") |>
  filter(f != 0)

# save  mass df
write.csv(kleb_mass_df, "kleb_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
write.csv(kleb_mass_df_non_exact, "kleb_bsi_kelborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#kleb_mass_df <- read.csv("kleb_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv")
#kleb_mass_df_non_exact <- read.csv("kleb_bsi_kelborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.1c Combined summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# combine dfs
ecoli_mass_df <- read.csv("ecoli_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv")
ecoli_mass_df_non_exact <- read.csv("ecoli_bsi_kelborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df <- read.csv("kleb_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df_non_exact <- read.csv("kleb_bsi_kelborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv")

mass_df_mlst <- rbind(ecoli_mass_df, kleb_mass_df)
mass_df_mlst_non_exact <- rbind(ecoli_mass_df_non_exact, kleb_mass_df_non_exact)
#View(mass_df_mlst)

# save merged mass df
write.csv(mass_df_mlst, "combined_bsi_kleborate_mlst_bayesboot_cumulative_frequency_mass_df.csv", row.names = FALSE)
write.csv(mass_df_mlst_non_exact, "combined_bsi_kleborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data
#mass_df_mlst <- read.csv("combined_bsi_kleborate_mlst_bayesboot_cumulative_frequency_mass_df.csv")
#mass_df_mlst_non_exact <- read.csv("combined_bsi_kleborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# summary tables 
# transform to long
mass_df_mlst_long <- mass_df_mlst |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_mlst_sample_coverage_summary_table <- mass_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "mean_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_mlst_sample_coverage_summary_table <- combined_bsi_mlst_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_mlst_sample_coverage_summary_table)
# save
write.csv(combined_bsi_mlst_sample_coverage_summary_table, "rarefaction/combined_bsi_mlst_sample_coverage_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * 1.1d Combined Plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# histogram of MLST frequencies 
ecoli_bsi_mlst_df <- ecoli_bsi_mlst_df |>
  mutate(frequency = count/sum(count))

posterior_mlst_freq_hist <- ggplot(data = ecoli_bsi_mlst_bayesboot_exact$summary_df, aes(x = mean) ) +
  geom_histogram(binwidth = 0.002, alpha = 0.6) +
  geom_histogram(aes(x = `q2.5.2.5%`), binwidth = 0.002, alpha = 0.3) +
  geom_histogram(aes(x = `q97.5.97.5%`), binwidth = 0.002, alpha = 0.1)  +
  # plot actual data
  # plot lower and upper CIs
  geom_histogram(data = ecoli_bsi_mlst_df, aes(x = frequency), fill = "seagreen3", alpha = 0.5, binwidth = 0.002) +
  theme_minimal()
posterior_mlst_freq_hist
ggsave("rarefaction/ecoli_posterior_mlst_freq_hist.png", posterior_mlst_freq_hist, width = 6, height = 4, units = "in", dpi = 300)
#~~~~~~~~~~~~~~~~~#
# histogram of MLST frequencies
kleb_bsi_mlst_df <- kleb_bsi_mlst_df |>
  mutate(frequency = count/sum(count))

posterior_mlst_freq_hist <- ggplot(data = kleb_bsi_mlst_bayesboot_exact$summary_df, aes(x = mean) ) +
  geom_histogram(binwidth = 0.002, alpha = 0.6) +
  geom_histogram(aes(x = `q2.5.2.5%`), binwidth = 0.002, alpha = 0.3) +
  geom_histogram(aes(x = `q97.5.97.5%`), binwidth = 0.002, alpha = 0.1)  +
  # plot actual data
  # plot lower and upper CIs
  geom_histogram(data = kleb_bsi_mlst_df, aes(x = frequency), fill = "darkorange", alpha = 0.5, binwidth = 0.002) +
  theme_minimal()
posterior_mlst_freq_hist
ggsave("rarefaction/kleb_posterior_mlst_freq_hist.png", posterior_mlst_freq_hist, width = 6, height = 4, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define colours
genus_colours <- c("Escherichia" = "seagreen3", "Klebsiella" = "darkorange")
# Plot (cumulative mass of population at frequency >= f)
cumulative_mass_plot <- ggplot(mass_df_mlst, aes(x = f, y = mean_mass, colour =  Genus, fill = Genus)) +
  geom_line(data = mass_df_mlst, aes(x = f, y = mean_mass, colour =  Genus)) +
  geom_ribbon(data = mass_df_mlst , aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Genus", values = genus_colours) +
  scale_colour_manual(name = "Genus", values = genus_colours) +
  scale_x_log10( breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(mass_df_mlst$f), 1)) +
  labs(x = "MLST frequency (f) (log scale)",
       y = "Proportion of population belonging to MLSTs of frequency ≥ f"
  ) +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
cumulative_mass_plot
# save
ggsave("rarefaction/combined_bsi_kleborate_mlst_bayesboot_cumulative_mass_plot.png", plot = cumulative_mass_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Cumulative fraction with sample size plot 
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
bootbayes_sample_coverage_plot <- ggplot(mass_df_mlst, aes(x = mean_mass, colour = Genus, fill = Genus)) +
  geom_line(aes(y = min_sample_90)) +
  geom_line(aes(y = min_sample_95) ) +
  geom_line(aes(y = min_sample_99)) +
  geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.35, colour = NA) +
  geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  #facet_wrap(~ Genus, ncol = 1) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(0,5000)) +
  labs(x = "Sample coverage",
       y = "Sample size",
       # title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bootbayes_sample_coverage_plot
ggsave("rarefaction/combined_bsi_kleborate_mlst_bootbayes_sample_coverage_plot_95.png", plot = bootbayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 1.2 fastBAPS clusters (level 3) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.2a E. coli fastbaps clusters ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prep data
ecoli_bsi_fastbaps_L3 <- ecoli_bsi_samples_metadata |>
  dplyr::group_by(Level.3) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(Level.3, count) |>
  dplyr::arrange(count, Level.3)
#View(ecoli_bsi_fastbaps_L3)

# check count distribution
table(ecoli_bsi_fastbaps_L3$count)
nrow(ecoli_bsi_fastbaps_L3) # 161

# estimates for alpha novel
# prior for alpha_novel
df <- ecoli_bsi_fastbaps_L3
N <- sum(df$count)
f1 <- sum(df$count == 1)  #       # number of singletons
q_hat <- f1 / N  # Good-Turing first-order - proportion of singletons
alpha_1 <- 1

# Preferred / exact approach (supports non-integer alphas):
ecoli_bsi_fastbaps_L3_bayesboot_exact <- bayesboot_mlst(ecoli_bsi_fastbaps_L3,
                                                        feature_col = "Level.3",
                                                        alpha_named = NULL, alpha_novel = alpha_1, 
                                                        B = 10000, use_exact_dirichlet = TRUE)
print(ecoli_bsi_fastbaps_L3_bayesboot_exact$summary_df)

# save
saveRDS(ecoli_bsi_fastbaps_L3_bayesboot_exact, "rarefaction/ecoli_bsi_fastbaps_L3_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#ecoli_bsi_fastbaps_L3_bayesboot_exact <- readRDS("rarefaction/ecoli_bsi_fastbaps_L3_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#

# compute mass >= f curve
ecoli_mass_df <- compute_mass_curve(ecoli_bsi_fastbaps_L3_bayesboot_exact$draws, f_grid = f_grid)

# add sample sizes
ecoli_mass_df <- ecoli_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(ecoli_mass_df, "ecoli_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#ecoli_mass_df <- read.csv("ecoli_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.2b Klebsiella fastbaps clusters ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare data - Kleb
kleb_bsi_fastbaps_L3_df <- kleb_bsi_samples_metadata |>
  dplyr::group_by(Level.3) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(Level.3, count) |>
  dplyr::arrange(count, Level.3)
#View(kleb_bsi_fastbaps_L3_df)

# check count distribution
table(kleb_bsi_fastbaps_L3_df$count)
nrow(kleb_bsi_fastbaps_L3_df) # 25
# estimates for alpha novel, ; prior for alpha_novel
df <- kleb_bsi_fastbaps_L3_df
N <- sum(df$count) # 468
alpha_1 < -1

# Preferred / exact approach (supports non-integer alphas):
kleb_bsi_fastbaps_L3_bayesboot_exact <- bayesboot_mlst(kleb_bsi_fastbaps_L3_df,
                                                       feature_col = "Level.3",
                                                       alpha_named = NULL, alpha_novel = alpha_1, 
                                                       B = 10000, use_exact_dirichlet = TRUE)

print(kleb_bsi_fastbaps_L3_bayesboot_exact$summary_df)

# save
saveRDS(kleb_bsi_fastbaps_L3_bayesboot_exact, "rarefaction/kleb_bsi_fastbaps_L3_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#kleb_bsi_fastbaps_L3_bayesboot_exact <- readRDS("rarefaction/kleb_bsi_fastbaps_L3_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
kleb_mass_df <- compute_mass_curve(kleb_bsi_fastbaps_L3_bayesboot_exact$draws, f_grid = f_grid)
#View(kleb_mass_df)

# add sample sizes
kleb_mass_df <- kleb_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Klebsiella",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(kleb_mass_df, "kleb_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Reload saved data 
#kleb_mass_df <- read.csv("kleb_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.2c Combined summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#
ecoli_mass_df <- read.csv("ecoli_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df <- read.csv("kleb_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
mass_df_fastbaps_L3 <- rbind(ecoli_mass_df, kleb_mass_df)
#View(mass_df_fastbaps_L3)
# save merged mass df
write.csv(mass_df_fastbaps_L3, "combined_bsi_fastbaps_L3_bayesboot_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#mass_df_fastbaps_L3 <- read.csv("combined_bsi_fastbaps_L3_bayesboot_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# summary tables 
# transform to long
mass_df_fastbaps_L3_long <- mass_df_fastbaps_L3 |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_fastbaps_L3_sample_coverage_summary_table <- mass_df_fastbaps_L3_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "mean_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_fastbaps_L3_sample_coverage_summary_table <- combined_bsi_fastbaps_L3_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_fastbaps_L3_sample_coverage_summary_table)
# save
write.csv(combined_bsi_fastbaps_L3_sample_coverage_summary_table, "rarefaction/combined_bsi_fastbaps_L3_sample_coverage_summary_table.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.2d Combined plots  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot (cumulative mass of population at frequency >= f)
cumulative_mass_plot <- ggplot(mass_df_fastbaps_L3, aes(x = f, y = mean_mass, colour =  Genus, fill = Genus)) +
  geom_line(aes(x = f, y = mean_mass, colour =  Genus)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Genus", values = genus_colours) +
  scale_colour_manual(name = "Genus", values = genus_colours) +
  scale_x_log10( breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(mass_df$f), 1)) +
  labs(x = "FastBAPS cluster frequency (f) (log scale)",
       y = "Proportion of population belonging to cluster of frequency ≥ f"
  ) +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
cumulative_mass_plot
# save
ggsave("rarefaction/combined_bsi_fastbaps_L3_bayesboot_cumulative_mass_plot.png", plot = cumulative_mass_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Cumulative arg fraction with sample size plot
bootbayes_sample_coverage_plot <- ggplot(mass_df_fastbaps_L3, aes(x = mean_mass, colour = Genus, fill = Genus)) +
  geom_line(aes(y = min_sample_90)) +
  geom_line(aes(y = min_sample_95) ) +
  geom_line(aes(y = min_sample_99)) +
  geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.35, colour = NA) +
  geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  #facet_wrap(~ Genus, ncol = 1) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(0,5000)) +
  labs(x = "Sample coverage",
       y = "Sample size",
       # title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bootbayes_sample_coverage_plot
ggsave("rarefaction/combined_bsi_fastbaps_L3_bootbayes_sample_coverage_plot_95.png", plot = bootbayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 1.3 SENSITIVITY ANALYSIS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# for varying alpha_initial (counts between 0.001), and alpha_novel
# initial starting values are uninformative uniform prions for alpha_inital (all 1s), and varying starting values of alpha_novel
# exact and non-exact dirichlet
# Lightweight sensitivity analysis: only keep mass_df, not draws

# install packages if required
#install.packages("future.apply")
library(future.apply)   # optional; will be used only if use_parallel = TRUE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.3a Exact ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set parameters
B <- 1000                        # posterior draws per run (reduce for speed)
f_grid <- c(0, seq(0.0001, 0.001, by = 0.00001),
            seq(0.0011, 0.01, by = 0.0001),
            seq(0.011, 1, by = 0.001))

alpha_named_values <- c(0.01, 0.1, 1)
use_parallel <- FALSE            # set TRUE to use parallel workers (requires future plan)
n_workers <- 12                  # only used when use_parallel = TRUE

#  helper that runs one scenario but only returns mass_df
run_one_scenario_mass_only <- function(df_counts, genus_name, alpha_named_val, alpha_novel_val, B_draws = B, use_exact = TRUE) {
  df_counts <- df_counts |> group_by(mlst_profile) |> summarise(count = sum(count)) |> ungroup()
  K <- nrow(df_counts)
  N <- sum(df_counts$count)
  f1 <- sum(df_counts$count == 1)
  alpha_named_vec <- rep(alpha_named_val, K)
  names(alpha_named_vec) <- df_counts$mlst_profile
  
  # Run the Dirichlet posterior draws but do not return them:
  # We rely on bayesboot_mlst returning a list with $draws — we will immediately compute the mass curve
  res <- bayesboot_mlst(df_counts,
                        alpha_named = alpha_named_vec,
                        alpha_novel = alpha_novel_val,
                        B = B_draws,
                        use_exact_dirichlet = use_exact)
  
  # compute mass curve from draws
  mass_df <- compute_mass_curve(res$draws, f_grid = f_grid)
  
  # free memory: remove draws and res quickly
  rm(res)
  gc()
  
  # add sample size columns and metadata; filter f>0
  mass_df <- mass_df |>
    filter(f != 0) |>
    mutate(
      min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
      min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
      min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f))),
      Genus = genus_name,
      alpha_named = alpha_named_val,
      alpha_novel = alpha_novel_val
    )
  
  return(mass_df)
}

# build scenario table 
make_genus_info <- function(df_counts, genus_name) {
  df_counts <- df_counts |> group_by(mlst_profile) |> summarise(count = sum(count)) |> ungroup()
  N <- sum(df_counts$count)
  f1 <- sum(df_counts$count == 1)
  K <- nrow(df_counts)
  q_hat <- if (N > 0) f1 / N else 0
  
  alpha_novel_null <- 0
  alpha_novel_0.5 <- 0.5
  alpha_novel_1 <- 1
  alpha_novel_2 <- 2
  alpha_novel_5 <- 5
  alpha_novel_7 <- 7
  alpha_novel_10 <- 10
  alpha_novel_20 <- 20
  alpha_novel_50 <- 50
  alpha_novel_75 <- 75
  alpha_novel_100 <- 100
  alpha_novel_200 <- 200
  alpha_novel_500 <- 500
  alpha_novel_750 <- 750
  alpha_novel_1000 <- 1000
  alpha_novel_2000 <- 2000
  alpha_novel_GT_anchored_basic <- (f1 / N) * K
  alpha_novel_f1 <- f1
  alpha_novel_2f1 <- 2 * f1
  
  tibble(
    Genus = genus_name,
    N = N, f1 = f1, q_hat = q_hat, K = K,
    alpha_novel_candidates = list(c(alpha_novel_null, alpha_novel_0.5, alpha_novel_1, 
                                    alpha_novel_2, alpha_novel_5, alpha_novel_7, 
                                    alpha_novel_10, alpha_novel_20, alpha_novel_50, alpha_novel_75, 
                                    alpha_novel_100,  alpha_novel_200,  alpha_novel_500,  alpha_novel_750, 
                                    alpha_novel_1000,   alpha_novel_2000, 
                                    alpha_novel_GT_anchored_basic, alpha_novel_f1, alpha_novel_2f1))
  )
}

ecoli_info <- make_genus_info(ecoli_bsi_mlst_df, "Escherichia")
kleb_info  <- make_genus_info(kleb_bsi_mlst_df, "Klebsiella")

# expand parameter grid
scenarios <- bind_rows(ecoli_info, kleb_info) |>
  unnest(alpha_novel_candidates) |>
  rename(alpha_novel_candidate = alpha_novel_candidates) |>
  crossing(tibble(alpha_named = alpha_named_values)) |>
  mutate(
    A_obs = alpha_named * K,
    alpha_novel_anchored_full = (q_hat / pmax(1 - q_hat, 1e-12)) * (N + A_obs),
    alpha_novel_final = map2(alpha_novel_candidate, alpha_novel_anchored_full, ~ c(.x, .y))
  ) |>
  select(Genus, N, f1, K, alpha_named, alpha_novel_final) |>
  unnest(alpha_novel_final) |>
  mutate(alpha_novel = as.numeric(alpha_novel_final)) |>
  select(-alpha_novel_final)

# optional: unique & sort
scenarios <- scenarios |>
  distinct(Genus, alpha_named, alpha_novel) |>
  arrange(Genus, alpha_named, alpha_novel)


#run scenarios (serial or parallel)
if (use_parallel) {
  plan(multisession, workers = n_workers)
  combined_mass_list <- future_pmap(
    list(
      df = list(ecoli_bsi_mlst_df, kleb_bsi_mlst_df)[match(scenarios$Genus, c("Escherichia", "Klebsiella"))],
      genus = scenarios$Genus,
      alpha_named = scenarios$alpha_named,
      alpha_novel = scenarios$alpha_novel
    ),
    run_one_scenario_mass_only,
    B_draws = B,
    use_exact = TRUE,
    .options = future_options(seed = TRUE)
  )
  plan(sequential)  # back to sequential
} else {
  combined_mass_list <- pmap(
    list(
      df = list(ecoli_bsi_mlst_df, kleb_bsi_mlst_df)[match(scenarios$Genus, c("Escherichia", "Klebsiella"))],
      genus = scenarios$Genus,
      alpha_named = scenarios$alpha_named,
      alpha_novel = scenarios$alpha_novel
    ),
    run_one_scenario_mass_only,
    B_draws = B,
    use_exact = TRUE
  )
}

# combine into single data.frame; attach scenario ids/labels
combined_mass_dfs <- map2_dfr(combined_mass_list, seq_along(combined_mass_list), function(mdf, i) {
  sc <- scenarios[i, ]
  mdf |>
    mutate(
      alpha_named = sc$alpha_named,
      alpha_novel = sc$alpha_novel,
      scenario_id = i
    )
})

# create nice factor labels for plotting
combined_mass_dfs <- combined_mass_dfs |>
  mutate(
    alpha_named_f = factor(alpha_named, levels = alpha_named_values, labels = paste0("alpha_named=", alpha_named_values)),
    alpha_novel_label = paste0("alpha_novel=", signif(alpha_novel, 3)),
    alpha_novel_f = factor(alpha_novel_label, levels = unique(paste0("alpha_novel=", signif(sort(unique(alpha_novel)), 3))))
  )

# Save compact results
saveRDS(combined_mass_dfs, "rarefaction/sensitivity_combined_mass_dfs_mass_only.rds")

# plot
genus_colours <- c("Escherichia" = "seagreen3", "Klebsiella" = "darkorange")

View(combined_mass_dfs)
combined_mass_df_filtered <- combined_mass_dfs |>
  filter(alpha_novel %in% c(0, 0.5, 1,2,5,7,10,20,50,75,100, 200, 500, 750, 1000))

cumulative_mass_plot_sensitivity <-  ggplot(combined_mass_df_filtered, aes(x = f, y = mean_mass, colour = Genus, fill = Genus,
                                                                           group = interaction(Genus, alpha_named_f, alpha_novel_f, scenario_id))) +
  geom_line(aes(linetype = alpha_named_f), alpha = 0.95) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.12, colour = NA) +
  #facet_grid(rows = vars(alpha_named_f), cols = vars(alpha_novel_f), labeller = label_value) +
  #facet_wrap( ~ alpha_novel, ncol = 5) +
  facet_wrap( Genus ~ alpha_named) +
  scale_colour_manual(values = genus_colours) +
  scale_fill_manual(values = genus_colours) +
  scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = label_number()) +
  labs(
    x = "Frequency (f) (log scale)",
    y = "Posterior mean of total mass with p >= f",
    title = "Sensitivity analysis: cumulative mass curves across priors" ,
    #subtitle = "Rows = alpha_named; Columns = alpha_novel"
    #subtitle = "Faceted by alpha_novel"
    subtitle = "Faceted by alpha_named and Genus"
    
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(face = "bold")
  )

print(cumulative_mass_plot_sensitivity)
#ggsave("rarefaction/cumulative_mass_sensitivity_mass_only_facet_grid.png", cumulative_mass_plot_sensitivity, width = 30, height = 7, dpi = 300)

# for facet alpha_novel, to see effect of alpha_named on same graph
#ggsave("rarefaction/cumulative_mass_sensitivity_mass_only_facet_alpha_novel.png", cumulative_mass_plot_sensitivity, width = 16, height = 8, dpi = 300)

# for facet by alpha_named, to see effect of increasing alpha_novel on same graph
ggsave("rarefaction/cumulative_mass_sensitivity_mass_only_facet_alpha_named.png", cumulative_mass_plot_sensitivity, width = 12, height = 6, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1. 3b Non-exact ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# parameter settings
B <- 1000                        # posterior draws per run (reduce for speed)
f_grid <- c(0, seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))

alpha_named_values <- c(1)
use_parallel <- FALSE            # set TRUE to use parallel workers (requires future plan)
n_workers <- 12                  # only used when use_parallel = TRUE


# build scenario table
make_genus_info_exact <- function(df_counts, genus_name) {
  df_counts <- df_counts |> group_by(mlst_profile) |> summarise(count = sum(count)) |> ungroup()
  N <- sum(df_counts$count)
  f1 <- sum(df_counts$count == 1)
  K <- nrow(df_counts)
  q_hat <- if (N > 0) f1 / N else 0
  
  alpha_novel_null <- 0
  alpha_novel_1 <- 1
  alpha_novel_2 <- 2
  alpha_novel_5 <- 5
  alpha_novel_7 <- 7
  alpha_novel_10 <- 10
  alpha_novel_20 <- 20
  alpha_novel_50 <- 50
  alpha_novel_75 <- 75
  alpha_novel_100 <- 100
  alpha_novel_200 <- 200
  alpha_novel_500 <- 500
  alpha_novel_750 <- 750
  alpha_novel_1000 <- 1000
  alpha_novel_2000 <- 2000
  
  tibble(
    Genus = genus_name,
    N = N, f1 = f1, q_hat = q_hat, K = K,
    alpha_novel_candidates = list(c(alpha_novel_null, alpha_novel_1, 
                                    alpha_novel_2, alpha_novel_5, alpha_novel_7, 
                                    alpha_novel_10, alpha_novel_20, alpha_novel_50, alpha_novel_75, 
                                    alpha_novel_100,  alpha_novel_200, alpha_novel_500,  alpha_novel_750, 
                                    alpha_novel_1000,   alpha_novel_2000))
  )
}

ecoli_info <- make_genus_info_exact(ecoli_bsi_mlst_df, "Escherichia")
kleb_info  <- make_genus_info_exact(kleb_bsi_mlst_df, "Klebsiella")

# expand parameter grid
scenarios <- bind_rows(ecoli_info, kleb_info) |>
  unnest(alpha_novel_candidates) |>
  rename(alpha_novel_candidate = alpha_novel_candidates) |>
  crossing(tibble(alpha_named = alpha_named_values)) |>
  mutate(
    A_obs = alpha_named * K,
    alpha_novel_anchored_full = (q_hat / pmax(1 - q_hat, 1e-12)) * (N + A_obs),
    alpha_novel_final = map2(alpha_novel_candidate, alpha_novel_anchored_full, ~ c(.x, .y))
  ) |>
  select(Genus, N, f1, K, alpha_named, alpha_novel_final) |>
  unnest(alpha_novel_final) |>
  mutate(alpha_novel = as.integer(alpha_novel_final)) |>
  select(-alpha_novel_final)

# optional: unique & sort
scenarios <- scenarios |>
  distinct(Genus, alpha_named, alpha_novel) |>
  arrange(Genus, alpha_named, alpha_novel) 
print(scenarios, n = 40)

message("Will run ", nrow(scenarios), " scenarios.")

# run scenarios (serial or parallel)
if (use_parallel) {
  plan(multisession, workers = n_workers)
  combined_mass_list <- future_pmap(
    list(
      df = list(ecoli_bsi_mlst_df, kleb_bsi_mlst_df)[match(scenarios$Genus, c("Escherichia", "Klebsiella"))],
      genus = scenarios$Genus,
      alpha_named = scenarios$alpha_named,
      alpha_novel = scenarios$alpha_novel
    ),
    run_one_scenario_mass_only,
    B_draws = B,
    use_exact = FALSE,
    .options = future_options(seed = TRUE)
  )
  plan(sequential)  # back to sequential
} else {
  combined_mass_list <- pmap(
    list(
      df = list(ecoli_bsi_mlst_df, kleb_bsi_mlst_df)[match(scenarios$Genus, c("Escherichia", "Klebsiella"))],
      genus = scenarios$Genus,
      alpha_named = scenarios$alpha_named,
      alpha_novel = scenarios$alpha_novel
    ),
    run_one_scenario_mass_only,
    B_draws = B,
    use_exact = FALSE
  )
}

# combine into single data.frame; attach scenario ids/labels
combined_mass_dfs_non_exact <- map2_dfr(combined_mass_list, seq_along(combined_mass_list), function(mdf, i) {
  sc <- scenarios[i, ]
  mdf |>
    mutate(
      alpha_named = sc$alpha_named,
      alpha_novel = sc$alpha_novel,
      scenario_id = i
    )
})

# create nice factor labels for plotting
combined_mass_dfs_non_exact <- combined_mass_dfs_non_exact |>
  mutate(
    alpha_named_f = factor(alpha_named, levels = alpha_named_values, labels = paste0("alpha_named=", alpha_named_values)),
    alpha_novel_label = paste0("alpha_novel=", signif(alpha_novel, 3)),
    alpha_novel_f = factor(alpha_novel_label, levels = unique(paste0("alpha_novel=", signif(sort(unique(alpha_novel)), 3))))
  )

# Save compact results
saveRDS(combined_mass_dfs_non_exact, "rarefaction/sensitivity_combined_mass_dfs_mass_only_non_exact.rds")

#~~~~~~~~~~~~~~~~~~#
# plot 
genus_colours <- c("Escherichia" = "seagreen3", "Klebsiella" = "darkorange")

View(combined_mass_dfs_non_exact)
combined_mass_df_non_exact_filtered <- combined_mass_dfs_non_exact |>
  filter(alpha_novel %in% c(0,1,2,5,7,10,20,50,75,100, 200, 500, 750, 1000)) |>
  mutate(Exact = "Non-exact")
combined_mass_df_filtered <- combined_mass_df_filtered |>
  mutate(Exact = "Exact")

# rbind
combined_mass_df_exact_non_exact <- rbind(combined_mass_df_filtered, combined_mass_df_non_exact_filtered)
combined_mass_df_exact_non_exact_filtered <- combined_mass_df_exact_non_exact |>
  filter(alpha_named == 1)


cumulative_mass_plot_sensitivity_non_exact <-  ggplot(combined_mass_df_exact_non_exact_filtered, aes(x = f, y = mean_mass, colour = Genus, fill = Genus,
                                                                                                     group = interaction(Genus, alpha_named_f, alpha_novel_f, scenario_id))) +
  geom_line(aes(linetype = Exact), alpha = 0.95) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.12, colour = NA) +
  #facet_grid(rows = vars(alpha_named_f), cols = vars(alpha_novel_f), labeller = label_value) +
  facet_wrap( ~ alpha_novel, ncol = 5) +
  #facet_wrap( ~ alpha_named) +
  scale_colour_manual(values = genus_colours) +
  scale_fill_manual(values = genus_colours) +
  scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = label_number()) +
  labs(
    x = "Frequency (f) (log scale)",
    y = "Posterior mean of total mass with p >= f",
    title = "Sensitivity analysis: cumulative mass curves across priors" ,
    #subtitle = "Rows = alpha_named; Columns = alpha_novel"
    subtitle = "Faceted by alpha_novel, Exact and Non-exact estimates"
    #subtitle = "Faceted by alpha_named and Genus"
    
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(face = "bold")
  )

print(cumulative_mass_plot_sensitivity_non_exact)
ggsave("cumulative_mass_sensitivity_mass_only_non_exact.png", cumulative_mass_plot_sensitivity_non_exact, width = 16, height = 8, dpi = 300)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. BAYESIAN BOOTSTRAPPING FOR GENES and PLASMIDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# account for within-isolate co-occurrence and distribution of features
# perform isolate-level Bayesian bootstrapping
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# write sub-isolate level bayesboot function
# Bayesian bootstrap on isolates: function that returns posterior draws of prevalences and other stats
isolate_bayesboot <- function(iso_df, gene_cols, B = 2000, alpha_isolate = 1, include_counts = FALSE, copy_cols = NULL) {
  # iso_df: rows = isolates, columns: isolate_id, species (optional), gene_X (0/1), gene_Y (0/1), (integer copy counts or binary presence/absence)
  # final column = r as sum of counts for each isolate
  N <- nrow(iso_df)
  
  # Precompute indicator matrix
  Z <- as.matrix(iso_df[, gene_cols, drop = FALSE])
  r_vec <- iso_df$r
  out <- vector("list", B)
  for (b in seq_len(B)) {
    # draw gamma weights (posterior Dirichlet with shapes = alpha_isolate) for number of genes per isolate
    Gdraws <- rgamma(N, shape = alpha_isolate, rate = 2)
    #normalised so add up to 1
    w <- Gdraws / sum(Gdraws)
    
    # marginal prevalence per gene:
    # multiply actual presence/absence (Z) with weights of number of genes per isolate samples from a gamma dist (w), then sum by columns to get per gene weight 
    p_g <- as.numeric(colSums(w * Z))
    names(p_g) <- gene_cols
    # co-occurrence matrix (pairwise), diagonal = p_g
    # compute pairwise prevalence p_{g1,g2}
    p_pair <- t(Z) %*% (w * Z)           # gives matrix of weighted co-occurrence counts
    p_pair <- as.matrix(p_pair)          # p_pair[g1,g2] = sum_i w_i * z_ig1 * z_ig2
    # per-isolate distribution: probability mass for each observed r value
    r_tab <- tapply(w, r_vec, sum)
    r_tab <- tibble(r = as.integer(names(r_tab)), prob = as.numeric(r_tab))
    # detection prob for sample size m: p_detect(m) = 1 - (1 - p_g)^m
    # we'll compute detection for a vector of m later outside loop if desired
    out[[b]] <- list(w = w, p_g = p_g, p_pair = p_pair, r_tab = r_tab)
  }
  out # return all draws
}


# Postprocess: extract posterior summaries for gene prevalences
extract_prevalence_df <- function(bbs, gene_cols) {
  p_mat <- sapply(bbs, function(x) x$p_g)  # genes x B
  # transpose to B x genes
  p_df <- as.data.frame(t(p_mat))
  colnames(p_df) <- gene_cols
  # compute summaries
  tibble(
    gene = gene_cols,
    mean = colMeans(p_df),
    sd = apply(p_df, 2, sd),
    q2.5 = apply(p_df, 2, quantile, 0.025),
    q97.5 = apply(p_df, 2, quantile, 0.975)
  )
}


# Cumulative population mass curves with increasing frequency 
##get posterior prevalence matrix (rows = draws, cols = genes/features)
#### assumes bbs is a list of length B, each element has $p_g (named vector)
get_p_matrix <- function(bbs) {
  p_list <- lapply(bbs, function(x) x$p_g)
  p_mat <- do.call(rbind, p_list)  # B x G
  colnames(p_mat) <- names(p_list[[1]])
  return(as.matrix(p_mat))
}


# Detection probability curves for a given gene 
compute_detection_curve <- function(bbs, gene, m_values = c(1,5,10,25,50,100,250,500)) {
  # p_g draws
  p_vec <- sapply(bbs, function(x) x$p_g[gene])
  tibble(m = m_values) |>
    mutate(
      detect_prob_mean = map_dbl(m, ~ mean(1 - (1 - p_vec)^.x)),
      detect_prob_lo   = map_dbl(m, ~ quantile(1 - (1 - p_vec)^.x, 0.025)),
      detect_prob_hi   = map_dbl(m, ~ quantile(1 - (1 - p_vec)^.x, 0.975))
    )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 2.1 AMR genes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.1a E. coli ARGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare E.coli data into format where 1 isolate per row, and genes are columns
ecoli_bsi_arg_df <- ecoli_bsi_amrfinder_metadata |>
  filter(Type =="AMR" | is.na(Type)) |>
  group_by(sample, Element.symbol) |> 
  summarise(count = n(),
            presence = case_when(count >0 ~ 1,
                                 count <=0 ~ 0,
                                 TRUE ~ 0)) 
ecoli_bsi_arg_presence <- ecoli_bsi_arg_df|>
  pivot_wider(id_cols = sample, names_from = Element.symbol, values_from = presence, values_fill =0) |>
  select(-c(`NA`)) |>
  ungroup()
#View(ecoli_bsi_arg_presence)
length(unique(ecoli_bsi_arg_presence$sample))

# Set params
iso_df <- ecoli_bsi_arg_presence
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))
#View(iso_df)
B <- 1000
# run
ecoli_bsi_bbs <- isolate_bayesboot(iso_df, gene_cols, B = B, alpha_isolate = 1)
# save
saveRDS(ecoli_bsi_bbs, "rarefaction/ecoli_bsi_ARG_bayesboot.rds")

#~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_bbs <- read_rds("rarefaction/ecoli_bsi_ARG_bayesboot.rds")
#~~~~~~~~~~~~~~~#

# Prevalence historgam
ecoli_bsi_prevalence_summary <- extract_prevalence_df(ecoli_bsi_bbs, gene_cols)
print(ecoli_bsi_prevalence_summary)

# plot histogram of detection probabilities (alpha varies by quantile)
arg_post_hist <- ggplot(data = ecoli_bsi_prevalence_summary) +
  geom_histogram(aes(x = mean), binwidth = 0.01, alpha= 0.5, fill = "seagreen3") +
  geom_histogram(aes(x = q2.5), binwidth = 0.01, alpha= 0.7, fill = "seagreen3") +
  geom_histogram(aes(x = q97.5), binwidth = 0.01, alpha= 0.2, fill = "seagreen3") +
  labs(title = "Histogram of posterior estimated AMR gene prevalence", 
       x = "AMR Gene prevalence\n(mean and 95% CIs)" ) +
  theme_minimal()
arg_post_hist
# save
ggsave("rarefaction/ecoli_bsi_arg_posterior_distribution_histogram.png", arg_post_hist, width = 6, height = 4, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative population mass curves
p_mat <- get_p_matrix(ecoli_bsi_bbs)   # B x G
B <- nrow(p_mat)
G <- ncol(p_mat)

# define frequency grid (avoid 0 due to log scales)
f_grid <- c(0, seq(0.00001, 0.0001, by = 0.000001), seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.0001))
f_grid <- unique(sort(f_grid))

# compute total mass per draw once 
total_mass_vec <- rowSums(p_mat)   # length B
# guard: if any total_mass == 0, set to NA to avoid division by zero later
total_mass_vec[total_mass_vec == 0] <- NA_real_

# normalised mass proportion for each f 
# returns matrix: rows = draws (B), cols = length(f_grid)
prop_mat <- sapply(f_grid, function(f) {
  mask <- p_mat >= f                 # B x G logical
  masked_mass <- rowSums(p_mat * mask) # B-length numeric
  prop <- masked_mass / total_mass_vec
  prop
}) # result is B x length(f_grid) matrix (as columns)
#View(prop_mat)

# summarise across draws for each f
mean_mass <- colMeans(prop_mat, na.rm = TRUE)
q2.5 <- apply(prop_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
q97.5 <- apply(prop_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

ecoli_bsi_mass_df <- tibble(
  f = f_grid,
  mean_mass = mean_mass,
  q2.5 = q2.5,
  q97.5 = q97.5
) |> filter(!is.na(mean_mass))

# min sample sizes using f as feature frequency
ecoli_bsi_mass_df <- ecoli_bsi_mass_df |>
  filter(f > 0) |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_mass_df)
# save
write.csv(ecoli_bsi_mass_df, "rarefaction/ecoli_bsi_ARG_bayesboot_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_mass_df <- read.csv("rarefaction/ecoli_bsi_ARG_bayesboot_mass_df.csv")
#~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.1b Klebsiella ARGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare Klebsiella data into format where 1 isolate per row, and genes are columns
kleb_bsi_arg_df <- kleb_bsi_amrfinder_metadata |>
  filter(Type =="AMR" | is.na(Type)) |>
  group_by(sample, Element.symbol) |> 
  summarise(count = n(),
            presence = case_when(count >0 ~ 1,
                                 count <=0 ~ 0,
                                 TRUE ~ 0)) 
kleb_bsi_arg_presence <- kleb_bsi_arg_df|>
  pivot_wider(id_cols = sample, names_from = Element.symbol, values_from = presence, values_fill =0) |>
  select(-c(`NA`)) |>
  ungroup()
#View(kleb_bsi_arg_presence)


# Set params
iso_df <- kleb_bsi_arg_presence
iso_df <- iso_df |> dplyr::ungroup()
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))
#View(iso_df)
#length(unique(iso_df$sample)) # 468, so all Klebs represented
B <- 1000
# run
kleb_bsi_bbs <- isolate_bayesboot(iso_df, gene_cols, B = B, alpha_isolate = 1)
# save
saveRDS(kleb_bsi_bbs, "rarefaction/kleb_bsi_ARG_bayesboot.rds")
#~~~~~~~~~~~~~#
# read in saved data
#kleb_bsi_bbs <- readRDS("rarefaction/kleb_bsi_ARG_bayesboot.rds")
#~~~~~~~~~~~~~#

#  Prevalence summary 
kleb_bsi_prevalence_summary <- extract_prevalence_df(kleb_bsi_bbs, gene_cols)
#print(kleb_bsi_prevalence_summary)

# plot histogram of detection probabilities:
arg_post_hist <- ggplot(data = kleb_bsi_prevalence_summary) +
  geom_histogram(aes(x = mean), binwidth = 0.01, alpha= 0.5, fill = "darkorange") +
  geom_histogram(aes(x = q2.5), binwidth = 0.01, alpha= 0.7, fill = "darkorange") +
  geom_histogram(aes(x = q97.5), binwidth = 0.01, alpha= 0.2, fill = "darkorange") +
  labs(title = "Histogram of posterior estimated AMR gene prevalence", 
       x = "AMR Gene prevalence\n(mean and 95% CIs)" ) +
  theme_minimal()
arg_post_hist
# save
ggsave("rarefaction/kleb_bsi_arg_posterior_distribution_histogram.png", arg_post_hist, width = 6, height = 4, units = "in", dpi = 300)

# cumulative population mass curves with increasing frequency 
p_mat <- get_p_matrix(kleb_bsi_bbs)   # B x G
B <- nrow(p_mat)
G <- ncol(p_mat)

f_grid <- c(0, seq(0.00001, 0.0001, by = 0.000001), seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
f_grid <- unique(sort(f_grid))

# compute total mass per draw once 
total_mass_vec <- rowSums(p_mat)   # length B
# guard: if any total_mass == 0, set to NA to avoid division by zero later
total_mass_vec[total_mass_vec == 0] <- NA_real_

# normalised mass proportion for each f 
prop_mat <- sapply(f_grid, function(f) {
  mask <- p_mat >= f                 # B x G logical
  masked_mass <- rowSums(p_mat * mask) # B-length numeric
  prop <- masked_mass / total_mass_vec
  prop
}) # result is B x length(f_grid) matrix (as columns)
#View(prop_mat)

# summarise across draws for each f 
mean_mass <- colMeans(prop_mat, na.rm = TRUE)
q2.5 <- apply(prop_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
q97.5 <- apply(prop_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

kleb_bsi_mass_df <- tibble(
  f = f_grid,
  mean_mass = mean_mass,
  q2.5 = q2.5,
  q97.5 = q97.5
) |> filter(!is.na(mean_mass))

# min sample sizes using f as feature frequency 
kleb_bsi_mass_df <- kleb_bsi_mass_df |>
  filter(f > 0) |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_mass_df)
# save
write.csv(kleb_bsi_mass_df, "rarefaction/kleb_bsi_arg_bayesboot_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~#
# read in saved df
#kleb_bsi_mass_df <- read.csv("rarefaction/kleb_bsi_arg_bayesboot_mass_df.csv")
#~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.1c Combined plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# combined Prevalence Histogram
ecoli_bsi_prevalence_summary <- ecoli_bsi_prevalence_summary |>
  mutate(Genus = "Escherichia")
kleb_bsi_prevalence_summary <- kleb_bsi_prevalence_summary |>
  mutate(Genus = "Klebsiella")
combined_bsi_prevalence_summary <- rbind(ecoli_bsi_prevalence_summary, kleb_bsi_prevalence_summary)
#View(combined_bsi_prevalence_summary)

# plot combined histogram of detection probabilities:
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
combined_bsi_arg_post_hist <- ggplot(data = combined_bsi_prevalence_summary, (aes(colour = Genus, fill = Genus))) +
  geom_histogram(aes(x = mean, y = after_stat(density)), binwidth = 0.01, alpha= 0.4, colour = NA, position = "identity") +
  geom_histogram(aes(x = q2.5, y = after_stat(density)), binwidth = 0.01, alpha= 0.6, colour = NA, position = "identity") +
  geom_histogram(aes(x = q97.5, y = after_stat(density)), binwidth = 0.01, alpha= 0.2, colour = NA, position = "identity") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "Histogram of posterior estimated AMR gene prevalence", 
    x = "AMR Gene prevalence\n(mean and 95% CIs)",
    y = "Density (%)") +
  theme_minimal()
combined_bsi_arg_post_hist
# save
ggsave("rarefaction/combined_bsi_arg_posterior_distribution_histogram.png", combined_bsi_arg_post_hist, width = 8, height = 4, units = "in", dpi = 300)

# combined E. coli and Klebsiella cumulative mass curve (normalised to 1) 
# rbind
combined_bsi_arg_mass_df <- rbind(ecoli_bsi_mass_df, kleb_bsi_mass_df)
# save
write.csv(combined_bsi_arg_mass_df, "rarefaction/combined_bsi_arg_bayesboot_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~#
# read in saved df
#combined_bsi_arg_mass_df <- read.csv("rarefaction/combined_bsi_arg_bayesboot_mass_df.csv")
#~~~~~~~~~~~~~#

# plot (mean + 95% CI) 
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
combined_cumulative_mass_plot <- ggplot(combined_bsi_arg_mass_df, aes(x = f, colour = Genus, fill = Genus)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, colour = NA) +
  geom_line(aes(y = mean_mass), size = 1) +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(min(mass_df$f), 1)) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(
    x = "ARG frequency (f) (log scale)",
    y = "Proportion of ARGs belonging to an allele of frequency ≥ f",
    #title = "Cumulative mass curve: proportion of population at frequency ≥ f",
    #subtitle = "Mean (line) and 95% posterior interval (ribbon)"
  ) +
  theme_minimal(base_size = 14)
combined_cumulative_mass_plot

#save
ggsave("rarefaction/combined_bsi_arg_cumulative_mass.png", combined_cumulative_mass_plot, width = 8, height = 6, dpi = 300)


# Combined E. coli and Klebsiella sample coverage plot
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
bootbayes_sample_coverage_plot <- ggplot(combined_bsi_arg_mass_df, aes(x = mean_mass, colour = Genus, fill = Genus)) +
  geom_line(aes(y = min_sample_90)) +
  geom_line(aes(y = min_sample_95) ) +
  geom_line(aes(y = min_sample_99)) +
  geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.15, colour = NA) +
  geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.3, colour = NA) +
  geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.45, colour = NA) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  geom_vline(xintercept = 0.80, linetype = "dashed") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(0,5000)) +
  labs(x = "Sample coverage",
       y = "Sample size" #,
       #title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bootbayes_sample_coverage_plot
ggsave("rarefaction/combined_bsi_arg_bootbayes_sample_coverage_plot_80.png", plot = bootbayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.1d Combined summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Summary table for E. coli and Klebsiella ARGs 
# transform to long
combined_bsi_arg_mass_df_long <- combined_bsi_arg_mass_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_arg_sample_coverage_summary_table <- combined_bsi_arg_mass_df_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "mean_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_arg_sample_coverage_summary_table <- combined_bsi_arg_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_arg_sample_coverage_summary_table)

write.csv(combined_bsi_arg_sample_coverage_summary_table, "rarefaction/combined_bsi_arg_sample_coverage_summary_table.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 2.2 Plasmids ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.2a E. coli plasmids ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare E. coli plasmids df
all_samples <- ecoli_bsi_samples_metadata |> distinct(isolateid)   

# detected pairs from AMRFinder (1 if detected)
detected <- ecoli_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)

# full grid of sample x feature using master sample list and the set of observed features
all_features <- detected |> pull(community_subcommunity) |> unique()
full_grid <- tidyr::expand_grid(sample = all_samples$isolateid,
                                community_subcommunity = all_features)
# left join detections onto the full grid and fill NAs with 0
presence_long <- full_grid |>
  left_join(detected, by = c("sample", "community_subcommunity")) |>
  mutate(presence = if_else(is.na(presence), 0L, presence))
# then pivot to wide
ecoli_bsi_pling_df <- presence_long |>
  pivot_wider(names_from = community_subcommunity,
              values_from = presence,
              values_fill = 0L) |>
  ungroup()
#View(ecoli_bsi_pling_df)
#length(unique(ecoli_bsi_pling_df$sample))# 1471

# Set params
iso_df <- ecoli_bsi_pling_df
iso_df <- iso_df |> dplyr::ungroup()
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))
#length(unique(iso_df$sample))# 1471
#View(iso_df)
B <- 1000

# Run isolate-level bayesian bootstrapping for E.coli pling plasmids
ecoli_bsi_pling_bbs <- isolate_bayesboot(iso_df, gene_cols, B = B, alpha_isolate = 1)

# save
saveRDS(ecoli_bsi_pling_bbs, "rarefaction/ecoli_bsi_PLING_bayesboot.rds")
#~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_pling_bbs <- read_rds("rarefaction/ecoli_bsi_PLING_bayesboot.rds")
#~~~~~~~~~~~~#

# Prevalence summary
ecoli_bsi_pling_prevalence_summary <- extract_prevalence_df(ecoli_bsi_pling_bbs, gene_cols)
#print(ecoli_bsi_pling_prevalence_summary)

# plot histogram of detection probabilities:
pling_post_hist <- ggplot(data = ecoli_bsi_pling_prevalence_summary) +
  geom_histogram(aes(x = mean), binwidth = 0.01, alpha= 0.5,  fill = "seagreen3") +
  geom_histogram(aes(x = q2.5), binwidth = 0.01, alpha= 0.7,  fill = "seagreen3") +
  geom_histogram(aes(x = q97.5), binwidth = 0.01, alpha= 0.2,  fill = "seagreen3") +
  labs(title = "Histogram of posterior estimated PLING plasmid subcommunity prevalence", 
       x = "AMR Gene prevalence\n(mean and 95% CIs)" ) +
  theme_minimal()
pling_post_hist
# save
ggsave("rarefaction/ecoli_bsi_pling_posterior_distribution_histogram.png", pling_post_hist, width = 8, height = 4, units = "in", dpi = 300)

# Cumulative population mass curves with increasing frequency 
p_mat <- get_p_matrix(ecoli_bsi_pling_bbs)   # B x G
B <- nrow(p_mat)
G <- ncol(p_mat)
min(p_mat)

#f_grid <- c(0, seq(0.00001, 0.0001, by = 0.000001), seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
#f_grid <- unique(sort(f_grid))

# compute total mass per draw once 
total_mass_vec <- rowSums(p_mat)   # length B
total_mass_vec[total_mass_vec == 0] <- NA_real_

# normalised mass proportion for each f
prop_mat <- sapply(f_grid, function(f) {
  mask <- p_mat >= f                 # B x G logical
  masked_mass <- rowSums(p_mat * mask) # B-length numeric
  prop <- masked_mass / total_mass_vec
  prop
}) # result is B x length(f_grid) matrix (as columns)
#View(prop_mat)

# summarise across draws for each f 
mean_mass <- colMeans(prop_mat, na.rm = TRUE)
q2.5 <- apply(prop_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
q97.5 <- apply(prop_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

ecoli_bsi_mass_df <- tibble(
  f = f_grid,
  mean_mass = mean_mass,
  q2.5 = q2.5,
  q97.5 = q97.5
) |> filter(!is.na(mean_mass))

# min sample sizes using f as feature frequency 
# (these are unchanged; they depend on f not on normalisation)
ecoli_bsi_mass_df <- ecoli_bsi_mass_df |>
  filter(f > 0) |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_mass_df)
# save
write.csv(ecoli_bsi_mass_df, "rarefaction/ecoli_bsi_pling_bayesboot_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_mass_df <-read.csv("rarefaction/ecoli_bsi_pling_bayesboot_mass_df.csv")
#~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.2b Klebsiella plasmid subcommunities ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prep data
all_samples <- kleb_bsi_samples_metadata |> distinct(isolateid)   

# detected pairs from AMRFinder (1 if detected)
detected <- kleb_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)

# full grid of sample x feature using master sample list and the set of observed features
all_features <- detected |> pull(community_subcommunity) |> unique()
full_grid <- tidyr::expand_grid(sample = all_samples$isolateid,
                                community_subcommunity = all_features)
# left join detections onto the full grid and fill NAs with 0
presence_long <- full_grid |>
  left_join(detected, by = c("sample", "community_subcommunity")) |>
  mutate(presence = if_else(is.na(presence), 0L, presence))
# then pivot to wide
kleb_bsi_pling_df <- presence_long |>
  pivot_wider(names_from = community_subcommunity,
              values_from = presence,
              values_fill = 0L) |>
  ungroup()
#View(kleb_bsi_pling_df)
#length(unique(kleb_bsi_pling_df$sample))# 468

# Set params
iso_df <- kleb_bsi_pling_df
iso_df <- iso_df |> dplyr::ungroup()
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))
#length(unique(iso_df$sample)) # 468
#View(iso_df)
B <- 1000

# Run isolate-level bayesian bootstrapping for E.coli pling plasmids 
kleb_bsi_pling_bbs <- isolate_bayesboot(iso_df, gene_cols, B = B, alpha_isolate = 1)
# save
saveRDS(kleb_bsi_pling_bbs, "rarefaction/kleb_bsi_PLING_bayesboot.rds")
#~~~~~~~~~~~~~#
# read in saved data
kleb_bsi_pling_bbs <- read_rds("rarefaction/kleb_bsi_PLING_bayesboot.rds")
#~~~~~~~~~~~~~#

# Prevalence summary 
kleb_bsi_pling_prevalence_summary <- extract_prevalence_df(kleb_bsi_pling_bbs, gene_cols)
#print(kleb_bsi_pling_prevalence_summary)

# plot histogram of detection probabilities:
pling_post_hist <- ggplot(data = kleb_bsi_pling_prevalence_summary) +
  geom_histogram(aes(x = mean), binwidth = 0.01, alpha= 0.5,  fill = "darkorange") +
  geom_histogram(aes(x = q2.5), binwidth = 0.01, alpha= 0.7,  fill = "darkorange") +
  geom_histogram(aes(x = q97.5), binwidth = 0.01, alpha= 0.2,  fill = "darkorange") +
  labs(title = "Histogram of posterior estimated PLING plasmid subcommunity prevalence", 
       x = "AMR Gene prevalence\n(mean and 95% CIs)" ) +
  theme_minimal()
pling_post_hist
# save
ggsave("rarefaction/kleb_bsi_pling_posterior_distribution_histogram.png", pling_post_hist, width = 8, height = 4, units = "in", dpi = 300)


# Cumulative population mass curves with increasing frequency
p_mat <- get_p_matrix(kleb_bsi_pling_bbs)   # B x G
B <- nrow(p_mat)
G <- ncol(p_mat)

f_grid <- c(0, seq(0.00001, 0.0001, by = 0.000001), seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
f_grid <- unique(sort(f_grid))

# compute total mass per draw once 
total_mass_vec <- rowSums(p_mat)   # length B
# guard: if any total_mass == 0, set to NA to avoid division by zero later
total_mass_vec[total_mass_vec == 0] <- NA_real_

# normalised mass proportion for each f
# returns matrix: rows = draws (B), cols = length(f_grid)
prop_mat <- sapply(f_grid, function(f) {
  mask <- p_mat >= f                 # B x G logical
  masked_mass <- rowSums(p_mat * mask) # B-length numeric
  prop <- masked_mass / total_mass_vec
  prop
}) # result is B x length(f_grid) matrix (as columns)
#View(prop_mat)

# summarise across draws for each f
mean_mass <- colMeans(prop_mat, na.rm = TRUE)
q2.5 <- apply(prop_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
q97.5 <- apply(prop_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

kleb_bsi_mass_df <- tibble(
  f = f_grid,
  mean_mass = mean_mass,
  q2.5 = q2.5,
  q97.5 = q97.5
) |> filter(!is.na(mean_mass))

# optional: min sample sizes using f as feature frequency 
# (these are unchanged; they depend on f not on normalisation)
kleb_bsi_mass_df <- kleb_bsi_mass_df |>
  filter(f > 0) |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_mass_df)
# save
write.csv(kleb_bsi_mass_df, "rarefaction/kleb_bsi_pling_bayesboot_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#kleb_bsi_mass_df <- read.csv("rarefaction/kleb_bsi_pling_bayesboot_mass_df.csv")
#~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.2c Combined plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# rbind
combined_bsi_pling_mass_df <- rbind(ecoli_bsi_mass_df, kleb_bsi_mass_df)
# save
write.csv(combined_bsi_pling_mass_df, "rarefaction/combined_bsi_pling_bayesboot_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~#
# read-in saved df 
#combined_bsi_pling_mass_df <- read.csv("rarefaction/combined_bsi_pling_bayesboot_mass_df.csv")
#~~~~~~~~~~~~~~~~~~#

# Combined Prevalence Histogram with E. coli
ecoli_bsi_pling_prevalence_summary <- ecoli_bsi_pling_prevalence_summary |>
  mutate(Genus = "Escherichia")
kleb_bsi_pling_prevalence_summary <- kleb_bsi_pling_prevalence_summary |>
  mutate(Genus = "Klebsiella")
combined_bsi_pling_prevalence_summary <- rbind(ecoli_bsi_pling_prevalence_summary, kleb_bsi_pling_prevalence_summary)
#View(combined_bsi_pling_prevalence_summary)
# plot combined histogram of detection probabilities:
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
combined_bsi_pling_post_hist <- ggplot(data = combined_bsi_pling_prevalence_summary, (aes(colour = Genus, fill = Genus))) +
  geom_histogram(aes(x = mean, y = after_stat(density)), binwidth = 0.01, alpha= 0.4, colour = NA, position = "identity") +
  geom_histogram(aes(x = q2.5, y = after_stat(density)), binwidth = 0.01, alpha= 0.6, colour = NA, position = "identity") +
  geom_histogram(aes(x = q97.5, y = after_stat(density)), binwidth = 0.01, alpha= 0.2, colour = NA, position = "identity") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "Histogram of posterior estimated AMR gene prevalence", 
    x = "AMR Gene prevalence\n(mean and 95% CIs)",
    y = "Density (%)") +
  theme_minimal()
combined_bsi_pling_post_hist
# save
ggsave("rarefaction/combined_bsi_pling_posterior_distribution_histogram.png", combined_bsi_pling_post_hist, width = 8, height = 4, units = "in", dpi = 300)


# combined cumulative mass plot (mean + 95% CI) 
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
combined_cumulative_mass_plot <- ggplot(combined_bsi_pling_mass_df, aes(x = f, colour = Genus, fill = Genus)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, colour = NA) +
  geom_line(aes(y = mean_mass), size = 1) +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(min(combined_bsi_pling_mass_df$f), 1)) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(
    x = "Plasmid subcommunity frequency (f) (log scale)",
    y = "Proportion of plasmids in a subcommunity of frequency ≥ f",
    #title = "Cumulative mass curve: proportion of population at frequency ≥ f",
    #subtitle = "Mean (line) and 95% posterior interval (ribbon)"
  ) +
  theme_minimal(base_size = 14)
combined_cumulative_mass_plot

#save
ggsave("rarefaction/combined_bsi_pling_cumulative_mass.png", combined_cumulative_mass_plot, width = 8, height = 6, dpi = 300)


# Combined sample coverage vs sample size plot
bootbayes_sample_coverage_plot <- ggplot(combined_bsi_pling_mass_df, aes(x = mean_mass, colour = Genus, fill = Genus)) +
  geom_line(aes(y = min_sample_90)) +
  geom_line(aes(y = min_sample_95) ) +
  geom_line(aes(y = min_sample_99)) +
  geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.15, colour = NA) +
  geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.3, colour = NA) +
  geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.45, colour = NA) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
 # facet_wrap(~ Genus, ncol = 1) +
  geom_vline(xintercept = 0.80, linetype = "dashed") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(0,5000)) +
  labs(x = "Sample coverage",
       y = "Sample size"#,
      # title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bootbayes_sample_coverage_plot
ggsave("rarefaction/combined_bsi_PLING_bootbayes_sample_coverage_plot_80.png", plot = bootbayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.2d Combined summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# transform to long
combined_bsi_pling_mass_df_long <- combined_bsi_pling_mass_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_pling_sample_coverage_summary_table <- combined_bsi_pling_mass_df_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "mean_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_pling_sample_coverage_summary_table <- combined_bsi_pling_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_pling_sample_coverage_summary_table)
# save
write.csv(combined_bsi_pling_sample_coverage_summary_table, "rarefaction/combined_bsi_pling_sample_coverage_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. PANEL PLOT FOR E. coli and Kleb posterior estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data
#mass_df_mlst <- read.csv("combined_bsi_kleborate_mlst_bayesboot_cumulative_frequency_mass_df.csv")
#mass_df_fastbaps_L3 <- read.csv("combined_bsi_fastbaps_L3_bayesboot_cumulative_frequency_mass_df.csv")
#combined_bsi_pling_mass_df <- read.csv("rarefaction/combined_bsi_pling_bayesboot_mass_df.csv")
#combined_bsi_arg_mass_df <- read.csv("rarefaction/combined_bsi_arg_bayesboot_mass_df.csv")

# define palate
genus_colours <- c("Escherichia" = "seagreen3", "Klebsiella"   = "darkorange")


# plot functions
big_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 13),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = grid::unit(0.28, "cm"),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(size = 15, face = "bold")
  )

make_mass_plot <- function(df, xlab, ylab, xlim_min = NULL) {
  p <- ggplot(df, aes(x = f, y = mean_mass, colour = Genus, fill = Genus)) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
    geom_line(linewidth = 1) +
    scale_fill_manual(values = genus_colours) +
    scale_colour_manual(values = genus_colours) +
    scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1),
      labels = c("0.0001", "0.001", "0.01", "0.1", "1")) +
    labs(x = xlab, y = ylab) +
    big_theme
  
  if (!is.null(xlim_min)) {
    p <- p + coord_cartesian(xlim = c(xlim_min, 1))
  }
  p
}

make_cov_plot <- function(df, xlab = "Sample coverage", ylab = "Sample size", y_limit = c(10, 10000)) {
    df_long <- df |>
    pivot_longer(cols = c(min_sample_90, min_sample_95, min_sample_99),
                 names_to = "confidence",
                 values_to = "sample_size" ) |>
    mutate(confidence = factor(confidence,
                               levels = c("min_sample_90", "min_sample_95", "min_sample_99"),
                               labels = c("90%", "95%", "99%"))) |>
      filter(mean_mass != 0)
    
    df_90 <- df_long |> filter(confidence == "90%")
    df_95 <- df_long |> filter(confidence == "95%")
    df_99 <- df_long |> filter(confidence == "99%") 

    ggplot(df_long, aes(x = mean_mass, colour = Genus)) +
      geom_ribbon(data = df_90,aes(xmin = q2.5, xmax = q97.5, y = sample_size, fill = Genus, group = Genus),
        inherit.aes = FALSE, alpha = 0.10, colour = NA) +
      geom_ribbon(data = df_95, aes(xmin = q2.5, xmax = q97.5, y = sample_size, fill = Genus,  group = Genus),
        inherit.aes = FALSE,alpha = 0.30, colour = NA) +
      geom_ribbon(data = df_99, aes(xmin = q2.5, xmax = q97.5, y = sample_size, fill = Genus,  group = Genus),
        inherit.aes = FALSE, alpha = 0.50, colour = NA) +
      geom_line(data = df_90, aes(y = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_line(data = df_95, aes(y = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_line(data = df_99, aes(y = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_vline(xintercept = 0.80, linetype = "dashed") +
      scale_colour_manual(values = genus_colours, name = "Genus") +
      scale_fill_manual(values = genus_colours, name = "Genus") +
      scale_alpha_manual( values = c("90%" = 0.3, "95%" = 0.55, "99%" = 0.85), name = "Certainty") +
      guides(colour = "none", fill = "none", alpha = guide_legend(title = "Certainty")) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_log10(limits = y_limit, 
                    breaks = c(10, 100, 1000, 10000),
                    labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      labs(x = xlab, y = ylab) +
      big_theme +
      theme(legend.position = "bottom")
}



# frequency vs cumulative population mass and sample coverage vs sample size panel plots
p1 <- make_mass_plot(
  mass_df_mlst,
  xlab = "MLST frequency (f)",
  ylab = "Proportion of bacterial population\nwith MLST of frequency ≥ f",
  xlim_min = min(mass_df_mlst$f[mass_df_mlst$f > 0], na.rm = TRUE)
)
p1

p2 <- make_cov_plot(mass_df_mlst)
p2

p3 <- make_mass_plot(
  mass_df_fastbaps_L3,
  xlab = "FastBAPS cluster frequency (f)",
  ylab = "Proportion of bacterial population\n in cluster of frequency ≥ f",
  xlim_min = min(mass_df_fastbaps_L3$f[mass_df_fastbaps_L3$f > 0], na.rm = TRUE)
)
p3

p4 <- make_cov_plot(mass_df_fastbaps_L3)
p4



p5 <- make_mass_plot(
  combined_bsi_pling_mass_df,
  xlab = "Plasmid subcommunity frequency (f)",
  ylab = "Proportion of plasmid population\n in subcommunity of frequency ≥ f",
  xlim_min = min(combined_bsi_pling_mass_df$f[combined_bsi_pling_mass_df$f > 0], na.rm = TRUE)
)
p5

p6 <- make_cov_plot(combined_bsi_pling_mass_df)
p6

p7 <- make_mass_plot(
  combined_bsi_arg_mass_df,
  xlab = "AMR gene frequency (f)",
  ylab = "Proportion of AMR gene population\nwith allele of frequency ≥ f",
  xlim_min = min(combined_bsi_arg_mass_df$f[combined_bsi_arg_mass_df$f > 0], na.rm = TRUE)
)
p7

p8 <- make_cov_plot(combined_bsi_arg_mass_df)
p8

# combine plots:
combined_figure <-
  (p1 | p2) /
  (p3 | p4) /
  (p5 | p6) /
  (p7 | p8) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_figure

# save
ggsave("combined_4x2_mass_and_coverage_panel_log_y.png",
  plot = combined_figure, width = 10, height = 15, units = "in", dpi = 300)



# * Summary plot for 80% sample coverage for various metrics
# make summary df (numbers taken from summary tables)
df <- tribble(
  ~metric, ~genus, ~mean, ~lower, ~upper,
  "MLST", "Escherichia", 1360, 1196, 1496,
  "MLST", "Klebsiella", 1196, 1150, 1300,
  "fastBAPS cluster", "Escherichia", 649, 574, 729, # L3 cluster
  "fastBAPS cluster", "Klebsiella", 109, 75, 141, # L3 cluster
  "AMR gene", "Escherichia", 24, 21, 27,
  "AMR gene", "Klebsiella", 67, 56, 80,
  "Plasmid", "Escherichia", 1497, 1301, 1761,
  "Plasmid", "Klebsiella", 998, 809, 1197
)

# add as percentage of annually occurring BSIs (42224 for E. coli, 13078 for Klebsiella)
df <- df |>
  mutate(total_pop = case_when(genus == "Escherichia" ~ 42224,
                               genus == "Klebsiella" ~ 13078),
         mean_pct = mean / total_pop * 100,
         lower_pct = lower / total_pop * 100,
         upper_pct = upper / total_pop * 100)

# order metrics by overall size (largest first)
metric_order <- df |>
  group_by(metric) |>
  summarise(m = max(mean_pct)) |>
  arrange(desc(m)) |>
  pull(metric)

df$metric <- factor(df$metric, levels = metric_order)
#View(df)

# plot
ss_80_barplot <- ggplot(df, aes(x = metric, y = mean, fill = genus)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.7),
    width = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Escherichia" = "seagreen3",
      "Klebsiella" = "darkorange"
    )
  ) +
  labs(x = NULL, y = "Count", fill = "Genus") +
  theme_classic(base_size = 14)

ss_80_barplot
ggsave("sample_size_80_cov_barplot.png", ss_80_barplot, width = 8, height = 4, dpi = 300)

# plot as percentages
ss_80_barplot <- ggplot(df, aes(x = metric, y = mean_pct, fill = genus)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_errorbar(
    aes(ymin = lower_pct, ymax = upper_pct),
    position = position_dodge(width = 0.7),
    width = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Escherichia" = "seagreen3",
      "Klebsiella" = "darkorange"
    )
  ) +
  labs(x = NULL, y = "Percent (%) of annual BSIs to be sampled", fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ss_80_barplot
ggsave("sample_size_80_cov_barplot_pct.png", ss_80_barplot, width = 5, height = 4, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. HIERARCHICAL BAYESIAN DIRICHLET MODEL - REGIONAL BAYESBOOT FOR MLSTS AND FASTBAPS CLUSTERS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# install and load packages if not already installed/loaded
#install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
#install.packages("loo")

#library(cmdstanr)
#library(loo)
#library(tibble)
#library(gtools)
#library(posterior)

# helper functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data prep
# - includes all MLSTs observed anywhere in the dataset
# - adds NOVEL as an explicit category
# - zero counts are kept
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
prep_mlst_data <- function(df,
                           region_col = "region",
                           feature_col = "mlst_profile",
                           count_col = "count",
                           novel_label = "NOVEL",
                           feature_universe = NULL,
                           alpha_novel = 1,
                           alpha_other = 1) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c(region_col, feature_col, count_col) %in% names(df)))
  
  df2 <- df |>
    group_by(.data[[region_col]], .data[[feature_col]]) |>
    summarise(count = sum(.data[[count_col]]), .groups = "drop")
  
  regions <- sort(unique(df2[[region_col]]))
  
  # All known MLSTs across the data, or supplied universe if you have one
  known_features <- sort(unique(df2[[feature_col]]))
  if (!is.null(feature_universe)) {
    known_features <- sort(unique(c(known_features, feature_universe)))
  }
  
  features <- c(setdiff(known_features, novel_label), novel_label)
  novel_idx <- match(novel_label, features)
  
  y <- matrix(0L,
              nrow = length(regions),
              ncol = length(features),
              dimnames = list(regions, features))
  
  rr <- match(df2[[region_col]], regions)
  cc <- match(df2[[feature_col]], features)
  y[cbind(rr, cc)] <- as.integer(df2$count)
  
  alpha <- rep(alpha_other, length(features))
  names(alpha) <- features
  alpha[novel_idx] <- alpha_novel
  
  list(
    y = y,
    regions = regions,
    features = features,
    novel_label = novel_label,
    novel_idx = novel_idx,
    alpha = alpha
  )
}


#~~~~~~~~~~~~~~~~~~~~~~~~#
# Stan code: shared tau
stan_shared_tau <- '
data {
  int<lower=1> R;
  int<lower=2> K;
  array[R, K] int<lower=0> y;
  vector<lower=0>[K] alpha;
}
parameters {
  simplex[K] pi;
  real<lower=1e-8> tau;
}
model {
  pi ~ dirichlet(alpha);
  tau ~ exponential(1);

  for (r in 1:R) {
    y[r] ~ dirichlet_multinomial(tau * pi);
  }
}
generated quantities {
  vector[R] log_lik;
  for (r in 1:R) {
    log_lik[r] = dirichlet_multinomial_lpmf(y[r] | tau * pi);
  }
}
'

#~~~~~~~~~~~~~~~~~~~~#
# Stan code: region-specific taus (hierarchical over log tau_r)
stan_region_tau <- '
data {
  int<lower=1> R;
  int<lower=2> K;
  array[R, K] int<lower=0> y;
  vector<lower=0>[K] alpha;
}
parameters {
  simplex[K] pi;
  real mu_log_tau;
  real<lower=0> sigma_log_tau;
  vector<lower=0>[R] tau_r;
}
model {
  pi ~ dirichlet(alpha);
  
  mu_log_tau ~ normal(log(20), 0.25);
  sigma_log_tau ~ exponential(2);
  tau_r ~ lognormal(mu_log_tau, sigma_log_tau);
  for (r in 1:R) {
    y[r] ~ dirichlet_multinomial(tau_r[r] * pi);
  }
}
generated quantities {
  vector[R] log_lik;
  for (r in 1:R) {
    log_lik[r] = dirichlet_multinomial_lpmf(y[r] | tau_r[r] * pi);
  }
}
'

write_stan_file_from_string <- function(code, file) {
  writeLines(code, con = file)
  file
}

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit wrappers
fit_shared_tau_model <- function(prep,
                                 stan_file = "shared_tau_model.stan",
                                 iter_warmup = 1000,
                                 iter_sampling = 1000,
                                 chains = 4,
                                 seed = 2026) {
  write_stan_file_from_string(stan_shared_tau, stan_file)
  mod <- cmdstan_model(stan_file)
  
  fit <- mod$sample(
    data = list(
      R = nrow(prep$y),
      K = ncol(prep$y),
      y = prep$y,
      alpha = prep$alpha
    ),
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  
  list(fit = fit, prep = prep, model = "shared_tau")
}

fit_region_tau_model <- function(prep,
                                 stan_file = "region_tau_model.stan",
                                 iter_warmup = 1000,
                                 iter_sampling = 1000,
                                 chains = 4,
                                 seed = 2026) {
  write_stan_file_from_string(stan_region_tau, stan_file)
  mod <- cmdstan_model(stan_file)
  
  fit <- mod$sample(
    data = list(
      R = nrow(prep$y),
      K = ncol(prep$y),
      y = prep$y,
      alpha = prep$alpha
    ),
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  
  list(fit = fit, prep = prep, model = "region_tau")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Helpers to extract posterior draws 
extract_matrix_from_draws_df <- function(draws_df, prefix, n, colnames = NULL) {
  cols <- sprintf("%s[%d]", prefix, seq_len(n))
  missing <- setdiff(cols, names(draws_df))
  if (length(missing) > 0) {
    stop("Missing columns in draws: ", paste(missing, collapse = ", "))
  }
  mat <- as.matrix(draws_df[, cols, drop = FALSE])
  if (!is.null(colnames)) colnames(mat) <- colnames
  mat
}

# extract the estimates of psterior frequencies of MLSTs/ clusters (features)
extract_draws <- function(fit_obj) {
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  
  pi_mat <- extract_matrix_from_draws_df(
    draws_df, "pi", n = length(fit_obj$prep$features), colnames = fit_obj$prep$features
  )
  
  out <- list(
    draws_df = draws_df,
    pi = pi_mat
  )
  
  if ("tau" %in% names(draws_df)) {
    out$tau <- draws_df$tau
  }
  
  tau_r_cols <- grep("^tau_r\\[", names(draws_df), value = TRUE)
  if (length(tau_r_cols) > 0) {
    out$tau_r <- as.matrix(draws_df[, tau_r_cols, drop = FALSE])
    colnames(out$tau_r) <- fit_obj$prep$regions
  }
  
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# K-fold cross validation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
log_mean_exp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(-Inf)
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# Dirichlet-multinomial log PMF
dmultinom_dirichlet_logpmf <- function(y, alpha) {
  y <- as.numeric(y)
  alpha <- as.numeric(alpha)
  
  if (any(!is.finite(alpha)) || any(alpha <= 0)) return(-Inf)
  if (any(!is.finite(y)) || any(y < 0)) return(-Inf)
  
  N <- sum(y)
  A <- sum(alpha)
  
  lgamma(N + 1) - sum(lgamma(y + 1)) +
    lgamma(A) - lgamma(N + A) +
    sum(lgamma(y + alpha) - lgamma(alpha))
}

extract_pi_draws <- function(draws_df, feature_names) {
  cols <- paste0("pi[", seq_along(feature_names), "]")
  if (!all(cols %in% names(draws_df))) {
    stop("Could not find all pi columns in draws.")
  }
  pi <- as.matrix(draws_df[, cols, drop = FALSE])
  colnames(pi) <- feature_names
  pi
}

subset_prep_rows <- function(prep, idx) {
  out <- prep
  out$y <- prep$y[idx, , drop = FALSE]
  out$regions <- prep$regions[idx]
  out
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Held-out log predictive density for one region
# Shared tau model: p(y_test | pi, tau) = Dirichlet-Multinomial(tau * pi)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

log_pred_shared_one_region <- function(y_test, draws_df, feature_names) {
  pi <- extract_pi_draws(draws_df, feature_names)
  tau <- draws_df$tau
  
  S <- nrow(pi)
  lps <- vapply(seq_len(S), function(s) {
    dmultinom_dirichlet_logpmf(y_test, tau[s] * pi[s, ])
  }, numeric(1))
  
  log_mean_exp(lps)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Held-out log predictive density for one region
# Hierarchical region-tau model:
#   pi ~ posterior draw
#   tau_new ~ LogNormal(mu_log_tau, sigma_log_tau)
#   y_test | tau_new, pi ~ Dirichlet-Multinomial(tau_new * pi)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

log_pred_regiontau_one_region <- function(y_test, draws_df, feature_names, n_tau_mc = 50) {
  pi <- extract_pi_draws(draws_df, feature_names)
  
  if (!all(c("mu_log_tau", "sigma_log_tau") %in% names(draws_df))) {
    stop("Expected mu_log_tau and sigma_log_tau in the region-tau fit.")
  }
  
  mu_log_tau <- draws_df$mu_log_tau
  sigma_log_tau <- draws_df$sigma_log_tau
  
  S <- nrow(pi)
  lps_s <- vapply(seq_len(S), function(s) {
    tau_new <- exp(rnorm(n_tau_mc, mean = mu_log_tau[s], sd = sigma_log_tau[s]))
    lps_tau <- vapply(tau_new, function(tn) {
      dmultinom_dirichlet_logpmf(y_test, tn * pi[s, ])
    }, numeric(1))
    log_mean_exp(lps_tau)
  }, numeric(1))
  
  log_mean_exp(lps_s)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# K-fold CV runner
# Uses region-level folds.
# You can set K_folds = nrow(prep$y) for leave-one-region-out.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

kfold_compare_tau_models <- function(prep,
                                     K_folds = min(5, nrow(prep$y)),
                                     iter_warmup = 1000,
                                     iter_sampling = 1000,
                                     chains = 4,
                                     seed = 2026,
                                     n_tau_mc = 50) {
  set.seed(seed)
  
  R <- nrow(prep$y)
  fold_id <- loo::kfold_split_random(K = K_folds, N = R)
  
  fold_results <- vector("list", K_folds)
  
  for (k in seq_len(K_folds)) {
    test_idx <- which(fold_id == k)
    train_idx <- which(fold_id != k)
    
    train_prep <- subset_prep_rows(prep, train_idx)
    test_prep  <- subset_prep_rows(prep, test_idx)
    
    # Fit shared-tau model on training regions
    fit_shared <- fit_shared_tau_model(
      train_prep,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = seed + k
    )
    draws_shared <- posterior::as_draws_df(fit_shared$fit$draws())
    
    # Fit region-tau model on training regions
    fit_region <- fit_region_tau_model(
      train_prep,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = seed + 1000 + k
    )
    draws_region <- posterior::as_draws_df(fit_region$fit$draws())
    
    # Compute held-out log predictive densities for each held-out region
    region_rows <- lapply(seq_len(nrow(test_prep$y)), function(j) {
      y_test <- as.numeric(test_prep$y[j, ])
      region_name <- test_prep$regions[j]
      
      elpd_shared <- log_pred_shared_one_region(
        y_test = y_test,
        draws_df = draws_shared,
        feature_names = prep$features
      )
      
      elpd_region <- log_pred_regiontau_one_region(
        y_test = y_test,
        draws_df = draws_region,
        feature_names = prep$features,
        n_tau_mc = n_tau_mc
      )
      
      data.frame(
        fold = k,
        region = region_name,
        elpd_shared = elpd_shared,
        elpd_region = elpd_region,
        delta = elpd_region - elpd_shared
      )
    })
    
    fold_results[[k]] <- bind_rows(region_rows)
  }
  
  per_region <- bind_rows(fold_results)
  
  summary <- per_region |>
    summarise(
      shared_elpd = sum(elpd_shared),
      region_elpd = sum(elpd_region),
      delta_elpd = sum(delta),
      se_delta = sqrt(n() * var(delta)),
      .groups = "drop"
    )
  
  list(
    fold_id = fold_id,
    per_region = per_region,
    summary = summary
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Compute mass curves 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
extract_pi_draws <- function(draws_df, feature_names) {
  cols <- paste0("pi[", seq_along(feature_names), "]")
  if (!all(cols %in% names(draws_df))) {
    stop("Missing pi columns in draws.")
  }
  pi <- as.matrix(draws_df[, cols, drop = FALSE])
  colnames(pi) <- feature_names
  pi
}

compute_region_global_mass_curves <- function(fit_obj,
                                              f_grid,
                                              n_post_draws = 400,
                                              seed = 2026) {
  set.seed(seed)
  
  prep <- fit_obj$prep
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  feature_names <- prep$features
  novel_idx <- prep$novel_idx
  
  pi_draws <- extract_pi_draws(draws_df, feature_names)
  B <- nrow(pi_draws)
  draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
  
  has_shared_tau <- "tau" %in% names(draws_df)
  has_region_tau <- any(grepl("^tau_r\\[", names(draws_df)))
  
  tau_shared <- if (has_shared_tau) draws_df$tau else NULL
  
  tau_r_mat <- NULL
  if (has_region_tau) {
    tau_r_cols <- grep("^tau_r\\[", names(draws_df), value = TRUE)
    tau_r_mat <- as.matrix(draws_df[, tau_r_cols, drop = FALSE])
    colnames(tau_r_mat) <- prep$regions
  }
  
  out_rows <- list()
  ii <- 0
  
  for (r in seq_along(prep$regions)) {
    region_name <- prep$regions[r]
    y_r <- as.numeric(prep$y[r, ])
    
    # storage for posterior masses across draws and thresholds
    region_mass_mat <- matrix(NA_real_, nrow = length(draw_idx), ncol = length(f_grid))
    global_mass_mat <- matrix(NA_real_, nrow = length(draw_idx), ncol = length(f_grid))
    global_all_mat <- matrix(NA_real_, nrow = length(draw_idx), ncol = length(f_grid))
    
    for (i in seq_along(draw_idx)) {
      b <- draw_idx[i]
      pi_b <- as.numeric(pi_draws[b, ])
      
      tau_b <- if (has_shared_tau) {
        tau_shared[b]
      } else if (has_region_tau) {
        tau_r_mat[b, region_name]
      } else {
        stop("No tau or tau_r found in fit.")
      }
      
      alpha_post <- y_r + tau_b * pi_b
      theta_b <- as.numeric(gtools::rdirichlet(1, alpha_post))
      
      for (j in seq_along(f_grid)) {
        f <- f_grid[j]
        
        sel_region <- theta_b >= f
        
        # region mass at threshold f
        region_mass_mat[i, j] <- sum(theta_b[sel_region])
        
        # global mass carried by those same MLSTs
        global_mass_mat[i, j] <- sum(pi_b[sel_region])
        
        # optional: global mass of all MLSTs whose global freq is >= f
        global_all_mat[i, j] <- sum(pi_b[pi_b >= f])
      }
    }
    
    region_summary <- tibble(
      region = region_name,
      f = f_grid,
      
      region_median = apply(region_mass_mat, 2, stats::median),
      region_q2.5 = apply(region_mass_mat, 2, stats::quantile, probs = 0.025, names = FALSE),
      region_q97.5 = apply(region_mass_mat, 2, stats::quantile, probs = 0.975, names = FALSE),
      
      global_median = apply(global_mass_mat, 2, stats::median),
      global_q2.5 = apply(global_mass_mat, 2, stats::quantile, probs = 0.025, names = FALSE),
      global_q97.5 = apply(global_mass_mat, 2, stats::quantile, probs = 0.975, names = FALSE),
      
      global_all_median = apply(global_all_mat, 2, stats::median),
      global_all_q2.5 = apply(global_all_mat, 2, stats::quantile, probs = 0.025, names = FALSE),
      global_all_q97.5 = apply(global_all_mat, 2, stats::quantile, probs = 0.975, names = FALSE)
    ) |>
      mutate(
        min_sample_90 = if_else(f > 0 & f < 1, ceiling(log(1 - 0.90) / log(1 - f)), NA_integer_),
        min_sample_95 = if_else(f > 0 & f < 1, ceiling(log(1 - 0.95) / log(1 - f)), NA_integer_),
        min_sample_99 = if_else(f > 0 & f < 1, ceiling(log(1 - 0.99) / log(1 - f)), NA_integer_)
      )
    
    out_rows[[ii <- ii + 1]] <- region_summary
  }
  
  bind_rows(out_rows)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 4.1 MLSTs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 4.1a E. coli MLSTs - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ecoli_bsi_count_region <- ecoli_bsi_samples_metadata |>
  group_by(region, escherichia__mlst_achtman__ST) |>
  summarise(count = n())

 prep <- prep_mlst_data(ecoli_bsi_count_region,
                        region_col = "region",
                        feature_col = "escherichia__mlst_achtman__ST",
                        count_col = "count",
                        novel_label = "NOVEL",
                        feature_universe = NULL,
                        alpha_novel = 1,
                        alpha_other = 1)
#View(prep)
str(prep) # list of 6

fit_u <- fit_shared_tau_model(prep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r <- fit_region_tau_model(prep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# save
saveRDS(fit_u, "ecoli_bsi_mlst_shared_tau_fit.rds")
saveRDS(fit_r, "ecoli_bsi_mlst_regional_tau_fit.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#fit_u <- readRDS("ecoli_bsi_mlst_shared_tau_fit.rds")
#fit_r <- readRDS("ecoli_bsi_mlst_regional_tau_fit.rds")
#~~~~~~~~~~~~~~~~#

# model summaries and diagnostics
str(fit_u)
str(fit_u$prep)
str(fit_u$model)

str(fit_u$fit)
print(fit_u$fit$summary(), n=300)
print(fit_r$fit$summary(), n=300)

fit_u$fit$cmdstan_summary()
fit_r$fit$cmdstan_summary()

fit_u$fit$cmdstan_diagnose()
fit_r$fit$cmdstan_diagnose()

fit_u$fit$diagnostic_summary()
fit_r$fit$diagnostic_summary()


fit_u$fit$sampler_diagnostics()

# k-fold cross-validation
kfold_compare <- kfold_compare_tau_models(prep,
                                      K_folds = min(5, nrow(prep$y)),
                                      iter_warmup = 1000,
                                      iter_sampling = 1000,
                                      chains = 4,
                                      seed = 2026,
                                      n_tau_mc = 50)
str(kfold_compare)
kfold_compare$summary
kfold_compare$per_region
kfold_compare$fold_id
# improved expected log predictive density (ELPD) for all regions, so region-specific tau improves all regions predictions
# save
saveRDS(kfold_compare, "ecoli_bsi_mlst_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#kfold_compare <- readRDS("ecoli_bsi_mlst_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Mass curves with global and regional masses  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
f_grid <- c(seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))

curve_df <- compute_region_global_mass_curves(
  fit_obj = fit_r,     # or fit_r
  f_grid = f_grid,
  n_post_draws = 1000)
#View(curve_df)
colnames(curve_df)
curves$global_curve
curves$region_curves

# save
write.csv(curve_df, "ecoli_bsi_mlst_mass_curves.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * * 4.1.ai Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# plot functions:
big_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 13),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = grid::unit(0.28, "cm"),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(size = 15, face = "bold")
  )

# plot regional mass curve - by global mass

make_mass_plot_region <- function(df, xlab, ylab, xlim_min = NULL, pal = region_pal) {
  p <- ggplot(df|> dplyr::filter(f > 0), aes(x = f, y = global_median, colour = region, fill = region)) +
    geom_ribbon(aes(ymin = global_q2.5, ymax = global_q97.5), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 1) +
    scale_fill_manual(values = region_pal) +
    scale_colour_manual(values = region_pal) +
    scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
    labs(x = xlab, y = ylab) +
    big_theme
  
  if (!is.null(xlim_min)) {
    p <- p + coord_cartesian(xlim = c(xlim_min, 1))
  }
  p
}


make_cov_plot_region <- function(df, xlab = "Sample coverage", ylab = "Sample size", y_limit = c(10, 10000), pal = region_pal) {
  ggplot(df |> dplyr::filter(f > 0), aes(x = global_median, colour = region, fill = region)) +
    geom_ribbon(aes(y = min_sample_95, xmin = global_q2.5, xmax = global_q97.5), alpha = 0.1, colour = NA) +
    geom_line(aes(y = min_sample_95), linewidth = 1.2) +
    geom_vline(xintercept = 0.80, linetype = "dashed") +
    scale_fill_manual(values = region_pal) +
    scale_colour_manual(values = region_pal) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_log10(limits = y_limit,
                  breaks = c(10, 100, 1000, 10000),
                  labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    labs(x = xlab, y = ylab) +
    big_theme
}

# colourBrewerSet2 + Dark2 themed
region_pal <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
)
region_pal <- rev(region_pal)

# make individual plots:
p1 <- make_mass_plot_region(
  curve_df,
  xlab = "MLST frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(curve_df$f[curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("ecoli_bsi_mlst_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_df, pal = region_pal)
p2
# save
ggsave("ecoli_bsi_mlst_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * * 4.1.aii Summary tables ####
# for E. coli MLSTs by region
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_df |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_df, national_df)

curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
table(curve_df_long$region_conf)

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# helper: for a single estimator df and a threshold, find mins for a given column
min_sample_at_or_above <- function(df, colname, thr) {
  # return NA if no rows meet the condition
  res <- df |>
    filter(!is.na(.data[[colname]])) |>
    filter(.data[[colname]] >= thr) |>
    summarise(min_ss = if (n() == 0) NA_real_ else min(sample_size, na.rm = TRUE)) |>
    pull(min_ss)
  if (length(res) == 0) NA_real_ else res
}


# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_mlst_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(ecoli_bsi_mlst_regional_bhm_summary)

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_mlst_regional_bhm_summary_cells_only <- ecoli_bsi_mlst_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_mlst_regional_bhm_summary)
print(ecoli_bsi_mlst_regional_bhm_summary_cells_only)
#View(ecoli_bsi_mlst_regional_bhm_summary)
#View(ecoli_bsi_mlst_regional_bhm_summary_cells_only)

# save
write.csv(ecoli_bsi_mlst_regional_bhm_summary_cells_only, "rarefaction/ecoli_bsi_mlst_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(ecoli_bsi_mlst_regional_bhm_summary, "rarefaction/ecoli_bsi_mlst_regional_bhm_summary.csv", row.names = FALSE)
#ecoli_bsi_mlst_regional_bhm_summary_cells_only <- read.csv("rarefaction/ecoli_bsi_mlst_regional_bhm_summary_cells_only.csv")
#ecoli_bsi_mlst_regional_bhm_summary <- read.csv("rarefaction/ecoli_bsi_mlst_regional_bhm_summary.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 4.1b Klebsiella MLSTs - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
colnames(kleb_bsi_samples_metadata)
kleb_bsi_count_region <- kleb_bsi_samples_metadata |>
  group_by(region, klebsiella_mlst_ST) |>
  summarise(count = n())

prep <- prep_mlst_data(kleb_bsi_count_region,
                       region_col = "region",
                       feature_col = "klebsiella_mlst_ST",
                       count_col = "count",
                       novel_label = "NOVEL",
                       feature_universe = NULL,
                       alpha_novel = 1,
                       alpha_other = 1)
#View(prep)
str(prep) # list of 6

fit_u <- fit_shared_tau_model(prep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r <- fit_region_tau_model(prep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# save
saveRDS(fit_u, "kleb_bsi_mlst_shared_tau_fit.rds")
saveRDS(fit_r, "kleb_bsi_mlst_regional_tau_fit.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#fit_u <- readRDS("kleb_bsi_mlst_shared_tau_fit.rds")
#fit_r <- readRDS("kleb_bsi_mlst_regional_tau_fit.rds")
#~~~~~~~~~~~~~~~~#
# model summaries and diagnostics
str(fit_u)
str(fit_u$prep)
str(fit_u$model)

str(fit_u$fit)
print(fit_u$fit$summary(), n=300)
print(fit_r$fit$summary(), n=300)

fit_u$fit$cmdstan_summary()
fit_r$fit$cmdstan_summary()

fit_u$fit$cmdstan_diagnose()
fit_r$fit$cmdstan_diagnose()

fit_u$fit$diagnostic_summary()
fit_r$fit$diagnostic_summary()


fit_u$fit$sampler_diagnostics()
# this is not really meaningful as only 10 regions. Better to do k-fold cross validation (see code below)
# this is the reason fo pareto-k-diagnostics being too high
#loo_shared <- fit_u$fit$loo()
#loo_region <- fit_r$fit$loo()
#cmp <- compare_models_loo(fit_u, fit_r)
#print(cmp$comparison)

# run k-fold cross-validation
kfold_compare <- kfold_compare_tau_models(prep,
                                          K_folds = min(5, nrow(prep$y)),
                                          iter_warmup = 1000,
                                          iter_sampling = 1000,
                                          chains = 4,
                                          seed = 2026,
                                          n_tau_mc = 50)
str(kfold_compare)
kfold_compare$summary
kfold_compare$per_region
kfold_compare$fold_id
# improved expected log predictive density (ELPD) for ALMOST all regions
# except North West and South East B (2 regions with fewest samples)
# overall delta 54 +/- 25 so still some evidence of overall improvement of using regional tau model for Klebs


# save
saveRDS(kfold_compare, "kleb_bsi_mlst_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#kfold_compare <- readRDS("kleb_bsi_mlst_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Mass curves with global and regional masses 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#f_grid <- c(seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))

curve_df <- compute_region_global_mass_curves(
  fit_obj = fit_r,     # or fit_u
  f_grid = f_grid,  
  n_post_draws = 1000)
#View(curve_df)
colnames(curve_df)
unique(curve_df$region)
# save
write.csv(curve_df, "kleb_bsi_mlst_mass_curves.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.1.bi Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# make individual plots:
p1 <- make_mass_plot_region(
  curve_df,
  xlab = "MLST frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(mass_df_mlst$f[mass_df_mlst$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("kleb_bsi_mlst_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_df, pal = region_pal)
p2
# save
ggsave("kleb_bsi_mlst_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.1.bii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_df |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_df, national_df)

curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
table(curve_df_long$region_conf)


#thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_mlst_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(kleb_bsi_mlst_regional_bhm_summary)

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_mlst_regional_bhm_summary_cells_only <- kleb_bsi_mlst_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_mlst_regional_bhm_summary)
print(kleb_bsi_mlst_regional_bhm_summary_cells_only)
#View(kleb_bsi_mlst_regional_bhm_summary)
#View(kleb_bsi_mlst_regional_bhm_summary_cells_only)

write.csv(kleb_bsi_mlst_regional_bhm_summary_cells_only, "rarefaction/kleb_bsi_mlst_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(kleb_bsi_mlst_regional_bhm_summary, "rarefaction/kleb_bsi_mlst_regional_bhm_summary.csv", row.names = FALSE)
#kleb_bsi_mlst_regional_bhm_summary_cells_only <- read.csv("rarefaction/kleb_bsi_mlst_regional_bhm_summary_cells_only.csv")
#kleb_bsi_mlst_regional_bhm_summary <- read.csv("rarefaction/kleb_bsi_mlst_regional_bhm_summary.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 4.2 FASTBAPS CLUSTERS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 4.2a E. coli fastBAPS - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
colnames(ecoli_bsi_samples_metadata)
ecoli_bsi_count_region <- ecoli_bsi_samples_metadata |>
  group_by(region, Level.3) |>
  summarise(count = n())

prep <- prep_mlst_data(ecoli_bsi_count_region,
                       region_col = "region",
                       feature_col = "Level.3",
                       count_col = "count",
                       novel_label = "NOVEL",
                       feature_universe = NULL,
                       alpha_novel = 1,
                       alpha_other = 1)
#View(prep)
str(prep) # list of 6

fit_u <- fit_shared_tau_model(prep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r <- fit_region_tau_model(prep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# save
saveRDS(fit_u, "ecoli_bsi_fastbaps_L3_shared_tau_fit.rds")
saveRDS(fit_r, "ecoli_bsi_fastbaps_L3_regional_tau_fit.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#fit_u <- readRDS("ecoli_bsi_fastbaps_L3_shared_tau_fit.rds")
#fit_r <- readRDS("ecoli_bsi_fastbaps_L3_regional_tau_fit.rds")
#~~~~~~~~~~~~~~~~#
# model summaries and diagnostics
str(fit_u)
str(fit_u$prep)
str(fit_u$model)

str(fit_u$fit)
print(fit_u$fit$summary(), n=300)
print(fit_r$fit$summary(), n=300)

fit_u$fit$cmdstan_summary()
fit_r$fit$cmdstan_summary()

fit_u$fit$cmdstan_diagnose()
fit_r$fit$cmdstan_diagnose()

fit_u$fit$diagnostic_summary()
fit_r$fit$diagnostic_summary()

fit_u$fit$sampler_diagnostics()
fit_r$fit$sampler_diagnostics()


# run k-fold cross-validation
kfold_compare <- kfold_compare_tau_models(prep,
                                          K_folds = min(5, nrow(prep$y)),
                                          iter_warmup = 1000,
                                          iter_sampling = 1000,
                                          chains = 4,
                                          seed = 2026,
                                          n_tau_mc = 50)
str(kfold_compare)
kfold_compare$summary
kfold_compare$per_region
kfold_compare$fold_id
# improved expected log predictive density (ELPD) for ALMOST all regions
# except North West and South East B (2 regions with fewest samples)
# overall delta 54 +/- 25 so still some evidence of overall improvement of using regional tau model for ecolis


# save
saveRDS(kfold_compare, "ecoli_bsi_fastbaps_L3_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#kfold_compare <- readRDS("ecoli_bsi_fastbaps_L3_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Mass curves with global and regional masses 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# call function
#f_grid <- c(seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
curve_df <- compute_region_global_mass_curves(
  fit_obj = fit_r,     # or fit_u
  f_grid = f_grid,  
  n_post_draws = 1000)
View(curve_df)
colnames(curve_df)
unique(curve_df$region)
# save
write.csv(curve_df, "ecoli_bsi_fastbaps_L3_mass_curves.csv", row.names = FALSE)
#curve_df <- read.csv("ecoli_bsi_fastbaps_L3_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2.ai Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# make individual plots:
p1 <- make_mass_plot_region(
  curve_df,
  xlab = "fastbaps_L3 frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(mass_df_fastbaps_L3$f[mass_df_fastbaps_L3$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("ecoli_bsi_fastbaps_L3_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_df, pal = region_pal)
p2
# save
ggsave("ecoli_bsi_fastbaps_L3_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2.aii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_df |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_df, national_df)

curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
table(curve_df_long$region_conf)


# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_fastbaps_L3_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(ecoli_bsi_fastbaps_L3_regional_bhm_summary)

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_fastbaps_L3_regional_bhm_summary_cells_only <- ecoli_bsi_fastbaps_L3_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_fastbaps_L3_regional_bhm_summary)
print(ecoli_bsi_fastbaps_L3_regional_bhm_summary_cells_only)
#View(ecoli_bsi_fastbaps_L3_regional_bhm_summary)
View(ecoli_bsi_fastbaps_L3_regional_bhm_summary_cells_only)

write.csv(ecoli_bsi_fastbaps_L3_regional_bhm_summary_cells_only, "rarefaction/ecoli_bsi_fastbaps_L3_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(ecoli_bsi_fastbaps_L3_regional_bhm_summary, "rarefaction/ecoli_bsi_fastbaps_L3_regional_bhm_summary.csv", row.names = FALSE)
#ecoli_bsi_fastbaps_L3_regional_bhm_summary <- read.csv("rarefaction/ecoli_bsi_fastbaps_L3_regional_bhm_summary.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 4.2b Klebsiella fastBAPS - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
colnames(kleb_bsi_samples_metadata)
kleb_bsi_count_region <- kleb_bsi_samples_metadata |>
  group_by(region, Level.3) |>
  summarise(count = n())

prep <- prep_mlst_data(kleb_bsi_count_region,
                       region_col = "region",
                       feature_col = "Level.3",
                       count_col = "count",
                       novel_label = "NOVEL",
                       feature_universe = NULL,
                       alpha_novel = 1,
                       alpha_other = 1)
#View(prep)
#str(prep) # list of 6

fit_u <- fit_shared_tau_model(prep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r <- fit_region_tau_model(prep, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# save
saveRDS(fit_u, "kleb_bsi_fastbaps_L3_shared_tau_fit.rds")
saveRDS(fit_r, "kleb_bsi_fastbaps_L3_regional_tau_fit.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#fit_u <- readRDS("kleb_bsi_fastbaps_L3_shared_tau_fit.rds")
#fit_r <- readRDS("kleb_bsi_fastbaps_L3_regional_tau_fit.rds")
#~~~~~~~~~~~~~~~~#
# model summaries and diagnostics
str(fit_u)
str(fit_u$prep)
str(fit_u$model)

str(fit_u$fit)
print(fit_u$fit$summary(), n=300)
print(fit_r$fit$summary(), n=300)

fit_u$fit$cmdstan_summary()
fit_r$fit$cmdstan_summary()

fit_u$fit$cmdstan_diagnose()
fit_r$fit$cmdstan_diagnose()

fit_u$fit$diagnostic_summary()
fit_r$fit$diagnostic_summary()


fit_u$fit$sampler_diagnostics()
fit_r$fit$sampler_diagnostics()


# run k-fold cross-validation
kfold_compare <- kfold_compare_tau_models(prep,
                                          K_folds = min(5, nrow(prep$y)),
                                          iter_warmup = 1000,
                                          iter_sampling = 1000,
                                          chains = 4,
                                          seed = 2026,
                                          n_tau_mc = 50)
str(kfold_compare)
kfold_compare$summary
kfold_compare$per_region
kfold_compare$fold_id
# improved expected log predictive density (ELPD) for ALMOST all regions
# except North West and South East B (2 regions with fewest samples)
# overall delta 54 +/- 25 so still some evidence of overall improvement of using regional tau model for Klebs


# save
saveRDS(kfold_compare, "kleb_bsi_fastbaps_L3_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#kfold_compare <- readRDS("kleb_bsi_fastbaps_L3_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Mass curves with global and regional masses 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# call function
#f_grid <- c(seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))

curve_df <- compute_region_global_mass_curves(
  fit_obj = fit_r,     # or fit_u
  f_grid = f_grid,  
  n_post_draws = 1000)
#View(curve_df)
colnames(curve_df)
unique(curve_df$region)
# save
write.csv(curve_df, "kleb_bsi_fastbaps_L3_mass_curves.csv", row.names = FALSE)
#curve_df <- read.csv("kleb_bsi_fastbaps_L3_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2.bi Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make individual plots:
p1 <- make_mass_plot_region(
  curve_df,
  xlab = "fastbaps_L3 frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(mass_df_fastbaps_L3$f[mass_df_fastbaps_L3$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("kleb_bsi_fastbaps_L3_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_df, pal = region_pal)
p2
# save
ggsave("kleb_bsi_fastbaos_L3_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2.bii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_df |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_df, national_df)

curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_fastbaps_L3_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(kleb_bsi_fastbaps_L3_regional_bhm_summary)

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_fastbaps_L3_regional_bhm_summary_cells_only <- kleb_bsi_fastbaps_L3_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_fastbaps_L3_regional_bhm_summary)
print(kleb_bsi_fastbaps_L3_regional_bhm_summary_cells_only)
#View(kleb_bsi_fastbaps_L3_regional_bhm_summary)
#View(kleb_bsi_fastbaps_L3_regional_bhm_summary_cells_only)

write.csv(kleb_bsi_fastbaps_L3_regional_bhm_summary_cells_only, "rarefaction/kleb_bsi_fastbaps_L3_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(kleb_bsi_fastbaps_L3_regional_bhm_summary, "rarefaction/kleb_bsi_fastbaps_L3_regional_bhm_summary.csv", row.names = FALSE)

#kleb_bsi_fastbaps_L3_regional_bhm_summary_cells_only <- read.csv("rarefaction/kleb_bsi_fastbaps_L3_regional_bhm_summary_cels_only.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. HIERARCHICAL BAYESIAN BETA-BINOMIAL MODEL FOR ARGS AND PLASMIDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# write feature-level bayesian bootstrapping functions, using Beta-Binomial hierarchical distributions
prep_subiso_data <- function(df,
                             region_col = "region",
                             isolate_id_col = "sample",
                             feature_cols,
                             presence_threshold = 0) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c(region_col, feature_cols) %in% names(df)))
  
  df2 <- df |>
    dplyr::mutate(
      dplyr::across(all_of(feature_cols), ~ as.integer(.x > presence_threshold))
    )
  
  regions <- sort(unique(df2[[region_col]]))
  features <- feature_cols
  
  # number of isolates per region
  n_by_region <- df2 |>
    dplyr::group_by(.data[[region_col]]) |>
    dplyr::summarise(n_iso = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(match(.data[[region_col]], regions))
  
  # region x feature counts of positive isolates
  y_mat <- matrix(0L, nrow = length(regions), ncol = length(features),
                  dimnames = list(regions, features))
  
  for (r in regions) {
    rows_r <- df2[[region_col]] == r
    y_mat[r, ] <- colSums(df2[rows_r, features, drop = FALSE])
  }
  
  # global positive counts and empirical frequencies (optional priors)
  global_pos <- colSums(y_mat)
  global_n <- sum(n_by_region$n_iso)
  global_freq <- global_pos / global_n
  
  list(
    y = y_mat,
    n_iso = setNames(n_by_region$n_iso, n_by_region[[region_col]]),
    regions = regions,
    features = features,
    global_pos = global_pos,
    global_n = global_n,
    global_freq = global_freq
  )
}

# shared tau model
stan_shared_tau_subiso <- '
data {
  int<lower=1> R;
  int<lower=1> K;
  array[R, K] int<lower=0> y;
  array[R] int<lower=1> n_iso;
  vector<lower=0>[K] alpha_prior;
  vector<lower=0>[K] beta_prior;
}
parameters {
  vector<lower=0, upper=1>[K] pi;
  real<lower=1e-8> tau;
}
model {
  for (k in 1:K) {
    pi[k] ~ beta(alpha_prior[k], beta_prior[k]);
  }

  tau ~ exponential(1);

  for (r in 1:R) {
    for (k in 1:K) {
      y[r, k] ~ beta_binomial(n_iso[r], tau * pi[k], tau * (1 - pi[k]));
    }
  }
}
generated quantities {
  vector[R * K] log_lik;
  {
    int idx = 1;
    for (r in 1:R) {
      for (k in 1:K) {
        log_lik[idx] = beta_binomial_lpmf(y[r, k] | n_iso[r], tau * pi[k], tau * (1 - pi[k]));
        idx += 1;
      }
    }
  }
}
'

# region-spcific tau model
stan_region_tau_subiso <- '
data {
  int<lower=1> R;
  int<lower=1> K;
  array[R, K] int<lower=0> y;
  array[R] int<lower=1> n_iso;
  vector<lower=0>[K] alpha_prior;
  vector<lower=0>[K] beta_prior;
}
parameters {
  vector<lower=0, upper=1>[K] pi;
  real mu_log_tau;
  real<lower=0> sigma_log_tau;
  vector<lower=0>[R] tau_r;
}
model {
  for (k in 1:K) {
    pi[k] ~ beta(alpha_prior[k], beta_prior[k]);
  }

  mu_log_tau ~ normal(log(20), 0.25);
  sigma_log_tau ~ exponential(2);
  tau_r ~ lognormal(mu_log_tau, sigma_log_tau);

  for (r in 1:R) {
    for (k in 1:K) {
      y[r, k] ~ beta_binomial(n_iso[r], tau_r[r] * pi[k], tau_r[r] * (1 - pi[k]));
    }
  }
}
generated quantities {
  vector[R * K] log_lik;
  {
    int idx = 1;
    for (r in 1:R) {
      for (k in 1:K) {
        log_lik[idx] = beta_binomial_lpmf(y[r, k] | n_iso[r], tau_r[r] * pi[k], tau_r[r] * (1 - pi[k]));
        idx += 1;
      }
    }
  }
}
'

# fit wrapers to write stan files
write_stan_file_from_string <- function(code, file) {
  writeLines(code, con = file)
  file
}


# wrappers to fit model
fit_shared_tau_subiso <- function(prep,
                                  alpha_prior = NULL,
                                  beta_prior = NULL,
                                  stan_file = "shared_tau_subiso.stan",
                                  iter_warmup = 1000,
                                  iter_sampling = 1000,
                                  chains = 4,
                                  seed = 2026) {
  write_stan_file_from_string(stan_shared_tau_subiso, stan_file)
  mod <- cmdstan_model(stan_file)
  
  K <- ncol(prep$y)
  if (is.null(alpha_prior)) alpha_prior <- rep(1, K)
  if (is.null(beta_prior))  beta_prior  <- rep(1, K)
  
  fit <- mod$sample(
    data = list(
      R = nrow(prep$y),
      K = K,
      y = prep$y,
      n_iso = as.integer(prep$n_iso[prep$regions]),
      alpha_prior = alpha_prior,
      beta_prior = beta_prior
    ),
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  
  list(fit = fit, prep = prep, model = "shared_tau_subiso")
}

fit_region_tau_subiso <- function(prep,
                                  alpha_prior = NULL,
                                  beta_prior = NULL,
                                  stan_file = "region_tau_subiso.stan",
                                  iter_warmup = 1000,
                                  iter_sampling = 1000,
                                  chains = 4,
                                  seed = 2026) {
  write_stan_file_from_string(stan_region_tau_subiso, stan_file)
  mod <- cmdstan_model(stan_file)
  
  K <- ncol(prep$y)
  if (is.null(alpha_prior)) alpha_prior <- rep(1, K)
  if (is.null(beta_prior))  beta_prior  <- rep(1, K)
  
  fit <- mod$sample(
    data = list(
      R = nrow(prep$y),
      K = K,
      y = prep$y,
      n_iso = as.integer(prep$n_iso[prep$regions]),
      alpha_prior = alpha_prior,
      beta_prior = beta_prior
    ),
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  
  list(fit = fit, prep = prep, model = "region_tau_subiso")
}

# rapper to extract posterior draws
extract_pi_draws_subiso <- function(draws_df, feature_names) {
  cols <- paste0("pi[", seq_along(feature_names), "]")
  if (!all(cols %in% names(draws_df))) {
    stop("Could not find all pi columns in draws.")
  }
  pi <- as.matrix(draws_df[, cols, drop = FALSE])
  colnames(pi) <- feature_names
  pi
}

extract_tau_draws_subiso <- function(draws_df, regions = NULL) {
  if ("tau" %in% names(draws_df)) {
    return(list(shared = draws_df$tau, region = NULL))
  }
  
  tau_r_cols <- grep("^tau_r\\[", names(draws_df), value = TRUE)
  if (length(tau_r_cols) > 0) {
    tau_r <- as.matrix(draws_df[, tau_r_cols, drop = FALSE])
    if (!is.null(regions)) colnames(tau_r) <- regions
    return(list(shared = NULL, region = tau_r))
  }
  
  stop("No tau or tau_r found in draws.")
}


# function to calculate cumulative mass curve (normalised)
compute_subiso_mass_curves <- function(fit_obj,
                                       f_grid,
                                       n_post_draws = 400,
                                       seed = 2026) {
  set.seed(seed)
  
  prep <- fit_obj$prep
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  pi_draws <- extract_pi_draws_subiso(draws_df, prep$features)
  tau_info <- extract_tau_draws_subiso(draws_df, regions = prep$regions)
  
  B <- nrow(pi_draws)
  draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
  
  out <- list()
  ii <- 0
  
  for (r in seq_along(prep$regions)) {
    region_name <- prep$regions[r]
    y_r <- as.numeric(prep$y[r, ])
    n_r <- as.integer(prep$n_iso[region_name])
    
    region_mass_mat <- matrix(NA_real_, nrow = length(draw_idx), ncol = length(f_grid))
    global_mass_mat  <- matrix(NA_real_, nrow = length(draw_idx), ncol = length(f_grid))
    global_all_mat   <- matrix(NA_real_, nrow = length(draw_idx), ncol = length(f_grid))
    
    for (i in seq_along(draw_idx)) {
      b <- draw_idx[i]
      pi_b <- as.numeric(pi_draws[b, ])
      
      tau_b <- if (!is.null(tau_info$shared)) {
        tau_info$shared[b]
      } else {
        tau_info$region[b, region_name]
      }
      
      # posterior draw of regional feature prevalences
      theta_b <- rbeta(
        n = length(pi_b),
        shape1 = y_r + tau_b * pi_b,
        shape2 = (n_r - y_r) + tau_b * (1 - pi_b)
      )
      
      # normalize within draw so the mass is on [0, 1]
      w_region <- theta_b / sum(theta_b)
      w_global <- pi_b / sum(pi_b)
      
      for (j in seq_along(f_grid)) {
        f <- f_grid[j]
        sel_region <- theta_b >= f
        
        region_mass_mat[i, j] <- sum(w_region[sel_region])
        global_mass_mat[i, j] <- sum(w_global[sel_region])
        global_all_mat[i, j] <- sum(w_global[w_global >= f])
      }
    }
    
    region_summary <- tibble(
      region = region_name,
      f = f_grid,
      
      region_median = apply(region_mass_mat, 2, median),
      region_q2.5 = apply(region_mass_mat, 2, quantile, probs = 0.025, names = FALSE),
      region_q97.5 = apply(region_mass_mat, 2, quantile, probs = 0.975, names = FALSE),
      
      global_median = apply(global_mass_mat, 2, median),
      global_q2.5 = apply(global_mass_mat, 2, quantile, probs = 0.025, names = FALSE),
      global_q97.5 = apply(global_mass_mat, 2, quantile, probs = 0.975, names = FALSE),
      
      global_all_median = apply(global_all_mat, 2, median),
      global_all_q2.5 = apply(global_all_mat, 2, quantile, probs = 0.025, names = FALSE),
      global_all_q97.5 = apply(global_all_mat, 2, quantile, probs = 0.975, names = FALSE)
    ) |>
      mutate(
        min_sample_90 = if_else(f > 0 & f < 1, ceiling(log(1 - 0.90) / log(1 - f)), NA_integer_),
        min_sample_95 = if_else(f > 0 & f < 1, ceiling(log(1 - 0.95) / log(1 - f)), NA_integer_),
        min_sample_99 = if_else(f > 0 & f < 1, ceiling(log(1 - 0.99) / log(1 - f)), NA_integer_)
      )
    
    out[[ii <- ii + 1]] <- region_summary
  }
  
  bind_rows(out)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  K-fold cross-validation 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Helpers
log_mean_exp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(-Inf)
  m <- max(x)
  m + log(mean(exp(x - m)))
}

dbetabinom_logpmf <- function(y, n, alpha, beta) {
  if (!is.finite(alpha) || !is.finite(beta) || alpha <= 0 || beta <= 0) return(-Inf)
  if (!is.finite(y) || !is.finite(n) || y < 0 || n < 0 || y > n) return(-Inf)
  
  lchoose(n, y) +
    lbeta(y + alpha, n - y + beta) -
    lbeta(alpha, beta)
}

extract_pi_draws_subiso <- function(draws_df, feature_names) {
  cols <- paste0("pi[", seq_along(feature_names), "]")
  if (!all(cols %in% names(draws_df))) {
    stop("Could not find all pi columns in draws.")
  }
  pi <- as.matrix(draws_df[, cols, drop = FALSE])
  colnames(pi) <- feature_names
  pi
}

subset_prep_rows_subiso <- function(prep, idx) {
  out <- prep
  out$y <- prep$y[idx, , drop = FALSE]
  out$regions <- prep$regions[idx]
  out$n_iso <- prep$n_iso[prep$regions[idx]]
  out
}

# Held-out log predictive density for one region
# shared tau model
log_pred_shared_subiso_one_region <- function(y_test_vec,
                                              n_test,
                                              draws_df,
                                              feature_names) {
  pi <- extract_pi_draws_subiso(draws_df, feature_names)
  
  if (!("tau" %in% names(draws_df))) {
    stop("Shared-tau model requires a tau column in draws.")
  }
  tau <- draws_df$tau
  
  S <- nrow(pi)
  
  lps <- vapply(seq_len(S), function(s) {
    sum(vapply(seq_along(feature_names), function(k) {
      dbetabinom_logpmf(
        y = y_test_vec[k],
        n = n_test,
        alpha = tau[s] * pi[s, k],
        beta  = tau[s] * (1 - pi[s, k])
      )
    }, numeric(1)))
  }, numeric(1))
  
  log_mean_exp(lps)
}

# Held-out log predictive density for one region
# region-specific tau model
log_pred_regiontau_subiso_one_region <- function(y_test_vec,
                                                 n_test,
                                                 draws_df,
                                                 feature_names,
                                                 n_tau_mc = 50) {
  pi <- extract_pi_draws_subiso(draws_df, feature_names)
  
  if (!all(c("mu_log_tau", "sigma_log_tau") %in% names(draws_df))) {
    stop("Region-tau model requires mu_log_tau and sigma_log_tau in draws.")
  }
  
  mu_log_tau <- draws_df$mu_log_tau
  sigma_log_tau <- draws_df$sigma_log_tau
  
  S <- nrow(pi)
  
  lps_s <- vapply(seq_len(S), function(s) {
    tau_new <- exp(rnorm(n_tau_mc, mean = mu_log_tau[s], sd = sigma_log_tau[s]))
    
    lps_tau <- vapply(tau_new, function(tn) {
      sum(vapply(seq_along(feature_names), function(k) {
        dbetabinom_logpmf(
          y = y_test_vec[k],
          n = n_test,
          alpha = tn * pi[s, k],
          beta  = tn * (1 - pi[s, k])
        )
      }, numeric(1)))
    }, numeric(1))
    
    log_mean_exp(lps_tau)
  }, numeric(1))
  
  log_mean_exp(lps_s)
}

# Main K-fold comparison function
kfold_compare_subiso_models <- function(prep,
                                        K_folds = min(5, length(prep$regions)),
                                        iter_warmup = 1000,
                                        iter_sampling = 1000,
                                        chains = 4,
                                        seed = 2026,
                                        n_tau_mc = 50) {
  stopifnot(is.list(prep))
  stopifnot(!is.null(prep$y), !is.null(prep$n_iso), !is.null(prep$regions), !is.null(prep$features))
  
  set.seed(seed)
  
  R <- nrow(prep$y)
  if (K_folds > R) K_folds <- R
  
  # region-level folds
  fold_id <- loo::kfold_split_random(K = K_folds, N = R)
  
  fold_results <- vector("list", K_folds)
  
  for (k in seq_len(K_folds)) {
    test_idx <- which(fold_id == k)
    train_idx <- which(fold_id != k)
    
    train_prep <- subset_prep_rows_subiso(prep, train_idx)
    test_prep  <- subset_prep_rows_subiso(prep, test_idx)
    
    # Fit shared-tau model on training regions
    fit_shared <- fit_shared_tau_subiso(
      train_prep,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = seed + k
    )
    draws_shared <- posterior::as_draws_df(fit_shared$fit$draws())
    
    # Fit region-tau model on training regions
    fit_region <- fit_region_tau_subiso(
      train_prep,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = seed + 1000 + k
    )
    draws_region <- posterior::as_draws_df(fit_region$fit$draws())
    
    # Predict each held-out region
    region_rows <- lapply(seq_len(nrow(test_prep$y)), function(j) {
      y_test_vec <- as.numeric(test_prep$y[j, ])
      region_name <- test_prep$regions[j]
      n_test <- as.integer(test_prep$n_iso[region_name])
      
      elpd_shared <- log_pred_shared_subiso_one_region(
        y_test_vec = y_test_vec,
        n_test = n_test,
        draws_df = draws_shared,
        feature_names = prep$features
      )
      
      elpd_region <- log_pred_regiontau_subiso_one_region(
        y_test_vec = y_test_vec,
        n_test = n_test,
        draws_df = draws_region,
        feature_names = prep$features,
        n_tau_mc = n_tau_mc
      )
      
      data.frame(
        fold = k,
        region = region_name,
        n_test = n_test,
        elpd_shared = elpd_shared,
        elpd_region = elpd_region,
        delta = elpd_region - elpd_shared
      )
    })
    
    fold_results[[k]] <- dplyr::bind_rows(region_rows)
  }
  
  per_region <- dplyr::bind_rows(fold_results)
  
  summary <- per_region |>
    summarise(
      shared_elpd = sum(elpd_shared),
      region_elpd = sum(elpd_region),
      delta_elpd = sum(delta),
      se_delta = sqrt(n() * var(delta)),
      .groups = "drop"
    )
  
  list(
    fold_id = fold_id,
    per_region = per_region,
    summary = summary
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 5.1 AMR genes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 5.1a E. coli AMR genes - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load and prepare data, if not already loaded
ecoli_bsi_amrfinder_metadata <- read.csv("rarefaction/ecoli_bsi_amrfinder_metadata.csv")

# prepare E.coli data into format where 1 isolate per row, and genes are columns
ecoli_bsi_arg_presence_region <- ecoli_bsi_amrfinder_metadata |>
  filter(region != "unknown") |>
  filter(Type =="AMR" | is.na(Type)) |>
  group_by(sample, region, Element.symbol) |> 
  summarise(count = n(),
            presence = case_when(count >0 ~ 1,
                                 count <=0 ~ 0,
                                 TRUE ~ 0)) |>
  pivot_wider(id_cols = c(sample, region), names_from = Element.symbol, values_from = presence, values_fill =0) |>
  select(-c(`NA`)) |>
  ungroup()
#View(ecoli_bsi_arg_presence)
#colnames(ecoli_bsi_arg_presence)
#length(unique(ecoli_bsi_arg_presence$sample)) # 1471
gene_cols <- setdiff(names(ecoli_bsi_arg_presence_region), c("sample", "region"))
prep_arg <- prep_subiso_data(
  df = ecoli_bsi_arg_presence_region,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols
)

#View(prep_arg$y)
fit_u_arg <- fit_shared_tau_subiso(prep_arg, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r_arg <- fit_region_tau_subiso(prep_arg, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# save
saveRDS(prep_arg, "ecoli_prep_arg_subiso.rds")
saveRDS(fit_u_arg, "ecoli_fit_u_arg_subiso.rds")
saveRDS(fit_r_arg, "ecoli_fit_r_arg_subiso.rds")
#read in
fit_u_arg_test <- readRDS("ecoli_fit_u_arg_subiso.rds")
fit_r_arg_test <- readRDS("ecoli_fit_r_arg_subiso.rds")

# check model diagnostics 


# call cross-validation function
kfold_compare_subiso <- kfold_compare_subiso_models(
  prep = prep_arg,
  K_folds = min(5, nrow(prep_arg$y)),
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  seed = 2026,
  n_tau_mc = 50
)

# check which better
kfold_compare_subiso$summary
kfold_compare_subiso$per_region
kfold_compare_subiso$fold_id
#save
saveRDS(kfold_compare_subiso, "ecoli_bsi_arg_by_region_kfold_compare.rds")

# cumulative mass grid
f_grid <- c(  seq(0.0001, 0.001, by = 0.00001),   seq(0.0011, 0.01, by = 0.0001),   seq(0.011, 1, by = 0.001))

curve_arg <- compute_subiso_mass_curves(
  fit_obj = fit_r_arg,
  f_grid = f_grid,
  n_post_draws = 1000)
#View(curve_arg)
colnames(curve_arg)
table(curve_arg$region)

write.csv(curve_arg, "ecoli_arg_subiso_mass_curves.csv", row.names = FALSE)
#curve_arg <- read.csv("ecoli_arg_subiso_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.ai Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# use region_pal defined above
# make individual plots:
p1 <- make_mass_plot_region(
  curve_arg,
  xlab = "AMR gene frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(curve_arg$f[curve_arg$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("ecoli_bsi_ARG_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_arg, pal = region_pal)
p2
# save
ggsave("ecoli_bsi_ARG_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.aii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_arg |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_arg, national_df)

# save
write.csv(curve_df2, "ecoli_arg_subiso_mass_curves_with_global.csv", row.names = FALSE)
#curve_df2 <- read.csv("ecoli_arg_subiso_mass_curves.csv")


curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)


# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_arg_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_arg_regional_bhm_summary_cells_only <- ecoli_bsi_arg_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_arg_regional_bhm_summary)
print(ecoli_bsi_arg_regional_bhm_summary_cells_only)
#View(ecoli_bsi_arg_regional_bhm_summary)
#View(ecoli_bsi_arg_regional_bhm_summary_cells_only)

# save
write.csv(ecoli_bsi_arg_regional_bhm_summary_cells_only, "rarefaction/ecoli_bsi_arg_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(ecoli_bsi_arg_regional_bhm_summary, "rarefaction/ecoli_bsi_arg_regional_bhm_summary.csv", row.names = FALSE)
#write.csv(ecoli_bsi_arg_regional_bhm_summary_cells_only, "rarefaction/ecoli_bsi_arg_regional_bhm_summary_cells_only.csv", row.names = FALSE)
#write.csv(ecoli_bsi_arg_regional_bhm_summary, "rarefaction/ecoli_bsi_arg_regional_bhm_summary.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 5.1b Klebsiella AMR genes - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare data
kleb_bsi_amrfinder_metadata <- read.csv("rarefaction/kleb_bsi_amrfinder_metadata.csv")

# prepare Klebsiella data into format where 1 isolate per row, and genes are columns
kleb_bsi_arg_presence_region <- kleb_bsi_amrfinder_metadata |>
  filter(region != "unknown") |>
  filter(Type =="AMR" | is.na(Type)) |>
  group_by(sample, region, Element.symbol) |> 
  summarise(count = n(),
            presence = case_when(count >0 ~ 1,
                                 count <=0 ~ 0,
                                 TRUE ~ 0))|>
  pivot_wider(id_cols = c(sample, region), names_from = Element.symbol, values_from = presence, values_fill =0) |>
  select(-c(`NA`)) |>
  ungroup()
#View(kleb_bsi_arg_presence)
#colnames(kleb_bsi_arg_presence)
#length(unique(kleb_bsi_arg_presence$sample)) # 468

gene_cols <- setdiff(names(kleb_bsi_arg_presence_region), c("sample", "region"))
prep_arg <- prep_subiso_data(
  df = kleb_bsi_arg_presence_region,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols)

fit_u_arg <- fit_shared_tau_subiso(prep_arg, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r_arg <- fit_region_tau_subiso(prep_arg, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# check diagnostics
# model summaries and diagnostics
str(fit_u_arg)
str(fit_u_arg$prep)
str(fit_u_arg$model)

str(fit_u_arg$fit)
print(fit_u_arg$fit$summary(), n=300)
print(fit_r_arg$fit$summary(), n=300)

fit_u_arg$fit$cmdstan_summary()
fit_r_arg$fit$cmdstan_summary()

fit_u_arg$fit$cmdstan_diagnose()
fit_r_arg$fit$cmdstan_diagnose()

fit_u_arg$fit$diagnostic_summary()
fit_r_arg$fit$diagnostic_summary()


fit_u_arg$fit$sampler_diagnostics()
fit_r_arg$fit$sampler_diagnostics()

# save
saveRDS(prep_arg, "kleb_arg_prep_arg_subiso.rds")
saveRDS(fit_u_arg, "kleb_arg_fit_u_arg_subiso.rds")
saveRDS(fit_r_arg, "kleb_arg_fit_r_arg_subiso.rds")

# call cross-validation function
kfold_compare_subiso <- kfold_compare_subiso_models(
  prep = prep_arg,
  K_folds = min(5, nrow(prep_arg$y)),
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  seed = 2026,
  n_tau_mc = 50
)

# check which better
kfold_compare_subiso$summary
kfold_compare_subiso$per_region
kfold_compare_subiso$fold_id

#save
saveRDS(kfold_compare_subiso, "kleb_bsi_arg_by_region_kfold_compare.rds")

# cumulative mass grid
f_grid <- c(  seq(0.0001, 0.001, by = 0.00001),   seq(0.0011, 0.01, by = 0.0001),   seq(0.011, 1, by = 0.001))

curve_arg <- compute_subiso_mass_curves(
  fit_obj = fit_r_arg,
  f_grid = f_grid,
  n_post_draws = 1000)

#View(curve_arg)
colnames(curve_arg)
table(curve_arg$region)
#save
write.csv(curve_arg, "kleb_arg_subiso_mass_curves.csv", row.names = FALSE)
#curve_arg <- read.csv("kleb_arg_subiso_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.bi Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
p1 <- make_mass_plot_region(
  curve_arg,
  xlab = "AMR gene frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(curve_arg$f[curve_arg$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("kleb_bsi_ARG_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_arg, pal = region_pal)
p2
# save
ggsave("kleb_bsi_ARG_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.bii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_arg |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_arg, national_df)

curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_arg_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_arg_regional_bhm_summary_cells_only <- kleb_bsi_arg_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_arg_regional_bhm_summary)
print(kleb_bsi_arg_regional_bhm_summary_cells_only)
#View(kleb_bsi_arg_regional_bhm_summary)
#View(kleb_bsi_arg_regional_bhm_summary_cells_only)

# save
write.csv(kleb_bsi_arg_regional_bhm_summary_cells_only, "rarefaction/kleb_bsi_arg_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(kleb_bsi_arg_regional_bhm_summary, "rarefaction/kleb_bsi_arg_regional_bhm_summary.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 5.2a E. coli Plasmids - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare E. coli plasmids df
# master sample list (one row per sample)
all_samples <- ecoli_bsi_samples_metadata |> distinct(sequencing_id)   

# detected pairs from PLING (1 if detected)
detected <- ecoli_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)
#View(detected)
nrow(detected)

# full grid of sample x feature using master sample list and the set of observed features
all_features <- detected |> pull(community_subcommunity) |> unique()
full_grid <- tidyr::expand_grid(sample = all_samples$sequencing_id,
                                community_subcommunity = all_features)
# left join detections onto the full grid and fill NAs with 0
presence_long <- full_grid |>
  left_join(detected, by = c("sample", "community_subcommunity")) |>
  mutate(presence = if_else(is.na(presence), 0L, presence))
# then pivot to wide
ecoli_bsi_pling_df <- presence_long |>
  pivot_wider(names_from = community_subcommunity,
              values_from = presence,
              values_fill = 0L) |>
  ungroup()

# add region data
ecoli_bsi_pling_df <- ecoli_bsi_pling_df |>
  left_join(ecoli_bsi_samples_metadata |> select(sequencing_id, region), by = c("sample" = "sequencing_id")) |>
  select(sample, region, everything())

#View(ecoli_bsi_pling_df)
length(unique(ecoli_bsi_pling_df$sample))# 1471
table(ecoli_bsi_pling_df$region)# 1471
sum(ecoli_bsi_pling_df[,3:827])# 3361

gene_cols <- setdiff(names(ecoli_bsi_pling_df), c("sample", "region"))
prep_pling <- prep_subiso_data(
  df = ecoli_bsi_pling_df,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols
)

fit_u_pling <- fit_shared_tau_subiso(prep_pling, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r_pling <- fit_region_tau_subiso(prep_pling, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# check diagnostics
# model summaries and diagnostics
str(fit_u_pling)
str(fit_u_pling$prep)
str(fit_u_pling$model)

str(fit_u_pling$fit)
print(fit_u_pling$fit$summary(), n=300)
print(fit_r_pling$fit$summary(), n=300)

fit_u_pling$fit$cmdstan_summary()
fit_r_pling$fit$cmdstan_summary()

fit_u_pling$fit$cmdstan_diagnose()
fit_r_pling$fit$cmdstan_diagnose()

fit_u_pling$fit$diagnostic_summary()
fit_r_pling$fit$diagnostic_summary()


fit_u_pling$fit$sampler_diagnostics()
fit_r_pling$fit$sampler_diagnostics()

# save
saveRDS(prep_pling, "ecoli_pling_prep_pling_subiso.rds")
saveRDS(fit_u_pling, "ecoli_pling_fit_u_pling_subiso.rds")
saveRDS(fit_r_pling, "ecoli_pling_fit_r_pling_subiso.rds")

# call cross-validation function
kfold_compare_subiso <- kfold_compare_subiso_models(
  prep = prep_pling,
  K_folds = min(5, nrow(prep_pling$y)),
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  seed = 2026,
  n_tau_mc = 50
)

# check which better
kfold_compare_subiso$summary
kfold_compare_subiso$per_region
kfold_compare_subiso$fold_id
#save
saveRDS(kfold_compare_subiso, "ecoli_bsi_pling_by_region_kfold_compare.rds")


# calculate cumulative mass grid
f_grid <- c(  seq(0.0001, 0.001, by = 0.00001),   seq(0.0011, 0.01, by = 0.0001),   seq(0.011, 1, by = 0.001))

curve_pling <- compute_subiso_mass_curves(
  fit_obj = fit_r_pling,
  f_grid = f_grid,
  n_post_draws = 1000)

#View(curve_pling)
colnames(curve_pling)
table(curve_pling$region)
#save
write.csv(curve_pling, "ecoli_pling_subiso_mass_curves.csv", row.names = FALSE)
#curve_pling <- read.csv("ecoli_pling_subiso_mass_curves.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2.ai Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
p1 <- make_mass_plot_region(
  curve_pling,
  xlab = "AMR gene frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(curve_pling$f[curve_pling$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("ecoli_bsi_pling_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_pling, pal = region_pal)
p2
# save
ggsave("ecoli_bsi_pling_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2.aii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_pling |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_pling, national_df)

curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.20, 0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_pling_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_pling_regional_bhm_summary_cells_only <- ecoli_bsi_pling_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_pling_regional_bhm_summary)
print(ecoli_bsi_pling_regional_bhm_summary_cells_only)
View(ecoli_bsi_pling_regional_bhm_summary)
View(ecoli_bsi_pling_regional_bhm_summary_cells_only)

# save
write.csv(ecoli_bsi_pling_regional_bhm_summary_cells_only, "rarefaction/ecoli_bsi_pling_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(ecoli_bsi_pling_regional_bhm_summary, "rarefaction/ecoli_bsi_pling_regional_bhm_summary.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 5.2b Klebsiella plasmids - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare Kleb plasmids df
# master sample list (one row per sample)
all_samples <- kleb_bsi_samples_metadata |> distinct(sequencing_id)   

# detected pairs from AMRFinder (1 if detected)
detected <- kleb_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)
nrow(detected)

# full grid of sample x feature using master sample list and the set of observed features
all_features <- detected |> pull(community_subcommunity) |> unique()
full_grid <- tidyr::expand_grid(sample = all_samples$sequencing_id,
                                community_subcommunity = all_features)
# left join detections onto the full grid and fill NAs with 0
presence_long <- full_grid |>
  left_join(detected, by = c("sample", "community_subcommunity")) |>
  mutate(presence = if_else(is.na(presence), 0L, presence))
# then pivot to wide
kleb_bsi_pling_df <- presence_long |>
  pivot_wider(names_from = community_subcommunity,
              values_from = presence,
              values_fill = 0L) |>
  ungroup()

# add region data
kleb_bsi_pling_df <- kleb_bsi_pling_df |>
  left_join(kleb_bsi_samples_metadata |> select(sequencing_id, region), by = c("sample" = "sequencing_id")) |>
  select(sample, region, everything())

#View(kleb_bsi_pling_df)
length(unique(kleb_bsi_pling_df$sample))# 468
table(kleb_bsi_pling_df$region)# 1471

gene_cols <- setdiff(names(kleb_bsi_pling_df), c("sample", "region"))
prep_pling <- prep_subiso_data(
  df = kleb_bsi_pling_df,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols
)

fit_u_pling <- fit_shared_tau_subiso(prep_pling, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r_pling <- fit_region_tau_subiso(prep_pling, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# check diagnostics
# model summaries and diagnostics
str(fit_u_pling)
str(fit_u_pling$prep)
str(fit_u_pling$model)

str(fit_u_pling$fit)
print(fit_u_pling$fit$summary(), n=300)
print(fit_r_pling$fit$summary(), n=300)

fit_u_pling$fit$cmdstan_summary()
fit_r_pling$fit$cmdstan_summary()

fit_u_pling$fit$cmdstan_diagnose()
fit_r_pling$fit$cmdstan_diagnose()

fit_u_pling$fit$diagnostic_summary()
fit_r_pling$fit$diagnostic_summary()


fit_u_pling$fit$sampler_diagnostics()
fit_r_pling$fit$sampler_diagnostics()

# save
saveRDS(prep_pling, "kleb_pling_prep_pling_subiso.rds")
saveRDS(fit_u_pling, "kleb_pling_fit_u_pling_subiso.rds")
saveRDS(fit_r_pling, "kleb_pling_fit_r_pling_subiso.rds")

# call cross-validation function
kfold_compare_subiso <- kfold_compare_subiso_models(
  prep = prep_pling,
  K_folds = min(5, nrow(prep_pling$y)),
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  seed = 2026,
  n_tau_mc = 50
)

# check which better
kfold_compare_subiso$summary
kfold_compare_subiso$per_region
kfold_compare_subiso$fold_id
#save
saveRDS(kfold_compare_subiso, "kleb_bsi_pling_by_region_kfold_compare.rds")


# cal cumulative mass grid
f_grid <- c(  seq(0.0001, 0.001, by = 0.00001),   seq(0.0011, 0.01, by = 0.0001),   seq(0.011, 1, by = 0.001))

curve_pling <- compute_subiso_mass_curves(
  fit_obj = fit_r_pling,
  f_grid = f_grid,
  n_post_draws = 1000)

#View(curve_pling)
colnames(curve_pling)
table(curve_pling$region)

write.csv(curve_pling, "kleb_pling_subiso_mass_curves.csv", row.names = FALSE)
#curve_pling <- read.csv("kleb_pling_subiso_mass_curves.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2.bi Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
p1 <- make_mass_plot_region(
  curve_pling,
  xlab = "AMR gene frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(curve_pling$f[curve_pling$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("kleb_bsi_pling_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_pling, pal = region_pal)
p2
# save
ggsave("kleb_bsi_pling_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2.bii Summary table ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_pling |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_pling, national_df)

curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.20, 0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_pling_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "mean (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(mean_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        mean_txt <- if (is.na(mean_ss)) "NA" else formatC(mean_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(mean_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             mean_ss = mean_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_pling_regional_bhm_summary_cells_only <- kleb_bsi_pling_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_pling_regional_bhm_summary)
print(kleb_bsi_pling_regional_bhm_summary_cells_only)
#View(kleb_bsi_pling_regional_bhm_summary)
#View(kleb_bsi_pling_regional_bhm_summary_cells_only)

# save
write.csv(kleb_bsi_pling_regional_bhm_summary_cells_only, "rarefaction/kleb_bsi_pling_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(kleb_bsi_pling_regional_bhm_summary, "rarefaction/kleb_bsi_pling_regional_bhm_summary.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. PANEL PLOT FOR ALL BY-REGION HIERARCHICAL MODELS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data
ecoli_mlst_curve_df <- read.csv("ecoli_bsi_mlst_mass_curves.csv")
kleb_mlst_curve_df <- read.csv("kleb_bsi_mlst_mass_curves.csv")
ecoli_fastbaps_L3_curve_df <- read.csv("ecoli_bsi_fastbaps_L3_mass_curves.csv")
kleb_fastbaps_L3_curve_df <- read.csv("kleb_bsi_fastbaps_L3_mass_curves.csv")
ecoli_pling_curve_df <- read.csv("ecoli_pling_subiso_mass_curves.csv")
kleb_pling_curve_df <- read.csv("kleb_pling_subiso_mass_curves.csv")
ecoli_arg_curve_df <- read.csv("ecoli_arg_subiso_mass_curves.csv")
kleb_arg_curve_df <- read.csv("kleb_arg_subiso_mass_curves.csv")

# plot functions and palette defined above
# colourBrewerSet2 + Dark2 themed
region_pal <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
)
region_pal <- rev(region_pal)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Cumulative mass plots for ecoli and kleb by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
p1 <- make_mass_plot_region(
  ecoli_mlst_curve_df,
  xlab = "MLST frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(ecoli_mlst_curve_df$f[ecoli_mlst_curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1

p2 <- make_mass_plot_region(
  kleb_mlst_curve_df,
  xlab = "MLST frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(kleb_mlst_curve_df$f[kleb_mlst_curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p2

p3 <- make_mass_plot_region(
  ecoli_fastbaps_L3_curve_df,
  xlab = "FastBAPS cluster frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(ecoli_fastbaps_L3_curve_df$f[ecoli_fastbaps_L3_curve_df$f > 0], na.rm = TRUE)
)
p3

p4 <- make_mass_plot_region(
  kleb_fastbaps_L3_curve_df,
  xlab = "FastBAPS cluster frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(kleb_fastbaps_L3_curve_df$f[kleb_fastbaps_L3_curve_df$f > 0], na.rm = TRUE)
)
p4

p5 <- make_mass_plot_region(
  ecoli_pling_curve_df,
  xlab = "Plasmid subcommunity frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(ecoli_pling_curve_df$f[ecoli_pling_curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p5

p6 <- make_mass_plot_region(
  kleb_pling_curve_df,
  xlab = "Plasmid subcommunity frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(kleb_pling_curve_df$f[kleb_pling_curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p6

p7 <- make_mass_plot_region(
  ecoli_arg_curve_df,
  xlab = "AMR gene frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(ecoli_arg_curve_df$f[ecoli_arg_curve_df$f > 0], na.rm = TRUE)
)
p7

p8 <- make_mass_plot_region(
  kleb_arg_curve_df,
  xlab = "AMR gene frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(kleb_arg_curve_df$f[kleb_arg_curve_df$f > 0], na.rm = TRUE)
)
p8

# combine plots:
combined_figure <-
  (p1 | p2) /
  (p3 | p4) /
  (p5 | p6) /
  (p7 | p8) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined_figure
# save
ggsave("combined_4x2_hierarchical_cumulative_mass_curves_by_region.png",
       plot = combined_figure, width = 10, height = 15, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Sample size plots for 95% confidence level for ecoli and kleb ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
p1 <- make_cov_plot_region(ecoli_mlst_curve_df, pal = region_pal)
p1
p2 <- make_cov_plot_region(kleb_mlst_curve_df, pal = region_pal)
p2
p3 <- make_cov_plot_region(ecoli_fastbaps_L3_curve_df, pal = region_pal)
p3
p4 <- make_cov_plot_region(kleb_fastbaps_L3_curve_df, pal = region_pal)
p4
p5 <- make_cov_plot_region(ecoli_pling_curve_df, pal = region_pal)
p5
p6 <- make_cov_plot_region(kleb_pling_curve_df, pal = region_pal)
p6
p7 <- make_cov_plot_region(ecoli_arg_curve_df, pal = region_pal)
p7
p8 <- make_cov_plot_region(kleb_arg_curve_df, pal = region_pal)
p8

# combine plots:
combined_figure <-
  (p1 | p2) /
  (p3 | p4) /
  (p5 | p6) /
  (p7 | p8) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_figure
# save
ggsave("combined_4x2_ss_vs_coverage_panel_by_region.png",
       plot = combined_figure, width = 10, height = 15, units = "in", dpi = 300)
