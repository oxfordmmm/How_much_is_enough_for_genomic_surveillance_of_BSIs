# app.R
# Bayesian sampling-frame explorer for genomic surveillance
# Statistical functions are aligned to the publication code supplied with this app.

suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(plotly)
  library(readr)
  library(htmltools)
  library(matrixStats)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

safe_num <- function(x, default = NA_real_) {
  out <- suppressWarnings(as.numeric(x))
  if (length(out) == 0 || is.na(out[1])) default else out[1]
}

validate_prob <- function(x, label) {
  x <- safe_num(x)
  if (!is.finite(x) || x <= 0 || x >= 1) {
    stop(sprintf("%s must be strictly between 0 and 1.", label))
  }
  x
}

parse_alpha_vector <- function(text, feature_names, label = "alpha_named prior") {
  txt <- trimws(text %||% "")
  if (txt == "" || is.na(txt)) {
    out <- rep(1, length(feature_names))
    names(out) <- feature_names
    return(out)
  }

  vals <- suppressWarnings(as.numeric(unlist(strsplit(txt, "[ ,;\\n\\t]+"))))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 1) {
    if (vals < 0) stop(sprintf("%s must be non-negative.", label))
    out <- rep(vals, length(feature_names))
    names(out) <- feature_names
    return(out)
  }
  if (length(vals) != length(feature_names)) {
    stop(sprintf("%s must be either one number or one value per observed feature.", label))
  }
  if (any(vals < 0)) stop(sprintf("%s must be non-negative.", label))
  stats::setNames(vals, feature_names)
}

read_csv_source <- function(file, text_input) {
  if (!is.null(file) && nzchar(file$name)) {
    readr::read_csv(file$datapath, show_col_types = FALSE, progress = FALSE)
  } else {
    txt <- trimws(text_input %||% "")
    if (!nzchar(txt)) stop("Provide either a CSV upload or pasted CSV text.")
    readr::read_csv(I(txt), show_col_types = FALSE, progress = FALSE)
  }
}

read_isolate_table <- function(file, text_input, feature_col, count_col) {
  df <- read_csv_source(file, text_input)
  if (!(feature_col %in% names(df))) {
    stop(sprintf("Feature column '%s' was not found in the CSV.", feature_col))
  }
  if (!(count_col %in% names(df))) {
    stop(sprintf("Count column '%s' was not found in the CSV.", count_col))
  }

  df <- df |>
    mutate(
      feature_value = as.character(.data[[feature_col]]),
      count_value = suppressWarnings(as.numeric(.data[[count_col]]))
    ) |>
    filter(!is.na(feature_value), nzchar(trimws(feature_value)),
           !is.na(count_value), is.finite(count_value))

  if (nrow(df) == 0) stop("No usable isolate-level rows were found.")
  if (any(df$count_value < 0)) stop("Counts must be non-negative.")

  df |>
    transmute(
      !!feature_col := feature_value,
      count = count_value
    ) |>
    group_by(.data[[feature_col]]) |>
    summarise(count = sum(count), .groups = "drop") |>
    arrange(.data[[feature_col]])
}

infer_id_column <- function(df) {
  if ("sample" %in% names(df)) return("sample")
  non_numeric <- names(df)[vapply(df, function(x) !is.numeric(x), logical(1))]
  if (length(non_numeric) >= 1) return(non_numeric[1])
  names(df)[1]
}

read_subiso_table <- function(file, text_input) {
  df <- read_csv_source(file, text_input)
  if (ncol(df) < 2) stop("The sub-isolate CSV must contain an isolate/ID column plus at least one gene/feature column.")

  id_col <- infer_id_column(df)
  gene_cols <- setdiff(names(df), id_col)
  if (length(gene_cols) < 1) stop("No gene/feature columns were found after excluding the isolate ID column.")

  numeric_df <- df |>
    mutate(across(all_of(gene_cols), ~ suppressWarnings(as.numeric(.x))))
  gene_matrix <- as.matrix(numeric_df[, gene_cols, drop = FALSE])
  if (any(!is.finite(gene_matrix))) {
    bad <- gene_cols[colSums(!is.finite(gene_matrix)) > 0]
    stop(sprintf("Gene/feature columns must contain numeric values. Problem columns: %s",
                 paste(bad, collapse = ", ")))
  }

  out <- numeric_df |>
    mutate(across(all_of(gene_cols), ~ as.integer(.x > 0)))

  list(data = out, id_col = id_col, gene_cols = gene_cols)
}

# -----------------------------------------------------------------------------
# Publication Framework 1: isolate-level features
# -----------------------------------------------------------------------------

fit_crp_theta_bayes <- function(counts,
                                theta_grid = exp(seq(log(1e-4), log(1e4), length.out = 5000)),
                                prior_shape = 1,
                                prior_rate = 1,
                                n_draws = 5000,
                                seed = 2026) {
  stopifnot(is.numeric(counts), all(counts >= 0), length(counts) >= 1)
  counts <- counts[counts > 0]
  N <- sum(counts)
  K <- length(counts)

  log_lik <- K * log(theta_grid) + lgamma(theta_grid) - lgamma(theta_grid + N)
  log_prior <- dgamma(theta_grid, shape = prior_shape, rate = prior_rate, log = TRUE)
  log_post <- log_lik + log_prior
  log_post <- log_post - max(log_post)
  w <- exp(log_post)
  w <- w / sum(w)

  set.seed(seed)
  theta_draws <- sample(theta_grid, size = n_draws, replace = TRUE, prob = w)
  novelty_prob_draws <- theta_draws / (N + theta_draws)

  summary_df <- data.frame(
    quantity = c("theta", "novelty_prob_next"),
    mean = c(mean(theta_draws), mean(novelty_prob_draws)),
    median = c(stats::median(theta_draws), stats::median(novelty_prob_draws)),
    q2.5 = c(stats::quantile(theta_draws, 0.025), stats::quantile(novelty_prob_draws, 0.025)),
    q97.5 = c(stats::quantile(theta_draws, 0.975), stats::quantile(novelty_prob_draws, 0.975))
  )

  list(
    N = N,
    K = K,
    theta_grid = theta_grid,
    posterior_weights = w,
    theta_draws = theta_draws,
    novelty_prob_draws = novelty_prob_draws,
    summary_df = summary_df
  )
}

set_alpha_novel_sum_from_crp <- function(alpha_named, novelty_prob_hat) {
  alpha_obs_total <- sum(alpha_named, na.rm = TRUE)
  alpha_obs_total * novelty_prob_hat / (1 - novelty_prob_hat)
}

isolate_bayes_dirichlet <- function(df,
                                   feature_col = "mlst_profile",
                                   alpha_named = NULL,
                                   alpha_novel_sum = 1,
                                   alpha_novel_num = 1,
                                   B = 1000,
                                   seed = 2026) {
  stopifnot(is.data.frame(df))
  stopifnot(is.character(feature_col) && length(feature_col) == 1)
  stopifnot(all(c(feature_col, "count") %in% colnames(df)))
  stopifnot(is.numeric(alpha_novel_sum), length(alpha_novel_sum) == 1, alpha_novel_sum >= 0)
  stopifnot(is.numeric(alpha_novel_num), length(alpha_novel_num) == 1, alpha_novel_num >= 0)
  stopifnot(B >= 1)
  set.seed(seed)

  alpha_novel_sum <- as.numeric(alpha_novel_sum)
  alpha_novel_num <- as.integer(alpha_novel_num)

  df <- df |>
    group_by(.data[[feature_col]]) |>
    summarise(count = sum(.data[["count"]]), .groups = "drop")

  K <- nrow(df)
  names_counts <- as.character(df[[feature_col]])
  n_k <- df$count
  N <- sum(n_k)

  if (is.null(alpha_named)) {
    alpha_named <- rep(1, K)
    names(alpha_named) <- names_counts
  } else {
    if (length(alpha_named) != K) stop("alpha_named must have the same length as unique features in df.")
    if (!is.null(names(alpha_named)) && all(names_counts %in% names(alpha_named))) {
      alpha_named <- alpha_named[names_counts]
    } else {
      names(alpha_named) <- names_counts
    }
  }
  alpha_named <- as.numeric(alpha_named)

  add_novel <- (alpha_novel_sum > 0) && (alpha_novel_num > 0)
  if (add_novel) {
    alpha_novel_each <- alpha_novel_sum / alpha_novel_num
    alpha_novel_vec <- rep(alpha_novel_each, alpha_novel_num)
    novel_names <- paste0("NOVEL_", seq_len(alpha_novel_num))
  } else {
    alpha_novel_each <- 0
    alpha_novel_vec <- numeric(0)
    novel_names <- character(0)
  }

  shapes_obs <- n_k + alpha_named
  shapes_all <- c(shapes_obs, alpha_novel_vec)
  out_names <- c(names_counts, novel_names)

  G <- matrix(
    rgamma(B * length(shapes_all), shape = shapes_all, rate = 1),
    nrow = B,
    ncol = length(shapes_all),
    byrow = TRUE
  )
  draws <- G / rowSums(G)
  colnames(draws) <- out_names

  summary_df <- data.frame(
    feature = out_names,
    median = matrixStats::colMedians(draws),
    sd = matrixStats::colSds(draws),
    q2.5 = matrixStats::colQuantiles(draws, probs = 0.025),
    q97.5 = matrixStats::colQuantiles(draws, probs = 0.975),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    draws = draws,
    summary_df = summary_df,
    alpha_novel_sum = alpha_novel_sum,
    alpha_novel_num = alpha_novel_num,
    alpha_novel_each = alpha_novel_each,
    N = N,
    K = K
  )
}

prep_bootstrap_draws <- function(draws) {
  x <- as.matrix(draws)
  x_sorted <- t(apply(x, 1, sort))
  cs <- t(apply(x_sorted, 1, cumsum))
  list(x_sorted = x_sorted, cumsums = cs, totals = rowSums(x_sorted))
}

default_f_grid <- function() {
  n_grid <- c(
    seq(0, 10000, by = 1),
    seq(10005, 50000, by = 5),
    seq(50010, 100000, by = 10)
  )
  1 - ((1 - 0.99)^(1 / n_grid))
}

f_grid_99 <- default_f_grid()

compute_species_richness_curve <- function(prep, f_grid = NULL) {
  x_sorted <- prep$x_sorted
  B <- nrow(x_sorted)
  K <- ncol(x_sorted)
  if (is.null(f_grid)) f_grid <- default_f_grid()
  F <- length(f_grid)
  counts <- matrix(0L, nrow = B, ncol = F)

  for (i in seq_len(B)) {
    idx <- findInterval(f_grid, x_sorted[i, ], left.open = TRUE, rightmost.closed = TRUE)
    counts[i, ] <- K - idx
  }

  prop <- counts / K
  q_counts <- matrixStats::colQuantiles(counts, probs = c(0.025, 0.975))
  q_prop_025 <- matrixStats::colQuantiles(prop, probs = 0.025)
  q_prop_975 <- matrixStats::colQuantiles(prop, probs = 0.975)

  data.frame(
    f = f_grid,
    median_species_count = matrixStats::colMedians(counts),
    sd_species_count = matrixStats::colSds(counts),
    q2.5 = q_counts[, 1],
    q97.5 = q_counts[, 2],
    median_species_proportion = matrixStats::colMedians(prop),
    sd_species_proportion = matrixStats::colSds(prop),
    q2.5_species_proportion = q_prop_025,
    q97.5_species_proportion = q_prop_975
  )
}

compute_mass_curve <- function(prep, f_grid = NULL) {
  x_sorted <- prep$x_sorted
  cs <- prep$cumsums
  totals <- prep$totals
  B <- nrow(x_sorted)
  if (is.null(f_grid)) f_grid <- default_f_grid()
  F <- length(f_grid)
  masses <- matrix(0, nrow = B, ncol = F)

  for (i in seq_len(B)) {
    idx <- findInterval(f_grid, x_sorted[i, ], left.open = TRUE, rightmost.closed = TRUE)
    cs_i <- c(0, cs[i, ])
    masses[i, ] <- totals[i] - cs_i[idx + 1]
  }

  q <- matrixStats::colQuantiles(masses, probs = c(0.025, 0.975))
  data.frame(
    f = f_grid,
    median_mass = matrixStats::colMedians(masses),
    sd_mass = matrixStats::colSds(masses),
    q2.5 = q[, 1],
    q97.5 = q[, 2]
  )
}

# -----------------------------------------------------------------------------
# Publication Framework 2: sub-isolate-level features
# -----------------------------------------------------------------------------

estimate_subiso_u_hat <- function(iso_df, gene_cols) {
  Z <- as.matrix(iso_df[, gene_cols, drop = FALSE])
  Z <- (Z > 0) * 1L
  feature_counts <- colSums(Z, na.rm = TRUE)
  total_feature_occurrences <- sum(feature_counts)
  if (!is.finite(total_feature_occurrences) || total_feature_occurrences <= 0) return(0)
  singleton_features <- sum(feature_counts == 1)
  singleton_features / total_feature_occurrences
}

subisolate_bayes_beta <- function(iso_df,
                                  gene_cols,
                                  B = 2000,
                                  beta_named = NULL,
                                  beta_novel_sum = 0,
                                  beta_novel_num = 0,
                                  seed = 2026) {
  stopifnot(is.data.frame(iso_df))
  stopifnot(length(gene_cols) >= 1)
  stopifnot(all(gene_cols %in% colnames(iso_df)))
  stopifnot(B >= 1)
  stopifnot(is.numeric(beta_novel_sum), length(beta_novel_sum) == 1, beta_novel_sum >= 0)
  stopifnot(is.numeric(beta_novel_num), length(beta_novel_num) == 1, beta_novel_num >= 0)
  set.seed(seed)

  beta_novel_num <- as.integer(beta_novel_num)
  N <- nrow(iso_df)
  Z <- as.matrix(iso_df[, gene_cols, drop = FALSE])
  Z <- (Z > 0) * 1L
  x_obs <- colSums(Z, na.rm = TRUE)
  K <- length(x_obs)
  obs_names <- colnames(Z)

  if (is.null(beta_named)) {
    beta_named <- rep(1, K)
    names(beta_named) <- obs_names
  } else {
    if (length(beta_named) != K) stop("beta_named must have the same length as gene_cols.")
    if (!is.null(names(beta_named)) && all(obs_names %in% names(beta_named))) {
      beta_named <- beta_named[obs_names]
    } else {
      names(beta_named) <- obs_names
    }
  }
  beta_named <- as.numeric(beta_named)

  add_novel <- (beta_novel_sum > 0) && (beta_novel_num > 0)
  if (add_novel) {
    beta_novel_each <- beta_novel_sum / beta_novel_num
    beta_novel_vec <- rep(beta_novel_each, beta_novel_num)
    novel_names <- paste0("NOVEL_", seq_len(beta_novel_num))
  } else {
    beta_novel_each <- 0
    beta_novel_vec <- numeric(0)
    novel_names <- character(0)
  }

  feature_names <- c(obs_names, novel_names)
  total_features <- length(feature_names)
  p_draws <- matrix(NA_real_, nrow = B, ncol = total_features)
  colnames(p_draws) <- feature_names

  for (b in seq_len(B)) {
    p_obs <- stats::rbeta(
      K,
      shape1 = x_obs + beta_named,
      shape2 = (N - x_obs) + beta_named
    )
    if (add_novel) {
      p_novel <- stats::rbeta(
        beta_novel_num,
        shape1 = beta_novel_vec,
        shape2 = N + beta_novel_vec
      )
      p_draws[b, ] <- c(p_obs, p_novel)
    } else {
      p_draws[b, ] <- p_obs
    }
  }

  summary_df <- data.frame(
    feature = feature_names,
    median_prevalence = matrixStats::colMedians(p_draws, na.rm = TRUE),
    sd_prevalence = matrixStats::colSds(p_draws, na.rm = TRUE),
    q2.5_prevalence = matrixStats::colQuantiles(p_draws, probs = 0.025, na.rm = TRUE),
    q97.5_prevalence = matrixStats::colQuantiles(p_draws, probs = 0.975, na.rm = TRUE),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    p_draws = p_draws,
    summary_df = summary_df,
    beta_novel_sum = beta_novel_sum,
    beta_novel_num = beta_novel_num,
    beta_novel_each = beta_novel_each
  )
}

get_p_matrix <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) return(as.matrix(x))
  if (is.list(x) && !is.null(x$p_draws)) return(as.matrix(x$p_draws))
  stop("Could not find p_draws in input.")
}

default_subiso_f_grid <- function() f_grid_99

prep_subiso_curve <- function(df) {
  stopifnot(is.data.frame(df) || is.matrix(df))
  p_mat <- as.matrix(df)
  p_sorted <- t(apply(p_mat, 1, sort))
  list(
    x_sorted = p_sorted,
    cumsums = matrixStats::rowCumsums(p_sorted),
    totals = rowSums(p_mat, na.rm = TRUE),
    K = ncol(p_mat)
  )
}

compute_species_richness_subiso <- function(x, f_grid = NULL) {
  p_mat <- get_p_matrix(x)
  B <- nrow(p_mat)
  K <- ncol(p_mat)
  if (is.null(f_grid)) f_grid <- default_subiso_f_grid()
  F <- length(f_grid)
  p_sorted <- t(apply(p_mat, 1, sort))
  rich_mat <- matrix(0L, nrow = B, ncol = F)

  for (i in seq_len(B)) {
    idx <- findInterval(f_grid, p_sorted[i, ], left.open = TRUE, rightmost.closed = TRUE)
    rich_mat[i, ] <- K - idx
  }

  prop_mat <- rich_mat / K
  q_count <- matrixStats::colQuantiles(rich_mat, probs = c(0.025, 0.975), na.rm = TRUE)
  q_prop <- matrixStats::colQuantiles(prop_mat, probs = c(0.025, 0.975), na.rm = TRUE)

  tibble::tibble(
    f = f_grid,
    median_species_count = matrixStats::colMedians(rich_mat, na.rm = TRUE),
    q2.5_species_count = q_count[, 1],
    q97.5_species_count = q_count[, 2],
    median_species_proportion = matrixStats::colMedians(prop_mat, na.rm = TRUE),
    q2.5_species_proportion = q_prop[, 1],
    q97.5_species_proportion = q_prop[, 2]
  ) |>
    filter(is.finite(f), !is.na(median_species_count))
}

compute_mass_curve_subiso <- function(x, f_grid = NULL) {
  p_mat <- get_p_matrix(x)
  B <- nrow(p_mat)
  if (is.null(f_grid)) f_grid <- default_subiso_f_grid()
  F <- length(f_grid)
  p_sorted <- t(apply(p_mat, 1, sort))
  cs <- matrixStats::rowCumsums(p_sorted)
  totals <- rowSums(p_mat, na.rm = TRUE)
  mass_mat <- matrix(NA_real_, nrow = B, ncol = F)

  for (i in seq_len(B)) {
    idx <- findInterval(f_grid, p_sorted[i, ], left.open = TRUE, rightmost.closed = TRUE)
    cs_i <- c(0, cs[i, ])
    if (!is.na(totals[i]) && totals[i] > 0) {
      mass_mat[i, ] <- (totals[i] - cs_i[idx + 1]) / totals[i]
    }
  }

  q_mass <- matrixStats::colQuantiles(mass_mat, probs = c(0.025, 0.975), na.rm = TRUE)
  tibble::tibble(
    f = f_grid,
    median_mass = matrixStats::colMedians(mass_mat, na.rm = TRUE),
    q2.5 = q_mass[, 1],
    q97.5 = q_mass[, 2]
  ) |>
    filter(is.finite(f), !is.na(median_mass))
}

summarise_subiso_posterior <- function(bbs, gene_cols, actual_vec, beta_named,
                                       beta_novel_sum, beta_novel_num,
                                       dataset_label = NA_character_) {
  p_mat <- get_p_matrix(bbs)
  if (!is.null(colnames(p_mat)) && all(gene_cols %in% colnames(p_mat))) {
    p_mat <- p_mat[, gene_cols, drop = FALSE]
  } else {
    colnames(p_mat) <- gene_cols
  }
  est_med <- apply(p_mat, 2, median, na.rm = TRUE)
  est_lo <- apply(p_mat, 2, quantile, probs = 0.025, na.rm = TRUE, names = FALSE)
  est_hi <- apply(p_mat, 2, quantile, probs = 0.975, na.rm = TRUE, names = FALSE)

  tibble::tibble(
    feature = gene_cols,
    actual = as.numeric(actual_vec[gene_cols]),
    estimate = as.numeric(est_med),
    lo = as.numeric(est_lo),
    hi = as.numeric(est_hi),
    beta_named = beta_named,
    beta_novel_sum = beta_novel_sum,
    beta_novel_num = beta_novel_num,
    dataset = dataset_label
  )
}

# -----------------------------------------------------------------------------
# Publication-style sample-size inversion / summaries
# -----------------------------------------------------------------------------

required_n_from_draw <- function(p, target, mode = c("coverage", "richness"), confidence = 0.95) {
  mode <- match.arg(mode)
  p <- p[is.finite(p) & p >= 0]
  if (length(p) == 0) return(NA_real_)
  total <- sum(p)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  p <- sort(p / total, decreasing = TRUE)

  cum <- if (mode == "coverage") cumsum(p) else seq_along(p) / length(p)
  idx <- which(cum >= target)[1]
  if (is.na(idx)) return(NA_real_)
  f_star <- p[idx]
  if (!is.finite(f_star) || f_star <= 0) return(NA_real_)
  if (f_star >= 1) return(1)
  n_star <- log(1 - confidence) / log(1 - f_star)
  ceiling(n_star)
}

summarise_required_n <- function(draws, target_coverage, target_richness, confidence) {
  p_mat <- as.matrix(draws)
  cov_n <- apply(p_mat, 1, required_n_from_draw,
                 target = target_coverage, mode = "coverage", confidence = confidence)
  rich_n <- apply(p_mat, 1, required_n_from_draw,
                  target = target_richness, mode = "richness", confidence = confidence)
  list(
    coverage = cov_n,
    richness = rich_n,
    coverage_summary = c(
      median = median(cov_n, na.rm = TRUE),
      q2.5 = quantile(cov_n, 0.025, na.rm = TRUE, names = FALSE),
      q97.5 = quantile(cov_n, 0.975, na.rm = TRUE, names = FALSE)
    ),
    richness_summary = c(
      median = median(rich_n, na.rm = TRUE),
      q2.5 = quantile(rich_n, 0.025, na.rm = TRUE, names = FALSE),
      q97.5 = quantile(rich_n, 0.975, na.rm = TRUE, names = FALSE)
    )
  )
}

format_n <- function(x) {
  if (!is.finite(x)) return("not estimable")
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------

plot_obs_vs_pred <- function(feature_summary, title) {
  ggplot(feature_summary, aes(
    x = actual, y = estimate,
    text = paste0(
      "Feature: ", feature,
      "<br>Observed: ", signif(actual, 4),
      "<br>Posterior median: ", signif(estimate, 4),
      "<br>95% CI: [", signif(lo, 4), ", ", signif(hi, 4), "]"
    )
  )) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, alpha = 0.45) +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Observed prevalence", y = "Posterior median prevalence", title = title) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

plot_mass_curve <- function(curve, title) {
  ggplot(curve, aes(x = f, y = median_mass, text = paste0(
    "Threshold f: ", signif(f, 4),
    "<br>Mass: ", round(median_mass, 4),
    "<br>95% CI: [", round(q2.5, 4), ", ", round(q97.5, 4), "]"
  ))) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.15) +
    geom_line(linewidth = 0.8) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Threshold f", y = "Population mass above f", title = title) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

plot_richness_curve <- function(curve, title) {
  ggplot(curve, aes(x = f, y = median_species_proportion, text = paste0(
    "Threshold f: ", signif(f, 4),
    "<br>Richness proportion: ", round(median_species_proportion, 4),
    "<br>95% CI: [", round(q2.5_species_proportion, 4), ", ", round(q97.5_species_proportion, 4), "]"
  ))) +
    geom_ribbon(aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.15) +
    geom_line(linewidth = 0.8) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Threshold f", y = "Feature richness proportion", title = title) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

plot_required_hist <- function(x, title) {
  df <- tibble::tibble(n_required = as.numeric(x)) |>
    filter(is.finite(n_required))
  ggplot(df, aes(n_required)) +
    geom_histogram(bins = 30) +
    labs(title = title, x = "Required sample size", y = "Posterior draw count") +
    theme_minimal(base_size = 11)
}

plotly_from_gg <- function(p) {
  plotly::ggplotly(p, tooltip = "text") |>
    plotly::layout(legend = list(orientation = "h", x = 0.1, y = -0.1))
}

# -----------------------------------------------------------------------------
# Templates
# -----------------------------------------------------------------------------

make_template_isolate <- function(path) {
  readr::write_csv(tibble::tibble(
    mlst_profile = c("A01", "A02", "A03"),
    count = c(10, 6, 1)
  ), path)
}

make_template_subiso <- function(path) {
  readr::write_csv(tibble::tibble(
    sample = c("S1", "S2", "S3"),
    geneA = c(1, 0, 1),
    geneB = c(0, 1, 1),
    geneC = c(1, 0, 0)
  ), path)
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

citation_text <- paste0(
  "Nagy et al., 2026. <strong>How much is enough? Optimising sampling frames for genomic surveillance of ",
  "<em>Escherichia coli</em> and <em>Klebsiella</em> spp. bloodstream infections &ndash; a retrospective study.</strong>"
)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("\
      .app-intro { color:#6f6f6f; font-style:italic; line-height:1.4; margin-bottom:14px; }\
      .param-help { color:#777; font-size:.85em; line-height:1.25; margin-top:-6px; margin-bottom:8px; }\
      .summary-box { padding:12px; border:1px solid #e5e5e5; border-radius:6px; background:#fff; line-height:1.45; }\
      .status-box { padding:8px 10px; background:#f7f7f7; border:1px solid #e5e5e5; border-radius:6px; margin-bottom:10px; }\
      .error-box { padding:8px 10px; background:#fff2f2; border:1px solid #f0c0c0; border-radius:6px; margin-bottom:10px; color:#9b1c1c; }\
      .detected-box { padding:8px 10px; background:#f7fafc; border:1px solid #dce6ef; border-radius:6px; margin-top:6px; margin-bottom:10px; }\
      .analysis-status { display:flex; align-items:center; gap:9px; padding:8px 10px; margin-bottom:10px; background:#f7f7f7; border:1px solid #e5e5e5; border-radius:6px; }\
      .analysis-spinner { width:16px; height:16px; border:2px solid #d9d9d9; border-top-color:#337ab7; border-radius:50%; animation:spin 0.8s linear infinite; flex:0 0 auto; }\
      @keyframes spin { to { transform:rotate(360deg); } }\
    "))
  ),
  titlePanel("How much is enough? Genomic surveillance sampling-frame estimator for bacterial pathogens"),
  div(class = "app-intro",
      "Bayesian analysis of isolate-level and sub-isolate-level genomic features to estimate genomic surveillance sample sizes for bacterial pathogens."),

  tabsetPanel(
    id = "feature_tabs",

    tabPanel(
      "Isolate-level features",
      sidebarLayout(
        sidebarPanel(
          radioButtons("input_mode", "Input mode",
                       choices = c("Upload CSV" = "upload", "Paste CSV text" = "paste"), selected = "upload"),
          div(class = "param-help",
              "CSV format: one feature column and one count column. The default feature column is mlst_profile."),
          fileInput("file", "Upload CSV", accept = ".csv"),
          textAreaInput("csv_text", "Paste CSV",
                        placeholder = "mlst_profile,count\nA01,10\nA02,6\nA03,1", rows = 8),
          downloadButton("download_template_isolate", "Download isolate template CSV"),
          textInput("feature_col", "Feature column", value = "mlst_profile"),
          textInput("count_col", "Count column", value = "count"),
          tags$hr(),
          numericInput("alpha_prior", "alpha prior", value = 1, min = 0, step = 0.1),
          div(class = "param-help",
              "A common prior is applied to every observed and novel feature. Novel prior mass is estimated from the CRP and then rounded upward before deriving the number of novel features."),
          numericInput("B", "Posterior draws", value = 2000, min = 100, step = 100),
          numericInput("seed", "Seed", value = 2026, step = 1),
          tags$hr(),
          numericInput("target_coverage", "Target sample coverage", value = 0.80, min = 0.001, max = 0.999, step = 0.01),
          numericInput("target_richness", "Target feature richness", value = 0.80, min = 0.001, max = 0.999, step = 0.01),
          numericInput("confidence_level", "Detection confidence", value = 0.95, min = 0.001, max = 0.999, step = 0.01),
          actionButton("run", "Run analysis", class = "btn-primary"),
          tags$hr(),
          downloadButton("download_csv", "Download summary CSV"),
          downloadButton("download_plots", "Download plots PDF")
        ),
        mainPanel(
          uiOutput("isolate_status"),
          fluidRow(
            column(5, htmlOutput("isolate_summary")),
            column(7, plotlyOutput("isolate_obs_pred", height = 420))
          ),
          fluidRow(
            column(6, plotlyOutput("isolate_mass", height = 360)),
            column(6, plotlyOutput("isolate_richness", height = 360))
          ),
          fluidRow(
            column(6, plotlyOutput("isolate_cov_req", height = 320)),
            column(6, plotlyOutput("isolate_rich_req", height = 320))
          ),
          tags$hr(),
          div(class = "app-intro", HTML(paste0("<strong>Citation:</strong> ", citation_text)))
        )
      )
    ),

    tabPanel(
      "Sub-isolate-level features",
      sidebarLayout(
        sidebarPanel(
          radioButtons("input_mode_sub", "Input mode",
                       choices = c("Upload CSV" = "upload", "Paste CSV text" = "paste"), selected = "upload"),
          div(class = "param-help",
              "CSV format: one isolate/ID column followed by gene columns. Gene names are extracted automatically from the remaining column headers; no feature-column field is required."),
          fileInput("file_sub", "Upload CSV", accept = ".csv"),
          textAreaInput("csv_text_sub", "Paste CSV",
                        placeholder = "sample,geneA,geneB,geneC\nS1,1,0,1\nS2,0,1,0\nS3,1,1,0", rows = 8),
          downloadButton("download_template_sub", "Download sub-isolate template CSV"),
          uiOutput("subiso_detected_ui"),
          tags$hr(),
          numericInput("beta_prior", "beta prior", value = 1, min = 0, step = 0.1),
          div(class = "param-help",
              "A common beta prior is applied to every observed and novel feature. The number of novel features is estimated from the supplied presence/absence data and the resulting novel prior mass is beta_novel_num × beta_prior."),
          numericInput("B_sub", "Posterior draws", value = 2000, min = 100, step = 100),
          numericInput("seed_sub", "Seed", value = 2026, step = 1),
          tags$hr(),
          numericInput("target_coverage_sub", "Target sample coverage", value = 0.80, min = 0.001, max = 0.999, step = 0.01),
          numericInput("target_richness_sub", "Target feature richness", value = 0.80, min = 0.001, max = 0.999, step = 0.01),
          numericInput("confidence_level_sub", "Detection confidence", value = 0.95, min = 0.001, max = 0.999, step = 0.01),
          actionButton("run_sub", "Run analysis", class = "btn-primary"),
          tags$hr(),
          downloadButton("download_csv_sub", "Download summary CSV"),
          downloadButton("download_plots_sub", "Download plots PDF")
        ),
        mainPanel(
          uiOutput("subiso_status"),
          fluidRow(
            column(5, htmlOutput("subiso_summary")),
            column(7, plotlyOutput("subiso_obs_pred", height = 420))
          ),
          fluidRow(
            column(6, plotlyOutput("subiso_mass", height = 360)),
            column(6, plotlyOutput("subiso_richness", height = 360))
          ),
          fluidRow(
            column(6, plotlyOutput("subiso_cov_req", height = 320)),
            column(6, plotlyOutput("subiso_rich_req", height = 320))
          ),
          tags$hr(),
          div(class = "app-intro", HTML(paste0("<strong>Citation:</strong> ", citation_text)))
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

server <- function(input, output, session) {
  isolate_results <- reactiveVal(NULL)
  subiso_results <- reactiveVal(NULL)
  isolate_error <- reactiveVal(NULL)
  subiso_error <- reactiveVal(NULL)
  isolate_running <- reactiveVal(FALSE)
  subiso_running <- reactiveVal(FALSE)

  output$isolate_status <- renderUI({
    if (isolate_running()) {
      return(div(class = "analysis-status",
                 div(class = "analysis-spinner"),
                 tags$span("Running analysis; this may take a little while...")))
    }
    err <- isolate_error()
    if (is.null(err)) return(NULL)
    div(class = "error-box", HTML(paste0("<strong>Error:</strong> ", htmlEscape(err))))
  })

  output$subiso_status <- renderUI({
    if (subiso_running()) {
      return(div(class = "analysis-status",
                 div(class = "analysis-spinner"),
                 tags$span("Running analysis; this may take a little while...")))
    }
    err <- subiso_error()
    if (is.null(err)) return(NULL)
    div(class = "error-box", HTML(paste0("<strong>Error:</strong> ", htmlEscape(err))))
  })

  output$download_template_isolate <- downloadHandler(
    filename = function() "isolate_level_template.csv",
    content = make_template_isolate
  )

  output$download_template_sub <- downloadHandler(
    filename = function() "sub_isolate_template.csv",
    content = make_template_subiso
  )

  output$subiso_detected_ui <- renderUI({
    df <- tryCatch(read_csv_source(input$file_sub, input$csv_text_sub), error = function(e) NULL)
    if (is.null(df) || ncol(df) < 2) return(NULL)
    id_col <- infer_id_column(df)
    genes <- setdiff(names(df), id_col)
    div(class = "detected-box",
        tags$strong("Automatically detected:"),
        tags$br(),
        paste0("ID column: ", id_col),
        tags$br(),
        paste0("Gene/feature columns (", length(genes), "): ", paste(genes, collapse = ", ")))
  })

  observeEvent(input$run, {
    isolate_error(NULL)
    isolate_running(TRUE)
    tryCatch({
      shiny::withProgress(message = "Running isolate-level analysis", value = 0, {
        shiny::incProgress(0.05, detail = "Reading input data")
      df <- read_isolate_table(input$file, input$csv_text, input$feature_col, input$count_col)
        shiny::incProgress(0.10, detail = "Preparing priors")
      alpha_prior <- safe_num(input$alpha_prior)
      if (!is.finite(alpha_prior) || alpha_prior < 0) stop("alpha prior must be non-negative.")
      alpha_named <- stats::setNames(rep(alpha_prior, nrow(df)), as.character(df[[input$feature_col]]))
      B <- as.integer(input$B)
      if (!is.finite(B) || B < 100) stop("Posterior draws must be at least 100.")
      seed <- as.integer(input$seed)
      if (!is.finite(seed)) stop("Seed must be an integer.")
      target_coverage <- validate_prob(input$target_coverage, "Target sample coverage")
      target_richness <- validate_prob(input$target_richness, "Target feature richness")
      confidence <- validate_prob(input$confidence_level, "Detection confidence")
      shiny::incProgress(0.10, detail = "Estimating novel-feature mass")
      crp <- fit_crp_theta_bayes(
        counts = df$count,
        n_draws = 5000,
        seed = seed
      )
      novelty_prob_hat <- mean(crp$novelty_prob_draws)
      alpha_novel_sum_crp <- set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat)
      if (!is.finite(alpha_novel_sum_crp) || alpha_novel_sum_crp < 0) {
        stop("CRP-derived novel prior mass must be non-negative.")
      }
      # Keep every novel feature on the same prior scale as the observed features.
      # First derive the number of novel features from the CRP-estimated novel mass,
      # then reconstruct the total novel prior mass from that count and alpha_prior.
      # This matches the manuscript calculation:
      #   alpha_novel_num <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
      #   alpha_novel_sum <- alpha_novel_num * alpha_prior
      alpha_novel_num <- as.integer(ceiling(alpha_novel_sum_crp))
      alpha_novel_sum <- alpha_novel_num * alpha_prior

      shiny::incProgress(0.20, detail = "Sampling posterior feature frequencies")
      fit <- isolate_bayes_dirichlet(
        df = df,
        feature_col = input$feature_col,
        alpha_named = alpha_named,
        alpha_novel_sum = alpha_novel_sum,
        alpha_novel_num = alpha_novel_num,
        B = B,
        seed = seed
      )

      observed <- df |>
        transmute(feature = as.character(.data[[input$feature_col]]), actual = count / sum(count))
      feature_summary <- fit$summary_df |>
        filter(!grepl("^NOVEL_", feature)) |>
        left_join(observed, by = "feature") |>
        mutate(estimate = median, lo = q2.5, hi = q97.5, dataset = "Uploaded data")

      shiny::incProgress(0.25, detail = "Calculating sampling curves")
      prep <- prep_bootstrap_draws(fit$draws)
      mass_curve <- compute_mass_curve(prep, f_grid = f_grid_99)
      richness_curve <- compute_species_richness_curve(prep, f_grid = f_grid_99)
      req <- summarise_required_n(fit$draws, target_coverage, target_richness, confidence)

      shiny::incProgress(0.20, detail = "Finalising results")
      isolate_results(list(
        raw = df,
        fit = fit,
        crp = crp,
        feature_summary = feature_summary,
        mass_curve = mass_curve,
        richness_curve = richness_curve,
        required = req,
        alpha_prior = alpha_prior,
        alpha_novel_sum = alpha_novel_sum,
        alpha_novel_sum_crp = alpha_novel_sum_crp,
        novelty_prob_hat = novelty_prob_hat,
        alpha_novel_num = alpha_novel_num,
        target_coverage = target_coverage,
        target_richness = target_richness,
        confidence = confidence,
        B = B,
        seed = seed
      ))
      })
    }, error = function(e) isolate_error(conditionMessage(e)), finally = isolate_running(FALSE))
  })

  observeEvent(input$run_sub, {
    subiso_error(NULL)
    subiso_running(TRUE)
    tryCatch({
      shiny::withProgress(message = "Running sub-isolate-level analysis", value = 0, {
        shiny::incProgress(0.05, detail = "Reading input data")
      parsed <- read_subiso_table(input$file_sub, input$csv_text_sub)
        shiny::incProgress(0.10, detail = "Detecting genes and preparing priors")
      iso_df <- parsed$data
      gene_cols <- parsed$gene_cols
      B <- as.integer(input$B_sub)
      if (!is.finite(B) || B < 100) stop("Posterior draws must be at least 100.")
      seed <- as.integer(input$seed_sub)
      if (!is.finite(seed)) stop("Seed must be an integer.")
      target_coverage <- validate_prob(input$target_coverage_sub, "Target sample coverage")
      target_richness <- validate_prob(input$target_richness_sub, "Target feature richness")
      confidence <- validate_prob(input$confidence_level_sub, "Detection confidence")
      beta_prior <- safe_num(input$beta_prior)
      if (!is.finite(beta_prior) || beta_prior < 0) stop("beta prior must be non-negative.")
      beta_named <- rep(beta_prior, length(gene_cols))
      names(beta_named) <- gene_cols

      u_hat <- estimate_subiso_u_hat(iso_df, gene_cols)
      K <- length(gene_cols)
      beta_novel_num <- if (u_hat > 0 && beta_prior > 0) ceiling((u_hat * K) / (1 - u_hat)) else 0L
      beta_novel_num <- as.integer(beta_novel_num)
      beta_novel_sum <- beta_novel_num * beta_prior

      shiny::incProgress(0.15, detail = "Sampling posterior feature prevalences")
      bbs <- subisolate_bayes_beta(
        iso_df = iso_df,
        gene_cols = gene_cols,
        B = B,
        beta_named = beta_named,
        beta_novel_sum = beta_novel_sum,
        beta_novel_num = beta_novel_num,
        seed = seed
      )

      p_mat <- get_p_matrix(bbs)
      actual_vec <- colMeans(iso_df[, gene_cols, drop = FALSE], na.rm = TRUE)
      actual_vec <- stats::setNames(actual_vec, gene_cols)
      feature_summary <- summarise_subiso_posterior(
        bbs = bbs,
        gene_cols = gene_cols,
        actual_vec = actual_vec,
        beta_named = beta_named,
        beta_novel_sum = beta_novel_sum,
        beta_novel_num = beta_novel_num,
        dataset_label = "Uploaded data"
      )

      shiny::incProgress(0.25, detail = "Calculating sampling curves")
      mass_curve <- compute_mass_curve_subiso(p_mat, f_grid = f_grid_99)
      richness_curve <- compute_species_richness_subiso(p_mat, f_grid = f_grid_99)
      req <- summarise_required_n(p_mat, target_coverage, target_richness, confidence)

      shiny::incProgress(0.20, detail = "Finalising results")
      subiso_results(list(
        raw = iso_df,
        id_col = parsed$id_col,
        gene_cols = gene_cols,
        fit = bbs,
        feature_summary = feature_summary,
        mass_curve = mass_curve,
        richness_curve = richness_curve,
        required = req,
        beta_prior = beta_prior,
        u_hat = u_hat,
        beta_named = beta_named,
        beta_novel_sum = beta_novel_sum,
        beta_novel_num = beta_novel_num,
        target_coverage = target_coverage,
        target_richness = target_richness,
        confidence = confidence,
        B = B,
        seed = seed
      ))
      })
    }, error = function(e) subiso_error(conditionMessage(e)), finally = subiso_running(FALSE))
  })

  output$isolate_summary <- renderUI({
    req(isolate_results())
    r <- isolate_results()
    tags$div(class = "summary-box",
      HTML(paste0(
        "<p><strong>Observed isolates:</strong> ", format(sum(r$raw$count), big.mark = ","), "</p>",
        "<p><strong>Observed feature categories:</strong> ", nrow(r$raw), "</p>",
        "<p><strong>CRP estimated next-feature probability:</strong> ", signif(r$novelty_prob_hat, 4), "</p>",
        "<p><strong>alpha prior:</strong> ", signif(r$alpha_prior, 5), "</p>",
        "<p><strong>CRP-derived novel prior mass (before ceiling):</strong> ", signif(r$alpha_novel_sum_crp, 5), "</p>",
        "<p><strong>Novel prior mass used:</strong> ", signif(r$alpha_novel_sum, 5), "</p>",
        "<p><strong>Novel feature categories derived:</strong> ", r$alpha_novel_num, "</p>",
        "<p><strong>Coverage sample size:</strong> ", format_n(r$required$coverage_summary["median"]),
        " (95% interval ", format_n(r$required$coverage_summary["q2.5"]), "-", format_n(r$required$coverage_summary["q97.5"]), ")</p>",
        "<p><strong>Richness sample size:</strong> ", format_n(r$required$richness_summary["median"]),
        " (95% interval ", format_n(r$required$richness_summary["q2.5"]), "-", format_n(r$required$richness_summary["q97.5"]), ")</p>",
        "<p><strong>99% detection f-grid:</strong> exact expression (no rounding; ", length(f_grid_99), " points; n = 0 through 100,000)</p>"
      ))
    )
  })

  output$subiso_summary <- renderUI({
    req(subiso_results())
    r <- subiso_results()
    tags$div(class = "summary-box",
      HTML(paste0(
        "<p><strong>Isolate ID column:</strong> ", htmlEscape(r$id_col), "</p>",
        "<p><strong>Isolates:</strong> ", nrow(r$raw), "</p>",
        "<p><strong>Automatically detected gene/feature columns:</strong> ", length(r$gene_cols), "</p>",
        "<p><strong>beta prior:</strong> ", signif(r$beta_prior, 4), "</p>",
        "<p><strong>Estimated u_hat:</strong> ", signif(r$u_hat, 5), "</p>",
        "<p><strong>Novel feature categories derived:</strong> ", r$beta_novel_num, "</p>",
        "<p><strong>beta_novel total mass:</strong> ", signif(r$beta_novel_sum, 5), "</p>",
        "<p><strong>Coverage sample size:</strong> ", format_n(r$required$coverage_summary["median"]),
        " (95% interval ", format_n(r$required$coverage_summary["q2.5"]), "-", format_n(r$required$coverage_summary["q97.5"]), ")</p>",
        "<p><strong>Richness sample size:</strong> ", format_n(r$required$richness_summary["median"]),
        " (95% interval ", format_n(r$required$richness_summary["q2.5"]), "-", format_n(r$required$richness_summary["q97.5"]), ")</p>",
        "<p><strong>99% detection f-grid:</strong> exact expression (no rounding; ", length(f_grid_99), " points; n = 0 through 100,000)</p>"
      ))
    )
  })

  output$isolate_obs_pred <- renderPlotly({ req(isolate_results()); plotly_from_gg(plot_obs_vs_pred(isolate_results()$feature_summary, "Observed vs posterior prevalence")) })
  output$isolate_mass <- renderPlotly({ req(isolate_results()); plotly_from_gg(plot_mass_curve(isolate_results()$mass_curve, "Posterior population mass above threshold")) })
  output$isolate_richness <- renderPlotly({ req(isolate_results()); plotly_from_gg(plot_richness_curve(isolate_results()$richness_curve, "Posterior feature richness above threshold")) })
  output$isolate_cov_req <- renderPlotly({ req(isolate_results()); plotly_from_gg(plot_required_hist(isolate_results()$required$coverage, "Required sample size: target coverage")) })
  output$isolate_rich_req <- renderPlotly({ req(isolate_results()); plotly_from_gg(plot_required_hist(isolate_results()$required$richness, "Required sample size: target richness")) })

  output$subiso_obs_pred <- renderPlotly({ req(subiso_results()); plotly_from_gg(plot_obs_vs_pred(subiso_results()$feature_summary, "Observed vs posterior prevalence")) })
  output$subiso_mass <- renderPlotly({ req(subiso_results()); plotly_from_gg(plot_mass_curve(subiso_results()$mass_curve, "Posterior population mass above threshold")) })
  output$subiso_richness <- renderPlotly({ req(subiso_results()); plotly_from_gg(plot_richness_curve(subiso_results()$richness_curve, "Posterior feature richness above threshold")) })
  output$subiso_cov_req <- renderPlotly({ req(subiso_results()); plotly_from_gg(plot_required_hist(subiso_results()$required$coverage, "Required sample size: target coverage")) })
  output$subiso_rich_req <- renderPlotly({ req(subiso_results()); plotly_from_gg(plot_required_hist(subiso_results()$required$richness, "Required sample size: target richness")) })

  output$download_csv <- downloadHandler(
    filename = function() paste0("nagy2026_isolate_summary_", Sys.Date(), ".csv"),
    content = function(file) {
      req(isolate_results())
      r <- isolate_results()
      out <- r$feature_summary |>
        select(feature, actual, estimate, lo, hi) |>
        mutate(alpha_prior = r$alpha_prior,
               alpha_novel_sum_crp = r$alpha_novel_sum_crp,
               alpha_novel_sum = r$alpha_novel_sum,
               alpha_novel_num = r$alpha_novel_num)
      readr::write_csv(out, file)
    }
  )

  output$download_csv_sub <- downloadHandler(
    filename = function() paste0("nagy2026_subisolate_summary_", Sys.Date(), ".csv"),
    content = function(file) {
      req(subiso_results())
      r <- subiso_results()
      out <- r$feature_summary |>
        mutate(beta_prior = r$beta_prior,
               u_hat = r$u_hat,
               beta_novel_sum = r$beta_novel_sum,
               beta_novel_num = r$beta_novel_num)
      readr::write_csv(out, file)
    }
  )

  output$download_plots <- downloadHandler(
    filename = function() paste0("nagy2026_isolate_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      req(isolate_results())
      r <- isolate_results()
      grDevices::pdf(file, width = 11, height = 8.5)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(plot_obs_vs_pred(r$feature_summary, "Observed vs posterior prevalence"))
      print(plot_mass_curve(r$mass_curve, "Posterior population mass above threshold"))
      print(plot_richness_curve(r$richness_curve, "Posterior feature richness above threshold"))
      print(plot_required_hist(r$required$coverage, "Required sample size: target coverage"))
      print(plot_required_hist(r$required$richness, "Required sample size: target richness"))
    }
  )

  output$download_plots_sub <- downloadHandler(
    filename = function() paste0("nagy2026_subisolate_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      req(subiso_results())
      r <- subiso_results()
      grDevices::pdf(file, width = 11, height = 8.5)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(plot_obs_vs_pred(r$feature_summary, "Observed vs posterior prevalence"))
      print(plot_mass_curve(r$mass_curve, "Posterior population mass above threshold"))
      print(plot_richness_curve(r$richness_curve, "Posterior feature richness above threshold"))
      print(plot_required_hist(r$required$coverage, "Required sample size: target coverage"))
      print(plot_required_hist(r$required$richness, "Required sample size: target richness"))
    }
  )
}

shinyApp(ui, server)
