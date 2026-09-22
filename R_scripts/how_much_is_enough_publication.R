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
install.packages("bayes")
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
install.packages("ggtree")
install.packages("vegan")
install.packages("cluster")


library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(purrr)
library(grafify)
library(tibble)
library(vegan)
library(cluster)
library(stringr)
library(patchwork)
library(gtsummary)
library(gt)
library(data.table)
library(viridis)
library(forcats)
library(scales)
library(gtools)
library(posterior)
library(cmdstanr)
library(loo)
library(future)
library(furrr)
library(ape)
library(ggtree)
library(phytools)
library(RColorBrewer)
library(ggnewscale)
library(matrixStats)
library(parallel)

# set workign directory
setwd("~/your_working_directory_here/")

# read-in cleaned/deduplicated data (under embargo until April 2027)
# load E.coli and Kleb MLST and fastbaps_L3 cluster data
ecoli_bsi_samples_metadata <- read.csv("neksus_ecoli_bsi_samples_metadata.csv")
kleb_bsi_samples_metadata <- read.csv("neksus_kleb_bsi_samples_metadata.csv")
colnames(ecoli_bsi_samples_metadata)
#View(ecoli_bsi_samples_metadata)
colnames(kleb_bsi_samples_metadata)


# load amrfinder and plasmid data 
## this df already had plasmid "community_subcommunity" (pling) and "rep_types_whole_plasmid" (mob-suite)
ecoli_bsi_amrfinder_metadata <- read.csv("neksus_ecoli_bsi_amrfinder_metadata.csv")
kleb_bsi_amrfinder_metadata <- read.csv("neksus_kleb_bsi_amrfinder_metadata.csv")
#colnames(ecoli_bsi_amrfinder_metadata)
#colnames(kleb_bsi_amrfinder_metadata)

# merge region
ecoli_regions <- ecoli_bsi_samples_metadata |>
  select(c(sequencing_id, run, region)) |>
  rename(sample = sequencing_id)
# 
kleb_regions <- kleb_bsi_samples_metadata |>
  select(c(sequencing_id, run, region)) |>
  rename(sample = sequencing_id)

# add region
ecoli_bsi_amrfinder_metadata <- ecoli_bsi_amrfinder_metadata |>
  left_join(ecoli_regions, by = c("sample" = "sample", "run" = "run"))

# add region to kleb
kleb_bsi_amrfinder_metadata <- kleb_bsi_amrfinder_metadata |>
  select(-c(region)) |>
  left_join(kleb_regions, by = c("sample" = "sample", "run" = "run"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# check
ecoli_arg_results <- readRDS("sens_results/ecoli_arg_beta_binomial_regional_checkpoint_novel.rds")
nrow(ecoli_arg_results)
kleb_arg_results <- readRDS("sens_results/kleb_arg_beta_binomial_regional_checkpoint_with_diagnostics_novel.rds")
nrow(kleb_arg_results)
ecoli_pling_results <- readRDS("sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel.rds")
nrow(ecoli_pling_results)
ecoli_pling_results_a0.1_b100 <- readRDS("sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a0.1_b100.rds")
nrow(ecoli_pling_results_a0.1_b100)
ecoli_pling_results_a0.1_b100_m1000 <- readRDS("sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a0.1_b100_m1000.rds")
nrow(ecoli_pling_results_a0.1_b100_m1000)
ecoli_pling_results_a1 <- readRDS("sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a1.rds")
nrow(ecoli_pling_results_a1)
ecoli_pling_results_a10 <- readRDS("sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a10.rds")
nrow(ecoli_pling_results_a10)
kleb_pling_results <- readRDS("sens_results/kleb_pling_beta_binomial_regional_checkpoint_with_diagnostics_novel.rds")
nrow(kleb_pling_results)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Descriptive summary of data: ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Frequency distribution histograms of MLSTs, fastbaps_L3 clusters, plasmids and ARGs 
ecoli_samples_df <- ecoli_bsi_samples_metadata
ecoli_amr_df     <- ecoli_bsi_amrfinder_metadata
klebsiella_samples_df <- kleb_bsi_samples_metadata   
klebsiella_amr_df     <- kleb_bsi_amrfinder_metadata
#colnames(ecoli_bsi_samples_metadata)
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

n_bins <- 40

# helper functions 
# given a samples df and the column name for MLST, compute counts per MLST
get_mlst_counts <- function(samples_df, mlst_col, genus_name) {
  samples_df |>
    filter(!is.na(.data[[mlst_col]])) |>
    group_by(!!sym(mlst_col)) |>
    summarise(count = n(), .groups = "drop") |>
    transmute(metric = "MLSTs", genus = genus_name, value = count)
}

get_mlst_counts_by_region <- function(samples_df, mlst_col, genus_name) {
  samples_df |>
    filter(!is.na(.data[[mlst_col]])) |>
    group_by(region, !!sym(mlst_col)) |>
    summarise(count = n(), .groups = "drop") |>
    transmute(metric = "MLSTs", region = region, genus = genus_name, value = count)
}

# given samples df and fastbaps_L3 level column name compute counts
get_fb_counts <- function(samples_df, level_col, genus_name, level_label) {
  samples_df |>
    filter(!is.na(.data[[level_col]])) |>
    group_by(!!sym(level_col)) |>
    summarise(count = n(), .groups = "drop") |>
    transmute(metric = paste0("fastBAPS clusters ", level_label), genus = genus_name, value = count)
}
get_fb_counts_by_region <- function(samples_df, level_col, genus_name, level_label) {
  samples_df |>
    filter(!is.na(.data[[level_col]])) |>
    group_by(region, !!sym(level_col)) |>
    summarise(count = n(), .groups = "drop") |>
    transmute(metric = paste0("fastBAPS clusters ", level_label),region = region, genus = genus_name, value = count)
}

# ARG counts (count occurrences of each gene across samples)
get_arg_counts <- function(amr_df, genus_name) {
  amr_df |>
    filter(!is.na(.data[[arg_symbol_col]]), .data[[arg_type_col]] == "AMR") |>
    group_by(!!sym(arg_symbol_col)) |>
    summarise(count = n(), .groups = "drop") |>
    transmute(metric = "ARGs", genus = genus_name, value = count)
}

get_arg_counts_by_region <- function(amr_df, genus_name) {
  amr_df |>
    filter(!is.na(.data[[arg_symbol_col]]), .data[[arg_type_col]] == "AMR") |>
    group_by(region, !!sym(arg_symbol_col)) |>
    summarise(count = n(), .groups = "drop") |>
    transmute(metric = "ARGs", region = region, genus = genus_name, value = count)
}

# plasmid subcommunity counts (deduplicate contigs per sample then count number of contigs per subcommunity)
get_plasmid_subcommunity_counts <- function(amr_df, genus_name) {
  amr_df |>
    filter(!is.na(.data[[subcommunity_col]])) |>
    # dedupe by sample/run/contig so contig not double counted
    group_by(!!!syms(sample_id_cols)) |>
    slice_head(n = 1) |>
    ungroup() |>
    group_by(!!sym(subcommunity_col)) |>
    summarise(count = n(), .groups = "drop") |>
    transmute(metric = "Plasmid subcommunities", genus = genus_name, value = count)
}

get_plasmid_subcommunity_counts_by_region <- function(amr_df, genus_name) {
  amr_df |>
    filter(!is.na(.data[[subcommunity_col]])) |>
    # dedupe by sample/run/contig so contig not double counted
    group_by(region, !!!syms(sample_id_cols)) |>
    slice_head(n = 1) |>
    ungroup() |>
    group_by(region, !!sym(subcommunity_col)) |>
    summarise(count = n(), .groups = "drop") |>
    transmute(metric = "Plasmid subcommunities", region = region, genus = genus_name, value = count)
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
  ecoli_mlst, ecoli_fb_l1, ecoli_fb_l2, ecoli_fb_l3,  ecoli_args, ecoli_plasm,
  kleb_mlst,  kleb_fb_l1,  kleb_fb_l2,  kleb_fb_l3,  kleb_args,  kleb_plasm
)
#View(all_counts)

# call functions - by region
# Escherichia
ecoli_mlst   <- get_mlst_counts_by_region(ecoli_samples_df, mlst_col_ecoli, "Escherichia")
ecoli_fb_l1  <- get_fb_counts_by_region(ecoli_samples_df, fb_cols[1], "Escherichia", "L1")
ecoli_fb_l2  <- get_fb_counts_by_region(ecoli_samples_df, fb_cols[2], "Escherichia", "L2")
ecoli_fb_l3  <- get_fb_counts_by_region(ecoli_samples_df, fb_cols[3], "Escherichia", "L3")
ecoli_args   <- get_arg_counts_by_region(ecoli_amr_df, "Escherichia")
ecoli_plasm  <- get_plasmid_subcommunity_counts_by_region(ecoli_amr_df, "Escherichia")

# Klebsiella
kleb_mlst   <- get_mlst_counts_by_region(klebsiella_samples_df, mlst_col_kleb, "Klebsiella")
kleb_fb_l1  <- get_fb_counts_by_region(klebsiella_samples_df, fb_cols[1], "Klebsiella", "L1")
kleb_fb_l2  <- get_fb_counts_by_region(klebsiella_samples_df, fb_cols[2], "Klebsiella", "L2")
kleb_fb_l3  <- get_fb_counts_by_region(klebsiella_samples_df, fb_cols[3], "Klebsiella", "L3")
kleb_args   <- get_arg_counts_by_region(klebsiella_amr_df, "Klebsiella")
kleb_plasm  <- get_plasmid_subcommunity_counts_by_region(klebsiella_amr_df, "Klebsiella")

# bind everything
all_counts_by_region <- bind_rows(
  ecoli_mlst, ecoli_fb_l1, ecoli_fb_l2, ecoli_fb_l3,  ecoli_args, ecoli_plasm,
  kleb_mlst,  kleb_fb_l1,  kleb_fb_l2,  kleb_fb_l3,  kleb_args,  kleb_plasm
)
#View(all_counts_by_region)

# add cumulative population mass
all_counts_cumulative_mass <- all_counts |>
  mutate(total_isolates = case_when(genus == "Escherichia" ~ 1471,
                                    genus == "Klebsiella" ~ 468,
                                    TRUE ~ NA_real_),
         frequency = value/total_isolates) |>
  group_by(metric, genus) |>
  mutate(total_features = sum(value)) |>
  mutate(frequency_among_features = value / total_features) |>
  ungroup () |>
  group_by(metric, genus, frequency, frequency_among_features) |>
  summarise(frequency_mass = sum(frequency),
            frequency_feature_mass = sum(frequency_among_features),
            .groups = "drop"
  ) |>
  group_by(metric, genus) |>
  arrange(desc(frequency_among_features), .by_group = TRUE) |>
  mutate(cumulative_frequency = cumsum(frequency_feature_mass))

#View(all_counts_cumulative_mass)


# label factor order for nice faceting
metric_levels <- c("MLSTs", "fastBAPS clusters L1", "fastBAPS clusters L2", "fastBAPS clusters L3",  "Plasmid subcommunities", "ARGs")
all_counts$metric <- factor(all_counts$metric, levels = metric_levels)
all_counts_by_region$metric <- factor(all_counts_by_region$metric, levels = metric_levels)
all_counts_cumulative_mass$metric <- factor(all_counts_cumulative_mass$metric, levels = metric_levels)

# keep only fastbaos level 3
plot_df <- all_counts |>  filter(!metric %in% c("fastBAPS clusters L1", "fastBAPS clusters L2"))
plot_df_region <- all_counts_by_region |>  filter(!metric %in% c("fastBAPS clusters L1", "fastBAPS clusters L2")) |>
  filter(region != "unknown")
cum_df <- all_counts_cumulative_mass |>  filter(!metric %in% c("fastBAPS clusters L1", "fastBAPS clusters L2"))


# make plot label df
label_df <- tibble::tibble(
  metric = factor(
    c("MLSTs",
      "fastBAPS clusters L3",
      "Plasmid subcommunities",
      "ARGs"),
    levels = metric_levels
  ),
  label_left  = c("a)", "b)", "c)", "d)"),
  label_right = c("e)", "f)", "g)", "h)")
)

# right-hand cumulative plot
cum_plot <- ggplot(cum_df, aes(x = frequency, y = cumulative_frequency, colour = genus)) +
  geom_step(direction = "hv", linewidth = 0.7) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = genus_cols, name = "Genus") +
  geom_text(data = label_df, aes(x = 0, y = Inf, label = label_right),
            inherit.aes = FALSE,  hjust = 1, vjust = -1, fontface = "bold", size = 4) +
  scale_x_log10() +
  labs(x = "Frequency", y = "Cumulative frequency") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")
cum_plot

counts_left <- ggplot(plot_df, aes(x = value, fill = genus)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.6, closed = "right") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = genus_cols, name = "Genus") +
  geom_text(data = label_df, aes(x = 0, y = Inf, label = label_left),
            inherit.aes = FALSE,  hjust = 1, vjust = -1, fontface = "bold", size = 4) +
  scale_x_log10() +
  labs(x = "Number of isolates per feature", y = "Number of unique features") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")
counts_left


counts_panel <- counts_left + cum_plot +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")
counts_panel
ggsave("model_results/counts_histogram_plus_cumulative_frequency.png",  counts_panel, width = 7, height = 9, dpi = 300)

# frequency historgram panle plot
freq_left <- ggplot(plot_df, aes(x = value, fill = genus, group = genus)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    bins = n_bins,
    position = "identity",
    alpha = 0.6,
    colour = NA
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(Escherichia = "seagreen3", Klebsiella = "darkorange")) +
  scale_x_log10() +
  labs(
    x = "Number of isolates per feature",
    y = "Proportion of unique features",
    fill = "Genus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

freq_panel <- freq_left + cum_plot +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")
freq_panel

ggsave("model_results/frequency_histogram_plus_cumulative_frequency.png",  freq_panel, width = 7,  height = 9, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make regionally-stratified histograms
hist_region <- ggplot(plot_df_region, aes(x = value, fill = genus)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.6, closed = "right") +
  facet_wrap(metric ~ region, scales = "free_y", ncol = 10) +
  scale_fill_manual(values = genus_cols, name = "Genus") +
  scale_x_log10() +
  labs(x = "Number of isolates per feature", y = "Number of unique features") +
  theme_minimal(base_size = 8) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
hist_region
ggsave("model_results/count_histogram_of_features_by_region.png",  hist_region, width = 14,  height = 7, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 0.a Phylogenetic trees ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 0.ai E. coli  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data <- ecoli_bsi_samples_metadata
# View(data)
# correct sample names and make matchable to tip labels 
data <- data |>
  mutate(seqid_corrected = case_when(run == "run0" ~ paste0("pilot_", sequencing_id),
                                     run != "run0" ~ gsub("AH", "AH-", sequencing_id)),
    tip_label = paste0(run, "_", seqid_corrected))
# check against tree tip labels
#tree <- read.tree("ecoli_bsi_verticall.newick")
#setdiff(data$tip_label, tree$tip.label)
#setdiff(tree$tip.label, data$tip_label)
nrow(data) #1471
colnames(data)
table(data$escherichia__mlst_achtman__ST, useNA = "ifany")

n_isolates <- length(unique(data$tip_label)) # 468

big_mlst <- data.frame(table(data$escherichia__mlst_achtman__ST)) 
single_mlst <- big_mlst |> filter(Freq == 1) 
big_mlst <- filter(big_mlst,!Var1=='-') |> 
  filter(Freq > (n_isolates * 0.02) ) |>
  arrange(desc(Freq)) # more than 5 isolates of same ST

data <- data |>
  mutate(mlst_cat = case_when(is.na(escherichia__mlst_achtman__ST) ~ "NA",
                              escherichia__mlst_achtman__ST %in% big_mlst$Var1 ~ escherichia__mlst_achtman__ST,
                              escherichia__mlst_achtman__ST %in% single_mlst$Var1 ~ "Single",
                              TRUE ~ "Minor"))
data <- data |>
  dplyr::select(tip_label, everything())

# order mlst cat
mlst_cat_order <- data |>
  group_by(mlst_cat) |>
  summarise(count = n()) |>
  arrange(count) |>
  pull(mlst_cat)
#extract brewer 2 colours
mlst_levels <- c(mlst_cat_order)
data$mlst_cat <- factor(data$mlst_cat, levels = mlst_levels)
table(data$mlst_cat, useNA = "ifany")


mlst_colours <- c(
  "Single" = "#BDBDBD",   # light grey
  "Minor"  = "#E6E6E6",   # very light grey
  
  "ST131"  = "#4477AA",   # blue
  "ST73"   = "#EE6677",   # coral red
  "ST69"   = "#228833",   # green
  "ST95"   = "#66CCEE",   # cyan
  "ST12"   = "#AA3377",   # purple
  "ST127"  = "#EE7733",   # orange
  
  "ST1193" = "#CCBB44",   # mustard yellow
  "ST10"   = "#44AA99",   # teal # removed when thresholds raised
  "ST404"  = "#882255",   # burgundy # removed when thresholds raised
  "ST141"  = "#999933"    # olive # removed when thresholds raised
)

#prep amrfinder data
amrfinder_data <- ecoli_bsi_amrfinder_metadata |>
  mutate(seqid_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                     run != "run0" ~ gsub("AH", "AH-", sample)),
         tip_label = paste0(run, "_", seqid_corrected)) |>
  filter(!is.na(Element.symbol)) |>
  filter(Type == "AMR")
data$region <- as.factor(data$region)
# check
#setdiff(amrfinder_data$tip_label, tree$tip.label)
#setdiff(tree$tip.label, amrfinder_data$tip_label)
length(unique(amrfinder_data$sample)) #1471
#colnames(amrfinder_data)

 
# create metadata panels data subsets
data$Level.3 <- as.factor(data$Level.3)
region_fastbaps <- data |>
  select(tip_label, region, Level.3) |>
  pivot_longer(cols = c(region, Level.3) , names_to = "var")
# set order
region_fastbaps$var <- factor(region_fastbaps$var, levels = c("region", "Level.3"))
#separate dfs
region_df <- region_fastbaps|>
  filter(var == "region")
region_df$value <- as.factor(region_df$value)
region_df <- region_df |>
  mutate(value = factor(value, levels = levels(value)),
         xlab = factor(""))

fastbaps_df <- region_fastbaps|>
  filter(var == "Level.3")
fastbaps_df$value <- as.factor(fastbaps_df$value)
fastbaps_df <- fastbaps_df |>
  mutate(xlab = factor(""))
table(fastbaps_df$value, useNA = "ifany")



# select genes occuring in > 5 isolates for plotting
gene_counts <- data.frame(table(amrfinder_data$Element.symbol)) |> 
  filter(Freq >= 5) #e.g. filter out if occur <5 times in population

#filter out gent res genes
carb <- filter(amrfinder_data,grepl('CARBAPENEM', Subclass)) |> 
  distinct(Element.symbol) #|> filter(Element.symbol %in% gene_counts$Var1) do not filter
carb_selected <- amrfinder_data |>
  filter(Element.symbol %in% carb$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(carb_selected)
carb_selected$present <- as.factor(carb_selected$present)

ceph <- filter(amrfinder_data,grepl('CEPHALOSPORIN',Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts$Var1) |>  filter(!Element.symbol=='blaEC-5')
ceph_selected <- amrfinder_data |>
  filter(Element.symbol %in% ceph$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(ceph_selected)
ceph_selected$present <- as.factor(ceph_selected$present)

beta_lactam <- filter(amrfinder_data, grepl('BETA-LACTAM', Subclass)) |> 
  distinct(Element.symbol)  |> filter(Element.symbol %in% gene_counts$Var1) |> 
  filter(!Element.symbol == 'blaEC')
beta_lactam_selected <- amrfinder_data |>
  filter(Element.symbol %in% beta_lactam$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(beta_lactam_selected)
beta_lactam_selected$present <- as.factor(beta_lactam_selected$present)

gent <- filter(amrfinder_data, grepl('GENTAMICIN', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts$Var1)
gent_selected <- amrfinder_data |>
  filter(Element.symbol %in% gent$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(gent_selected)
gent_selected$present <- as.factor(gent_selected$present)

quin <- filter(amrfinder_data, grepl('QUINOLONE', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts$Var1)
quin_selected <- amrfinder_data |>
  filter(Element.symbol %in% quin$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(quin_selected)
quin_selected$present <- as.factor(quin_selected$present)

# plasmids
unitig_info_mobsuite_combined <- read.csv("autocycer_unitig_info_mobsuite_combined.csv")
unitig_info_mobsuite_combined <- unitig_info_mobsuite_combined |>
  mutate(seqid_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                     run != "run0" ~ gsub("AH", "AH-", sample)),
         tip_label = paste0(run, "_", seqid_corrected))

# keep just ecolisiella isolates
plasmid_data <- unitig_info_mobsuite_combined |>
  filter(tip_label %in% data$tip_label)
nrow(plasmid_data) # 4037
colnames(plasmid_data) # 4037
length(unique(plasmid_data$tip_label))
#View(plasmid_data)


# get top contig reps (accounting for multi-replicon combos)
contig_reps <- plasmid_data |>
  select(tip_label, autocycler_contig_id, rep_type.s.) |>
  filter(rep_type.s.!= "-") |>
  arrange(rep_type.s.) |>
  group_by(tip_label, autocycler_contig_id) |>
  summarise(contig_rep_types = paste(rep_type.s., collapse = ",")) |> # makes no difference to collapse rep types by contig, as 991 annotations both before and after. 
  ungroup() |>
  group_by(contig_rep_types) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  ungroup()
#View(contig_reps)
sum(contig_reps$count) # 3636
nrow(contig_reps) # 211

# get top 20
top10_plasmids <- contig_reps |>
  slice_head(n=10) |>
  pull(contig_rep_types)
top10_plasmids

# make plasmid summary
plasmid_plot_df <- plasmid_data |>
  select(tip_label, autocycler_contig_id, rep_type.s.) |>
  filter(rep_type.s.!= "-") |>
  arrange(rep_type.s.) |>
  group_by(tip_label, autocycler_contig_id) |>
  summarise(contig_rep_types = paste(rep_type.s., collapse = ",")) |> # makes no difference to collapse rep types by contig, as 991 annotations both before and after. 
  ungroup() |>
  mutate(rep_types_cat = case_when(contig_rep_types %in% top10_plasmids ~ contig_rep_types,
                                   grepl("rep_cluster", contig_rep_types) & !grepl(",", contig_rep_types) ~ "other rep_clusters",
                                   grepl("Inc", contig_rep_types) & !grepl(",", contig_rep_types) ~ "other Inc types",
                                   grepl("Col", contig_rep_types) & !grepl(",", contig_rep_types) ~ "other Col-like",
                                   grepl(",", contig_rep_types) ~ "other multi-replicon",
                                   )) |>
  select(-c(autocycler_contig_id, contig_rep_types)) |>
  mutate(present = "present")
#table(plasmid_plot_df$rep_types_cat)
# set order of plasmids
other_levels <- plasmid_plot_df |>
  count(rep_types_cat, sort = TRUE) |>
  filter(grepl("^other", rep_types_cat)) |>
  pull(rep_types_cat)
plasmid_levels <- c(top10_plasmids, other_levels)
plasmid_plot_df <- plasmid_plot_df |>
  mutate(rep_types_cat = factor(rep_types_cat, levels = plasmid_levels))

#View(plasmid_plot_df)
#nrow(plasmids_data)
#View(plasmids_data)
#colnames(amrfinder_data)


# define colours
colors <- c("#0d0c3b", "#ffffff")
region_pal <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
)
region_pal <- rev(region_pal)

# set E. coli fastBPAS cluster colour palette
# order E. coli clusters numerically
num_suffix <- function(x) as.integer(sub(".*_", "", x))

ecoli_levels <- unique(fastbaps_df$value)
ecoli_levels <- ecoli_levels[order(num_suffix(ecoli_levels))]

# muted, colourblind-friendly-ish HCL palette for 161 categories
# low chroma + high luminance = softer colours
ecoli_cols <- grDevices::hcl(
  h = seq(15, 375, length.out = length(ecoli_levels) + 1)[-1],
  c = rep(c(26, 30, 28, 32), length.out = length(ecoli_levels)),
  l = rep(c(84, 80, 82, 78), length.out = length(ecoli_levels))
)

ecoli_pal <- setNames(ecoli_cols, ecoli_levels)
if FALSE {
ecoli_2   ecoli_3   ecoli_4   ecoli_6   ecoli_7   ecoli_8   ecoli_9  ecoli_11  ecoli_12  ecoli_13  ecoli_14  ecoli_15 
"#F0C8C5" "#E9BCB7" "#ECC3BC" "#E4B7AD" "#EFC9C1" "#E6BDB1" "#EAC4B7" "#E1B8A8" "#ECCBBC" "#E3BFAD" "#E6C5B3" "#DEBAA3" 
ecoli_16  ecoli_17  ecoli_18  ecoli_19  ecoli_20  ecoli_21  ecoli_22  ecoli_23  ecoli_24  ecoli_25  ecoli_26  ecoli_27 
"#E9CCB9" "#E0C1A9" "#E3C7AF" "#D9BC9F" "#E5CEB6" "#DBC3A5" "#DFC9AC" "#D4BE9C" "#E1CFB3" "#D6C4A3" "#DACBAA" "#CFC099" 
ecoli_28  ecoli_29  ecoli_30  ecoli_31  ecoli_32  ecoli_33  ecoli_34  ecoli_35  ecoli_36  ecoli_37  ecoli_38  ecoli_39 
"#DCD1B2" "#D1C6A1" "#D5CCA9" "#C9C298" "#D7D3B1" "#CBC8A0" "#CFCEA8" "#C3C498" "#D2D4B1" "#C5CAA0" "#CAD0A9" "#BDC699" 
ecoli_40  ecoli_41  ecoli_42  ecoli_43  ecoli_44  ecoli_45  ecoli_46  ecoli_47  ecoli_48  ecoli_49  ecoli_50  ecoli_51 
"#CDD6B2" "#BFCCA2" "#C4D1AA" "#B6C79B" "#C7D7B3" "#B9CDA4" "#BED3AD" "#AFC99E" "#C1D8B6" "#B2CFA7" "#B8D4B0" "#A8CAA1" 
ecoli_52  ecoli_53  ecoli_54  ecoli_55  ecoli_56  ecoli_57  ecoli_58  ecoli_59  ecoli_60  ecoli_61  ecoli_62  ecoli_63 
"#BCD9B9" "#ACD0AB" "#B2D5B3" "#A1CBA6" "#B7DABD" "#A6D1AF" "#ADD6B8" "#9BCCAB" "#B2DBC1" "#A1D1B4" "#A8D6BC" "#96CDB0" 
ecoli_64  ecoli_65  ecoli_66  ecoli_67  ecoli_68  ecoli_69  ecoli_70  ecoli_71  ecoli_72  ecoli_73  ecoli_74  ecoli_75 
"#AFDCC5" "#9CD2BA" "#A4D7C1" "#91CDB6" "#ABDCCA" "#98D2BF" "#A1D7C6" "#8DCDBC" "#A9DCCF" "#96D2C5" "#9FD7CC" "#8BCDC2" 
ecoli_76  ecoli_77  ecoli_78  ecoli_79  ecoli_80  ecoli_81  ecoli_82  ecoli_83  ecoli_84  ecoli_85  ecoli_86  ecoli_87 
"#A8DCD4" "#95D2CA" "#9ED7D1" "#8BCDC8" "#A8DBD8" "#95D1D0" "#9FD6D6" "#8CCCCD" "#A9DADD" "#97D0D5" "#A1D5DA" "#8FCBD2" 
ecoli_88  ecoli_89  ecoli_90  ecoli_91  ecoli_92  ecoli_93  ecoli_94  ecoli_95  ecoli_96  ecoli_97  ecoli_98  ecoli_99 
"#ACDAE1" "#9ACFD9" "#A4D4DF" "#93CAD7" "#AFD8E5" "#9ECEDD" "#A8D3E2" "#98C8DB" "#B3D7E8" "#A4CCE1" "#AED1E6" "#9FC6DF" 
ecoli_100 ecoli_101 ecoli_102 ecoli_103 ecoli_104 ecoli_105 ecoli_106 ecoli_107 ecoli_108 ecoli_109 ecoli_110 ecoli_111 
"#B8D6EB" "#AACBE4" "#B4CFE8" "#A6C4E1" "#BED4ED" "#B1C9E6" "#BACEEA" "#AEC2E3" "#C4D2EE" "#B8C7E8" "#C0CCEB" "#B6C0E5" 
ecoli_112 ecoli_113 ecoli_114 ecoli_115 ecoli_117 ecoli_118 ecoli_119 ecoli_120 ecoli_121 ecoli_122 ecoli_123 ecoli_124 
"#CAD0EF" "#BFC5E8" "#C7CAEC" "#BDBEE5" "#D0CFEF" "#C7C2E8" "#CEC8EB" "#C5BCE4" "#D6CDEF" "#CDC0E7" "#D4C6EA" "#CCB9E3" 
ecoli_125 ecoli_126 ecoli_127 ecoli_128 ecoli_129 ecoli_130 ecoli_131 ecoli_132 ecoli_133 ecoli_134 ecoli_135 ecoli_136 
"#DCCBED" "#D4BFE6" "#DAC4E9" "#D2B7E1" "#E1CAEB" "#D9BDE3" "#DFC3E6" "#D8B6DE" "#E5C8E9" "#DEBBE0" "#E3C1E3" "#DDB4DA" 
ecoli_137 ecoli_138 ecoli_139 ecoli_140 ecoli_141 ecoli_142 ecoli_143 ecoli_144 ecoli_145 ecoli_146 ecoli_147 ecoli_148 
"#E9C7E6" "#E3BADC" "#E7C0E0" "#E1B3D6" "#EDC7E2" "#E6B9D8" "#EAC0DB" "#E4B3D1" "#EFC6DE" "#E9B9D3" "#EDBFD7" "#E6B2CC" 
ecoli_149 ecoli_150 ecoli_151 ecoli_152 ecoli_153 ecoli_154 ecoli_155 ecoli_156 ecoli_157 ecoli_158 ecoli_159 ecoli_160 
"#F1C6D9" "#EAB9CE" "#EEBFD2" "#E8B2C6" "#F2C6D5" "#EBB9C9" "#EFC0CD" "#E8B3C0" "#F2C6D0" "#EBBAC3" "#EFC0C8" "#E8B4BA" 
ecoli_161 ecoli_162 ecoli_163 ecoli_164 ecoli_165 
"#F2C7CB" "#EBBBBD" "#EEC1C2" "#E7B5B4" "#F1C8C6"
}

# build tree
tree <- read.tree("ecoli_bsi_verticall.newick")
tree$edge.length[tree$edge.length < 0] <- 0
tree <- midpoint_root(tree) #so does this
tree$edge.length <- (tree$edge.length)^(1/3) #this just makes the tree look better
tree_plot <- ggtree(tree, options(ignore.negative.edge=TRUE))

p <- tree_plot %<+% data + 
  geom_tippoint(aes(color = mlst_cat), size =1.2) +  
  scale_color_manual(values = mlst_colours) + 
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_treescale(x = 0, y = length(tree$tip.label) * 0.90, offset = length(tree$tip.label) * 0.01, width = 0.1, color = "black") +  
  theme_tree() 
p

p2 <- p + new_scale_fill() +
  geom_facet(panel = "Region", data = region_df, geom = geom_tile,
      aes( x = xlab, fill = value), width = 1) + 
  scale_fill_manual(
    values = region_pal,
    guide = "none"
  ) +
  
  new_scale_fill() +
  geom_facet(panel = "fast-\nBAPS", data = fastbaps_df, geom = geom_tile,
      aes( x = xlab, fill = value), width = 1) + 
  scale_fill_manual(
    values = ecoli_pal,
    guide = guide_legend(title = "fastBAPS\ncluster", ncol = 7, byrow = TRUE)
  ) +
  
  new_scale_fill() +
  geom_facet(panel = "CARB", data = carb_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) + 
  geom_facet(panel = "CEPH", data = ceph_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) +
  geom_facet(panel = "B-L", data = beta_lactam_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) +
  geom_facet(panel = "GENT", data = gent_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) + 
  geom_facet(panel = "QUINOLONE", data = quin_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) +
  geom_facet(panel = "Plasmids", data = plasmid_plot_df, geom = geom_tile, 
             aes(x = rep_types_cat, fill = present),  width = 1) +
  scale_fill_manual(values = colors, guide = "none") +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), 
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "white", color = "black"),
        legend.position =  "right",
        legend.box = "vertical") +
  scale_x_discrete() +
  labs(color="MLST")
p2


# adjust panel size
gt <- ggplotGrob(p2)
panel_widths <- c(3, 0.4, 0.4, 0.6, 0.8, 0.5, 0.3, 1.7, 1.5)  # same order as geom_facet() calls
panel_idx <- which(grepl("^panel", gt$layout$name))
gt$widths[gt$layout$l[panel_idx]] <- unit(panel_widths, "null")

# preview
grid.newpage()
grid.draw(gt)

# save
ggsave("ecoli_bsi_verticall_tree.png", plot = gt, width = 17, height = 11, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 0.aii Klebsiella  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data <- kleb_bsi_samples_metadata
# View(data)
# correct sample names and make matchable to tip labels 
data <- data |>
  mutate(seqid_corrected = case_when(run == "run0" ~ paste0("pilot_", sequencing_id),
                                     run != "run0" ~ gsub("AH", "AH-", sequencing_id)),
    tip_label = paste0(run, "_", seqid_corrected))
# check
#setdiff(data$tip_label, tree$tip.label)
#setdiff(tree$tip.label, data$tip_label)
nrow(data) #468
colnames(data)

table(data$klebsiella_mlst_ST, useNA = "ifany")

n_isolates <- length(unique(data$tip_label)) # 468

big_mlst <- data.frame(table(data$klebsiella_mlst_ST)) 
single_mlst <- big_mlst |> filter(Freq == 1) 
big_mlst <- filter(big_mlst,!Var1=='-') |> 
  filter(Freq > (n_isolates * 0.015) ) |>
  arrange(desc(Freq)) # more than 5 isolates of same ST

data <- data |>
  mutate(mlst_cat = case_when(is.na(klebsiella_mlst_ST) ~ "NA",
                              klebsiella_mlst_ST %in% big_mlst$Var1 ~ klebsiella_mlst_ST,
                              klebsiella_mlst_ST %in% single_mlst$Var1 ~ "Single",
                              TRUE ~ "Minor"))
table(data$mlst_cat, useNA = "ifany")
data <- data |>
  dplyr::select(tip_label, everything())

# order mlst cat
mlst_cat_order <- data |>
  group_by(mlst_cat) |>
  summarise(count = n()) |>
  arrange(count) |>
  pull(mlst_cat)
#extract brewer 2 colours
mlst_levels <- c(mlst_cat_order)
data$mlst_cat <- factor(data$mlst_cat, levels = mlst_levels)


mlst_colours <- c(
  "Single" = "#BBBBBB",
  "Minor" = "#CCBB44",
  "kpsc ST45" = "#4477AA",
  "kpsc ST14"  =  "#EE6677",
  "kpsc ST29" = "#228833",
  "kpsc ST17" = "#66CCEE",
  "kpsc ST23"  = "#AA3377",
  "kpsc ST307"=  "#EE7733",
  "kpsc ST54" =  "#000000"
)


#prep amrfinder data
amrfinder_data <- kleb_bsi_amrfinder_metadata |>
  mutate(seqid_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                     run != "run0" ~ gsub("AH", "AH-", sample)),
         tip_label = paste0(run, "_", seqid_corrected)) |>
  filter(!is.na(Element.symbol)) |>
  filter(Type == "AMR")
data$region <- as.factor(data$region)
# check
#setdiff(amrfinder_data$tip_label, tree$tip.label)
#setdiff(tree$tip.label, amrfinder_data$tip_label)
#length(unique(amrfinder_data$sample)) #468
#colnames(amrfinder_data)

 
# create metadata panels data subsets
region_fastbaps <- data |>
  select(tip_label, region, Level.3) |>
  pivot_longer(cols = c(region, Level.3) , names_to = "var")
# set order
region_fastbaps$var <- factor(region_fastbaps$var, levels = c("region", "Level.3"))
#separate dfs
region_df <- region_fastbaps|>
  filter(var == "region")
region_df$value <- as.factor(region_df$value)

fastbaps_df <- region_fastbaps|>
  filter(var == "Level.3")
table(fastbaps_df$value, useNA = "ifany")
fastbaps_df$value <- as.factor(fastbaps_df$value)
table(fastbaps_df$value)


# select genes occuring in > 5 isolates for plotting
gene_counts <- data.frame(table(amrfinder_data$Element.symbol)) |> 
  filter(Freq >= 5) #e.g. filter out if occur <5 times in population

#filter out gent res genes
carb <- filter(amrfinder_data,grepl('CARBAPENEM', Subclass)) |> 
  distinct(Element.symbol) #|> filter(Element.symbol %in% gene_counts$Var1) do not filter
carb_selected <- amrfinder_data |>
  filter(Element.symbol %in% carb$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(carb_selected)
carb_selected$present <- as.factor(carb_selected$present)

ceph <- filter(amrfinder_data,grepl('CEPHALOSPORIN',Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts$Var1) |>  filter(!Element.symbol=='blaEC-5')
ceph_selected <- amrfinder_data |>
  filter(Element.symbol %in% ceph$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(ceph_selected)
ceph_selected$present <- as.factor(ceph_selected$present)

beta_lactam <- filter(amrfinder_data, grepl('BETA-LACTAM', Subclass)) |> 
  distinct(Element.symbol)  |> filter(Element.symbol %in% gene_counts$Var1) |> 
  filter(!Element.symbol == 'blaEC')
beta_lactam_selected <- amrfinder_data |>
  filter(Element.symbol %in% beta_lactam$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(beta_lactam_selected)
beta_lactam_selected$present <- as.factor(beta_lactam_selected$present)

gent <- filter(amrfinder_data, grepl('GENTAMICIN', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts$Var1)
gent_selected <- amrfinder_data |>
  filter(Element.symbol %in% gent$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(gent_selected)
gent_selected$present <- as.factor(gent_selected$present)

quin <- filter(amrfinder_data, grepl('QUINOLONE', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts$Var1)
quin_selected <- amrfinder_data |>
  filter(Element.symbol %in% quin$Element.symbol) |>
  select(tip_label, Element.symbol) |>
  mutate(present = "present")
#View(quin_selected)
quin_selected$present <- as.factor(quin_selected$present)

# plasmids
unitig_info_mobsuite_combined <- read.csv("autocycer_unitig_info_mobsuite_combined.csv")
unitig_info_mobsuite_combined <- unitig_info_mobsuite_combined |>
  mutate(seqid_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                     run != "run0" ~ gsub("AH", "AH-", sample)),
         tip_label = paste0(run, "_", seqid_corrected))

# keep just Klebsiella isolates
plasmid_data <- unitig_info_mobsuite_combined |>
  filter(tip_label %in% data$tip_label)
nrow(plasmid_data) # 4037
colnames(plasmid_data) # 4037
length(unique(plasmid_data$tip_label))
#View(plasmid_data)


# get top contig reps (accounting for multi-replicon combos)
contig_reps <- plasmid_data |>
  select(tip_label, autocycler_contig_id, rep_type.s.) |>
  filter(rep_type.s.!= "-") |>
  arrange(rep_type.s.) |>
  group_by(tip_label, autocycler_contig_id) |>
  summarise(contig_rep_types = paste(rep_type.s., collapse = ",")) |> # makes no difference to collapse rep types by contig, as 991 annotations both before and after. 
  ungroup() |>
  group_by(contig_rep_types) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  ungroup()
#View(contig_reps)
sum(contig_reps$count) # 991
nrow(contig_reps)

# get top 20
top10_plasmids <- contig_reps |>
  slice_head(n=10) |>
  pull(contig_rep_types)
top10_plasmids

# make plasmid summary
plasmid_plot_df <- plasmid_data |>
  select(tip_label, autocycler_contig_id, rep_type.s.) |>
  filter(rep_type.s.!= "-") |>
  arrange(rep_type.s.) |>
  group_by(tip_label, autocycler_contig_id) |>
  summarise(contig_rep_types = paste(rep_type.s., collapse = ",")) |> # makes no difference to collapse rep types by contig, as 991 annotations both before and after. 
  ungroup() |>
  mutate(rep_types_cat = case_when(contig_rep_types %in% top10_plasmids ~ contig_rep_types,
                                   grepl("rep_cluster", contig_rep_types) & !grepl(",", contig_rep_types) ~ "other rep_clusters",
                                   grepl("Inc", contig_rep_types) & !grepl(",", contig_rep_types) ~ "other Inc types",
                                   grepl("Col", contig_rep_types) & !grepl(",", contig_rep_types) ~ "other Col-like",
                                   grepl(",", contig_rep_types) ~ "other multi-replicon",
                                   )) |>
  select(-c(autocycler_contig_id, contig_rep_types)) |>
  mutate(present = "present")
#table(plasmid_plot_df$rep_types_cat)
# set order of plasmids
other_levels <- plasmid_plot_df |>
  count(rep_types_cat, sort = TRUE) |>
  filter(grepl("^other", rep_types_cat)) |>
  pull(rep_types_cat)
plasmid_levels <- c(top10_plasmids, other_levels)
plasmid_plot_df <- plasmid_plot_df |>
  mutate(rep_types_cat = factor(rep_types_cat, levels = plasmid_levels))

#View(plasmid_plot_df)
#nrow(plasmids_data)
#View(plasmids_data)
#colnames(amrfinder_data)


# define colours
colors <- c("#0d0c3b", "#ffffff")
region_pal <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
)
region_pal <- rev(region_pal)

# fastBAPS labels
fastbaps_levels <- c(
  "kaerogenes_1", "kaerogenes_2",
  "kosc_10", "kosc_11", "kosc_12", "kosc_13", "kosc_14", "kosc_15",
  "kosc_2", "kosc_3", "kosc_5", "kosc_6", "kosc_7", "kosc_8", "kosc_9",
  "kpsc_1", "kpsc_11", "kpsc_12", "kpsc_13", "kpsc_14", "kpsc_15",
  "kpsc_2", "kpsc_3", "kpsc_4", "kpsc_5", "kpsc_6", "kpsc_7", "kpsc_8", "kpsc_9"
)

num_suffix <- function(x) as.integer(sub(".*_", "", x))

kaero_levels <- fastbaps_levels[grepl("^kaerogenes_", fastbaps_levels)]
kaero_levels <- kaero_levels[order(num_suffix(kaero_levels))]

kosc_levels <- fastbaps_levels[grepl("^kosc_", fastbaps_levels)]
kosc_levels <- kosc_levels[order(num_suffix(kosc_levels))]

kpsc_levels <- fastbaps_levels[grepl("^kpsc_", fastbaps_levels)]
kpsc_levels <- kpsc_levels[order(num_suffix(kpsc_levels))]

# kaerogenes: pale straw -> muted ochre
kaero_cols <- colorRampPalette(
  c("#E6DDAA", "#C7AA72")
)(length(kaero_levels))

# kosc: pale blue -> dusty teal -> muted green
kosc_cols <- colorRampPalette(
  c("#B7D7E8", "#8DBCB6", "#6F9F7F")
)(length(kosc_levels))

# kpsc: pale lavender -> dusty mauve -> muted rose
kpsc_cols <- colorRampPalette(
  c("#D5BDD6", "#C7A2BE", "#C18A96")
)(length(kpsc_levels))

fastbaps_pal <- c(
  setNames(kaero_cols, kaero_levels),
  setNames(kosc_cols, kosc_levels),
  setNames(kpsc_cols, kpsc_levels)
)

if FALSE {
kaerogenes_1 kaerogenes_2       kosc_2       kosc_3       kosc_5       kosc_6       kosc_7 
"#E6DDAA"    "#C7AA72"    "#B7D7E8"    "#B0D2DF"    "#A9CED7"    "#A2C9CF"    "#9BC5C6" 
kosc_8       kosc_9      kosc_10      kosc_11      kosc_12      kosc_13      kosc_14 
"#94C0BE"    "#8DBCB6"    "#88B7AC"    "#83B2A3"    "#7EAD9A"    "#79A891"    "#74A388" 
kosc_15       kpsc_1       kpsc_2       kpsc_3       kpsc_4       kpsc_5       kpsc_6 
"#6F9F7F"    "#D5BDD6"    "#D2B8D2"    "#D0B4CE"    "#CEB0CA"    "#CCACC7"    "#CAA8C3" 
kpsc_7       kpsc_8       kpsc_9      kpsc_11      kpsc_12      kpsc_13      kpsc_14 
"#C8A4BF"    "#C6A0BA"    "#C59CB4"    "#C498AE"    "#C395A8"    "#C291A2"    "#C18D9C" 
kpsc_15 
"#C18A96" 
}

region_df <- region_df |>
  mutate(value = factor(value, levels = levels(value)),
         xlab = factor(""))
fastbaps_df <- fastbaps_df |>
  mutate(value = factor(value, levels = levels(value)),
         xlab = factor(""))

# build tree
tree <- read.tree("kleb_bsi_verticall.newick")
tree <- midpoint_root(tree) #so does this
tree$edge.length[tree$edge.length < 0] <- 0
tree$edge.length <- (tree$edge.length)^(1/3) #this just makes the tree look better
tree_plot <- ggtree(tree, options(ignore.negative.edge=TRUE))

p <- tree_plot %<+% data + 
  geom_tippoint(aes(color = mlst_cat), size =1.2) +  
  scale_color_manual(values = mlst_colours) + 
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_treescale(x = 0, y = length(tree$tip.label) * 0.90, offset = length(tree$tip.label) * 0.01, width = 0.1, color = "black") +  
  theme_tree() 
p

p2 <- p + new_scale_fill() +
  geom_facet(panel = "Region", data = region_df, geom = geom_tile,
      aes( x = xlab, fill = value), width = 1) + 
  scale_fill_manual(
    values = region_pal,
    guide = guide_legend(title = "Region")
  ) +
  
  new_scale_fill() +
  geom_facet(panel = "fast-\nBAPS", data = fastbaps_df, geom = geom_tile,
      aes( x = xlab, fill = value), width = 1) + 
  scale_fill_manual(
    values = fastbaps_pal,
    guide = guide_legend(title = "fastBAPS\ncluster", ncol = 3, byrow = TRUE)
  ) +
  
  new_scale_fill() +
    geom_facet(panel = "CARB", data = carb_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) + 
  geom_facet(panel = "CEPH", data = ceph_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) +
  geom_facet(panel = "B-L", data = beta_lactam_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) +
  geom_facet(panel = "GENT", data = gent_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) + 
  geom_facet(panel = "QUINOLONE", data = quin_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = present),  width = 1) +
  geom_facet(panel = "Plasmids", data = plasmid_plot_df, geom = geom_tile, 
             aes(x = rep_types_cat, fill = present),  width = 1) +
  scale_fill_manual(values = colors, guide = "none") +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), 
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "white", color = "black"),
        legend.position =  "right",
        legend.box = "vertical") +
  scale_x_discrete() +
  labs(color="MLST")
p2

# adjust panel size
gt <- ggplotGrob(p2)
panel_widths <- c(3, 0.4, 0.4, 1.0, 1.1, 1.2, 0.4, 1.5, 1.4)  # same order as geom_facet() calls
panel_idx <- which(grepl("^panel", gt$layout$name))
gt$widths[gt$layout$l[panel_idx]] <- unit(panel_widths, "null")

# preview
grid.newpage()
grid.draw(gt)

# save
ggsave("kleb_bsi_verticall_tree.png", plot = gt, width = 16, height = 10, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 0.b AMR Profiles ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 0.bi E. coli  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data
ecoli_bsi_amrfinder_metadata <- read.csv("neksus_ecoli_bsi_amrfinder_metadata.csv")
colnames(ecoli_bsi_amrfinder_metadata)
ecoli_amrfinder <- ecoli_bsi_amrfinder_metadata |>
  filter(Type == "AMR")

# make istance matrix
amr_mat <- ecoli_amrfinder |>
  filter(!is.na(`Element.symbol`), `Element.symbol` != "") |>
  distinct(sample, feature = `Element.symbol`) |>
  mutate(present = 1L) |>
  pivot_wider(
    names_from = feature,
    values_from = present,
    values_fill = 0L
  ) |>
  column_to_rownames("sample")

# Jaccard distance for binary data
d  <- vegdist(amr_mat, method = "jaccard", binary = TRUE)
hc <- hclust(d, method = "average")
n  <- nrow(amr_mat)

# Function to compute mean silhouette width for a given k
sil_for_k <- function(k, dist_obj) {
  cl <- cutree(hclust(dist_obj, method = "average"), k = k)
  if (length(unique(cl)) < 2 ) {
    return(NA_real_)
  }
  mean(silhouette(cl, dist_obj)[, 3])
}


cut_height_for_k <- function(hc, k) {
  n <- length(hc$order)
  stopifnot(k >= 2, k <= n - 1)
  
  lo_idx <- n - k
  hi_idx <- n - k + 1
  
  lo <- if (lo_idx >= 1) hc$height[lo_idx] else 0
  hi <- hc$height[hi_idx]
  
  (lo + hi) / 2
}

# compute silhouette across k, mapped to tree height
k_vals <- 2:(n - 1)

sil_df <- tibble(
  k = k_vals,
  silhouette = sapply(k_vals, sil_for_k, dist_obj = d)
) %>%
  mutate(height = sapply(k, function(x) cut_height_for_k(hc, x))) %>%
  arrange(height)

best_k <- sil_df %>%
  filter(!is.na(silhouette)) %>%
  slice_max(silhouette, n = 1)

best_h <- best_k$height[[1]]
max_h  <- max(hc$height)

x_lim <- c(0, max_h * 1.05)


# plot top panel: silhouette vs tree height
p_sil <- ggplot(sil_df, aes(x = height, y = silhouette)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.8) +
  geom_vline(xintercept = best_h, linetype = "dashed", color = "red", linewidth = 1.0) +
  annotate(
    "text",
    x = best_h,
    y = max(sil_df$silhouette, na.rm = TRUE),
    label = paste0("best k = ", best_k$k),
    color = "red",
    vjust = -0.6,
    hjust = 0,
    size = 3.5
  ) +
  scale_x_continuous(limits = x_lim, expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "Dendrogram height",
    y = "Mean silhouette width"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 0, 5)
  )

# bottom panel: dendrogram, no tip labels
p_tree <- ggtree(as.phylo(hc), layout = "rectangular") +
  geom_vline(xintercept = best_h, linetype = "dashed", color = "red", linewidth = 1.0) +
  scale_x_continuous(limits = x_lim, expand = expansion(mult = c(0, 0.02))) +
  theme_tree2() +
  theme(
    text = element_text(size = 11),
    plot.margin = margin(0, 5, 5, 5)
  )


# combine
final_plot <- p_sil / p_tree + plot_layout(heights = c(1, 4))
final_plot














# Compute silhouette over a range of k
k_vals <- 2:(n - 1)

sil_df <- tibble(
  k = k_vals,
  silhouette = sapply(k_vals, sil_for_k, dist_obj = d)
)
#View(sil_df)

# Best k by silhouette
best_k <- sil_df %>%
  filter(!is.na(silhouette)) %>%
  slice_max(silhouette, n = 1)

best_k # 698


# Plot dendrogram
p_sil <- ggplot(sil_df, aes(x = k, y = silhouette)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.8) +
  geom_vline(xintercept = best_k$k, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    x = "Number of clusters (k)",
    y = "Mean silhouette width"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 0, 5)
  )

# save
ggsave("ecoli_bsi_silhouette_plot.png", p_sil, units = "in" , width = 6, height = 6, dpi = 300)

p_tree <- ggtree(as.phylo(hc), layout = "rectangular") +
  geom_tiplab(size = 2, align = TRUE, linesize = 0.2, linetype = "dotted") +
  theme_tree2() +
  theme(
    text = element_text(size = 11),
    plot.margin = margin(0, 5, 5, 5)
  ) 

p_tree

# If the tree is still upside down in your setup, swap scale_x_reverse()
# for scale_y_reverse().

final_plot <- p_sil / p_tree + plot_layout(heights = c(1, 4))
final_plot

# # profiles cut at best silhouette level
groups <- cutree(hc, k = best_k$k)

profile_df <- tibble(
  sample = names(groups),
  AMR_profile = paste0("Profile_", groups)
)

#View(profile_df)
table(profile_df$AMR_profile)

# plot denndrogram with line indicating cut level
# Number of desired clusters
k <- best_k$k

# Height that gives k clusters
h <- mean(hc$height[c(length(hc$height) - k + 1,
                      length(hc$height) - k + 2)])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 0.bii Klebsiella  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
kleb_bsi_amrfinder_metadata <- read.csv("neksus_kleb_bsi_amrfinder_metadata.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 0.c Verticall distance summary ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 0.bi E. coli  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# read in .tsv files of pairwise distances
# generate genome sizes
generate_filenames <- function(df, sample_column, run_column) {
  sample_col <- sym(sample_column)
  run_col <- sym(run_column)
  df <- df |>
    mutate(
      dir = !!run_col,
      #dir = gsub("run0", "assembly_main_pipeline_no_db", dir),
      #dir = gsub("microbesng_failed", "microbesng_failed_assembly", dir),
      filename = case_when(!!run_col == "run0" ~ paste0("pilot_", !!sample_col),
                           TRUE ~ !!sample_col)) |>
    mutate(filename = case_when(run != "run0" ~ gsub("AHA", "AH-A", filename),
                                run == "run0" ~ filename),
           filename = case_when(run != "run0" ~ gsub("AHB", "AH-B", filename),
                                run == "run0" ~ filename)
    ) |>
    mutate(filename = gsub("AB90", "AB90-A", filename)) |>
    mutate(filename = gsub("AF210", "AF210-A", filename)) |>
    mutate(filename = paste0(dir, "_" , filename)) |>
    select(filename, everything())
  return(df)
  
}


# ecoli
ecoli_bsi_genome_sizes <- generate_filenames(ecoli_bsi_samples_metadata, "sequencing_id", "run") |>
  select(c(filename, Genome_Size))
setDT(ecoli_bsi_genome_sizes)
#View(ecoli_bsi_genome_sizes)

# kleb
kleb_bsi_genome_sizes <- generate_filenames(kleb_bsi_samples_metadata, "sequencing_id", "run") |>
  select(c(filename, Genome_Size))
setDT(kleb_bsi_genome_sizes)
#View(kleb_bsi_genome_sizes)

# fread for fast read-in of big files
ecoli_dt <- fread("ecoli_bsi_verticall.tsv", sep = "\t")

# create row-wise canonical pair id and absolute genomic distance
ecoli_dt[, pair_id := paste(pmin(assembly_a, assembly_b), pmax(assembly_a, assembly_b), sep = "|")]
ECOLI_dt[, absolute_genomic_distance := mean_window_distance * window_size * window_count]

# reorder columns so pair_id is first
setcolorder(dt, c("pair_id", setdiff(names(dt), "pair_id")))

# sort by pair_id
setorder(dt, pair_id)

# ecoli
ecoli_bsi_verticall_distances <- read_delim("ecoli_bsi_verticall.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
nrow(ecoli_bsi_verticall_distances)
colnames(ecoli_bsi_verticall_distances)
View(ecoli_bsi_verticall_distances)
rm(ecoli_bsi_verticall)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 0.bii Klebsiella  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# fread is much faster than read.csv/read_tsv for big files
dt <- fread("kleb_bsi_verticall.tsv", sep = "\t")
#setDT(dt)
#setDT(kleb_bsi_genome_sizes)

# create row-wise canonical pair id and absolute genomic distance
dt[, pair_id := paste(pmin(assembly_a, assembly_b), pmax(assembly_a, assembly_b), sep = "|")]
dt[, absolute_genomic_distance := mean_window_distance * window_size * window_count]
# reorder columns so pair_id is first
setcolorder(dt, c("pair_id", setdiff(names(dt), "pair_id")))
# sort by pair_id
setorder(dt, pair_id)
names(dt)

# compare between pairs in either direction
# get vertical genome legths

# helper functions
pct_to_num <- function(x) {
  as.numeric(sub("%$", "", trimws(x)))
}

rm_to_num <- function(x) {
  x <- trimws(x)
  x[x == "undefined" | x == ""] <- NA_character_
  as.numeric(x)
}

region_length <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  parts <- trimws(parts[nzchar(parts)])
  
  vals <- vapply(parts, function(p) {
    m <- regexec("^[^:]+:(\\d+)-(\\d+)$", p)
    mm <- regmatches(p, m)[[1]]
    if (length(mm) == 0) return(NA_real_)
    as.numeric(mm[3]) - as.numeric(mm[2])
  }, numeric(1))
  
  if (all(is.na(vals))) NA_real_ else sum(vals, na.rm = TRUE)
}


#1. convert character columns to numeric
dt[, alignments_vertical_fraction := pct_to_num(alignments_vertical_fraction)]
dt[, alignments_horizontal_fraction := pct_to_num(alignments_horizontal_fraction)]
dt[, `r/m` := rm_to_num(`r/m`)]

# 2) add vertical region lengths from the region strings
dt[, assembly_a_vertical_region_length := vapply(assembly_a_vertical_regions, region_length, numeric(1))]
dt[, assembly_b_vertical_region_length := vapply(assembly_b_vertical_regions, region_length, numeric(1))]

#  3) join genome sizes to assembly_a and assembly_b 
genome_lengths[, genome_size := as.numeric(genome_size)]

dt[genome_lengths, on = .(assembly_a = sample), assembly_a_genome_size := i.genome_size]
dt[genome_lengths, on = .(assembly_b = sample), assembly_b_genome_size := i.genome_size]

#  4) compute vertical lengths from genome size * vertical fraction 
# fraction column is stored as percent, so divide by 100
dt[, assembly_a_vertical_length := assembly_a_genome_size * get(vert_frac_col) / 100]
dt[, assembly_b_vertical_length := assembly_b_genome_size * get(vert_frac_col) / 100]

#  5) compare the two vertical-length methods 
dt[, assembly_a_vertical_length_diff :=
     assembly_a_vertical_region_length - assembly_a_vertical_length]

dt[, assembly_b_vertical_length_diff :=
     assembly_b_vertical_region_length - assembly_b_vertical_length]

#  6) make an unordered pair key so (a,b) == (b,a) 
dt[, pair_id := paste(pmin(assembly_a, assembly_b),
                      pmax(assembly_a, assembly_b),
                      sep = "|")]

#  7) summary for primary rows only 
primary_dt <- dt[result_level == "primary"]

pair_summary <- primary_dt[, .(
  n_rows = .N,
  
  mean_mean_vertical_window_distance =
    mean(mean_vertical_window_distance, na.rm = TRUE),
  diff_mean_vertical_window_distance =
    diff(range(mean_vertical_window_distance, na.rm = TRUE)),
  
  mean_median_vertical_window_distance =
    mean(median_vertical_window_distance, na.rm = TRUE),
  diff_median_vertical_window_distance =
    diff(range(median_vertical_window_distance, na.rm = TRUE)),
  
  mean_mean_vertical_distance =
    mean(mean_vertical_distance, na.rm = TRUE),
  diff_mean_vertical_distance =
    diff(range(mean_vertical_distance, na.rm = TRUE))
), by = pair_id]

# optional: keep only pairs that appear more than once in primary rows
pair_summary <- pair_summary[n_rows > 1]

# 8) histogram of the within-pair differences
diff_long <- melt(
  pair_summary,
  id.vars = c("pair_id", "n_rows"),
  measure.vars = patterns("^diff_"),
  variable.name = "metric",
  value.name = "difference"
)

ggplot(diff_long, aes(x = difference)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ metric, scales = "free") +
  theme_bw() +
  labs(x = "Within-pair difference", y = "Count")

# 9) write outputs 
fwrite(dt, "ecoli_bsi_verticall_augmented.tsv", sep = "\t", na = "NA")
fwrite(pair_summary, "primary_pair_summary.tsv", sep = "\t", na = "NA")


# plot
ggplot(dt) +  geom_histogram(aes(x = mean_window_distance), bins = 100) # trimodal distribution
ggplot(dt) +  geom_histogram(aes(x = absolute_genomic_distance), bins = 100)
ggplot(dt) +  geom_histogram(aes(x = peak_mass), bins = 100) +
  scale_x_log10()
ggplot(dt) +  geom_histogram(aes(x = peak_window_distance), bins = 100) # trimodal distribution
ggplot(dt) +  geom_histogram(aes(x = mean_vertical_distance), bins = 100) # trimodal distribution
ggplot(dt) +  geom_histogram(aes(x = mean_horizontal_distance), bins = 100) # trimodal distribution


ggplot(dt) +  geom_histogram(aes(x = mean_vertical_window_distance), bins = 100) # trimodal distribution
ggplot(dt) +  geom_histogram(aes(x = median_vertical_window_distance), bins = 100) # trimodal distribution
ggplot(dt) +  geom_histogram(aes(x = mean_vertical_distance*5000000), bins = 100) +
  scale_x_log10()# trimodal distribution





kleb_bsi_verticall_distances <- read_delim("kleb_bsi_verticall.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
nrow(kleb_bsi_verticall_distances)
colnames(kleb_bsi_verticall_distances)
View(kleb_bsi_verticall_distances)
# mutate to be able to identif same pair by other way around
kleb_bsi_verticall_distances <- kleb_bsi_verticall_distances |>
  mutate(pair_id = paste(sort(c(assembly_a, assembly_b)), collapse = "|"))|>
  select(pair_id, everything()) |>
  arrange(pair_id)

# add absolute genomic distnace
kleb_bsi_verticall_distances <- kleb_bsi_verticall_distances |>
   mutate(absolute_genomic_distance = mean_window_distance * window_size * window_count)


# plot histogram of lengths
ggplot(data = kleb_bsi_verticall_distances) +
  geom_histogram(aes(x = median_genomic_distance))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# load ncbi db data, and filter only E. coli and Kleb relevant genes
library(readr)
refgenes <- read_delim("C:/Users/dnagy/Downloads/refgenes.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
View(refgenes)
dim(refgenes) # 11357 x 16
str(refgenes)
colnames(refgenes)
table(refgenes$`Whitelisted taxa`, useNA = "ifany")
table(refgenes$Type, useNA = "ifany") # 10064 ARGs

# filter only those relevant to E. coli
ecoli_refgenes <- refgenes |>
  filter(Type == "AMR") |>
  filter(!grepl("Escherichia", `Blacklisted taxa` )) |>
  filter((grepl("Escherichia", `Whitelisted taxa`) | is.na(`Whitelisted taxa`)))
nrow(ecoli_refgenes)# 825
# filter only those relevant to Klebsiella
kleb_refgenes <- refgenes |>
  filter(Type == "AMR") |>
  filter(!grepl("Klebsiella", `Blacklisted taxa` )) |>
  filter((grepl("Klebsiella", `Whitelisted taxa`) | is.na(`Whitelisted taxa`)))
nrow(kleb_refgenes)# 8556



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Framework 1: Overall isolate-level features (MLST AND fastbaps_L3 clusters) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Chinese restaurant process function to estimate novel mass
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
  
  # CRP / Ewens likelihood up to a theta-free constant:
  # p(counts | theta) ∝ theta^K * Gamma(theta) / Gamma(theta + N)
  log_lik <- K * log(theta_grid) + lgamma(theta_grid) - lgamma(theta_grid + N)
  
  # Gamma prior on theta
  log_prior <- dgamma(theta_grid, shape = prior_shape, rate = prior_rate, log = TRUE)
  
  log_post <- log_lik + log_prior
  
  # normalize on the grid
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

# function to get novel mass from posterior estimate of CPR
set_alpha_novel_sum_from_crp <- function(alpha_named, novelty_prob_hat) {
  alpha_obs_total <- sum(alpha_named, na.rm = TRUE)
  alpha_obs_total * novelty_prob_hat / (1 - novelty_prob_hat)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTION: isolate_bayes_dirichlet
# Inputs:
#   df: data.frame with columns 'mlst_profile' and 'count'
#   alpha_named: named numeric vector of alpha pseudo-counts for each observed MLST
#                (names must match df$mlst_profile). May be 0 for some types.
#   alpha_novel: numeric pseudo-count for the grouped "novel" category (can be non-integer)
#   B: number of posterior draws
#   use_dirichlet: if TRUE, use direct Dirichlet sampling (supports non-integer alphas)
#                         if FALSE, attempt to use bayes by adding integer pseudo-observations
# Returns:
#   A list with:
#     draws: matrix B x (K+1) of posterior frequency draws (columns named by mlst + "NOVEL")
#     summary_df: data.frame summarizing median, sd, quantiles per category
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# functions
isolate_bayes_dirichlet <- function(df,
                           feature_col = "kleborate_mlst",
                           alpha_named = NULL,
                           alpha_novel_sum = 1,
                           alpha_novel_num = 1,
                           B = 1000,
                           seed = 2026) {
  # dependencies: dplyr, tidyr, bayes
  stopifnot(requireNamespace("dplyr", quietly = TRUE),
            requireNamespace("tidyr", quietly = TRUE))
  set.seed(seed)
  
  # basic checks
  stopifnot(is.data.frame(df))
  stopifnot(is.character(feature_col) && length(feature_col) == 1)
  stopifnot(all(c(feature_col, "count") %in% colnames(df)))
  stopifnot(is.numeric(alpha_novel_sum), length(alpha_novel_sum) == 1, alpha_novel_sum >= 0)
  stopifnot(is.numeric(alpha_novel_num), length(alpha_novel_num) == 1, alpha_novel_num >= 0)
  stopifnot(B >= 1)
  
  alpha_novel_sum <- as.numeric(alpha_novel_sum)
  alpha_novel_num <- as.integer(alpha_novel_num)
  
  # collapse to unique feature rows (group by supplied column name)
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
    if(length(alpha_named) != K){
      stop("alpha_named must have the same length as unique features in df[[feature_col]].")
    }
    # if names are present, reorder to match observed features; otherwise assume supplied in the right order
    if (!is.null(names(alpha_named)) && all(names_counts %in% names(alpha_named))) {
      alpha_named <- alpha_named[names_counts]
    } else {
      names(alpha_named) <- names_counts
    }
  }
  
  alpha_named <- as.numeric(alpha_named)
  
  # add novel categories if sum and num are not 0
  add_novel <- (alpha_novel_sum > 0) && (alpha_novel_num > 0)
  
  if (add_novel) {
    alpha_novel_each <- alpha_novel_sum / alpha_novel_num
    alpha_novel_vec <- rep(alpha_novel_each, alpha_novel_num)
    novel_names <- paste0("NOVEL_", seq_len(alpha_novel_num))
    
    shapes_all <- c(n_k + alpha_named, alpha_novel_vec)
    out_names <- c(names_counts, novel_names)
    
  } else {
    alpha_novel_each <- 0
    alpha_novel_vec <- numeric(0)
    novel_names <- character(0)
    
    shapes_all <- n_k + alpha_named
    out_names <- names_counts
    
  }
  
  
  # Exact Dirichlet sampling: draw directly with rgamma:
  shapes_obs <- n_k + alpha_named
  shapes_all <- c(shapes_obs, alpha_novel_vec)
  out_names <- c(names_counts, novel_names)
  
  # exact Dirichlet via independent gamma draws
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
  
  return(list(
    draws = draws,
    summary_df = summary_df,
    alpha_novel_sum = alpha_novel_sum,
    alpha_novel_num = alpha_novel_num,
    alpha_novel_each = alpha_novel_each
  ))
  
}
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# helper for computing richness curves
prep_bootstrap_draws <- function(draws) {
  x <- as.matrix(draws)
    # sort each row once (ascending)
  x_sorted <- t(apply(x, 1, sort))
    # cumulative sums for the mass curve
  cs <- t(apply(x_sorted, 1, cumsum))
    list(
    x_sorted = x_sorted,
    cumsums  = cs,
    totals   = rowSums(x_sorted)
  )
}

default_f_grid <- function() {
  sort(unique(c(0, round(10^seq(-5, 0, length.out = 200), 8))))
}

compute_species_richness_curve <- function(prep, f_grid = NULL) {
  x_sorted <- prep$x_sorted
  B <- nrow(x_sorted)
  K <- ncol(x_sorted)
  if (is.null(f_grid)) f_grid <- default_f_grid()
  F <- length(f_grid)
  counts <- matrix(0L, nrow = B, ncol = F)
  
  for (i in seq_len(B)) {
    # number of values <= f for each threshold
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
    sd = matrixStats::colSds(counts),
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
    # idx = number of values <= f
    idx <- findInterval(f_grid, x_sorted[i, ], left.open = TRUE, rightmost.closed = TRUE)
    
    # cumulative sum up to idx, with 0 prepended
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# sensitivity analysis wrappers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
alpha_vals <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)

# Fit metric: RMSE on log10 prevalence scale
compute_log_rmse <- function(actual, estimate, eps = NULL, weights = NULL) {
  stopifnot(length(actual) == length(estimate))
  
  if (is.null(eps)) {
    eps <- 1e-12
  }
  
  err <- log10(actual + eps) - log10(estimate + eps)
  sq <- err^2
  
  if (is.null(weights)) {
    return(sqrt(mean(sq, na.rm = TRUE)))
  }
  
  weights <- weights / sum(weights, na.rm = TRUE)
  sqrt(sum(weights * sq, na.rm = TRUE))
}

# Optional helper: draw a darker shade sequence
make_green_shades <- function(n) {
  grDevices::colorRampPalette(c("#d9f0d3", "#7fc97f", "#238b45"))(n)
}

make_orange_shades <- function(n) {
  grDevices::colorRampPalette(c("#fee6ce", "#fdae6b", "#e6550d"))(n)
}

# Main wrapper
run_mlst_sensitivity_grid <- function(df,
                                      feature_col = "kleborate_mlst",
                                      alpha_named_vals = c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000),
                                      alpha_novel_sum_vals = c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000),
                                      alpha_novel_num_vals = c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000),
                                      B = 1000,
                                      seed = 2026,
                                      out_dir = ".",
                                      prefix = "mlst_sensitivity",
                                      dataset_label = NA_character_,
                                      save_csv = TRUE,
                                      keep_draws = FALSE) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c(feature_col, "count") %in% names(df)))
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Collapse to unique MLST rows
  df <- df |>
    dplyr::group_by(.data[[feature_col]]) |>
    dplyr::summarise(count = sum(.data[["count"]]), .groups = "drop")
  
  N <- sum(df$count)
  obs_lookup <- df |>
    dplyr::transmute(feature = .data[[feature_col]], actual = count / N)
  
  grid <- tidyr::crossing(alpha_named = alpha_named_vals, alpha_novel_sum = alpha_novel_sum_vals, alpha_novel_num = alpha_novel_num_vals) |>
    dplyr::mutate(grid_id = dplyr::row_number())
  
  feature_list <- vector("list", nrow(grid))
  fit_list <- vector("list", nrow(grid))
  mass_list <- vector("list", nrow(grid))
  richness_list <- vector("list", nrow(grid))
  draws_list <- if (keep_draws) vector("list", nrow(grid)) else NULL
  
  for (i in seq_len(nrow(grid))) {
    a_named <- grid$alpha_named[i]
    a_novel_sum <- grid$alpha_novel_sum[i]
    a_novel_num <- grid$alpha_novel_num[i]
    
    alpha_named_vec <- stats::setNames(rep(a_named, nrow(df)), as.character(df[[feature_col]]))
    
    fit <- isolate_bayes_dirichlet(
      df = df,
      feature_col = feature_col,
      alpha_named = alpha_named_vec,
      alpha_novel_sum = a_novel_sum,
      alpha_novel_num = a_novel_num,
      B = B,
      seed = seed + i
    )
    
    feat <- fit$summary_df |>
      #dplyr::rename(feature = dplyr::all_of(feature_col)) |>
      dplyr::filter(!grepl("NOVEL", feature)) |>
      dplyr::left_join(obs_lookup, by = "feature") |>
      dplyr::mutate(
        alpha_named = a_named,
        alpha_novel_sum = a_novel_sum,
        alpha_novel_num = a_novel_num,
        dataset = dataset_label,
        estimate = median,
        lo = q2.5,
        hi = q97.5
      )
    
    metric <- tibble::tibble(
      alpha_named = a_named,
      alpha_novel_sum = a_novel_sum,
      alpha_novel_num = a_novel_num,
      dataset = dataset_label,
      rmse_log10 = compute_log_rmse(
        actual = feat$actual,
        estimate = feat$estimate,
        eps = 0.5 / N
      ),
      mae_log10 = mean(abs(log10(feat$actual + 0.5 / N) - log10(feat$estimate + 0.5 / N)), na.rm = TRUE)
    )
    
    # prep draws
    prep <- prep_bootstrap_draws(fit$draws)
    
    mass <- compute_mass_curve(prep) |>
      dplyr::mutate(
        alpha_named = a_named,
        alpha_novel_sum = a_novel_sum,
        alpha_novel_num = a_novel_num,
        dataset = dataset_label
      )
    
    richness <- compute_species_richness_curve(prep) |>
      dplyr::mutate(
        alpha_named = a_named,
        alpha_novel_sum = a_novel_sum,
        alpha_novel_num = a_novel_num,
        dataset = dataset_label
      )
    
    feature_list[[i]] <- feat
    fit_list[[i]] <- metric
    mass_list[[i]] <- mass
    richness_list[[i]] <- richness
    if (keep_draws) draws_list[[i]] <- fit$draws
  }
  
  feature_summary <- dplyr::bind_rows(feature_list)
  fit_metrics <- dplyr::bind_rows(fit_list)
  mass_curves <- dplyr::bind_rows(mass_list)
  richness_curves <- dplyr::bind_rows(richness_list)
  
  if (save_csv) {
    readr::write_csv(feature_summary, file.path(out_dir, paste0(prefix, "_feature_summary.csv")))
    readr::write_csv(fit_metrics, file.path(out_dir, paste0(prefix, "_fit_metrics.csv")))
    readr::write_csv(mass_curves, file.path(out_dir, paste0(prefix, "_mass_curves.csv")))
    readr::write_csv(richness_curves, file.path(out_dir, paste0(prefix, "_richness_curves.csv")))
  }
  
  list(
    feature_summary = feature_summary,
    fit_metrics = fit_metrics,
    mass_curves = mass_curves,
    richness_curves = richness_curves,
    draws = draws_list,
    alpha_vals = alpha_vals
  )
}

# Heatmap: lower RMSE = better fit
plot_fit_heatmap <- function(fit_metrics, alpha_named_vals, alpha_novel_sum_vals, alpha_novel_num_vals, fill_val = "rmse_log10", fill_label = "log-RMSE") {
  fill_sym <- rlang::sym(fill_val)
  
  # make facet labels
  facet_labs <- setNames(
    paste0("alpha[1:k] * ': ' * ", format(alpha_named_vals, trim = TRUE, scientific = FALSE)),
    as.character(alpha_named_vals)
  )
  
  fit_metrics |>
    dplyr::mutate(
      alpha_named = factor(alpha_named, levels = alpha_named_vals),
      alpha_novel_sum = factor(alpha_novel_sum, levels = alpha_novel_sum_vals),
      alpha_novel_num = factor(alpha_novel_num, levels = alpha_novel_num_vals)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = alpha_novel_sum, y = alpha_novel_num, fill = !!fill_sym)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_viridis_c(direction = -1, option = "C") +
    ggplot2::facet_wrap(
      ~ alpha_named,
      labeller = ggplot2::as_labeller(facet_labs, label_parsed)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = expression(alpha[novel]~Mass),
      y = expression(alpha[novel]~Number),
      fill = fill_label
    ) +
    ggplot2::theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
}

# Observed vs posterior median panel plot
plot_obs_vs_pred <- function(feature_summary, alpha_named_vals, alpha_novel_sum_vals, alpha_novel_num_vals) {
  
  # discrete steel blue palette
  n_cols <- length(alpha_novel_num_vals)
  steelblue_pal <- grDevices::colorRampPalette(c("steelblue4", "steelblue2", "lightsteelblue1"))(n_cols)
  
  # add labels
  named_labeller <- ggplot2::as_labeller(
    setNames(paste0("alpha[1:k]: ", alpha_named_vals),
             alpha_named_vals),label_parsed)
  
  novel_sum_labeller <- ggplot2::as_labeller(
    setNames(paste0("alpha[novel]~Mass: ", alpha_novel_sum_vals),
      alpha_novel_sum_vals),label_parsed)
  
  feature_summary |>
    dplyr::mutate(
      alpha_named = factor(alpha_named, levels = alpha_named_vals),
      alpha_novel_sum = factor(alpha_novel_sum, levels = rev(alpha_novel_sum_vals)),
      alpha_novel_num = factor(alpha_novel_num, levels = rev(alpha_novel_num_vals))
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = actual, y = estimate, colour = alpha_novel_num)) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", linewidth = 0.4, color = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo, ymax = hi),
      width = 0, alpha = 0.7, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.7, alpha = 0.9) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(
      values = setNames(steelblue_pal, alpha_novel_num_vals),
      name = expression(alpha[novel]~Number)
    ) +
    ggplot2::facet_grid(
      alpha_novel_sum ~ alpha_named,
      labeller = ggplot2::labeller(
        alpha_named = named_labeller,
        alpha_novel_sum = novel_sum_labeller
      )
    ) +
    ggplot2::labs(
      x = "Observed prevalence",
      y = "Posterior median prevalence"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

plot_curve_facet <- function(curve_df,
                             alpha_named_vals, 
                             alpha_novel_sum_vals, 
                             alpha_novel_num_vals, 
                             value_col, 
                             lo_col, hi_col, 
                             y_lab, 
                             title, 
                             curve_type = c("mass", "richness")) {
  
  curve_type <- match.arg(curve_type) 
  datasets <- unique(curve_df$dataset) 
  datasets <- datasets[!is.na(datasets)] 
  alpha_named_levels <- sort(unique(curve_df$alpha_named)) 
  alpha_named_labels <- as.character(alpha_named_levels) 
  alpha_novel_sum_levels <- sort(unique(curve_df$alpha_novel_sum)) 
  alpha_novel_sum_labels <- as.character(alpha_novel_sum_levels) 
  alpha_novel_num_levels <- sort(unique(curve_df$alpha_novel_num)) 
  alpha_novel_num_labels <- as.character(alpha_novel_num_levels) 
  
  pal <- c() 
  if ("Escherichia" %in% datasets) { 
    ecoli_cols <- make_green_shades(length(alpha_novel_levels)) 
    names(ecoli_cols) <- paste("Escherichia", alpha_novel_num_labels, sep = "__") 
    pal <- c(pal, ecoli_cols) 
  } 
  if ("Klebsiella" %in% datasets) { 
    kleb_cols <- make_orange_shades(length(alpha_novel_levels)) 
    names(kleb_cols) <- paste("Klebsiella", alpha_novel_num_labels, sep = "__") 
    pal <- c(pal, kleb_cols) 
  } 
  
  curve_df <- curve_df |> 
    dplyr::mutate( dataset = factor(dataset, 
                                    levels = c("Escherichia", "Klebsiella")), 
                   alpha_named = factor(alpha_named, levels = alpha_named_vals), 
                   alpha_novel_sum = factor(alpha_novel_sum, levels = alpha_novel_sum_vals), 
                   alpha_novel_num = factor(alpha_novel_num, levels = alpha_novel_num_vals), 
                   line_id = paste(dataset, alpha_novel_num, sep = "__") ) 
  ecoli_levels <- paste("Escherichia", alpha_named_vals , sep = "__") 
  kleb_levels <- paste("Klebsiella", alpha_named_vals , sep = "__") 
  line_levels <- c(ecoli_levels, kleb_levels) 
  
  curve_df <- curve_df |> 
    dplyr::mutate(line_id = factor(line_id, levels = line_levels)) 
  legend_labels <- c(paste("E. coli α =", alpha_named_vals), paste("Klebsiella α =", alpha_named_vals)) 
  
  ggplot2::ggplot(curve_df, ggplot2::aes(x = f, y = .data[[value_col]], group = line_id, color = line_id, fill = line_id ) ) + 
    ggplot2::geom_ribbon( ggplot2::aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]), alpha = 0.12, colour = NA ) + 
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::facet_wrap( ggplot2::vars(dataset, alpha_novel), ncol = 7, labeller = ggplot2::label_both, scales = "free_y" ) + 
    ggplot2::scale_color_manual(values = pal, breaks = line_levels, labels = legend_labels, name = expression(alpha[novel]~Number)) + 
    ggplot2::scale_fill_manual(values = pal, guide = "none") + 
    ggplot2::scale_x_log10() + 
    ggplot2::facet_grid(
      alpha_novel_sum ~ alpha_named,
      labeller = ggplot2::labeller(
        alpha_named = named_labeller,
        alpha_novel_sum = novel_sum_labeller
      )
    ) +
    ggplot2::labs(x = "Threshold f", y = y_lab, title = title, color = "Dataset / Number of novel features" ) + 
    ggplot2::theme_minimal(base_size = 12) + 
    ggplot2::theme( strip.text = ggplot2::element_text(face = "bold"), legend.position = "right" ) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.1a E. coli MLSTs #### 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare data - E. coli
ecoli_bsi_mlst_df <- ecoli_bsi_samples_metadata |>
  dplyr::rename(kleborate_mlst = escherichia__mlst_achtman__ST) |>
  dplyr::group_by(kleborate_mlst) |>
  dplyr::summarise(count = n(), .groups = "drop") |>
  dplyr::select(kleborate_mlst, count) |>
  dplyr::arrange(count, kleborate_mlst)
#View(ecoli_bsi_mlst_df)
#table(ecoli_bsi_mlst_df$count)
#nrow(ecoli_bsi_mlst_df) # 263
#write.csv(ecoli_bsi_mlst_df, "NEKSUS_ecoli_bsi_mlst_df.csv", row.names = FALSE, quote = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 1.1ai Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
alpha_named_vals <- c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
alpha_novel_sum_vals <- c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
alpha_novel_num_vals <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)

res_ecoli <- run_mlst_sensitivity_grid(
   ecoli_bsi_mlst_df,
   feature_col = "kleborate_mlst",
   alpha_named_vals = alpha_named_vals,
   alpha_novel_sum_vals = alpha_novel_sum_vals,
   alpha_novel_num_vals = alpha_novel_num_vals,
   B = 1000,
   seed = 2026,
   out_dir = "sens_results/ecoli_mlst",
   prefix = "ecoli",
   dataset_label = "Escherichia",
   save_csv = TRUE,
   keep_draws = FALSE
 )
#View(res_ecoli$fit_metrics)
#View(res_ecoli$feature_summary)
res_ecoli <- c()
res_ecoli$fit_metrics <- read.csv("sens_results/ecoli_mlst/ecoli_fit_metrics.csv")
res_ecoli$feature_summary <- read.csv("sens_results/ecoli_mlst/ecoli_feature_summary.csv")

# rmse heatmap
heatmap_ecoli <- plot_fit_heatmap(res_ecoli$fit_metrics, 
                                  alpha_named_vals = alpha_named_vals,
                                  alpha_novel_sum_vals = alpha_novel_sum_vals,
                                  alpha_novel_num_vals = alpha_novel_num_vals,
                                  fill_val = "rmse_log10" , fill_label = "log-RMSE")
heatmap_ecoli
ggsave("sens_results/ecoli_mlst_simpleBB_sensitivity_heatmap.png", heatmap_ecoli, width = 10, height = 8, dpi = 300)

# sparse panel plot
alpha_named_min <- c(0.5, 1, 10, 100)
alpha_novel_sum_min <- c(0, 1, 10, 100, 250, 500)
alpha_novel_num_min <- c(0, 1, 10, 100, 250, 500)

thinned <- res_ecoli$feature_summary |>
  filter(alpha_named %in% alpha_named_min 
         & alpha_novel_sum %in% alpha_novel_sum_min
         & alpha_novel_num %in% alpha_novel_num_min)
panel_ecoli <- plot_obs_vs_pred(thinned, alpha_named_min, alpha_novel_sum_min, alpha_novel_num_min)
panel_ecoli
ggsave("sens_results/ecoli_mlst_simpleBB_sensitivity_panel_thinned.png", panel_ecoli, width = 10, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 1.1aii Run single prior parameter set ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
df <- ecoli_bsi_mlst_df
counts <- df$count
N <- sum(df$count)
K <- length(df$count) # 263  # num unique features
f1 <- sum(df$count == 1)  # 161       # number of singletons
q_hat <- f1 / N  # 11% Good-Turing first-order - proportion of singletons

# Estimate novel/ unseen mass using Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(counts)
crp_fit$summary_df # mean: theta = 62.5 -> 4% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.0407

# set priors
alpha_prior <- 1
alpha_named <- rep(alpha_prior, K)  # set uninformative priors
names(alpha_named) <- df$kleborate_mlst
alpha_novel_num <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum <- alpha_novel_num * alpha_prior

# apply bayes function to data:
ecoli_bsi_mlst_bayes <- isolate_bayes_dirichlet(ecoli_bsi_mlst_df,
                                           feature_col = "kleborate_mlst",
                                           alpha_named = alpha_named,
                                           alpha_novel_sum = alpha_novel_sum,
                                           alpha_novel_num = alpha_novel_sum,
                                           B = 10000)

#View(ecoli_bsi_mlst_bayes$summary_df)
#print(res$draws)

# save
saveRDS(ecoli_bsi_mlst_bayes, "model_results/ecoli_bsi_mlst_bayes.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#ecoli_bsi_mlst_bayes <- readRDS("model_results/ecoli_bsi_mlst_bayes.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# define f-grid based on sample sizes:
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
#f_grid_90 <- 1-((1-0.90)^(1/n_grid))
#f_grid_95 <- 1-((1-0.95)^(1/n_grid))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))
#f_grid <- unique(c(f_grid_90, f_grid_95, f_grid_99))
#f_grid <- c(seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))

# cumulative species richness curve (>=f; number and proportion of MLSTs at least as frequent as f)
prep_draws <- prep_bootstrap_draws(ecoli_bsi_mlst_bayes$draws)
ecoli_species_richness_df <- compute_species_richness_curve(prep_draws, f_grid = f_grid_99)
#View(ecoli_species_richness_df)

# add sample sizes
ecoli_species_richness_df <- ecoli_species_richness_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save species_richness df
write.csv(ecoli_species_richness_df, "model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#ecoli_species_richness_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
ecoli_mass_df <- compute_mass_curve(prep_draws, f_grid = f_grid_99)
#View(ecoli_mass_df)

# add sample sizes
ecoli_mass_df <- ecoli_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
         min_sample_95 = log(1-0.95)/log(1-f),
         min_sample_99 = log(1-0.99)/log(1-f),
         Genus = "Escherichia",
         Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(ecoli_mass_df, "model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#ecoli_mass_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# illustrative results plot for figure
post_density_plot <- ggplot(ecoli_bsi_mlst_bayes$draws) +
  geom_density(
    aes(x = ST10),
    fill = "lightsteelblue4",
    colour = "lightsteelblue4",
    alpha = 0.5,
    linewidth = 0.8
  ) +
  geom_density(
    aes(x = ST131),
    fill = "#4292c6",
    colour = "#4292c6",
    alpha = 0.5,
    linewidth = 0.8
  ) +
  geom_density(
    aes(x = ST88),
    fill = "#08519c",
    colour = "#08519c",
    alpha = 0.5,
    linewidth = 0.8
  ) +
  geom_density(
    aes(x = ST69),
    fill = "lightsteelblue1",
    colour = "lightsteelblue1",
    alpha = 0.5,
    linewidth = 0.8
  ) +
  labs(
    x = "Posterior frequency",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) 
post_density_plot
ggsave("post_density_plot.png", post_density_plot, units = "in", width = 4, height = 3, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.1b Klebsiella MLSTs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare data - Kleb
kleb_bsi_mlst_df <- kleb_bsi_samples_metadata |>
  dplyr::rename(kleborate_mlst = klebsiella_mlst_ST) |>
  dplyr::group_by(kleborate_mlst) |>
  dplyr::summarise(count = n(), .groups = "drop") |>
  dplyr::select(kleborate_mlst, count) |>
  dplyr::arrange(count, kleborate_mlst)
#View(kleb_bsi_mlst_df)
# check MLST count distribution
table(kleb_bsi_mlst_df$count)
nrow(kleb_bsi_mlst_df) # 295
write.csv(kleb_bsi_mlst_df, "NEKSUS_kleb_bsi_mlst_df.csv", row.names = FALSE, quote = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 1.1bi Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
alpha_named_vals <- c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
alpha_novel_sum_vals <- c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
alpha_novel_num_vals <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)


res_kleb <- run_mlst_sensitivity_grid(
  kleb_bsi_mlst_df,
  feature_col = "kleborate_mlst",
  alpha_named_vals = alpha_named_vals,
  alpha_novel_sum_vals = alpha_novel_sum_vals,
  alpha_novel_num_vals = alpha_novel_num_vals,
  B = 1000,
  seed = 2026,
  out_dir = "sens_results/kleb_mlst",
  prefix = "kleb",
  dataset_label = "Klebsiella",
  save_csv = TRUE,
  keep_draws = FALSE
)
#View(res_kleb$fit_metrics)
res_kleb <- c()
res_kleb$fit_metrics <- read.csv("sens_results/kleb_mlst/kleb_fit_metrics.csv")
res_kleb$feature_summary <- read.csv("sens_results/kleb_mlst/kleb_feature_summary.csv")

# rmse heatmap
heatmap_kleb <- plot_fit_heatmap(res_kleb$fit_metrics, 
                                  alpha_named_vals = alpha_named_vals,
                                  alpha_novel_sum_vals = alpha_novel_sum_vals,
                                  alpha_novel_num_vals = alpha_novel_num_vals,
                                  fill_val = "rmse_log10" , fill_label = "log-RMSE")
heatmap_kleb
ggsave("sens_results/kleb_mlst_simpleBB_sensitivity_heatmap.png", heatmap_kleb, width = 10, height = 8, dpi = 300)

# sparse panel plot
alpha_named_min <- c(0.5, 1, 10, 100)
alpha_novel_sum_min <- c(0, 1, 10, 100, 250, 500)
alpha_novel_num_min <- c(0, 1, 10, 100, 250, 500)

thinned <- res_kleb$feature_summary |>
  filter(alpha_named %in% alpha_named_min 
         & alpha_novel_sum %in% alpha_novel_sum_min
         & alpha_novel_num %in% alpha_novel_num_min)
panel_kleb <- plot_obs_vs_pred(thinned, alpha_named_min, alpha_novel_sum_min, alpha_novel_num_min)
panel_kleb
ggsave("sens_results/kleb_mlst_simpleBB_sensitivity_panel_thinned.png", panel_kleb, width = 10, height = 8, dpi = 300)


#
mass_plot <- plot_curve_facet(
  dplyr::bind_rows(res_ecoli$mass_curves, res_kleb$mass_curves),
  alpha_named_vals = alpha_named_vals,
  alpha_novel_sum_vals = alpha_novel_sum_vals,
  alpha_novel_num_vals = alpha_novel_num_vals,
  value_col = "median_mass",
  lo_col = "q2.5",
  hi_col = "q97.5",
  y_lab = "Total mass above threshold",
  title = "Cumulative population mass curves",
  curve_type = "mass"
)
mass_plot
ggsave("sens_results/ecoli_kleb_mlst_simpleBB_sensitivity_mass_plot.png", mass_plot, width = 15, height = 10, dpi = 300)

richness_plot <- plot_curve_facet(
  dplyr::bind_rows(res_ecoli$richness_curves, res_kleb$richness_curves),
  alpha_named_vals = alpha_named_vals,
  alpha_novel_sum_vals = alpha_novel_sum_vals,
  alpha_novel_num_vals = alpha_novel_num_vals,
  value_col = "median_species_count",
  lo_col = "q2.5",
  hi_col = "q97.5",
  y_lab = "Number of species above threshold",
  title = "Cumulative species richness curves",
  curve_type = "richness"
)

richness_plot
ggsave("sens_results/ecoli_kleb_mlst_simpleBB_sensitivity_richness_plot.png", richness_plot, width = 15, height = 10, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 1.1bii Run single prior parameter set ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
df <- kleb_bsi_mlst_df
counts <- df$count
N <- sum(df$count)
K <- nrow(df)
f1 <- sum(df$count == 1)  # 223       # number of singletons
q_hat <- f1 / N  # 47% Good-Turing first-order - proportion of singletons

# fit data to estimate Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(counts)
crp_fit$summary_df # mean: theta = 111.2 -> 19% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.192

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$kleborate_mlst
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #71

# Preferred / exact approach (supports non-integer alphas):
kleb_bsi_mlst_bayes <- isolate_bayes_dirichlet(kleb_bsi_mlst_df,
                                          feature_col = "kleborate_mlst",
                                          alpha_named = alpha_named, 
                                          alpha_novel_sum =  alpha_novel_sum, 
                                          alpha_novel_num =  alpha_novel_sum, 
                                          B = 10000)

#print(kleb_bsi_mlst_bayes$summary_df)

# save
saveRDS(kleb_bsi_mlst_bayes, "model_results/kleb_bsi_kleborate_mlst_bayes.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#kleb_bsi_mlst_bayes <- readRDS("model_results/kleb_bsi_kleborate_mlst_bayes.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# define f-grid based on sample sizes:
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

# cumulative species richness curve (>=f; number and proportion of MLSTs at least as frequent as f)
prep_draws <- prep_bootstrap_draws(kleb_bsi_mlst_bayes$draws)
kleb_species_richness_df <- compute_species_richness_curve(prep_draws, f_grid = f_grid_99)
#View(kleb_species_richness_df)

# add sample sizes
kleb_species_richness_df <- kleb_species_richness_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Klebsiella",
                Exact = "Exact") |>
  filter(f != 0)
# save  species_richness df
write.csv(kleb_species_richness_df, "model_results/kleb_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#kleb_species_richness_df <- read.csv("model_results/kleb_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
kleb_mass_df <- compute_mass_curve(prep_draws, f_grid = f_grid_99)
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
write.csv(kleb_mass_df, "model_results/kleb_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#kleb_mass_df <- read.csv("model_results/kleb_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.1c Combined summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# combine dfs species_richness
ecoli_species_richness_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
kleb_species_richness_df <- read.csv("model_results/kleb_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
species_richness_df_mlst <- rbind(ecoli_species_richness_df, kleb_species_richness_df)
#View(species_richness_df_mlst)

# save merged species_richness df
write.csv(species_richness_df_mlst, "model_results/combined_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data
#species_richness_df_mlst <- read.csv("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# combine dfs mass
ecoli_mass_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
kleb_mass_df <- read.csv("model_results/kleb_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
mass_df_mlst <- rbind(ecoli_mass_df, kleb_mass_df)
#View(mass_df_mlst)

# save merged mass df
write.csv(mass_df_mlst, "model_results/combined_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data
#mass_df_mlst <- read.csv("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
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

# summary tables - population species_richness
# transform to long
species_richness_df_mlst_long <- species_richness_df_mlst |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))
#View(species_richness_df_mlst)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_mlst_species_richness_summary_table <- species_richness_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_species_proportion", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5_species_proportion", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5_species_proportion", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(combined_bsi_mlst_species_richness_summary_table)

# tidy table - extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_mlst_species_richness_summary_table <- combined_bsi_mlst_species_richness_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_mlst_species_richness_summary_table)
#View(combined_bsi_mlst_species_richness_summary_table)
# save
write.csv(combined_bsi_mlst_species_richness_summary_table, "model_results/combined_bsi_mlst_species_richness_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~#
# summary tables 
# transform to long
mass_df_mlst_long <- mass_df_mlst |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))
#View(mass_df_mlst_long)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_mlst_sample_coverage_summary_table <- mass_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_mlst_sample_coverage_summary_table <- combined_bsi_mlst_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_mlst_sample_coverage_summary_table)
#View(combined_bsi_mlst_sample_coverage_summary_table)
# save
write.csv(combined_bsi_mlst_sample_coverage_summary_table, "model_results/combined_bsi_mlst_sample_coverage_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * 1.1d Combined Plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# histogram of MLST frequencies 
ecoli_bsi_mlst_df <- ecoli_bsi_mlst_df |>
  mutate(frequency = count/sum(count))

posterior_mlst_freq_hist <- ggplot(data = ecoli_bsi_mlst_bayes$summary_df, aes(x = median) ) +
  geom_histogram(binwidth = 0.002, alpha = 0.6) +
  geom_histogram(aes(x = `q2.5`), binwidth = 0.002, alpha = 0.3) +
  geom_histogram(aes(x = `q97.5`), binwidth = 0.002, alpha = 0.1)  +
  scale_y_log10() +
 # scale_x_log10() +
  # plot actual data
  # plot lower and upper CIs
  geom_histogram(data = ecoli_bsi_mlst_df, aes(x = frequency), fill = "seagreen3", alpha = 0.5, binwidth = 0.002) +
  theme_minimal()
posterior_mlst_freq_hist
ggsave("model_results/ecoli_posterior_mlst_freq_hist.png", posterior_mlst_freq_hist, width = 6, height = 4, units = "in", dpi = 300)
#~~~~~~~~~~~~~~~~~#
# histogram of MLST frequencies
kleb_bsi_mlst_df <- kleb_bsi_mlst_df |>
  mutate(frequency = count/sum(count))

posterior_mlst_freq_hist <- ggplot(data = kleb_bsi_mlst_bayes$summary_df, aes(x = median) ) +
  geom_histogram(binwidth = 0.002, alpha = 0.6) +
  geom_histogram(aes(x = `q2.5`), binwidth = 0.002, alpha = 0.3) +
  geom_histogram(aes(x = `q97.5`), binwidth = 0.002, alpha = 0.1)  +
  # plot actual data
  # plot lower and upper CIs
  geom_histogram(data = kleb_bsi_mlst_df, aes(x = frequency), fill = "darkorange", alpha = 0.5, binwidth = 0.002) +
  theme_minimal()
posterior_mlst_freq_hist
ggsave("model_results/kleb_posterior_mlst_freq_hist.png", posterior_mlst_freq_hist, width = 6, height = 4, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot obs vs pred
plot_df <- ecoli_bsi_mlst_bayes$summary_df |>
  left_join(df, by = c("feature" = "kleborate_mlst")) |>
  mutate(count = case_when(is.na(count) ~ 0,
                           TRUE ~ count)) |>
  mutate(actual = count / sum(count),
         estimate = median, 
         lo = q2.5,
         hi = q97.5) |>
  filter(count !=0)


plot <- ggplot(data = plot_df, aes(x = actual, y = estimate)) +
  geom_abline(intercept = 0, slope = 1,
    linetype = "dashed", linewidth = 0.4, color = "grey50"
  ) +
  geom_errorbar(
    ggplot2::aes(ymin = lo, ymax = hi),
    width = 0, alpha = 0.45, color = "steelblue"
  ) +
  ggplot2::geom_point(size = 1.6, alpha = 0.8, color = "steelblue") +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  #ggplot2::facet_grid(alpha_novel ~ alpha_named, labeller = ggplot2::label_both) +
  ggplot2::labs(
    x = "Observed prevalence",
    y = "Posterior median prevalence"
  ) +
  ggplot2::theme_minimal(base_size = 11)
plot




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define colours
genus_colours <- c("Escherichia" = "seagreen3", "Klebsiella" = "darkorange")

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot (cumulative species_richness of MLSTs at frequency >= f)
cumulative_species_richness_plot <- ggplot(species_richness_df_mlst, aes(x = f, y = median_species_proportion, colour =  Genus, fill = Genus)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Genus", values = genus_colours) +
  scale_colour_manual(name = "Genus", values = genus_colours) +
  scale_x_log10( breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(species_richness_df_mlst$f), 1)) +
  labs(x = "MLST frequency (f) (log scale)",
       y = "Proportion of MLSTs of frequency ≥ f"
  ) +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
cumulative_species_richness_plot
# save
ggsave("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_species_richness_plot.png", plot = cumulative_species_richness_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Cumulative MLST fraction with sample size plot 
bayes_species_richness_plot <- ggplot(species_richness_df_mlst, aes(y = median_species_proportion, colour = Genus, fill = Genus)) +
  geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  geom_line(aes(x = min_sample_99)) +
  geom_ribbon(aes(x = min_sample_90, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.35, colour = NA) +
  geom_ribbon(aes(x = min_sample_99, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  #facet_wrap(~ Genus, ncol = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,5000)) +
  labs(y = "Proportion of MLSTs detected",
       x = "Sample size",
       ) +
  theme_minimal()
bayes_species_richness_plot
ggsave("model_results/combined_bsi_kleborate_mlst_bayes_ss_vs_species_richness_plot_80.png", plot = bayes_species_richness_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Plot (cumulative mass of population at frequency >= f)
cumulative_mass_plot <- ggplot(mass_df_mlst, aes(x = f, y = median_mass, colour =  Genus, fill = Genus)) +
  geom_line(data = mass_df_mlst, aes(x = f, y = median_mass, colour =  Genus)) +
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
ggsave("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_mass_plot.png", plot = cumulative_mass_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Cumulative fraction with sample size plot 
bayes_sample_coverage_plot <- ggplot(mass_df_mlst, aes(y = median_mass, colour = Genus, fill = Genus)) +
  geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  geom_line(aes(x = min_sample_99)) +
  geom_ribbon(aes(x = min_sample_90, ymin = q2.5, ymax = q97.5), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5, ymax = q97.5), alpha = 0.35, colour = NA) +
  geom_ribbon(aes(x = min_sample_99, ymin = q2.5, ymax = q97.5), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  #facet_wrap(~ Genus, ncol = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,5000)) +
  labs(y = "coverage",
       x = "Sample size",
       # title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/combined_bsi_kleborate_mlst_bayes_sample_coverage_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 1.2 fastbaps_L3 clusters (level 3) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.2a E. coli fastbaps_L3 clusters ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prep data
ecoli_bsi_fastbaps_L3 <- ecoli_bsi_samples_metadata |>
  dplyr::group_by(Level.3) |>
  dplyr::summarise(count = n(), .groups = "drop") |>
  dplyr::select(Level.3, count) |>
  dplyr::arrange(count, Level.3) |>
  dplyr::mutate(Level.3 = paste0("cluster_", Level.3)) 
#View(ecoli_bsi_fastbaps_L3)

# check count distribution
table(ecoli_bsi_fastbaps_L3$count)
nrow(ecoli_bsi_fastbaps_L3) # 161

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 1.2ai Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
alpha_named_vals <- c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
alpha_novel_sum_vals <- c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
alpha_novel_num_vals <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)

res_ecoli_fastbaps_L3 <- run_mlst_sensitivity_grid(
  ecoli_bsi_fastbaps_L3,
  feature_col = "Level.3",
  alpha_named_vals = alpha_named_vals,
  alpha_novel_sum_vals = alpha_novel_sum_vals,
  alpha_novel_num_vals = alpha_novel_num_vals,
  B = 1000,
  seed = 2026,
  out_dir = "sens_results/ecoli_fastbaps_L3",
  prefix = "ecoli",
  dataset_label = "Escherichia",
  save_csv = TRUE,
  keep_draws = FALSE
)
#View(res_ecoli_fastbaps_L3$fit_metrics)
#View(res_ecoli_fastbaps_L3$feature_summary)
res_ecoli_fastbaps_L3 <- c()
res_ecoli_fastbaps_L3$fit_metrics <- read.csv("sens_results/ecoli_fastbaps_L3/ecoli_fit_metrics.csv")
res_ecoli_fastbaps_L3$feature_summary <- read.csv("sens_results/ecoli_fastbaps_L3/ecoli_feature_summary.csv")


# rmse heatmap
heatmap_ecoli <- plot_fit_heatmap(res_ecoli_fastbaps_L3$fit_metrics, 
                                  alpha_named_vals = alpha_named_vals,
                                  alpha_novel_sum_vals = alpha_novel_sum_vals,
                                  alpha_novel_num_vals = alpha_novel_num_vals,
                                  fill_val = "rmse_log10" , fill_label = "log-RMSE")
heatmap_ecoli
ggsave("sens_results/ecoli_fastbaps_L3_simpleBB_sensitivity_heatmap.png", heatmap_ecoli, width = 10, height = 8, dpi = 300)

# sparse panel plot
alpha_named_min <- c(0.5, 1, 10, 100)
alpha_novel_sum_min <- c(0, 1, 10, 100, 250, 500)
alpha_novel_num_min <- c(0, 1, 10, 100, 250, 500)

thinned <- res_ecoli_fastbaps_L3$feature_summary |>
  filter(alpha_named %in% alpha_named_min 
         & alpha_novel_sum %in% alpha_novel_sum_min
         & alpha_novel_num %in% alpha_novel_num_min)
panel_ecoli <- plot_obs_vs_pred(thinned, alpha_named_min, alpha_novel_sum_min, alpha_novel_num_min)
panel_ecoli
ggsave("sens_results/ecoli_fastbaps_L3_simpleBB_sensitivity_panel_thinned.png", panel_ecoli, width = 10, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 1.2aii Run single prior parameter set ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
df <- ecoli_bsi_fastbaps_L3
counts <- df$count
N <- sum(df$count)
K <- nrow(df)
f1 <- sum(df$count == 1)  # 41      # number of singletons
q_hat <- f1 / N  # 2.8% Good-Turing first-order - proportion of singletons

# fit data to estimate Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(counts)
crp_fit$summary_df # mean: theta = 33.5-> 2% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.022

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$Level.3
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #4

# Preferred / exact approach (supports non-integer alphas):
ecoli_bsi_fastbaps_L3_bayes <- isolate_bayes_dirichlet(ecoli_bsi_fastbaps_L3,
                                                  feature_col = "Level.3",
                                                  alpha_named = alpha_named, 
                                                  alpha_novel_sum =  alpha_novel_sum, 
                                                  alpha_novel_num =  alpha_novel_sum, 
                                                  B = 10000)
#print(ecoli_bsi_fastbaps_L3_bayes$summary_df)

# save
saveRDS(ecoli_bsi_fastbaps_L3_bayes, "model_results/ecoli_bsi_fastbaps_L3_bayes.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#ecoli_bsi_fastbaps_L3_bayes <- readRDS("model_results/ecoli_bsi_fastbaps_L3_bayes.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# define f-grid based on sample sizes:
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_90 <- 1-((1-0.90)^(1/n_grid))
f_grid_95 <- 1-((1-0.95)^(1/n_grid))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))
#f_grid <- unique(c(f_grid_90, f_grid_95, f_grid_99))
#f_grid <- c(seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))

# cumulative species richness curve (>=f; number and proportion of fastbaps_L3s at least as frequent as f)
prep_draws <- prep_bootstrap_draws(ecoli_bsi_fastbaps_L3_bayes$draws)
ecoli_fastbaps_L3_species_richness_df <- compute_species_richness_curve(prep_draws, f_grid = f_grid_99)
#View(ecoli_fastbaps_L3_species_richness_df)

# add sample sizes
ecoli_fastbaps_L3_species_richness_df <- ecoli_fastbaps_L3_species_richness_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)
# save  species_richness df
write.csv(ecoli_fastbaps_L3_species_richness_df, "model_results/ecoli_bsi_fastbaps_L3_bayes_cumulative_species_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#ecoli_fastbaps_L3_species_richness_df <- read.csv("model_results/ecoli_bsi_fastbaps_L3_bayes_cumulative_species_richness_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
ecoli_mass_df <- compute_mass_curve(prep_draws, f_grid = f_grid_99)

# add sample sizes
ecoli_mass_df <- ecoli_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(ecoli_mass_df, "model_results/ecoli_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#ecoli_mass_df <- read.csv("model_results/ecoli_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.2b Klebsiella fastbaps_L3 clusters ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare data - Kleb
kleb_bsi_fastbaps_L3_df <- kleb_bsi_samples_metadata |>
  dplyr::group_by(Level.3) |>
  dplyr::summarise(count = n(), .groups = "drop") |>
  dplyr::select(Level.3, count) |>
  dplyr::arrange(count, Level.3) |>
  dplyr::mutate(Level.3 = paste0("cluster_", Level.3))
#View(kleb_bsi_fastbaps_L3_df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 1.2bi Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
alpha_named_vals <- c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
alpha_novel_sum_vals <- c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
alpha_novel_num_vals <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)

res_kleb_fastbaps_L3 <- run_mlst_sensitivity_grid(
  kleb_bsi_fastbaps_L3_df,
  feature_col = "Level.3",
  alpha_named_vals = alpha_named_vals,
  alpha_novel_sum_vals = alpha_novel_sum_vals,
  alpha_novel_num_vals = alpha_novel_num_vals,
  B = 1000,
  seed = 2026,
  out_dir = "sens_results/kleb_fastbaps_L3",
  prefix = "kleb",
  dataset_label = "Klebsiella",
  save_csv = TRUE,
  keep_draws = FALSE
)
#View(res_kleb_fastbaps_L3$fit_metrics)
res_kleb_fastbaps_L3 <- c()
res_kleb_fastbaps_L3$fit_metrics <- read.csv("sens_results/kleb_fastbaps_L3/kleb_fit_metrics.csv")
res_kleb_fastbaps_L3$feature_summary <- read.csv("sens_results/kleb_fastbaps_L3/kleb_feature_summary.csv")



# rmse heatmap
heatmap_kleb <- plot_fit_heatmap(res_kleb_fastbaps_L3$fit_metrics, 
                                  alpha_named_vals = alpha_named_vals,
                                  alpha_novel_sum_vals = alpha_novel_sum_vals,
                                  alpha_novel_num_vals = alpha_novel_num_vals,
                                  fill_val = "rmse_log10" , fill_label = "log-RMSE")
heatmap_kleb
ggsave("sens_results/kleb_fastbaps_L3_simpleBB_sensitivity_heatmap.png", heatmap_kleb, width = 10, height = 8, dpi = 300)

# sparse panel plot
alpha_named_min <- c(0.5, 1, 10, 100)
alpha_novel_sum_min <- c(0, 1, 10, 100, 250, 500)
alpha_novel_num_min <- c(0, 1, 10, 100, 250, 500)

thinned <- res_kleb_fastbaps_L3$feature_summary |>
  filter(alpha_named %in% alpha_named_min 
         & alpha_novel_sum %in% alpha_novel_sum_min
         & alpha_novel_num %in% alpha_novel_num_min)
panel_kleb <- plot_obs_vs_pred(thinned, alpha_named_min, alpha_novel_sum_min, alpha_novel_num_min)
panel_kleb
ggsave("sens_results/kleb_fastbaps_L3_simpleBB_sensitivity_panel_thinned.png", panel_kleb, width = 10, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~#
# plot combined curves
mass_plot <- plot_curve_facet(
  dplyr::bind_rows(res_ecoli$mass_curves, res_kleb_fastbaps_L3$mass_curves),
  alpha_vals = alpha_vals,
  value_col = "median_mass",
  lo_col = "q2.5",
  hi_col = "q97.5",
  y_lab = "Total mass above threshold",
  title = "Cumulative population mass curves",
  curve_type = "mass"
)
mass_plot
ggsave("sens_results/ecoli_kleb_fastbaps_L3_simpleBB_sensitivity_mass_plot.png", mass_plot, width = 15, height = 10, dpi = 300)

richness_plot <- plot_curve_facet(
  dplyr::bind_rows(res_ecoli$richness_curves, res_kleb_fastbaps_L3$richness_curves),
  alpha_vals = alpha_vals,
  value_col = "median_species_count",
  lo_col = "q2.5",
  hi_col = "q97.5",
  y_lab = "Number of species above threshold",
  title = "Cumulative species richness curves",
  curve_type = "richness"
)

richness_plot
ggsave("sens_results/ecoli_kleb_fastbaps_L3_simpleBB_sensitivity_richness_plot.png", richness_plot, width = 15, height = 10, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 1.2bii Run single prior parameter set ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
df <- kleb_bsi_fastbaps_L3_df
counts <- df$count
N <- sum(df$count)
K <- nrow(df)
f1 <- sum(df$count == 1)  # 11      # number of singletons
q_hat <- f1 / N  # 2.3% Good-Turing first-order - proportion of singletons

# fit data to estimate Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(counts)
crp_fit$summary_df # mean: theta = 5.2 -> 1% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.0108

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$Level.3
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat)) # round up
alpha_novel_sum #1

# Preferred / exact approach (supports non-integer alphas):
kleb_bsi_fastbaps_L3_bayes <- isolate_bayes_dirichlet(kleb_bsi_fastbaps_L3_df,
                                                 feature_col = "Level.3",
                                                 alpha_named = alpha_named, 
                                                 alpha_novel_sum =  alpha_novel_sum, 
                                                 alpha_novel_num =  alpha_novel_sum, 
                                                 B = 10000)
#print(kleb_bsi_fastbaps_L3_bayes$summary_df)
# save
saveRDS(kleb_bsi_fastbaps_L3_bayes, "model_results/kleb_bsi_fastbaps_L3_bayes.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in 
#kleb_bsi_fastbaps_L3_bayes <- readRDS("model_results/kleb_bsi_fastbaps_L3_bayes.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# define f-grid based on sample sizes:
n_grid_course <- c(seq(0, 10000, by = 10), seq(10010, 50000, by = 20))
f_grid_99 <- 1-((1-0.99)^(1/n_grid_course))

# cumulative species richness curve (>=f; number and proportion of fastbaps_L3s at least as frequent as f)
prep_draws <- prep_bootstrap_draws(kleb_bsi_fastbaps_L3_bayes$draws)
kleb_fastbaps_L3_species_richness_df <- compute_species_richness_curve(prep_draws, f_grid = f_grid_99)
#View(kleb_fastbaps_L3_species_richness_df)


# add sample sizes
kleb_fastbaps_L3_species_richness_df <- kleb_fastbaps_L3_species_richness_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Klebsiella",
                Exact = "Exact") |>
  filter(f != 0)
# save  species_richness df
write.csv(kleb_fastbaps_L3_species_richness_df, "model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_coarse.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#kleb_fastbaps_L3_species_richness_df <- read.csv("model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_species_richness_df.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# define f-grid based on sample sizes:
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

# cumulative species richness curve (>=f; number and proportion of fastbaps_L3s at least as frequent as f)
kleb_fastbaps_L3_species_richness_df_fine <- compute_species_richness_curve(prep_draws, f_grid = f_grid_99)
#View(kleb_fastbaps_L3_species_richness_df)


# add sample sizes
kleb_fastbaps_L3_species_richness_df_fine <- kleb_fastbaps_L3_species_richness_df_fine |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Klebsiella",
                Exact = "Exact") |>
  filter(f != 0)
# save  species_richness df
write.csv(kleb_fastbaps_L3_species_richness_df_fine, "model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_fine.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#kleb_fastbaps_L3_species_richness_df_fine <- read.csv("model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_fine.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# redifeine finder grain f-grid
#n_grid_course <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5))
#f_grid_99 <- 1-((1-0.99)^(1/n_grid))

# compute mass >= f curve
kleb_mass_df <- compute_mass_curve(prep_draws, f_grid = f_grid_99)
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
write.csv(kleb_mass_df, "model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Reload saved data 
#kleb_mass_df <- read.csv("model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df_fine.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.2c Combined summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#
ecoli_species_richness_df <- read.csv("model_results/ecoli_bsi_fastbaps_L3_bayes_cumulative_species_richness_df.csv")
kleb_species_richness_df_fine <- read.csv("model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_fine.csv")
species_richness_df_fastbaps_L3_fine <- rbind(ecoli_species_richness_df, kleb_species_richness_df_fine)
#View(species_richness_df_fastbaps_L3)
# save merged species_richness df
write.csv(species_richness_df_fastbaps_L3_fine, "model_results/combined_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_fine.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#species_richness_df_fastbaps_L3_fine <- read.csv("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_fine.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# load coarse grid version
kleb_species_richness_df_coarse <- read.csv("model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_coarse.csv")
species_richness_df_fastbaps_L3_coarse <- rbind(ecoli_species_richness_df, kleb_species_richness_df_coarse)
#View(species_richness_df_fastbaps_L3_coarse)
# save merged species_richness df
write.csv(species_richness_df_fastbaps_L3_coarse, "model_results/combined_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_coarse.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#species_richness_df_fastbaps_L3_coarse <- read.csv("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_coarse.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
ecoli_mass_df <- read.csv("model_results/ecoli_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv")
kleb_mass_df <- read.csv("model_results/kleb_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv")
mass_df_fastbaps_L3 <- rbind(ecoli_mass_df, kleb_mass_df)
#View(mass_df_fastbaps_L3)
# save merged mass df
write.csv(mass_df_fastbaps_L3, "model_results/combined_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#mass_df_fastbaps_L3 <- read.csv("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# summary tables - population species_richness
# transform to long
species_richness_df_fastbaps_L3_long <- species_richness_df_fastbaps_L3_fine |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))
#View(species_richness_df_fastbaps_L3_long)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_fastbaps_L3_species_richness_summary_table <- species_richness_df_fastbaps_L3_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_species_proportion", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5_species_proportion", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5_species_proportion", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(combined_bsi_fastbaps_L3_species_richness_summary_table)

# tidy table - extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_fastbaps_L3_species_richness_summary_table <- combined_bsi_fastbaps_L3_species_richness_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_fastbaps_L3_species_richness_summary_table)
#View(combined_bsi_fastbaps_L3_species_richness_summary_table)
# save
write.csv(combined_bsi_fastbaps_L3_species_richness_summary_table, "model_results/combined_bsi_fastbaps_L3_species_richness_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# summary tables - population mass
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
      median_ss  <- min_sample_at_or_above(df, "median_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
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
write.csv(combined_bsi_fastbaps_L3_sample_coverage_summary_table, "model_results/combined_bsi_fastbaps_L3_sample_coverage_summary_table.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 1.2d Combined plots  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define colours
genus_colours <- c("Escherichia" = "seagreen3", "Klebsiella" = "darkorange")

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot (cumulative species_richness of MLSTs at frequency >= f)
cumulative_species_richness_plot <- ggplot(species_richness_df_fastbaps_L3_coarse, aes(x = f, y = median_species_proportion, colour =  Genus, fill = Genus)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Genus", values = genus_colours) +
  scale_colour_manual(name = "Genus", values = genus_colours) +
  scale_x_log10( breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(species_richness_df_fastbaps_L3_coarse$f), 1)) +
  labs(x = "fastBAPS cluster frequency (f) (log scale)",
       y = "Proportion of fastbaps_L3 clusters of frequency ≥ f"
  ) +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
cumulative_species_richness_plot
# save
ggsave("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_species_richness_plot.png", plot = cumulative_species_richness_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Cumulative fastBAS clusters fraction with sample size plot 
bayes_sample_coverage_plot <- ggplot(species_richness_df_fastbaps_L3_coarse, aes(y = median_species_proportion, colour = Genus, fill = Genus)) +
  geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  geom_line(aes(x = min_sample_99)) +
  geom_ribbon(aes(x = min_sample_90, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.35, colour = NA) +
  geom_ribbon(aes(x = min_sample_99, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  #facet_wrap(~ Genus, ncol = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,5000)) +
  labs(y = "Proportion of fastbaps_L3 clusters detected",
       x = "Sample size",
  ) +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/combined_bsi_kleborate_mlst_bayes_ss_vs_species_richness_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)
# jagged edge from klebsiella comes from rounding error of f-grib being too fine-scaled. 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot (cumulative mass of population at frequency >= f)
cumulative_mass_plot <- ggplot(mass_df_fastbaps_L3, aes(x = f, y = median_mass, colour =  Genus, fill = Genus)) +
  geom_line(aes(x = f, y = median_mass, colour =  Genus)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Genus", values = genus_colours) +
  scale_colour_manual(name = "Genus", values = genus_colours) +
  scale_x_log10( breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(mass_df_fastbaps_L3$f), 1)) +
  labs(x = "fastBAPS cluster frequency (f) (log scale)",
       y = "Proportion of population belonging to cluster of frequency ≥ f"
  ) +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
cumulative_mass_plot
# save
ggsave("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_mass_plot.png", plot = cumulative_mass_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Cumulative arg fraction with sample size plot
bayes_sample_coverage_plot <- ggplot(mass_df_fastbaps_L3, aes(y = median_mass, colour = Genus, fill = Genus)) +
  geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  geom_line(aes(x = min_sample_99)) +
  geom_ribbon(aes(x = min_sample_90, ymin = q2.5, ymax = q97.5), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5, ymax = q97.5), alpha = 0.35, colour = NA) +
  geom_ribbon(aes(x = min_sample_99, ymin = q2.5, ymax = q97.5), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  #facet_wrap(~ Genus, ncol = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,5000)) +
  labs(y = "coverage",
       x = "Sample size",
       # title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/combined_bsi_fastbaps_L3_bayes_sample_coverage_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)
#~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Framework 2: Overall subisolate-level features (ARGs and plasmids) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# sub-isolate-level (feature-level) beta-binomial bootstrap (accounting for unseen features)
# - observed features: posterior Beta(count + beta, N - count + beta)
# - novel features: Beta(beta_each, beta_each), with total novel mass split evenly, across beta_novel_num novel features
# Returns:
# - p_draws:     B x (K + novel_num) sampled prevalences
# - summary_df:   posterior summaries of prevalence
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
subisolate_bayes_beta <- function(iso_df,
                                    gene_cols,
                                    B = 2000,
                                    beta_named = NULL,
                                    beta_novel_sum = 0,
                                    beta_novel_num = 0,
                                    seed = 2026) {
  
  stopifnot(requireNamespace("matrixStats", quietly = TRUE))
  set.seed(seed)
  
  # basic checks
  stopifnot(is.data.frame(iso_df))
  stopifnot(length(gene_cols) >= 1)
  stopifnot(all(gene_cols %in% colnames(iso_df)))
  stopifnot(B >= 1)
  stopifnot(is.numeric(beta_novel_sum), length(beta_novel_sum) == 1, beta_novel_sum >= 0)
  stopifnot(is.numeric(beta_novel_num), length(beta_novel_num) == 1, beta_novel_num >= 0)
  
  beta_novel_num <- as.integer(beta_novel_num)
  N <- nrow(iso_df)
  
  # presence/absence matrix for observed features
  Z <- as.matrix(iso_df[, gene_cols, drop = FALSE])
  Z <- (Z > 0) * 1L
  
  # observed feature counts across isolates
  x_obs <- colSums(Z, na.rm = TRUE)
  K <- length(x_obs)
  obs_names <- colnames(Z)
  
  # default beta pseudocounts for observed features
  if (is.null(beta_named)) {
    beta_named <- rep(1, K)
    names(beta_named) <- obs_names
  } else {
    if (length(beta_named) != K) {
      stop("beta_named must have the same length as gene_cols.")
    }
    
    # if names are present and match, reorder to align with gene_cols
    if (!is.null(names(beta_named)) && all(obs_names %in% names(beta_named))) {
      beta_named <- beta_named[obs_names]
    } else {
      names(beta_named) <- obs_names
    }
  }
  
  beta_named <- as.numeric(beta_named)
  
  # decide whether to add novel features
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
  
  # posterior parameters for observed features

    # storage
  p_draws <- matrix(NA_real_, nrow = B, ncol = total_features)
  colnames(p_draws) <- feature_names
  
  for (b in seq_len(B)) {
    # observed features
    p_obs <- stats::rbeta(K, shape1 = x_obs + beta_named, shape2 = (N - x_obs) + beta_named)
    
    if (add_novel) {
      # novel features: symmetric Beta(beta_each, beta_each) prior
      # then Binomial(N, p) sampling
      p_novel <- stats::rbeta(beta_novel_num, shape1 = beta_novel_vec, shape2 = N + beta_novel_vec)
      
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
  
  out <- list(
    p_draws = p_draws,
    summary_df = summary_df,
    beta_novel_sum = beta_novel_sum,
    beta_novel_num = beta_novel_num,
    beta_novel_each = beta_novel_each
  )
  
  out
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Postprocess: extract posterior summaries for gene prevalences
extract_prevalence_df <- function(bbs, gene_cols) {
  p_mat <- sapply(bbs, function(x) x$p_draws)  # genes x B
  # transpose to B x genes
  p_df <- as.data.frame(t(p_mat))
  colnames(p_df) <- gene_cols
  # compute summaries
  tibble(
    gene = gene_cols,
    median = apply(p_df, 2, quantile, 0.5),
    sd = apply(p_df, 2, sd),
    q2.5 = apply(p_df, 2, quantile, 0.025),
    q97.5 = apply(p_df, 2, quantile, 0.975)
  )
}


# helper function to extract posterior prevalence matrix  from model output (rows = draws, cols = genes/features)
get_p_matrix <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(as.matrix(x))
  }
  if (is.list(x) && !is.null(x$p_draws)) {
    return(as.matrix(x$p_draws))
  }
  stop("Could not find p_draws in input.")
}


# function to make sumultive species richness table for given frequency threshold, f
default_subiso_f_grid <- function() {
  n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5))
  f_grid <- 1-((1-0.99)^(1/n_grid)) # 99% confidence of detection
  return(f_grid)
}

prep_subiso_curve <- function(df) {
  stopifnot(is.data.frame(df) || is.matrix(df))
  p_mat <- as.matrix(df)
  
  # sort each posterior draw once
  p_sorted <- t(apply(p_mat, 1, sort))
  
  list(
    x_sorted = p_sorted,
    cumsums  = matrixStats::rowCumsums(p_sorted),
    totals   = rowSums(p_mat, na.rm = TRUE),
    K        = ncol(p_mat)
  )
}

compute_species_richness_subiso <- function(x, f_grid = NULL) {
  stopifnot(requireNamespace("matrixStats", quietly = TRUE))
  
  p_mat <- get_p_matrix(x)
  B <- nrow(p_mat)
  K <- ncol(p_mat)
  
  if (is.null(f_grid)) f_grid <- default_subiso_f_grid()
  F <- length(f_grid)
  
  # sort each posterior draw once
  p_sorted <- t(apply(p_mat, 1, sort))
  
  rich_mat <- matrix(0L, nrow = B, ncol = F)
  
  for (i in seq_len(B)) {
    idx <- findInterval(
      f_grid,
      p_sorted[i, ],
      left.open = TRUE,
      rightmost.closed = TRUE
    )
    rich_mat[i, ] <- K - idx
  }
  
  prop_mat <- rich_mat / K
  
  q_count <- matrixStats::colQuantiles(rich_mat, probs = c(0.025, 0.975), na.rm = TRUE)
  q_prop  <- matrixStats::colQuantiles(prop_mat, probs = c(0.025, 0.975), na.rm = TRUE)
  
  out <- tibble::tibble(
    f = f_grid,
    median_species_count = matrixStats::colMedians(rich_mat, na.rm = TRUE),
    q2.5_species_count = q_count[, 1],
    q97.5_species_count = q_count[, 2],
    median_species_proportion = matrixStats::colMedians(prop_mat, na.rm = TRUE),
    q2.5_species_proportion = q_prop[, 1],
    q97.5_species_proportion = q_prop[, 2]
  )
  
  dplyr::filter(out, !is.na(median_species_count))
}


compute_mass_curve_subiso <- function(x, f_grid = NULL) {
  stopifnot(requireNamespace("matrixStats", quietly = TRUE))
  
  p_mat <- get_p_matrix(x)
  B <- nrow(p_mat)
  
  if (is.null(f_grid)) f_grid <- default_subiso_f_grid()
  F <- length(f_grid)
  
  p_sorted <- t(apply(p_mat, 1, sort))
  cs <- matrixStats::rowCumsums(p_sorted)
  totals <- rowSums(p_mat, na.rm = TRUE)
  
  mass_mat <- matrix(NA_real_, nrow = B, ncol = F)
  
  for (i in seq_len(B)) {
    idx <- findInterval(
      f_grid,
      p_sorted[i, ],
      left.open = TRUE,
      rightmost.closed = TRUE
    )
    
    cs_i <- c(0, cs[i, ])
    
    if (!is.na(totals[i]) && totals[i] > 0) {
      mass_mat[i, ] <- (totals[i] - cs_i[idx + 1]) / totals[i]
    }
  }
  
  q_mass <- matrixStats::colQuantiles(mass_mat, probs = c(0.025, 0.975), na.rm = TRUE)
  
  out <- tibble::tibble(
    f = f_grid,
    median_mass = matrixStats::colMedians(mass_mat, na.rm = TRUE),
    q2.5 = q_mass[, 1],
    q97.5 = q_mass[, 2]
  )
  
  dplyr::filter(out, !is.na(median_mass))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# sensitivity analysis wrappers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit metric: RMSE on log10 prevalence scale
compute_log_rmse_subiso <- function(actual, estimate, eps = NULL, weights = NULL) {
  stopifnot(length(actual) == length(estimate))
  if (is.null(eps)) eps <- 0.5 / max(1, length(actual))
  err <- log10(actual + eps) - log10(estimate + eps)
  sq <- err^2
  if (is.null(weights)) {
    return(sqrt(mean(sq, na.rm = TRUE)))
  }
  weights <- weights / sum(weights, na.rm = TRUE)
  sqrt(sum(weights * sq, na.rm = TRUE))
}

compute_log_mae_subiso <- function(actual, estimate, eps = NULL, weights = NULL) {
  stopifnot(length(actual) == length(estimate))
  if (is.null(eps)) eps <- 0.5 / max(1, length(actual))
  err <- abs(log10(actual + eps) - log10(estimate + eps))
  if (is.null(weights)) {
    return(mean(err, na.rm = TRUE))
  }
  weights <- weights / sum(weights, na.rm = TRUE)
  sum(weights * err, na.rm = TRUE)
}

compute_coverage95_subiso <- function(actual, lo, hi) {
  mean(actual >= lo & actual <= hi, na.rm = TRUE)
}

summarise_subiso_posterior <- function(bbs, gene_cols, actual_vec, beta_named, beta_novel_sum, beta_novel_num, dataset_label = NA_character_) {
  p_mat <- get_p_matrix(bbs)
  # ensure column order
  if (!is.null(colnames(p_mat)) && all(gene_cols %in% colnames(p_mat))) {
    p_mat <- p_mat[, gene_cols, drop = FALSE]
  } else {
    colnames(p_mat) <- gene_cols
  }
  est_med <- apply(p_mat, 2, stats::median, na.rm = TRUE)
  est_lo  <- apply(p_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE, names = FALSE)
  est_hi  <- apply(p_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE, names = FALSE)
  
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


# Optional helper: draw a darker shade sequence
make_green_shades <- function(n) {
  grDevices::colorRampPalette(c("#d9f0d3", "#7fc97f", "#238b45"))(n)
}

make_orange_shades <- function(n) {
  grDevices::colorRampPalette(c("#fee6ce", "#fdae6b", "#e6550d"))(n)
}

# Main wrapper
run_subiso_sensitivity_grid <- function(iso_df,
                                        gene_cols = gene_cols,
                                        beta_named = c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000),
                                        beta_novel_sum = c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000),
                                        beta_novel_num = c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000),
                                        B = 1000,
                                        seed = 2026,
                                        out_dir = ".",
                                        prefix = "subiso",
                                        dataset_label = NA_character_,
                                        f_grid = NULL,
                                        save_csv = TRUE,
                                        keep_draws = FALSE,
                                        compute_curves = TRUE) {
  stopifnot(is.data.frame(iso_df))
  stopifnot(all(gene_cols %in% names(iso_df)))
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  beta_named <- as.numeric(beta_named)
  beta_novel_sum <- as.numeric(beta_novel_sum)
  beta_novel_num <- as.numeric(beta_novel_num)
  
  beta_grid <- tidyr::crossing(beta_named = beta_named, beta_novel_sum = beta_novel_sum, beta_novel_num = beta_novel_num) |>
    dplyr::mutate(grid_id = dplyr::row_number())
    
  actual_vec <- colMeans(iso_df[, gene_cols, drop = FALSE], na.rm = TRUE)
  actual_vec <- stats::setNames(actual_vec, gene_cols)
  
  
  n_grid <- nrow(beta_grid)
  feature_list <- vector("list", n_grid)
  metric_list <- vector("list", n_grid)
  
  mass_list <- if (compute_curves) vector("list", n_grid) else NULL
  rich_list <- if (compute_curves) vector("list", n_grid) else NULL
  draws_list <- if (keep_draws) vector("list", n_grid) else NULL
  
  K <- length(gene_cols)
  
  for (i in seq_len(n_grid)) {
    b_named <- beta_grid$beta_named[[i]]
    b_novel_sum <- beta_grid$beta_novel_sum[[i]]
    b_novel_num <- beta_grid$beta_novel_num[[i]]
    
    b_named_vec <- rep.int(b_named, K)
    
    bbs <- subisolate_bayes_beta(
      iso_df = iso_df,
      gene_cols = gene_cols,
      beta_named = b_named_vec,
      beta_novel_sum = b_novel_sum,
      beta_novel_num = b_novel_num,
      B = B,
      seed = seed + i
    )
    
    p_mat <- bbs$p_draws
    
    feat <- summarise_subiso_posterior(
      bbs = bbs,
      gene_cols = gene_cols,
      actual_vec = actual_vec,
      beta_named = b_named,
      beta_novel_sum = b_novel_sum,
      beta_novel_num = b_novel_num,
      dataset_label = dataset_label
    )
    
    eps <- 0.5 / nrow(iso_df)
    
    metric <- tibble::tibble(
      beta_named = b_named,
      beta_novel_sum = b_novel_sum,
      beta_novel_num = b_novel_num,
      dataset = dataset_label,
      rmse_log10 = compute_log_rmse_subiso(feat$actual, feat$estimate, eps = eps),
      mae_log10  = compute_log_mae_subiso(feat$actual, feat$estimate, eps = eps),
      coverage95  = compute_coverage95_subiso(feat$actual, feat$lo, feat$hi)
    )
    
    if (compute_curves) {
      mass_list[[i]] <- compute_mass_curve_subiso(p_mat, f_grid = f_grid) |>
        dplyr::mutate(
          beta_named = b_named,
          beta_novel_sum = b_novel_sum,
          beta_novel_num = b_novel_num,
          dataset = dataset_label
        )
      
      rich_list[[i]] <- compute_species_richness_subiso(p_mat, f_grid = f_grid) |>
        dplyr::mutate(
          beta_named = b_named,
          beta_novel_sum = b_novel_sum,
          beta_novel_num = b_novel_num,
          dataset = dataset_label
        )
    }
    
    feature_list[[i]] <- feat
    metric_list[[i]] <- metric
    if (keep_draws) draws_list[[i]] <- bbs
    
  }
  feature_summary <- dplyr::bind_rows(feature_list)
  fit_metrics <- dplyr::bind_rows(metric_list)
  
  mass_curves <- if (compute_curves) dplyr::bind_rows(mass_list) else NULL
  richness_curves <- if (compute_curves) dplyr::bind_rows(rich_list) else NULL
  
  if (save_csv) {
    readr::write_csv(feature_summary, file.path(out_dir, paste0(prefix, "_feature_summary.csv")))
    readr::write_csv(fit_metrics, file.path(out_dir, paste0(prefix, "_fit_metrics.csv")))
    
    if (compute_curves) {
      readr::write_csv(mass_curves, file.path(out_dir, paste0(prefix, "_mass_curves.csv")))
      readr::write_csv(richness_curves, file.path(out_dir, paste0(prefix, "_richness_curves.csv")))
    }
  }
  
  list(
    feature_summary = feature_summary,
    fit_metrics = fit_metrics,
    mass_curves = mass_curves,
    richness_curves = richness_curves,
    draws = draws_list,
    beta_named_grid = beta_grid$beta_named,
    beta_novel_sum_grid = beta_grid$beta_novel_sum,
    beta_novel_num_grid = beta_grid$beta_novel_num
  )
}
  

# Heatmap: lower RMSE = better fit
plot_fit_heatmap_subiso <- function(fit_metrics, beta_named_vals, beta_novel_sum_vals, beta_novel_num_vals, fill_val = "rmse_log10", fill_label = "log-RMSE") {
  fill_sym <- rlang::sym(fill_val)
  
  # make facet labels
  facet_labs <- setNames(
    paste0("beta[1:k] * ': ' * ", format(beta_named_vals, trim = TRUE, scientific = FALSE)),
    as.character(beta_named_vals)
  )
  
  ggplot2::ggplot(
    fit_metrics |>
      dplyr::mutate(beta_named = factor(beta_named, levels = beta_named_vals),
                    beta_novel_sum = factor(beta_novel_sum, levels = beta_novel_sum_vals),
                    beta_novel_num = factor(beta_novel_num, levels = beta_novel_num_vals),
                  ),
    ggplot2::aes(x = beta_novel_sum, y = beta_novel_num, fill = !!fill_sym)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_y_discrete(breaks = NULL, labels = NULL) +
    ggplot2::scale_fill_viridis_c(direction = -1, option = "C") +
    ggplot2::facet_wrap(~ beta_named,
                        labeller = ggplot2::as_labeller(facet_labs, label_parsed)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = expression(beta[novel]~Mass),
      y = expression(beta[novel]~Number),
      fill = fill_label
    ) +
    ggplot2::theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
}


# Observed vs posterior median panel plot
plot_obs_vs_pred_subiso <- function(feature_summary, beta_named_vals, beta_novel_sum_vals, beta_novel_num_vals) {
  
  # discrete steel blue palette
  n_cols <- length(beta_novel_num_vals)
  steelblue_pal <- grDevices::colorRampPalette(c("steelblue4", "steelblue2", "lightsteelblue1"))(n_cols)
  
  # add labels
  named_labeller <- ggplot2::as_labeller(
    setNames(paste0("beta[1:k]: ", beta_named_vals),
             beta_named_vals),label_parsed)
  
  novel_sum_labeller <- ggplot2::as_labeller(
    setNames(paste0("beta[novel]~Mass: ", beta_novel_sum_vals),
             beta_novel_sum_vals),label_parsed)
  
  feature_summary |>
    dplyr::mutate(
      beta_named = factor(beta_named, levels = beta_named_vals),
      beta_novel_sum = factor(beta_novel_sum, levels = rev(beta_novel_sum_vals)),
      beta_novel_num = factor(beta_novel_num, levels = rev(beta_novel_num_vals))
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = actual, y = estimate, colour = beta_novel_num)) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", linewidth = 0.4, color = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo, ymax = hi),
      width = 0, beta = 0.7, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.7, beta = 0.9, alpha = 0.5) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(
      values = setNames(steelblue_pal, beta_novel_num_vals),
      name = expression(beta[novel]~Number)
    ) +
    ggplot2::facet_grid(
      beta_novel_sum ~ beta_named,
      labeller = ggplot2::labeller(
        beta_named = named_labeller,
        beta_novel_sum = novel_sum_labeller
      )
    ) +
    ggplot2::labs(
      x = "Observed prevalence",
      y = "Posterior median prevalence"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

# Curve plot: colour lines by beta_novel_num, facet by beta_novel_sum, and beta_named
plot_curve_facet_subiso <- function(curve_df,
                             beta_named_vals, 
                             beta_novel_sum_vals, 
                             beta_novel_num_vals, 
                             value_col, 
                             lo_col, hi_col, 
                             y_lab, 
                             title, 
                             curve_type = c("mass", "richness")) {
  
  curve_type <- match.arg(curve_type) 
  datasets <- unique(curve_df$dataset) 
  datasets <- datasets[!is.na(datasets)] 
  beta_named_levels <- sort(unique(curve_df$beta_named)) 
  beta_named_labels <- as.character(beta_named_levels) 
  beta_novel_sum_levels <- sort(unique(curve_df$beta_novel_sum)) 
  beta_novel_sum_labels <- as.character(beta_novel_sum_levels) 
  beta_novel_num_levels <- sort(unique(curve_df$beta_novel_num)) 
  beta_novel_num_labels <- as.character(beta_novel_num_levels) 
  
  pal <- c() 
  
  if ("Escherichia" %in% datasets) { 
    ecoli_cols <- make_green_shades(length(beta_novel_levels)) 
    names(ecoli_cols) <- paste("Escherichia", beta_novel_num_labels, sep = "__") 
    pal <- c(pal, ecoli_cols) 
  } 
  if ("Klebsiella" %in% datasets) { 
    kleb_cols <- make_orange_shades(length(beta_novel_levels)) 
    names(kleb_cols) <- paste("Klebsiella", beta_novel_num_labels, sep = "__") 
    pal <- c(pal, kleb_cols) 
  } 
  
  curve_df <- curve_df |> 
    dplyr::mutate( dataset = factor(dataset, 
                                    levels = c("Escherichia", "Klebsiella")), 
                   beta_named = factor(beta_named, levels = beta_named_vals), 
                   beta_novel_sum = factor(beta_novel_sum, levels = beta_novel_sum_vals), 
                   beta_novel_num = factor(beta_novel_num, levels = beta_novel_num_vals), 
                   line_id = paste(dataset, beta_novel_num, sep = "__") ) 
  ecoli_levels <- paste("Escherichia", beta_named_vals , sep = "__") 
  kleb_levels <- paste("Klebsiella", beta_named_vals , sep = "__") 
  line_levels <- c(ecoli_levels, kleb_levels) 
  
  curve_df <- curve_df |> 
    dplyr::mutate(line_id = factor(line_id, levels = line_levels)) 
  legend_labels <- c(paste("E. coli α =", beta_named_vals), paste("Klebsiella α =", beta_named_vals)) 
  
  ggplot2::ggplot(curve_df, ggplot2::aes(x = f, y = .data[[value_col]], group = line_id, color = line_id, fill = line_id ) ) + 
    ggplot2::geom_ribbon( ggplot2::aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]), beta = 0.12, colour = NA ) + 
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::facet_wrap( ggplot2::vars(dataset, beta_novel), ncol = 7, labeller = ggplot2::label_both, scales = "free_y" ) + 
    ggplot2::scale_color_manual(values = pal, breaks = line_levels, labels = legend_labels, name = expression(beta[novel]~Number)) + 
    ggplot2::scale_fill_manual(values = pal, guide = "none") + 
    ggplot2::scale_x_log10() + 
    ggplot2::facet_grid(
      beta_novel_sum ~ beta_named,
      labeller = ggplot2::labeller(
        beta_named = named_labeller,
        beta_novel_sum = novel_sum_labeller
      )
    ) +
    ggplot2::labs(x = "Threshold f", y = y_lab, title = title, color = "Dataset / Number of novel features" ) + 
    ggplot2::theme_minimal(base_size = 12) + 
    ggplot2::theme( strip.text = ggplot2::element_text(face = "bold"), legend.position = "right" ) 
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 2.1 ARGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.1a E. coli ARGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare E.coli data into format where 1 isolate per row, and genes are columns
ecoli_bsi_arg_df <- ecoli_bsi_amrfinder_metadata |>
  filter(Type =="AMR" & !is.na(Type)) |>
  group_by(sample, Element.symbol) |> 
  summarise(count = n(),
            presence = case_when(count >0 ~ 1,
                                 count <=0 ~ 0,
                                 TRUE ~ 0))

ecoli_bsi_arg_presence <- ecoli_bsi_arg_df|>
  pivot_wider(id_cols = sample, names_from = Element.symbol, values_from = presence, values_fill =0) |>
  ungroup()
#View(ecoli_bsi_arg_presence)
length(unique(ecoli_bsi_arg_presence$sample))
write.csv(ecoli_bsi_arg_presence, "NEKSUS_ecoli_bsi_arg_presence_matrix.csv", row.names = FALSE, quote = FALSE)


# Set params and add 'r' as row counts
iso_df <- ecoli_bsi_arg_presence
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 2.1ai Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
beta_named <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
beta_novel_sum <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
beta_novel_num <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)


res_ecoli_arg <- run_subiso_sensitivity_grid(
  iso_df,
  gene_cols = gene_cols,
  beta_named = beta_named, 
  beta_novel_sum = beta_novel_sum, 
  beta_novel_num = beta_novel_num, 
  B = 1000,
  seed = 2026,
  out_dir = "sens_results/ecoli_arg",
  prefix = "ecoli",
  dataset_label = "Escherichia",
  save_csv = TRUE,
  keep_draws = FALSE,
  compute_curves = FALSE
)
#View(res_ecoli_arg$fit_metrics)
#View(res_ecoli_arg$feature_summary)
res_ecoli_arg <- c()
res_ecoli_arg$fit_metrics <- read.csv("sens_results/ecoli_arg/ecoli_fit_metrics.csv")
res_ecoli_arg$feature_summary <- read.csv("sens_results/ecoli_arg/ecoli_feature_summary.csv")

# rmse heatmap 
heatmap_df <- res_ecoli_arg$fit_metrics |>
  filter(beta_named %in% c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000))
heatmap_ecoli <- plot_fit_heatmap_subiso(heatmap_df, 
                                         beta_named_vals = beta_named, 
                                         beta_novel_sum_vals  = beta_novel_sum, 
                                         beta_novel_num_vals = beta_novel_num,
                                         fill_val = "rmse_log10" , 
                                         fill_label = "log-RMSE")
heatmap_ecoli
ggsave("sens_results/ecoli_arg_simpleBB_sensitivity_heatmap.png", heatmap_ecoli, width = 10, height = 8, dpi = 300)

# full panel plot
panel_df <- res_ecoli_arg$feature_summary |>
  filter(beta_named %in% c(0.5, 1, 10, 100, 1000),
         beta_novel_sum %in% c(0, 1, 10, 100, 1000),
         beta_novel_num %in% c(0, 1, 10, 100, 1000))

panel_ecoli <- plot_obs_vs_pred_subiso(panel_df, 
                                       beta_named_vals = c(0.5, 1, 10, 100, 1000),
                                       beta_novel_sum_vals = c(0, 1, 10, 100, 1000),
                                       beta_novel_num_vals = c(0, 1, 10, 100, 1000))
panel_ecoli
ggsave("sens_results/ecoli_arg_simpleBB_sensitivity_panel.png", panel_ecoli, width = 10, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 2.1aii Run single prior parameter set ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
df <- ecoli_bsi_arg_df |>
  group_by(Element.symbol) |>
  summarise(count = n()) |>
  arrange(count)
N_isolate <- length(unique(ecoli_bsi_arg_df$sample)) # number of genes
N_gene <- sum(df$count) # number of genes
K <- length(unique(ecoli_bsi_arg_df$Element.symbol)) # 189 # num unique features
f1 <- sum(df$count == 1) # 61 - no. singletons
u_hat <- f1 / N_isolate  # 0.005 - frequency of singletons
u_hat

# set priors
beta_prior <- 1 
beta_named <- rep(beta_prior, K)  # set uninformative priors
names(beta_named) <- df$Element.symbol
# set novel priors to preserve 'novel' mass, without and with including singletons 
beta_novel_num = ceiling((u_hat * K)/ (1 - u_hat)) # 188 - derived to make the prior novel frequency the same as in the data
beta_novel_sum = beta_novel_num * beta_prior 
beta_novel_num
beta_novel_sum

# run
ecoli_bsi_arg_bbs <- subisolate_bayes_beta(
  iso_df = iso_df,
  gene_cols = gene_cols,
  beta_named = beta_named,
  beta_novel_sum = beta_novel_sum,
  beta_novel_num = beta_novel_sum, # many uncommon novel features, as opposed to 1 common novel block
  B = 10000)

# save
saveRDS(ecoli_bsi_arg_bbs, "model_results/ecoli_bsi_ARG_bayes.rds")
#~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_arg_bbs <- read_rds("model_results/ecoli_bsi_ARG_bayes.rds")
#~~~~~~~~~~~~~~~#
# Prevalence historgam
ecoli_bsi_prevalence_summary <- ecoli_bsi_arg_bbs$summary_df
#View(ecoli_bsi_arg_bbs$summary_df)
# plot histogram of detection probabilities (alpha varies by quantile)
arg_post_hist <- ggplot(data = ecoli_bsi_prevalence_summary) +
  geom_histogram(aes(x = median_prevalence), binwidth = 0.01, alpha= 0.5, fill = "seagreen3") +
  geom_histogram(aes(x = q2.5_prevalence), binwidth = 0.01, alpha= 0.7, fill = "seagreen3") +
  geom_histogram(aes(x = q97.5_prevalence), binwidth = 0.01, alpha= 0.2, fill = "seagreen3") +
  labs(title = "Histogram of posterior estimated ARG prevalence", 
       x = "ARG prevalence\n(median and 95% CIs)" ) +
  theme_minimal()
arg_post_hist
# save
ggsave("model_results/ecoli_bsi_arg_posterior_distribution_histogram.png", arg_post_hist, width = 6, height = 4, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define frequency matrix
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

#extract posterior prevalence matrix 
p_mat <- get_p_matrix(ecoli_bsi_arg_bbs)

# cumulative species richness curves
ecoli_bsi_arg_richness_df <- compute_species_richness_subiso(p_mat, f_grid = f_grid_99)
# cumulative population mass curves
ecoli_bsi_arg_mass_df  <- compute_mass_curve_subiso(p_mat, f_grid = f_grid_99)

# add sample sizes 
ecoli_bsi_arg_richness_df <- ecoli_bsi_arg_richness_df |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_arg_richness_df)
# save
write.csv(ecoli_bsi_arg_richness_df, "model_results/ecoli_bsi_ARG_bayes_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_arg_richness_df <- read.csv("model_results/ecoli_bsi_ARG_bayes_richness_df.csv")
#~~~~~~~~~~~~~~~#

# min sample sizes using f as feature frequency
ecoli_bsi_arg_mass_df <- ecoli_bsi_arg_mass_df |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_arg_mass_df)
# save
write.csv(ecoli_bsi_arg_mass_df, "model_results/ecoli_bsi_ARG_bayes_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_arg_mass_df <- read.csv("model_results/ecoli_bsi_ARG_bayes_mass_df.csv")
#~~~~~~~~~~~~~~~#

# quick code to plot illustrative results for figure
set.seed(2026)

p_mat <- ecoli_bsi_arg_bbs$p_draws

# sample 4 random feature columns
sel_cols <- sample(colnames(p_mat), 4)

# long format for plotting
plot_df <- as.data.frame(p_mat[, sel_cols, drop = FALSE]) |>
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Feature",
    values_to = "Posterior frequency"
  )
plot_df$Feature <- factor(plot_df$Feature, levels = c("tet(B)","NOVEL_36","ompF_Q88STOP", "NOVEL_32"))

# blue shades, one per selected feature
blue_cols <- c("#c6dbef", "#9ecae1", "#4292c6", "#08519c")
names(blue_cols) <- sel_cols

post_beta_binomial_dist <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = `Posterior frequency`, fill = Feature, colour = Feature)
) +
  ggplot2::geom_density(alpha = 0.45, linewidth = 0.8) +
  ggplot2::facet_wrap(~ Feature, nrow = 1, scales = "free") +
  ggplot2::scale_fill_manual(values = blue_cols) +
  ggplot2::scale_colour_manual(values = blue_cols) +
  ggplot2::labs(
    x = "Posterior frequency",
    y = "Density"
  ) +
  ggplot2::theme_minimal(base_size = 8) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(face = "bold")
  )
post_beta_binomial_dist
ggsave("post_beta_binomial_dist.png", post_beta_binomial_dist, units = "in", width = 8, height = 2, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.1b Klebsiella ARGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare Klebsiella data into format where 1 isolate per row, and genes are columns
kleb_bsi_arg_df <- kleb_bsi_amrfinder_metadata |>
  filter(Type =="AMR" & !is.na(Type)) |>
  group_by(sample, Element.symbol) |> 
  summarise(count = n(),
            presence = case_when(count >0 ~ 1,
                                 count <=0 ~ 0,
                                 TRUE ~ 0)) 
kleb_bsi_arg_presence <- kleb_bsi_arg_df|>
  pivot_wider(id_cols = sample, names_from = Element.symbol, values_from = presence, values_fill =0) |>
  ungroup()
#View(kleb_bsi_arg_presence)
write.csv(kleb_bsi_arg_presence, "NEKSUS_kleb_bsi_arg_presence_matrix.csv", row.names = FALSE, quote = FALSE)


# Set params
iso_df <- kleb_bsi_arg_presence
iso_df <- iso_df |> dplyr::ungroup()
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))
#View(iso_df)
#length(unique(iso_df$sample)) # 468

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 2.1bi Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
beta_named <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
beta_novel_sum <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
beta_novel_num <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)

res_kleb_arg <- run_subiso_sensitivity_grid(
  iso_df,
  gene_cols = gene_cols,
  beta_named = beta_named, 
  beta_novel_sum = beta_novel_sum, 
  beta_novel_num = beta_novel_num, 
  B = 1000,
  seed = 2026,
  out_dir = "sens_results/kleb_arg",
  prefix = "kleb",
  dataset_label = "Klebsiella",
  save_csv = TRUE,
  keep_draws = FALSE,
  compute_curves = FALSE
)
#View(res_kleb_arg$fit_metrics)
#View(res_kleb_arg$feature_summary)
res_kleb_arg <- c()
res_kleb_arg$fit_metrics <- read.csv("sens_results/kleb_arg/kleb_fit_metrics.csv")
res_kleb_arg$feature_summary <- read.csv("sens_results/kleb_arg/kleb_feature_summary.csv")


# rmse heatmap 
heatmap_df <- res_kleb_arg$fit_metrics |>
  filter(beta_named %in% c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000))
heatmap_kleb <- plot_fit_heatmap_subiso(heatmap_df, 
                                         beta_named_vals = beta_named, 
                                         beta_novel_sum_vals  = beta_novel_sum, 
                                         beta_novel_num_vals = beta_novel_num,
                                         fill_val = "rmse_log10" , 
                                         fill_label = "log-RMSE")
heatmap_kleb
ggsave("sens_results/kleb_arg_simpleBB_sensitivity_heatmap.png", heatmap_kleb, width = 10, height = 8, dpi = 300)

# full panel plot
panel_df <- res_kleb_arg$feature_summary |>
  filter(beta_named %in% c(0.5, 1, 10, 100, 1000),
         beta_novel_sum %in% c(0, 1, 10, 100, 1000),
         beta_novel_num %in% c(0, 1, 10, 100, 1000))

panel_kleb <- plot_obs_vs_pred_subiso(panel_df, 
                                       beta_named_vals = c(0.5, 1, 10, 100, 1000),
                                       beta_novel_sum_vals = c(0, 1, 10, 100, 1000),
                                       beta_novel_num_vals = c(0, 1, 10, 100, 1000))
panel_kleb
ggsave("sens_results/kleb_arg_simpleBB_sensitivity_panel.png", panel_kleb, width = 10, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
combined_mass_df <-  dplyr::bind_rows(res_ecoli_arg$mass_curves, res_kleb_arg$mass_curves) |>
  filter(alpha %in% alpha_vals)
mass_plot  <- plot_curve_subiso_facet(
  combined_mass_df,
  alpha_vals,
  value_col = "median_mass", lo_col = "q2.5", hi_col = "q97.5",
  y_lab = "Proportion of mass above threshold",
  title = "Cumulative population mass curves",
  ncol = 6
)
mass_plot
ggsave("sens_results/ecoli_kleb_arg_simpleBB_sensitivity_mass_plot.png", mass_plot, width = 15, height = 10, dpi = 300)


combined_richness_df <-  dplyr::bind_rows(res_ecoli_arg$richness_curves, res_kleb_arg$richness_curves) |>
  filter(alpha %in% alpha_vals)

richness_plot  <- plot_curve_subiso_facet(
  combined_richness_df,
  alpha_vals,
  value_col = "median_species_proportion",
  lo_col = "q2.5_species_proportion",
  hi_col = "q97.5_species_proportion",
  y_lab = "Proportion of genes above threshold",
  title = "Cumulative species richness curves",
  ncol = 6
)
richness_plot
ggsave("sens_results/ecoli_kleb_arg_simpleBB_sensitivity_richness_plot.png", richness_plot, width = 15, height = 10, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 2.1bii Run with single prior parameter values ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
df <- kleb_bsi_arg_df |>
  group_by(Element.symbol) |>
  summarise(count = n()) |>
  arrange(count)
N_isolate <- length(unique(kleb_bsi_arg_df$sample)) # 468
N_gene <- sum(df$count) # number of genes 3258
K <- length(unique(kleb_bsi_arg_df$Element.symbol)) # num unique features 199
f1 <- sum(df$count == 1) # 64 - no. singletons
u_hat <- f1 / N_isolate  # 0.137 - frequency of singletons
u_hat

# set priors
beta_named <- rep(1, K)  # set uninformative priors
names(beta_named) <- df$Element.symbol
# set novel priors to preserve 'novel' mass 
beta_novel_num = ceiling((u_hat * K)/ (1 - u_hat)) # derived to make the prior novel frequency the same as in the data
beta_novel_sum = beta_novel_num * beta_prior #

# run
kleb_bsi_arg_bbs <- subisolate_bayes_beta(
  iso_df = iso_df,
  gene_cols = gene_cols,
  beta_named = beta_named,
  beta_novel_sum = beta_novel_sum,
  beta_novel_num = beta_novel_sum,
  B = 10000)

# save
saveRDS(kleb_bsi_arg_bbs, "model_results/kleb_bsi_ARG_bayes.rds")
#~~~~~~~~~~~~~#
# read in saved data
#kleb_bsi_arg_bbs <- readRDS("model_results/kleb_bsi_ARG_bayes.rds")
#~~~~~~~~~~~~~#
#  Prevalence summary 
kleb_bsi_prevalence_summary <- kleb_bsi_arg_bbs$summary_df
#View(kleb_bsi_prevalence_summary)

# plot histogram of detection probabilities:
arg_post_hist <- ggplot(data = kleb_bsi_prevalence_summary) +
  geom_histogram(aes(x = median_prevalence), binwidth = 0.01, alpha= 0.5, fill = "darkorange") +
  geom_histogram(aes(x = q2.5_prevalence), binwidth = 0.01, alpha= 0.7, fill = "darkorange") +
  geom_histogram(aes(x = q97.5_prevalence), binwidth = 0.01, alpha= 0.2, fill = "darkorange") +
  labs(title = "Histogram of posterior estimated ARG prevalence", 
       x = "ARG prevalence\n(median and 95% CIs)" ) +
  theme_minimal()
arg_post_hist
# save
ggsave("model_results/kleb_bsi_arg_posterior_distribution_histogram.png", arg_post_hist, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define frequency matrix
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

#extract posterior prevalence matrix 
p_mat <- get_p_matrix(kleb_bsi_arg_bbs)

# cumulative species richness curves
kleb_bsi_arg_richness_df <- compute_species_richness_subiso(p_mat, f_grid = f_grid_99)
# cumulative population mass curves
kleb_bsi_arg_mass_df  <- compute_mass_curve_subiso(p_mat, f_grid = f_grid_99)

# add sample sizes 
kleb_bsi_arg_richness_df <- kleb_bsi_arg_richness_df |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_arg_richness_df)
# save
write.csv(kleb_bsi_arg_richness_df, "model_results/kleb_bsi_ARG_bayes_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#kleb_bsi_arg_richness_df <- read.csv("model_results/kleb_bsi_ARG_bayes_richness_df.csv")
#~~~~~~~~~~~~~~~#

# min sample sizes using f as feature frequency
kleb_bsi_arg_mass_df <- kleb_bsi_arg_mass_df |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_arg_mass_df)
# save
write.csv(kleb_bsi_arg_mass_df, "model_results/kleb_bsi_ARG_bayes_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#kleb_bsi_arg_mass_df <- read.csv("model_results/kleb_bsi_ARG_bayes_mass_df.csv")
#~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.1c Combined plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# combined Prevalence Histogram
ecoli_bsi_prevalence_summary <- ecoli_bsi_prevalence_summary |>
  mutate(Genus = "Escherichia")
kleb_bsi_prevalence_summary <- kleb_bsi_prevalence_summary |>
  mutate(Genus = "Klebsiella")
combined_bsi_prevalence_summary <- rbind(ecoli_bsi_prevalence_summary, kleb_bsi_prevalence_summary)
combined_bsi_prevalence_summary <- combined_bsi_prevalence_summary |>
  rename(median = median_prevalence, 
         q2.5 = q2.5_prevalence,
         q97.5 = q97.5_prevalence
         )
#View(combined_bsi_prevalence_summary)

# plot combined histogram of detection probabilities:
genus_colours <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
combined_bsi_arg_post_hist <- ggplot(data = combined_bsi_prevalence_summary, (aes(colour = Genus, fill = Genus))) +
  geom_histogram(aes(x = median, y = after_stat(density)), binwidth = 0.01, alpha= 0.4, colour = NA, position = "identity") +
  geom_histogram(aes(x = q2.5, y = after_stat(density)), binwidth = 0.01, alpha= 0.6, colour = NA, position = "identity") +
  geom_histogram(aes(x = q97.5, y = after_stat(density)), binwidth = 0.01, alpha= 0.2, colour = NA, position = "identity") +
  scale_colour_manual(values = genus_colours) +
  scale_fill_manual(values = genus_colours) +
  labs(
    title = "Histogram of posterior estimated ARG prevalence", 
    x = "ARG prevalence\n(median and 95% CIs)",
    y = "Density (%)") +
  theme_minimal()
combined_bsi_arg_post_hist
# save
ggsave("model_results/combined_bsi_arg_posterior_distribution_histogram.png", combined_bsi_arg_post_hist, width = 8, height = 4, units = "in", dpi = 300)

# combined E. coli and Klebsiella cumulative species richness curve 
combined_bsi_arg_richness_df <- rbind(ecoli_bsi_arg_richness_df , kleb_bsi_arg_richness_df )
# save
write.csv(combined_bsi_arg_richness_df, "model_results/combined_bsi_arg_bayes_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~#
# read in saved df
#combined_bsi_arg_richness_df <- read.csv("model_results/combined_bsi_arg_bayes_richness_df.csv")
#~~~~~~~~~~~~~#

# combined E. coli and Klebsiella cumulative mass curve (normalised to 1) 
combined_bsi_arg_mass_df <- rbind(ecoli_bsi_arg_mass_df , kleb_bsi_arg_mass_df )
#View(combined_bsi_arg_mass_df)
# save
write.csv(combined_bsi_arg_mass_df, "model_results/combined_bsi_arg_bayes_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~#
# read in saved df
#combined_bsi_arg_mass_df <- read.csv("model_results/combined_bsi_arg_bayes_mass_df.csv")
#~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~#
# plot cumulative spcies richness count
# define colours
genus_colours <- c("Escherichia" = "seagreen3", "Klebsiella" = "darkorange")
#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot (cumulative species_richness of MLSTs at frequency >= f)
cumulative_species_richness_plot <- ggplot(combined_bsi_arg_richness_df, aes(x = f, y = median_species_proportion, colour =  Genus, fill = Genus)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Genus", values = genus_colours) +
  scale_colour_manual(name = "Genus", values = genus_colours) +
  scale_x_log10( breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(combined_bsi_arg_richness_df$f), 1)) +
  labs(x = "ARG freqeuncy (f)",
       y = "Proportion of ARGs of frequency ≥ f"
  ) +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
cumulative_species_richness_plot
# save
ggsave("model_results/combined_bsi_arg_bayes_cumulative_species_richness_plot.png", plot = cumulative_species_richness_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Cumulative fastBAS clusters fraction with sample size plot 
bayes_sample_coverage_plot <- ggplot(combined_bsi_arg_richness_df, aes(y = median_species_proportion, colour = Genus, fill = Genus)) +
  geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  geom_line(aes(x = min_sample_99)) +
  geom_ribbon(aes(x = min_sample_90, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.35, colour = NA) +
  geom_ribbon(aes(x = min_sample_99, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  #facet_wrap(~ Genus, ncol = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,5000)) +
  labs(y = "Proportion of unique ARGs detected",
       x = "Sample size",
  ) +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/combined_bsi_arg_bayes_ss_vs_species_richness_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~#
# plot (median + 95% CI) 
combined_cumulative_mass_plot <- ggplot(combined_bsi_arg_mass_df, aes(x = f, colour = Genus, fill = Genus)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, colour = NA) +
  geom_line(aes(y = median_mass), size = 1) +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(min(combined_bsi_arg_mass_df$f), 1)) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  labs(
    x = "ARG frequency (f) (log scale)",
    y = "Proportion of ARGs belonging to an allele of frequency ≥ f",
    #title = "Cumulative mass curve: proportion of population at frequency ≥ f",
    #subtitle = "median (line) and 95% posterior interval (ribbon)"
  ) +
  theme_minimal(base_size = 14)
combined_cumulative_mass_plot
#save
ggsave("model_results/combined_bsi_arg_cumulative_mass.png", combined_cumulative_mass_plot, width = 8, height = 6, dpi = 300)

#~~~~~~~~~~~~~~~~~~~#
# Combined E. coli and Klebsiella coverage plot
bayes_sample_coverage_plot <- ggplot(combined_bsi_arg_mass_df, aes(y = median_mass, colour = Genus, fill = Genus)) +
  geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  geom_line(aes(x = min_sample_99)) +
  geom_ribbon(aes(x = min_sample_90, ymin = q2.5, ymax = q97.5), alpha = 0.15, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5, ymax = q97.5), alpha = 0.3, colour = NA) +
  geom_ribbon(aes(x = min_sample_99, ymin = q2.5, ymax = q97.5), alpha = 0.45, colour = NA) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,5000)) +
  labs(y = "Coverage",
       x = "Sample size") +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/combined_bsi_arg_bayes_sample_coverage_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.1d Combined summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

# Summary tables for E. coli and Klebsiella ARGs - species rich tables
# transform to long
combined_bsi_arg_richness_df_long <- combined_bsi_arg_richness_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_arg_richness_summary_table <- combined_bsi_arg_richness_df_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_species_proportion", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5_species_proportion", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5_species_proportion", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# tidy
combined_bsi_arg_richness_summary_table <- combined_bsi_arg_richness_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_arg_richness_summary_table)
#View(combined_bsi_arg_richness_summary_table)

write.csv(combined_bsi_arg_richness_summary_table, "model_results/combined_bsi_arg_species_richness_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Summary tables for E. coli and Klebsiella ARGs - cumulative population mass 
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
      median_ss  <- min_sample_at_or_above(df, "median_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_arg_sample_coverage_summary_table <- combined_bsi_arg_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_arg_sample_coverage_summary_table)

write.csv(combined_bsi_arg_sample_coverage_summary_table, "model_results/combined_bsi_arg_sample_coverage_summary_table.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 2.2 Plasmids ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.2a E. coli plasmids ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare E. coli plasmids df
all_samples <- ecoli_bsi_amrfinder_metadata |> distinct(sample)   

# detected pairs from AMRFinder (1 if detected)
detected <- ecoli_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)
length(unique(ecoli_bsi_amrfinder_metadata$community_subcommunity))

# full grid of sample x feature using master sample list and the set of observed features
all_features <- detected |> pull(community_subcommunity) |> unique()
full_grid <- tidyr::expand_grid(sample = all_samples$sample,
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
counts <- presence_long |>
  group_by(community_subcommunity) |>
  summarise(count = sum(presence), .groups = "drop")
#table(counts$count)
sum(counts[counts$count ==1,]$count)

# Set params
iso_df <- ecoli_bsi_pling_df
iso_df <- iso_df |> dplyr::ungroup()
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))
#length(unique(iso_df$sample))# 1471
#View(iso_df)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 2.2ai Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
beta_named <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
beta_novel_sum <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
beta_novel_num <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)


res_ecoli_pling <- run_subiso_sensitivity_grid(
  iso_df,
  gene_cols = gene_cols,
  beta_named = beta_named, 
  beta_novel_sum = beta_novel_sum, 
  beta_novel_num = beta_novel_num, 
  B = 1000,
  seed = 2026,
  out_dir = "sens_results/ecoli_pling",
  prefix = "ecoli",
  dataset_label = "Escherichia",
  save_csv = TRUE,
  keep_draws = FALSE,
  compute_curves = FALSE
)
#View(res_ecoli_pling$fit_metrics)
#View(res_ecoli_pling$feature_summary)
res_ecoli_pling <- c()
res_ecoli_pling$fit_metrics <- read.csv("sens_results/ecoli_pling/ecoli_fit_metrics.csv")
res_ecoli_pling$feature_summary <- read.csv("sens_results/ecoli_pling/ecoli_feature_summary.csv")


# rmse heatmap 
heatmap_df <- res_ecoli_pling$fit_metrics |>
  filter(beta_named %in% c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000))
heatmap_ecoli <- plot_fit_heatmap_subiso(heatmap_df, 
                                         beta_named_vals = beta_named, 
                                         beta_novel_sum_vals  = beta_novel_sum, 
                                         beta_novel_num_vals = beta_novel_num,
                                         fill_val = "rmse_log10" , 
                                         fill_label = "log-RMSE")
heatmap_ecoli
ggsave("sens_results/ecoli_pling_simpleBB_sensitivity_heatmap.png", heatmap_ecoli, width = 10, height = 8, dpi = 300)

# full panel plot
panel_df <- res_ecoli_pling$feature_summary |>
  filter(beta_named %in% c(0.5, 1, 10, 100, 1000),
         beta_novel_sum %in% c(0, 1, 10, 100, 1000),
         beta_novel_num %in% c(0, 1, 10, 100, 1000))

panel_ecoli <- plot_obs_vs_pred_subiso(panel_df, 
                                       beta_named_vals = c(0.5, 1, 10, 100, 1000),
                                       beta_novel_sum_vals = c(0, 1, 10, 100, 1000),
                                       beta_novel_num_vals = c(0, 1, 10, 100, 1000))
panel_ecoli
ggsave("sens_results/ecoli_pling_simpleBB_sensitivity_panel.png", panel_ecoli, width = 10, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 2.2aii Run single prior parameter set ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
N_isolate <- length(unique(ecoli_bsi_pling_df$sample)) # number of genes
N_gene <- sum(presence_long$presence) # 3361 number of plasmids
K <- length(unique(presence_long$community_subcommunity)) # 825 num unique features
f1 <- sum(counts$count == 1) #  622- no. singletons
u_hat <- f1 / N_isolate  # frequency of singletons
u_hat

# set priors
beta_prior <- 1
beta_named <- rep(beta_prior, K)  # set uninformative priors
names(beta_named) <- gene_cols
# set novel priors to preserve 'novel' mass, without and with including singletons 
beta_novel_num = ceiling((u_hat * K)/ (1 - u_hat)) # derived to make the prior novel frequency the same as in the data
beta_novel_sum = beta_novel_num * beta_prior
beta_novel_num
beta_novel_sum

# run
ecoli_bsi_pling_bbs <- subisolate_bayes_beta(
  iso_df = iso_df,
  gene_cols = gene_cols,
  beta_named = beta_named,
  beta_novel_sum = beta_novel_sum,
  beta_novel_num = beta_novel_sum, # many uncommon novel features, as opposed to 1 common novel block
  B = 10000)

# save
saveRDS(ecoli_bsi_pling_bbs, "model_results/ecoli_bsi_PLING_bayes.rds")
#~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_pling_bbs <- read_rds("model_results/ecoli_bsi_PLING_bayes.rds")
#~~~~~~~~~~~~#
# Prevalence summary
ecoli_bsi_pling_prevalence_summary <- ecoli_bsi_pling_bbs$summary_df
#print(ecoli_bsi_pling_prevalence_summary)

# plot histogram of detection probabilities:
pling_post_hist <- ggplot(data = ecoli_bsi_pling_prevalence_summary) +
  geom_histogram(aes(x = median_prevalence), binwidth = 0.01, alpha= 0.5,  fill = "seagreen3") +
  geom_histogram(aes(x = q2.5_prevalence), binwidth = 0.01, alpha= 0.7,  fill = "seagreen3") +
  geom_histogram(aes(x = q97.5_prevalence), binwidth = 0.01, alpha= 0.2,  fill = "seagreen3") +
  labs(title = "Histogram of posterior estimated PLING plasmid subcommunity prevalence", 
       x = "ARG prevalence\n(median and 95% CIs)" ) +
  theme_minimal()
pling_post_hist
# save
ggsave("model_results/ecoli_bsi_pling_posterior_distribution_histogram.png", pling_post_hist, width = 8, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define frequency matrix
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

#extract posterior prevalence matrix 
p_mat <- get_p_matrix(ecoli_bsi_pling_bbs)

# cumulative species richness curves
ecoli_bsi_pling_richness_df <- compute_species_richness_subiso(p_mat, f_grid = f_grid_99)
# cumulative population mass curves
ecoli_bsi_pling_mass_df  <- compute_mass_curve_subiso(p_mat, f_grid = f_grid_99)

# add sample sizes 
ecoli_bsi_pling_richness_df <- ecoli_bsi_pling_richness_df |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_pling_richness_df)
# save
write.csv(ecoli_bsi_pling_richness_df, "model_results/ecoli_bsi_pling_bayes_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_pling_richness_df <- read.csv("model_results/ecoli_bsi_pling_bayes_richness_df.csv")
#~~~~~~~~~~~~~~~#

# min sample sizes using f as feature frequency
ecoli_bsi_pling_mass_df <- ecoli_bsi_pling_mass_df |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_pling_mass_df)
# save
write.csv(ecoli_bsi_pling_mass_df, "model_results/ecoli_bsi_pling_bayes_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_pling_mass_df <- read.csv("model_results/ecoli_bsi_pling_bayes_mass_df.csv")
#~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.2b Klebsiella plasmid subcommunities ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prep data
all_samples <- kleb_bsi_amrfinder_metadata |> distinct(sample)   

# detected pairs from AMRFinder (1 if detected)
detected <- kleb_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)

# full grid of sample x feature using master sample list and the set of observed features
all_features <- detected |> pull(community_subcommunity) |> unique()
full_grid <- tidyr::expand_grid(sample = all_samples$sample,
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

counts <- presence_long |>
  group_by(community_subcommunity) |>
  summarise(count = sum(presence), .groups = "drop")
#table(counts$count)
sum(counts[counts$count ==1,]$count)

# Set params
iso_df <- kleb_bsi_pling_df
iso_df <- iso_df |> dplyr::ungroup()
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))
#length(unique(iso_df$sample)) # 468
#View(iso_df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 2.2ai Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
beta_named <- c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
beta_novel_sum <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)
beta_novel_num <- c(0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000)


res_kleb_pling <- run_subiso_sensitivity_grid(
  iso_df,
  gene_cols = gene_cols,
  beta_named = beta_named, 
  beta_novel_sum = beta_novel_sum, 
  beta_novel_num = beta_novel_num, 
  B = 1000,
  seed = 2026,
  out_dir = "sens_results/kleb_pling",
  prefix = "kleb",
  dataset_label = "Escherichia",
  save_csv = TRUE,
  keep_draws = FALSE,
  compute_curves = FALSE
)
#View(res_kleb_pling$fit_metrics)
#View(res_kleb_pling$feature_summary)
res_kleb_pling <- c()
res_kleb_pling$fit_metrics <- read.csv("sens_results/kleb_pling/kleb_fit_metrics.csv")
res_kleb_pling$feature_summary <- read.csv("sens_results/kleb_pling/kleb_feature_summary.csv")

# rmse heatmap 
heatmap_df <- res_kleb_arg$fit_metrics |>
  filter(beta_named %in% c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000))
heatmap_kleb <- plot_fit_heatmap_subiso(heatmap_df, 
                                         beta_named_vals = beta_named, 
                                         beta_novel_sum_vals  = beta_novel_sum, 
                                         beta_novel_num_vals = beta_novel_num,
                                         fill_val = "rmse_log10" , 
                                         fill_label = "log-RMSE")
heatmap_kleb
ggsave("sens_results/kleb_pling_simpleBB_sensitivity_heatmap.png", heatmap_kleb, width = 10, height = 8, dpi = 300)

# full panel plot
panel_df <- res_kleb_arg$feature_summary |>
  filter(beta_named %in% c(0.5, 1, 10, 100, 1000),
         beta_novel_sum %in% c(0, 1, 10, 100, 1000),
         beta_novel_num %in% c(0, 1, 10, 100, 1000))

panel_kleb <- plot_obs_vs_pred_subiso(panel_df, 
                                       beta_named_vals = c(0.5, 1, 10, 100, 1000),
                                       beta_novel_sum_vals = c(0, 1, 10, 100, 1000),
                                       beta_novel_num_vals = c(0, 1, 10, 100, 1000))
panel_kleb
ggsave("sens_results/kleb_pling_simpleBB_sensitivity_panel.png", panel_kleb, width = 10, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
combined_mass_df <-  dplyr::bind_rows(res_ecoli_pling$mass_curves, res_kleb_pling$mass_curves) |>
  filter(alpha %in% alpha_vals)
mass_plot  <- plot_curve_subiso_facet(
  combined_mass_df,
  alpha_vals,
  value_col = "median_mass", lo_col = "q2.5", hi_col = "q97.5",
  y_lab = "Proportion of mass above threshold",
  title = "Cumulative population mass curves",
  ncol = 6
)
mass_plot
ggsave("sens_results/ecoli_kleb_pling_simpleBB_sensitivity_mass_plot.png", mass_plot, width = 15, height = 10, dpi = 300)


combined_richness_df <-  dplyr::bind_rows(res_ecoli_pling$richness_curves, res_kleb_pling$richness_curves) |>
  filter(alpha %in% alpha_vals)

richness_plot  <- plot_curve_subiso_facet(
  combined_richness_df,
  alpha_vals,
  value_col = "median_species_proportion",
  lo_col = "q2.5_species_proportion",
  hi_col = "q97.5_species_proportion",
  y_lab = "Proportion of genes above threshold",
  title = "Cumulative species richness curves",
  ncol = 6
)
richness_plot
ggsave("sens_results/ecoli_kleb_pling_simpleBB_sensitivity_richness_plot.png", richness_plot, width = 15, height = 10, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 2.2bii Run single prior parameter set ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
N_isolate <- length(unique(kleb_bsi_pling_df$sample)) # 468
N_gene <- sum(presence_long$presence) # 847 number of plasmids
K <- length(unique(presence_long$community_subcommunity)) # 417 num unique features
f1 <- sum(counts$count == 1) #  330- no. singletons
u_hat <- f1 / N_isolate # frequency of singletons
u_hat

# set priors
beta_prior <- 1
beta_named <- rep(beta_prior, K)  # set uninformative priors
names(beta_named) <- gene_cols
# set novel priors to preserve 'novel' mass, without and with including singletons 
beta_novel_num = ceiling((u_hat * K)/ (1 - u_hat)) # derived to make the prior novel frequency the same as in the data
beta_novel_sum = beta_novel_num * beta_prior # 998
beta_novel_num
beta_novel_sum

# run
kleb_bsi_pling_bbs <- subisolate_bayes_beta(
  iso_df = iso_df,
  gene_cols = gene_cols,
  beta_named = beta_named,
  beta_novel_sum = beta_novel_sum,
  beta_novel_num = beta_novel_num, # many uncommon novel features, as opposed to 1 common novel block
  B = 10000)

# save
saveRDS(kleb_bsi_pling_bbs, "model_results/kleb_bsi_PLING_bayes.rds")
#~~~~~~~~~~~~~#
# read in saved data
#kleb_bsi_pling_bbs <- read_rds("model_results/kleb_bsi_PLING_bayes.rds")
#~~~~~~~~~~~~~#

# Prevalence summary 
kleb_bsi_pling_prevalence_summary <- kleb_bsi_pling_bbs$summary_df
#print(kleb_bsi_pling_prevalence_summary)

# plot histogram of detection probabilities:
pling_post_hist <- ggplot(data = kleb_bsi_pling_prevalence_summary) +
  geom_histogram(aes(x = median_prevalence), binwidth = 0.01, alpha= 0.5,  fill = "darkorange") +
  geom_histogram(aes(x = q2.5_prevalence), binwidth = 0.01, alpha= 0.7,  fill = "darkorange") +
  geom_histogram(aes(x = q97.5_prevalence), binwidth = 0.01, alpha= 0.2,  fill = "darkorange") +
  labs(title = "Histogram of posterior estimated PLING plasmid subcommunity prevalence", 
       x = "ARG prevalence\n(median and 95% CIs)" ) +
  theme_minimal()
pling_post_hist
# save
ggsave("model_results/kleb_bsi_pling_posterior_distribution_histogram.png", pling_post_hist, width = 8, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define frequency matrix
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

#extract posterior prevalence matrix 
p_mat <- get_p_matrix(kleb_bsi_pling_bbs)

# cumulative species richness curves
kleb_bsi_pling_richness_df <- compute_species_richness_subiso(p_mat, f_grid = f_grid_99)
# cumulative population mass curves
kleb_bsi_pling_mass_df <- compute_mass_curve_subiso(p_mat, f_grid = f_grid_99)

# add sample sizes 
kleb_bsi_pling_richness_df <- kleb_bsi_pling_richness_df |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_pling_richness_df)
# save
write.csv(kleb_bsi_pling_richness_df, "model_results/kleb_bsi_pling_bayes_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#kleb_bsi_pling_richness_df <- read.csv("model_results/kleb_bsi_pling_bayes_richness_df.csv")
#~~~~~~~~~~~~~~~#

# min sample sizes using f as feature frequency
kleb_bsi_pling_mass_df <- kleb_bsi_pling_mass_df |>
  mutate(
    min_sample_90 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.90) / log(1 - f))),
    min_sample_95 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.95) / log(1 - f))),
    min_sample_99 = ifelse(f >= 1, NA_real_, ceiling(log(1 - 0.99) / log(1 - f)))
  ) |>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_pling_mass_df)
# save
write.csv(kleb_bsi_pling_mass_df, "model_results/kleb_bsi_pling_bayes_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in saved data
#kleb_bsi_pling_mass_df <- read.csv("model_results/kleb_bsi_pling_bayes_mass_df.csv")
#~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.2c Combined plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# combined population richness df
combined_bsi_pling_richness_df <- rbind(ecoli_bsi_pling_richness_df, kleb_bsi_pling_richness_df)
# save
write.csv(combined_bsi_pling_richness_df, "model_results/combined_bsi_pling_bayes_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~#
# read-in saved df 
#combined_bsi_pling_richness_df <- read.csv("model_results/combined_bsi_pling_bayes_richness_df.csv")
#~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~#
# combined population mass df
combined_bsi_pling_mass_df <- rbind(ecoli_bsi_pling_mass_df, kleb_bsi_pling_mass_df)
# save
write.csv(combined_bsi_pling_mass_df, "model_results/combined_bsi_pling_bayes_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~#
# read-in saved df 
#combined_bsi_pling_mass_df <- read.csv("model_results/combined_bsi_pling_bayes_mass_df.csv")
#View(combined_bsi_pling_mass_df)
#~~~~~~~~~~~~~~~~~~#

# Combined Prevalence Histogram with E. coli
ecoli_bsi_pling_prevalence_summary <- ecoli_bsi_pling_prevalence_summary |>
  mutate(Genus = "Escherichia")
kleb_bsi_pling_prevalence_summary <- kleb_bsi_pling_prevalence_summary |>
  mutate(Genus = "Klebsiella")
combined_bsi_pling_prevalence_summary <- rbind(ecoli_bsi_pling_prevalence_summary, kleb_bsi_pling_prevalence_summary)
#View(combined_bsi_pling_prevalence_summary)

# define colours
genus_colours <- c("Escherichia" = "seagreen3", "Klebsiella" = "darkorange")

# plot combined histogram of detection probabilities:
combined_bsi_pling_post_hist <- ggplot(data = combined_bsi_pling_prevalence_summary, (aes(colour = Genus, fill = Genus))) +
  geom_histogram(aes(x = median_prevalence, y = after_stat(density)), binwidth = 0.01, alpha= 0.4, colour = NA, position = "identity") +
  geom_histogram(aes(x = q2.5_prevalence, y = after_stat(density)), binwidth = 0.01, alpha= 0.6, colour = NA, position = "identity") +
  geom_histogram(aes(x = q97.5_prevalence, y = after_stat(density)), binwidth = 0.01, alpha= 0.2, colour = NA, position = "identity") +
  scale_colour_manual(values = genus_colours) +
  scale_fill_manual(values = genus_colours) +
  labs(
    title = "Histogram of posterior estimated plasmid subcommunity prevalence", 
    x = "Plasmid subcommunity prevalence\n(median and 95% CIs)",
    y = "Density (%)") +
  theme_minimal()
combined_bsi_pling_post_hist
# save
ggsave("model_results/combined_bsi_pling_posterior_distribution_histogram.png", combined_bsi_pling_post_hist, width = 8, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~#
# plot cumulative spcies richness count
#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot (cumulative species_richness of MLSTs at frequency >= f)
cumulative_species_richness_plot <- ggplot(combined_bsi_pling_richness_df, aes(x = f, y = median_species_proportion, colour =  Genus, fill = Genus)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Genus", values = genus_colours) +
  scale_colour_manual(name = "Genus", values = genus_colours) +
  scale_x_log10( breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(combined_bsi_pling_richness_df$f), 1)) +
  labs(x = "Plasmid subcommunity frequency (f) (log scale)",
       y = "Proportion of Plasmid subcommunities of frequency ≥ f"
  ) +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
cumulative_species_richness_plot
# save
ggsave("model_results/combined_bsi_pling_bayes_cumulative_species_richness_plot.png", plot = cumulative_species_richness_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Cumulative feature proportion with sample size plot 
bayes_sample_coverage_plot <- ggplot(combined_bsi_pling_richness_df, aes(y = median_species_proportion, colour = Genus, fill = Genus)) +
  geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  geom_line(aes(x = min_sample_99)) +
  geom_ribbon(aes(x = min_sample_90, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.35, colour = NA) +
  geom_ribbon(aes(x = min_sample_99, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  #facet_wrap(~ Genus, ncol = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,5000)) +
  labs(y = "Proportion of unique Plasmid subcommunities detected",
       x = "Sample size",
  ) +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/combined_bsi_pling_bayes_ss_vs_species_richness_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~#
# combined cumulative mass plot (median + 95% CI) 
combined_cumulative_mass_plot <- ggplot(combined_bsi_pling_mass_df, aes(x = f, colour = Genus, fill = Genus)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, colour = NA) +
  geom_line(aes(y = median_mass), size = 1) +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(min(combined_bsi_pling_mass_df$f), 1)) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
  labs(
    x = "Plasmid subcommunity frequency (f) (log scale)",
    y = "Proportion of plasmids in a subcommunity of frequency ≥ f",
    #title = "Cumulative mass curve: proportion of population at frequency ≥ f",
    #subtitle = "median (line) and 95% posterior interval (ribbon)"
  ) +
  theme_minimal(base_size = 14)
combined_cumulative_mass_plot
#save
ggsave("model_results/combined_bsi_pling_cumulative_mass.png", combined_cumulative_mass_plot, width = 8, height = 6, dpi = 300)


# Combined coverage vs sample size plot
bayes_sample_coverage_plot <- ggplot(combined_bsi_pling_mass_df, aes(y = median_mass, colour = Genus, fill = Genus)) +
  geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  geom_line(aes(x = min_sample_99)) +
  geom_ribbon(aes(x = min_sample_90, ymin = q2.5, ymax = q97.5), alpha = 0.15, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5, ymax = q97.5), alpha = 0.3, colour = NA) +
  geom_ribbon(aes(x = min_sample_99, ymin = q2.5, ymax = q97.5), alpha = 0.45, colour = NA) +
  scale_fill_manual(values = genus_colours) +
  scale_colour_manual(values = genus_colours) +
 # facet_wrap(~ Genus, ncol = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_y_continuous() +
  scale_x_continuous(limits = c(0,5000)) +
  labs(y = "coverage",
       x = "Sample size") +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/combined_bsi_PLING_bayes_sample_coverage_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 2.2d Combined summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# transform to long
combined_bsi_pling_richness_df_long <- combined_bsi_pling_richness_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
combined_bsi_pling_richness_summary_table <- combined_bsi_pling_richness_df_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_species_proportion", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5_species_proportion", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5_species_proportion", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
combined_bsi_pling_richness_summary_table <- combined_bsi_pling_richness_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(combined_bsi_pling_richness_summary_table)
#View(combined_bsi_pling_richness_summary_table)
# save
write.csv(combined_bsi_pling_richness_summary_table, "model_results/combined_bsi_pling_richness_summary_table.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative mass df
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
      median_ss  <- min_sample_at_or_above(df, "median_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
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
write.csv(combined_bsi_pling_sample_coverage_summary_table, "model_results/combined_bsi_pling_sample_coverage_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. PANEL PLOT FOR E. coli and Kleb posterior estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load cumulative species richness data 
richness_df_mlst <- read.csv("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
richness_df_fastbaps_L3 <- read.csv("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_fine.csv")
combined_bsi_pling_richness_df <- read.csv("model_results/combined_bsi_pling_bayes_richness_df.csv")
combined_bsi_arg_richness_df <- read.csv("model_results/combined_bsi_arg_bayes_richness_df.csv")

# load mass curve data
mass_df_mlst <- read.csv("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
mass_df_fastbaps_L3 <- read.csv("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv")
combined_bsi_pling_mass_df <- read.csv("model_results/combined_bsi_pling_bayes_mass_df.csv")
combined_bsi_arg_mass_df <- read.csv("model_results/combined_bsi_arg_bayes_mass_df.csv")

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

make_richness_plot <- function(df, xlab, ylab, xlim_min = NULL) {
  p <- ggplot(df, aes(x = f, y = median_species_proportion , colour = Genus, fill = Genus)) +
    geom_ribbon(aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.4, colour = NA) +
    geom_line(linewidth = 1) +
    scale_fill_manual(values = genus_colours) +
    scale_colour_manual(values = genus_colours) +
    scale_x_log10(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1),
      labels = c("0.00001","0.0001", "0.001", "0.01", "0.1", "1")) +
    labs(x = xlab, y = ylab) +
    big_theme
  
  if (!is.null(xlim_min)) {
    p <- p + coord_cartesian(xlim = c(xlim_min, 1))
  }
  p
}

make_ss_vs_richness_plot <- function(df, ylab = "Proportion of features detected", xlab = "Sample size", x_limit = c(10, 10000)) {
    df_long <- df |>
    pivot_longer(cols = c(min_sample_90, min_sample_95, min_sample_99),
                 names_to = "confidence",
                 values_to = "sample_size" ) |>
    mutate(confidence = factor(confidence,
                               levels = c("min_sample_90", "min_sample_95", "min_sample_99"),
                               labels = c("90%", "95%", "99%"))) |>
      filter(median_species_proportion != 0)
    
    df_90 <- df_long |> filter(confidence == "90%")
    df_95 <- df_long |> filter(confidence == "95%")
    df_99 <- df_long |> filter(confidence == "99%") 

    ggplot(df_long, aes(y = median_species_proportion, colour = Genus)) +
      geom_ribbon(data = df_90,aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion, x = sample_size, fill = Genus, group = Genus),
        inherit.aes = FALSE, alpha = 0.10, colour = NA) +
      geom_ribbon(data = df_95, aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion, x = sample_size, fill = Genus,  group = Genus),
        inherit.aes = FALSE,alpha = 0.30, colour = NA) +
      geom_ribbon(data = df_99, aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion, x = sample_size, fill = Genus,  group = Genus),
        inherit.aes = FALSE, alpha = 0.50, colour = NA) +
      geom_line(data = df_90, aes(x = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_line(data = df_95, aes(x = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_line(data = df_99, aes(x = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_hline(yintercept = 0.80, linetype = "dashed") +
      scale_colour_manual(values = genus_colours, name = "Genus") +
      scale_fill_manual(values = genus_colours, name = "Genus") +
      scale_alpha_manual( values = c("90%" = 0.3, "95%" = 0.55, "99%" = 0.85), name = "Certainty") +
      guides(colour = "none", fill = "none", alpha = guide_legend(title = "Certainty")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_log10(limits = x_limit, 
                    breaks = c(10, 100, 1000, 10000),
                    labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      labs(x = xlab, y = ylab) +
      big_theme +
      theme(legend.position = "bottom")
}

# cumulative mass plot function
make_mass_plot <- function(df, xlab, ylab, xlim_min = NULL) {
  p <- ggplot(df, aes(x = f, y = median_mass, colour = Genus, fill = Genus)) +
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

# sample cove vs sample size function
make_cov_plot <- function(df, ylab = "coverage", xlab = "Sample size", x_limit = c(10, 10000)) {
    df_long <- df |>
    pivot_longer(cols = c(min_sample_90, min_sample_95, min_sample_99),
                 names_to = "confidence",
                 values_to = "sample_size" ) |>
    mutate(confidence = factor(confidence,
                               levels = c("min_sample_90", "min_sample_95", "min_sample_99"),
                               labels = c("90%", "95%", "99%"))) |>
      filter(median_mass != 0)
    
    df_90 <- df_long |> filter(confidence == "90%")
    df_95 <- df_long |> filter(confidence == "95%")
    df_99 <- df_long |> filter(confidence == "99%") 

    ggplot(df_long, aes(y = median_mass, colour = Genus)) +
      geom_ribbon(data = df_90,aes(ymin = q2.5, ymax = q97.5, x = sample_size, fill = Genus, group = Genus),
        inherit.aes = FALSE, alpha = 0.10, colour = NA) +
      geom_ribbon(data = df_95, aes(ymin = q2.5, ymax = q97.5, x = sample_size, fill = Genus,  group = Genus),
        inherit.aes = FALSE,alpha = 0.30, colour = NA) +
      geom_ribbon(data = df_99, aes(ymin = q2.5, ymax = q97.5, x = sample_size, fill = Genus,  group = Genus),
        inherit.aes = FALSE, alpha = 0.50, colour = NA) +
      geom_line(data = df_90, aes(x = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_line(data = df_95, aes(x = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_line(data = df_99, aes(x = sample_size, alpha = confidence, group = Genus), linewidth = 0.8) +
      geom_hline(yintercept = 0.80, linetype = "dashed") +
      scale_colour_manual(values = genus_colours, name = "Genus") +
      scale_fill_manual(values = genus_colours, name = "Genus") +
      scale_alpha_manual( values = c("90%" = 0.3, "95%" = 0.55, "99%" = 0.85), name = "Certainty") +
      guides(colour = "none", fill = "none", alpha = guide_legend(title = "Certainty")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_log10(limits = x_limit, 
                    breaks = c(10, 100, 1000, 10000),
                    labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      labs(x = xlab, y = ylab) +
      big_theme +
      theme(legend.position = "bottom")
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# frequency vs cumulative species richness and coverage vs sample size panel plots
p1 <- make_richness_plot(
  richness_df_mlst,
  xlab = "MLST frequency (f)",
  ylab = "Proportion of MLSTs\nof frequency ≥ f",
  xlim_min = min(richness_df_mlst$f[richness_df_mlst$f > 0], na.rm = TRUE)
)
p1

p2 <- make_ss_vs_richness_plot(richness_df_mlst, 
                               x_limit = c(10, 20000),
                               ylab = "Proportion of MLSTs\ndetected")
p2

p3 <- make_richness_plot(
  richness_df_fastbaps_L3,
  xlab = "fastBAPS cluster frequency (f)",
  ylab = "Proportion of fastBAPS clusters\nof frequency ≥ f",
  xlim_min = min(richness_df_fastbaps_L3$f[richness_df_fastbaps_L3$f > 0], na.rm = TRUE)
)
p3

p4 <- make_ss_vs_richness_plot(richness_df_fastbaps_L3, 
                               x_limit = c(10, 20000),
                               ylab = "Proportion of fastBAPS clusters\ndetected")
p4


p5 <- make_richness_plot(
  combined_bsi_pling_richness_df,
  xlab = "Plasmid subcommunity frequency (f)",
  ylab = "Proportion of plasmid subcommunities\nof frequency ≥ f",
  xlim_min = min(combined_bsi_pling_richness_df$f[combined_bsi_pling_richness_df$f > 0], na.rm = TRUE)
)
p5

p6 <- make_ss_vs_richness_plot(combined_bsi_pling_richness_df,
                               x_limit = c(10, 20000),
                               ylab = "Proportion of plasmid subcommunities\ndetected")
p6

p7 <- make_richness_plot(
  combined_bsi_arg_richness_df,
  xlab = "ARG frequency (f)",
  ylab = "Proportion of ARGs\nof frequency ≥ f",
  xlim_min = min(combined_bsi_arg_richness_df$f[combined_bsi_arg_richness_df$f > 0], na.rm = TRUE)
)
p7

p8 <- make_ss_vs_richness_plot(combined_bsi_arg_richness_df, 
                               x_limit = c(10, 20000),
                               ylab = "Proportion of ARGs\ndetected")
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
ggsave("model_results/combined_4x2_species_richness_and_ss_panel_log_y.png",
  plot = combined_figure, width = 11, height = 16.5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# frequency vs cumulative population mass and coverage vs sample size panel plots
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
  xlab = "fastBAPS cluster frequency (f)",
  ylab = "Proportion of bacterial population\n in fastBAPS cluster of frequency ≥ f",
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
  xlab = "ARG frequency (f)",
  ylab = "Proportion of ARG population\nwith allele of frequency ≥ f",
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
ggsave("model_results/combined_4x2_mass_and_coverage_panel_log_y.png",
  plot = combined_figure, width = 11, height = 16.5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. HIERARCHICAL BAYESIAN DIRICHLET MODEL - REGIONAL bayes FOR MLSTS AND fastbaps_L3 CLUSTERS ####
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
                           feature_universe = NULL,
                           region_universe = NULL,
                           alpha_novel_sum = 1,
                           alpha_novel_num = 1,
                           alpha_other = 1,
                           novel_prefix = "NOVEL_") {
  stopifnot(is.data.frame(df))
  stopifnot(all(c(region_col, feature_col, count_col) %in% names(df)))
  stopifnot(is.numeric(alpha_novel_sum), length(alpha_novel_sum) == 1, alpha_novel_sum >= 0)
  stopifnot(is.numeric(alpha_novel_num), length(alpha_novel_num) == 1, alpha_novel_num >= 0)
  
  alpha_novel_num <- as.integer(alpha_novel_num)
  
  df2 <- df |>
    dplyr::group_by(.data[[region_col]], .data[[feature_col]]) |>
    dplyr::summarise(count = sum(.data[[count_col]]), .groups = "drop")
  
  regions <- if (is.null(region_universe)) {
    sort(unique(df2[[region_col]]))
  } else {
    region_universe
  }
  
  # All known features across the data, or supplied universe if you have one
  known_features <- sort(unique(df2[[feature_col]]))
  if (!is.null(feature_universe)) {
    known_features <- sort(unique(c(known_features, feature_universe)))
  }
  
  # Build the novel placeholder labels
  add_novel <- (alpha_novel_sum > 0) && (alpha_novel_num > 0)
  
  if (add_novel) {
    novel_labels <- paste0(novel_prefix, seq_len(alpha_novel_num))
    alpha_novel_each <- alpha_novel_sum / alpha_novel_num
  } else {
    novel_labels <- character(0)
    alpha_novel_each <- 0
  }
  
  # Keep known features first, then append novel placeholders
  features <- c(setdiff(known_features, novel_labels), novel_labels)
  
  y <- matrix(
    0L,
    nrow = length(regions),
    ncol = length(features),
    dimnames = list(regions, features)
  )
  
  rr <- match(df2[[region_col]], regions)
  cc <- match(df2[[feature_col]], features)
  
  # Only fill values that matched a known feature column
  keep <- !is.na(rr) & !is.na(cc)
  y[cbind(rr[keep], cc[keep])] <- as.integer(df2$count[keep])
  
  # Prior vector alpha
  alpha <- rep(alpha_other, length(features))
  names(alpha) <- features
  
  if (length(novel_labels) > 0) {
    alpha[novel_labels] <- alpha_novel_each
  }
  
  novel_idx <- match(novel_labels, features)
  
  list(
    y = y,
    regions = regions,
    features = features,
    novel_labels = novel_labels,
    novel_idx = novel_idx,
    alpha = alpha,
    alpha_novel_sum = alpha_novel_sum,
    alpha_novel_num = alpha_novel_num,
    alpha_novel_each = alpha_novel_each
  )
}

set_alpha_on_prep <- function(prep, alpha_val) {
  prep$alpha <- rep(alpha_val, length(prep$features))
  names(prep$alpha) <- prep$features
  
  prep$alpha_other <- alpha_val
  if (!is.null(prep$alpha_novel_num)) {
    prep$alpha_novel_sum <- alpha_val * prep$alpha_novel_num
    prep$alpha_novel_each <- alpha_val
  }
  
  prep
}

# stratified hold-out split by region 
split_holdout_stratified_by_region <- function(df,
                                               region_col = "region",
                                               feature_col = "escherichia__mlst_achtman__ST",
                                               count_col = "count",
                                               train_prop = 2/3,
                                               seed = 2026) {
  stopifnot(all(c(region_col, feature_col, count_col) %in% names(df)))
  set.seed(seed)
  
  # Expand to isolate-level rows
  long_df <- tidyr::uncount(df, weights = .data[[count_col]], .remove = FALSE)
  
  # Faster than group_modify for many groups
  split_list <- split(long_df, long_df[[region_col]], drop = TRUE)
  
  split_list <- lapply(split_list, function(.x) {
    n <- nrow(.x)
    
    if (n <= 1L) {
      .x$split <- "train"
      return(.x)
    }
    
    n_train <- floor(train_prop * n)
    n_train <- max(1L, min(n - 1L, n_train))
    
    idx <- sample.int(n)
    .x$split <- "test"
    .x$split[idx[seq_len(n_train)]] <- "train"
    .x
  })
  
  split_long <- bind_rows(split_list)
  
  train_df <- split_long |>
    filter(split == "train") |>
    group_by(.data[[region_col]], .data[[feature_col]]) |>
    summarise(count = n(), .groups = "drop")
  
  test_df <- split_long |>
    filter(split == "test") |>
    group_by(.data[[region_col]], .data[[feature_col]]) |>
    summarise(count = n(), .groups = "drop")
  
  list(train = train_df, test = test_df, long = split_long)
}


build_holdout_preps_from_split <- function(split,
                                           df,
                                           region_col = "region",
                                           feature_col = "mlst",
                                           count_col = "count",
                                           alpha_val = 1,
                                           alpha_novel_num = 1,
                                           alpha_other = alpha_val,
                                           novel_prefix = "NOVEL_") {
  regions_universe <- sort(unique(df[[region_col]]))
  features_universe <- sort(unique(df[[feature_col]]))
  
  alpha_novel_sum <- alpha_val * alpha_novel_num
  
  train_prep <- prep_mlst_data(
    split$train,
    region_col = region_col,
    feature_col = feature_col,
    count_col = "count",
    region_universe = regions_universe,
    feature_universe = features_universe,
    alpha_novel_sum = alpha_novel_sum,
    alpha_novel_num = alpha_novel_num,
    alpha_other = alpha_other,
    novel_prefix = novel_prefix
  )
  
  test_prep <- prep_mlst_data(
    split$test,
    region_col = region_col,
    feature_col = feature_col,
    count_col = "count",
    region_universe = regions_universe,
    feature_universe = features_universe,
    alpha_novel_sum = alpha_novel_sum,
    alpha_novel_num = alpha_novel_num,
    alpha_other = alpha_other,
    novel_prefix = novel_prefix
  )
  
  # Set the actual alpha vector used for fitting
  train_prep <- set_alpha_on_prep(train_prep, alpha_val)
  test_prep  <- set_alpha_on_prep(test_prep, alpha_val)
  
  list(train = train_prep, test = test_prep, split = split)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~#
# Stan code: region-specific taus (hierarchical over log tau_r)
stan_region_tau <- '
data {
  int<lower=1> R;
  int<lower=2> K;
  array[R, K] int<lower=0> y;
  vector<lower=0>[K] alpha;
  real<lower=0> mu_tau_mean;
  real<lower=0> mu_tau_sd;
  real<lower=0> sigma_tau_rate;
}
parameters {
  simplex[K] pi;
  real mu_log_tau;
  real<lower=0> sigma_log_tau;
  vector<lower=0>[R] tau_r;
}
model {
  pi ~ dirichlet(alpha);
  
  mu_log_tau ~ normal(log(mu_tau_mean), mu_tau_sd);
  sigma_log_tau ~ exponential(sigma_tau_rate);
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# checking compiled models

# (Optional) make CmdStanR write to a persistent directory instead of temp
options(cmdstanr_output_dir = "cmdstan_output")

# Write Stan code only if needed
write_stan_file_if_needed <- function(code, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  code_lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  
  if (!file.exists(file) || !identical(readLines(file), code_lines)) {
    writeLines(code_lines, con = file)
  }
  normalizePath(file, mustWork = TRUE)
}

# Compile once and reuse in-session
.stan_model_cache <- new.env(parent = emptyenv())

get_compiled_model <- function(model_name,
                               code,
                               stan_file,
                               exe_dir = "stan_executables",
                               force_recompile = FALSE) {
  if (exists(model_name, envir = .stan_model_cache, inherits = FALSE) && !force_recompile) {
    cached <- get(model_name, envir = .stan_model_cache, inherits = FALSE)
    if (!is.null(cached$exe_file) && file.exists(cached$exe_file)) {
      return(cached)
    }
  }
  
  stan_file <- write_stan_file_if_needed(code, stan_file)
  dir.create(exe_dir, recursive = TRUE, showWarnings = FALSE)
  
  mod <- cmdstan_model(stan_file, compile = FALSE)
  mod$compile(dir = exe_dir, force_recompile = force_recompile)
  
  out <- list(
    model = mod,
    stan_file = stan_file,
    exe_file = mod$exe_file()
  )
  
  assign(model_name, out, envir = .stan_model_cache)
  out
}

fit_region_tau_model <- function(prep,
                                 mu_tau_mean = 20, 
                                 mu_tau_sd = 0.25,
                                 sigma_tau_rate = 2,
                                 iter_warmup = 1000,
                                 iter_sampling = 1000,
                                 chains = 4,
                                 seed = 2026,
                                 refresh = 50,
                                 model_cache_dir = "stan_cache",
                                 exe_dir = "stan_executables",
                                 force_recompile = FALSE) {
  compiled <- get_compiled_model(
    model_name = "region_tau_model",
    code = stan_region_tau,
    stan_file = file.path(model_cache_dir, "region_tau_model.stan"),
    exe_dir = exe_dir,
    force_recompile = force_recompile
  )
  
  fit <- compiled$model$sample(
    data = list(
      R = nrow(prep$y),
      K = ncol(prep$y),
      y = prep$y,
      alpha = prep$alpha,
      mu_tau_mean = mu_tau_mean,
      mu_tau_sd = mu_tau_sd,
      sigma_tau_rate = sigma_tau_rate
    ),
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh
  )
  
  list(
    fit = fit,
    prep = prep,
    model = "region_tau",
    mu_tau_mean = mu_tau_mean,
    mu_tau_sd = mu_tau_sd,
    sigm_tau_rate = sigma_tau_rate,
    stan_file = compiled$stan_file,
    exe_file = compiled$exe_file
  )
}

# Save / load helpers
save_fit_bundle <- function(fit_obj, dir, basename) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  fit_file <- file.path(dir, paste0(basename, "_fit.rds"))
  prep_file <- file.path(dir, paste0(basename, "_prep.rds"))
  meta_file <- file.path(dir, paste0(basename, "_meta.rds"))
  #csv_dir <- file.path(dir, paste0(basename, "_cmdstan_files"))
  
  # optional create csv_dir
  #dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  
  # safest save for CmdStanR fits
  fit_obj$fit$save_object(fit_file)
  # optional
  #fit_obj$fit$save_output_files(dir = csv_dir)
  
  saveRDS(fit_obj$prep, prep_file)
  saveRDS(
    fit_obj[setdiff(names(fit_obj), c("fit", "prep"))],
    meta_file
  )
  
  invisible(list(
    fit_file = fit_file,
    prep_file = prep_file,
    meta_file = meta_file #,
    #csv_dir = csv_dir
  ))
}

load_fit_bundle <- function(dir, basename) {
  fit_file <- file.path(dir, paste0(basename, "_fit.rds"))
  prep_file <- file.path(dir, paste0(basename, "_prep.rds"))
  meta_file <- file.path(dir, paste0(basename, "_meta.rds"))
  
  list(
    fit = readRDS(fit_file),
    prep = readRDS(prep_file),
    meta = readRDS(meta_file)
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
extract_pi_draws <- function(draws_df, feature_names) {
  cols <- paste0("pi[", seq_along(feature_names), "]")
  if (!all(cols %in% names(draws_df))) {
    stop("Could not find all pi columns in draws.")
  }
  pi <- as.matrix(draws_df[, cols, drop = FALSE])
  colnames(pi) <- feature_names
  pi
}

extract_tau_draws <- function(draws_df, regions = NULL) {
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# helpers for log held-out log predictive score
log_mean_exp <- function(x) {
  if (length(x) == 0L) return(NA_real_)
  if (all(!is.finite(x))) return(-Inf)
  m <- max(x, na.rm = TRUE)
  m + log(mean(exp(x - m)))
}

dmultinom_logpmf <- function(y, alpha) {
  if (any(!is.finite(alpha)) || any(alpha <= 0)) return(-Inf)
  if (any(!is.finite(y)) || any(y < 0)) return(-Inf)
  
  n <- sum(y)
  if (n == 0) return(0)
  
  lgamma(n + 1) - sum(lgamma(y + 1)) +
    lgamma(sum(alpha)) - lgamma(sum(alpha) + n) +
    sum(lgamma(y + alpha) - lgamma(alpha))
}

extract_param_matrix <- function(draws_df, base_name, n) {
  cols <- paste0(base_name, "[", seq_len(n), "]")
  missing <- setdiff(cols, names(draws_df))
  if (length(missing) > 0) {
    stop("Missing expected draw columns: ", paste(missing, collapse = ", "))
  }
  
  as.matrix(draws_df[, cols, drop = FALSE])
}

score_holdout_regional <- function(fit_obj, test_prep, exclude_novel = TRUE) {
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  
  K <- ncol(test_prep$y)
  R <- nrow(test_prep$y)
  
  pi_mat  <- extract_param_matrix(draws_df, "pi", K)
  tau_mat <- extract_param_matrix(draws_df, "tau_r", R)
  
  keep <- if (exclude_novel && length(test_prep$novel_labels) > 0) {
    !(test_prep$features %in% test_prep$novel_labels)
  } else {
    rep(TRUE, K)
  }
  
  if (!any(keep)) {
    stop("All features removed when excluding novel category — check novel_labels.")
  }
  
  pi_mat <- pi_mat[, keep, drop = FALSE]
  y_kept  <- test_prep$y[, keep, drop = FALSE]
  n_test  <- sum(y_kept, na.rm = TRUE)
  
  S <- nrow(pi_mat)
  lp_draws <- numeric(S)
  
  for (s in seq_len(S)) {
    lp_r <- vapply(seq_len(R), function(r) {
      alpha_sr <- as.numeric(pi_mat[s, ]) * tau_mat[s, r]
      dmultinom_logpmf(y_kept[r, ], alpha_sr)
    }, numeric(1))
    
    lp_draws[s] <- sum(lp_r)
  }
  
  total_score <- log_mean_exp(lp_draws)
  
  tibble(
    holdout_log_score_total = total_score,
    holdout_log_score_per_isolate = if (n_test > 0) total_score / n_test else NA_real_,
    n_test = n_test
  )
}


# function to extract model run diagnostics
extract_diagnostics <- function(fit_obj, iter_sampling = NULL, chains = NULL) {
  
  diag <- fit_obj$fit$diagnostic_summary()
  
  # aggregate chain-level diagnostics 
  num_divergent <- sum(diag$num_divergent)
  num_max_treedepth <- sum(diag$num_max_treedepth)
  ebfmi <- mean(diag$ebfmi, na.rm = TRUE)
  
  # compute fraction of treedepth hits
  if (!is.null(iter_sampling) && !is.null(chains)) {
    treedepth_frac <- num_max_treedepth / (iter_sampling * chains)
  } else {
    treedepth_frac <- NA_real_
  }
  
  # parameter summaries 
  summ <- fit_obj$fit$summary()
  
  rhat <- summ$rhat
  ess  <- summ$ess_bulk
  
  summarise_vec <- function(x, prefix) {
    tibble::tibble(
      !!paste0(prefix, "_min")   := min(x, na.rm = TRUE),
      !!paste0(prefix, "_q2.5")  := quantile(x, 0.025, na.rm = TRUE),
      !!paste0(prefix, "_q25")   := quantile(x, 0.25,  na.rm = TRUE),
      !!paste0(prefix, "_median"):= median(x, na.rm = TRUE),
      !!paste0(prefix, "_q75")   := quantile(x, 0.75,  na.rm = TRUE),
      !!paste0(prefix, "_q97.5") := quantile(x, 0.975, na.rm = TRUE),
      !!paste0(prefix, "_max")   := max(x, na.rm = TRUE)
    )
  }
  
  dplyr::bind_cols(
    tibble::tibble(
      num_divergent = num_divergent,
      num_max_treedepth = num_max_treedepth,
      treedepth_frac = treedepth_frac,
      ebfmi = ebfmi
    ),
    summarise_vec(rhat, "rhat"),
    summarise_vec(ess, "ess_bulk")
  )
}


empty_diagnostics <- function() {
  tibble::tibble(
    num_divergent = NA_real_,
    num_max_treedepth = NA_real_,
    treedepth_frac = NA_real_,
    ebfmi = NA_real_,
    
    rhat_min = NA_real_,
    rhat_q2.5 = NA_real_,
    rhat_q25 = NA_real_,
    rhat_median = NA_real_,
    rhat_q75 = NA_real_,
    rhat_q97.5 = NA_real_,
    rhat_max = NA_real_,
    
    ess_bulk_min = NA_real_,
    ess_bulk_q2.5 = NA_real_,
    ess_bulk_q25 = NA_real_,
    ess_bulk_median = NA_real_,
    ess_bulk_q75 = NA_real_,
    ess_bulk_q97.5 = NA_real_,
    ess_bulk_max = NA_real_
  )
}

score_holdout_regional <- function(fit_obj, test_prep, exclude_novel = TRUE) {
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  K <- ncol(test_prep$y)
  R <- nrow(test_prep$y)
  pi_mat <- extract_param_matrix(draws_df, "pi", K)
  tau_mat <- extract_param_matrix(draws_df, "tau_r", R)
  
  # identify which features to keep
  if (exclude_novel) {
    keep <- test_prep$features != test_prep$novel_label
  } else {
    keep <- rep(TRUE, K)
  }
  
  # filter posterior draws
  pi_mat <- pi_mat[, keep, drop = FALSE]
  S <- nrow(pi_mat)
  lp_draws <- numeric(S)
  
  for (s in seq_len(S)) {
    lp_r <- vapply(seq_len(R), function(r) {
        # filter alpha
      alpha_sr <- as.numeric(pi_mat[s, ]) * tau_mat[s, r]
        # filter observed counts
      y_r <- test_prep$y[r, keep]
        dmultinom_logpmf(y_r, alpha_sr)
      }, numeric(1))
    lp_draws[s] <- sum(lp_r)
  }
  total_score <- log_mean_exp(lp_draws)
  # filter out novel category 
  # ensure we have features left
  if (!any(keep)) {
    stop("All features removed when excluding novel category — check novel_label.")
  }
  # safe sum
  n_test <- sum(test_prep$y[, keep, drop = FALSE], na.rm = TRUE)
  
  tibble(
    holdout_log_score_total = total_score,
    holdout_log_score_per_isolate = if (n_test > 0) total_score / n_test else NA_real_,
    n_test = n_test
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# sensitivity analysis functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
run_one_regional_combo_holdout <- function(train_prep,
                                           test_prep,
                                           alpha_val,
                                           alpha_novel_num,
                                           mu_tau_mean,
                                           mu_tau_sd,
                                           sigma_tau_rate,
                                           iter_warmup = 1000,
                                           iter_sampling = 1000,
                                           chains = 4,
                                           seed = 2026,
                                           refresh = 50,
                                           force_recompile = FALSE) {
  tryCatch({
    train2 <- set_alpha_on_prep(train_prep, alpha_val)
    
    fit_r <- fit_region_tau_model(
      prep = train2,
      mu_tau_mean = mu_tau_mean,
      mu_tau_sd = mu_tau_sd,
      sigma_tau_rate = sigma_tau_rate,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = seed,
      refresh = refresh,
      force_recompile = force_recompile
    )
    
    score_r <- score_holdout_regional(fit_r, test_prep)
    diag_r  <- extract_diagnostics(fit_r, iter_sampling = iter_sampling, chains = chains)
    
    bind_cols(
      tibble(
        alpha_val = alpha_val,
        alpha_novel_num = alpha_novel_num,
        alpha_novel_sum = alpha_val * alpha_novel_num,
        mu_tau_mean = mu_tau_mean,
        mu_tau_sd = mu_tau_sd,
        sigma_tau_rate = sigma_tau_rate,
        status = "ok",
        error = NA_character_,
        holdout_log_score_total = score_r$holdout_log_score_total,
        holdout_log_score_per_isolate = score_r$holdout_log_score_per_isolate,
        n_test = score_r$n_test
      ),
      diag_r
    )
  }, error = function(e) {
    diag_r <- empty_diagnostics()
    
    bind_cols(
      tibble(
        alpha_val = alpha_val,
        alpha_novel_num = alpha_novel_num,
        alpha_novel_sum = alpha_val * alpha_novel_num,
        mu_tau_mean = mu_tau_mean,
        mu_tau_sd = mu_tau_sd,
        sigma_tau_rate = sigma_tau_rate,
        status = "failed",
        error = conditionMessage(e),
        holdout_log_score_total = NA_real_,
        holdout_log_score_per_isolate = NA_real_,
        n_test = sum(test_prep$y)
      ),
      diag_r
    )
  })
}

make_combo_id_regional <- function(alpha_val, alpha_novel_num, mu_tau_mean, mu_tau_sd, sigma_tau_rate) {
  paste0(
    "a", alpha_val,
    "_n", alpha_novel_num,
    "_mu", mu_tau_mean,
    "_sd", mu_tau_sd,
    "_sr", sigma_tau_rate
  )
}


run_sensitivity_grid_regional_resumable <- function(
    df,
    region_col = "region",
    feature_col = "mlst",
    count_col = "count",
    alpha_val_grid,
    alpha_novel_num_grid,
    mu_tau_mean_grid,
    mu_tau_sd_grid,
    sigma_tau_rate_grid,
    train_prop = 2/3,
    seed_split = 2026,
    iter_warmup = 500,
    iter_sampling = 500,
    chains = 1,
    seed_fit = 2026,
    refresh = 100,
    save_dir = "sens_results",
    checkpoint_file = "sens_results/regional_checkpoint.rds",
    force_recompile = FALSE) {
  
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load checkpoint if it exists
  if (file.exists(checkpoint_file)) {
    results <- readRDS(checkpoint_file)
    cat("Loaded checkpoint with", nrow(results), "rows\n")
  } else {
    results <- tibble()
  }
  
  done_ids <- if (nrow(results) > 0 && "combo_id" %in% names(results)) results$combo_id else character(0)
  
  # Full grid
  grid <- tidyr::crossing(
    alpha_val = alpha_val_grid,
    alpha_novel_num = alpha_novel_num_grid,
    mu_tau_mean = mu_tau_mean_grid,
    mu_tau_sd = mu_tau_sd_grid,
    sigma_tau_rate = sigma_tau_rate_grid
  ) |>
    mutate(combo_id = make_combo_id_regional(alpha_val, alpha_novel_num, mu_tau_mean, mu_tau_sd, sigma_tau_rate))
  
  grid <- grid |>
    filter(!combo_id %in% done_ids)
  
  cat("Remaining models:", nrow(grid), "\n")
  
  if (nrow(grid) == 0) return(results)
  
  # Split once
  split <- split_holdout_stratified_by_region(
    df = df,
    region_col = region_col,
    feature_col = feature_col,
    count_col = count_col,
    train_prop = train_prop,
    seed = seed_split
  )
  
  # Speed-up: build preps once per alpha_novel_num
  alpha_groups <- split(grid, grid$alpha_novel_num)
  
  for (alpha_novel_num_chr in names(alpha_groups)) {
    subgrid <- alpha_groups[[alpha_novel_num_chr]]
    alpha_novel_num <- subgrid$alpha_novel_num[[1]]
    
    cat("\nBuilding preps for alpha_novel_num =", alpha_novel_num, "\n")
    
    preps <- build_holdout_preps_from_split(
      split = split,
      df = df,
      region_col = region_col,
      feature_col = feature_col,
      count_col = count_col,
      alpha_val = 1,
      alpha_novel_num = alpha_novel_num
    )
    
    for (i in seq_len(nrow(subgrid))) {
      row <- subgrid[i, ]
      
      cat("Running", i, "of", nrow(subgrid), "|", row$combo_id, "\n")
      
      res <- tryCatch({
        run_one_regional_combo_holdout(
          train_prep = preps$train,
          test_prep = preps$test,
          alpha_val = row$alpha_val,
          alpha_novel_num = row$alpha_novel_num,
          mu_tau_mean = row$mu_tau_mean,
          mu_tau_sd = row$mu_tau_sd,
          sigma_tau_rate = row$sigma_tau_rate,
          iter_warmup = iter_warmup,
          iter_sampling = iter_sampling,
          chains = chains,
          seed = seed_fit,
          refresh = refresh,
          force_recompile = force_recompile
        ) |>
          mutate(combo_id = row$combo_id)
      }, error = function(e) {
        tibble(
          alpha_val = row$alpha_val,
          alpha_novel_num = row$alpha_novel_num,
          alpha_novel_sum = row$alpha_val * row$alpha_novel_num,
          mu_tau_mean = row$mu_tau_mean,
          mu_tau_sd = row$mu_tau_sd,
          sigma_tau_rate = row$sigma_tau_rate,
          combo_id = row$combo_id,
          status = "failed",
          error = conditionMessage(e),
          holdout_log_score_total = NA_real_,
          holdout_log_score_per_isolate = NA_real_
        )
      })
      
      results <- bind_rows(results, res)
      
      saveRDS(results, checkpoint_file)
      readr::write_csv(results, sub("\\.rds$", ".csv", checkpoint_file))
    }
  }
  
  results
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
flag_valid_fit <- function(df, iterations = 250) {
  df |>
    mutate(num_divergent_frac = num_divergent / iterations) |>
    mutate(
      valid_fit = status == "ok" &
        is.finite(num_divergent) &
        is.finite(treedepth_frac) &
        is.finite(rhat_max) &
        is.finite(ess_bulk_min) &
        num_divergent_frac <0.05 &
        treedepth_frac < 0.05 &
        rhat_q75 < 1.05 &
        ess_bulk_q25 > 100 
    )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Heatmap plots
plot_regional_holdout_heatmap <- function(regional_results,
                                          value_col = "holdout_log_score_per_isolate") {
  
  stopifnot(value_col %in% names(regional_results))
  stopifnot(all(c("alpha_val", "alpha_novel_num", "alpha_novel_sum",
                  "mu_tau_mean", "mu_tau_sd", "sigma_tau_rate",
                  "valid_fit") %in% names(regional_results)))
  
  sigma_levels <- sort(unique(regional_results$sigma_tau_rate))
  
  fill_limits <- range(regional_results[[value_col]], na.rm = TRUE)
  fill_breaks <- pretty(fill_limits, n = 5)
  
  make_one_plot <- function(sig_val) {
    df <- regional_results %>%
      filter(sigma_tau_rate == sig_val) %>%
      mutate(
        valid_fit = dplyr::coalesce(valid_fit, FALSE),
        fill_val = .data[[value_col]],
        tile_alpha = if_else(valid_fit, 1, 0.3),
        alpha_val = factor(alpha_val, levels = sort(unique(alpha_val))),
        alpha_novel_num = factor(alpha_novel_num, levels = sort(unique(alpha_novel_num))),
        mu_tau_mean = factor(
          paste0("mu[mu*log(tau)]==", mu_tau_mean),
          levels = paste0("mu[mu*log(tau)]==", sort(unique(regional_results$mu_tau_mean)))
        ),
        mu_tau_sd = factor(
          paste0("sigma[mu*log(tau)]==", mu_tau_sd),
          levels = paste0("sigma[mu*log(tau)]==", sort(unique(regional_results$mu_tau_sd)))
        )
      )
    
    ggplot(df, aes(x = alpha_val, y = alpha_novel_num, fill = fill_val, alpha = tile_alpha)) +
      geom_tile(width = 0.95, height = 0.95, colour = "white", linewidth = 0.25) +
      facet_grid(
        rows = vars(mu_tau_sd),
        cols = vars(mu_tau_mean),
        labeller = labeller(
          mu_tau_mean = label_parsed,
          mu_tau_sd = label_parsed
        )
      ) +
      scale_fill_viridis_c(
        name = "Hold-out score",
        limits = fill_limits,
        breaks = fill_breaks,
        oob = scales::squish,
        na.value = "grey90"
      ) +
      scale_alpha_identity(guide = "none") +
      labs(
        x = expression(alpha[1:k]),
        y = expression(m),
        title = bquote(lambda[sigma*log(tau)] == .(sig_val))
      ) +
      theme_minimal(base_size = 11) +
      theme(
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom"
      )
  }
  
  plots <- lapply(sigma_levels, make_one_plot)
  
  if (length(plots) == 1L) {
    return(plots[[1]])
  }
  
  wrap_plots(plots, nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

subset_prep_rows <- function(prep, idx) {
  out <- prep
  out$y <- prep$y[idx, , drop = FALSE]
  out$regions <- prep$regions[idx]
  out
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Compute richness and mass curves 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fast helper: summarize each column of a draw-by-threshold matrix
summarise_draw_matrix_ms <- function(mat) {
  q <- matrixStats::colQuantiles(
    mat,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    drop = TRUE
  )
  
  tibble::tibble(
    region_median = matrixStats::colMedians(mat, na.rm = TRUE),
    region_q2.5   = as.numeric(q[1, ]),
    region_q97.5  = as.numeric(q[2, ])
  )
}

# Fast helper for threshold tail counts / sums on descending-sorted vectors
tail_count_desc <- function(x_sorted_desc, f_grid) {
  findInterval(-f_grid, -x_sorted_desc)
}

tail_sum_desc <- function(x_sorted_desc, f_grid) {
  idx <- findInterval(-f_grid, -x_sorted_desc)
  cs <- cumsum(x_sorted_desc)
  out <- numeric(length(f_grid))
  keep <- idx > 0L
  out[keep] <- cs[idx[keep]]
  out
}

# Core worker for one region
region_curve_worker <- function(
    region_name,
    y_r,
    pi_draws,
    tau_vec,
    f_grid,
    total_k,
    mode = c("richness", "mass")) {
  
  mode <- match.arg(mode)
  
  nd <- nrow(pi_draws)
  nf <- length(f_grid)
  y_r_non_zero <- sum(y_r > 0)
  
  region_mat     <- matrix(NA_real_, nrow = nd, ncol = nf)
  global_mat     <- matrix(NA_real_, nrow = nd, ncol = nf)
  global_all_mat <- matrix(NA_real_, nrow = nd, ncol = nf)
  
  for (i in seq_len(nd)) {
    pi_b <- as.numeric(pi_draws[i, ])
    tau_b <- tau_vec[i]
    
    alpha_post <- y_r + tau_b * pi_b
    
    # Faster than gtools::rdirichlet(1, alpha_post)
    theta_b <- rgamma(length(alpha_post), shape = alpha_post, rate = 1)
    theta_b <- theta_b / sum(theta_b)
    
    ord_theta <- order(theta_b, decreasing = TRUE)
    theta_sorted <- theta_b[ord_theta]
    pi_in_theta_order <- pi_b[ord_theta]
    
    k_theta <- tail_count_desc(theta_sorted, f_grid)
    
    if (mode == "richness") {
      # Count of features above each threshold
      # k_theta[j] = number of MLSTs with theta >= f_grid[j]
      k_theta <- tail_count_desc(theta_sorted, f_grid)
      
      # Richness as a proportion of total features
      region_mat[i, ] <- k_theta / total_k
      
      # Same selected feature set, also as proportion of total features
      global_mat[i, ] <- k_theta / total_k
      
      # Global richness of all MLSTs whose global freq is >= f
      ord_pi <- order(pi_b, decreasing = TRUE)
      pi_sorted <- pi_b[ord_pi]
      k_pi <- tail_count_desc(pi_sorted, f_grid)
      global_all_mat[i, ] <- k_pi / total_k
    } else {
      cs_theta <- cumsum(theta_sorted)
      tmp_region <- numeric(nf)
      keep <- k_theta > 0L
      tmp_region[keep] <- cs_theta[k_theta[keep]]
      region_mat[i, ] <- tmp_region
      
      cs_pi_theta <- cumsum(pi_in_theta_order)
      tmp_global <- numeric(nf)
      tmp_global[keep] <- cs_pi_theta[k_theta[keep]]
      global_mat[i, ] <- tmp_global
      
      ord_pi <- order(pi_b, decreasing = TRUE)
      pi_sorted <- pi_b[ord_pi]
      k_pi <- tail_count_desc(pi_sorted, f_grid)
      
      cs_pi <- cumsum(pi_sorted)
      tmp_all <- numeric(nf)
      keep_pi <- k_pi > 0L
      tmp_all[keep_pi] <- cs_pi[k_pi[keep_pi]]
      global_all_mat[i, ] <- tmp_all
    }
  }
  
  region_stats <- tibble::tibble(
    region_median = matrixStats::colMedians(region_mat, na.rm = TRUE),
    region_q2.5   = as.numeric(matrixStats::colQuantiles(region_mat, probs = 0.025, na.rm = TRUE, drop = TRUE)),
    region_q97.5  = as.numeric(matrixStats::colQuantiles(region_mat, probs = 0.975, na.rm = TRUE, drop = TRUE))
  )
  
  global_stats <- tibble::tibble(
    global_median = matrixStats::colMedians(global_mat, na.rm = TRUE),
    global_q2.5   = as.numeric(matrixStats::colQuantiles(global_mat, probs = 0.025, na.rm = TRUE, drop = TRUE)),
    global_q97.5  = as.numeric(matrixStats::colQuantiles(global_mat, probs = 0.975, na.rm = TRUE, drop = TRUE))
  )
  
  global_all_stats <- tibble::tibble(
    global_all_median = matrixStats::colMedians(global_all_mat, na.rm = TRUE),
    global_all_q2.5   = as.numeric(matrixStats::colQuantiles(global_all_mat, probs = 0.025, na.rm = TRUE, drop = TRUE)),
    global_all_q97.5  = as.numeric(matrixStats::colQuantiles(global_all_mat, probs = 0.975, na.rm = TRUE, drop = TRUE))
  )
  
  tibble::tibble(
    region = region_name,
    f = f_grid
  ) |>
    dplyr::bind_cols(region_stats, global_stats, global_all_stats) |>
    dplyr::mutate(
      min_sample_90 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.90) / log(1 - f)), NA_integer_),
      min_sample_95 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.95) / log(1 - f)), NA_integer_),
      min_sample_99 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.99) / log(1 - f)), NA_integer_)
    )
}

compute_region_global_richness_curves_ms <- function(
    fit_obj,
    f_grid,
    n_post_draws = 1000,
    seed = 2026,
    cores = max(1L, parallel::detectCores(logical = FALSE) - 1L)) {
  
  set.seed(seed)
  
  prep <- fit_obj$prep
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  feature_names <- prep$features
  
  pi_draws <- extract_pi_draws(draws_df, feature_names)
  B <- nrow(pi_draws)
  draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
  
  pi_draws <- pi_draws[draw_idx, , drop = FALSE]
  
  has_shared_tau <- "tau" %in% names(draws_df)
  has_region_tau <- any(grepl("^tau_r\\[", names(draws_df)))
  
  tau_shared <- if (has_shared_tau) draws_df$tau[draw_idx] else NULL
  
  tau_r_mat <- NULL
  if (has_region_tau) {
    tau_r_cols <- grep("^tau_r\\[", names(draws_df), value = TRUE)
    tau_r_mat <- as.matrix(draws_df[draw_idx, tau_r_cols, drop = FALSE])
    colnames(tau_r_mat) <- prep$regions
  }
  
  if (!has_shared_tau && !has_region_tau) {
    stop("No tau or tau_r found in fit.")
  }
  
  region_names <- prep$regions
  
  worker_args <- lapply(seq_along(region_names), function(r) {
    region_name <- region_names[r]
    y_r <- as.numeric(prep$y[r, ])
    
    tau_vec <- if (has_shared_tau) {
      tau_shared
    } else {
      tau_r_mat[, region_name]
    }
    
    list(
      region_name = region_name,
      y_r = y_r,
      pi_draws = pi_draws,
      tau_vec = tau_vec,
      f_grid = f_grid,
      total_k = length(feature_names),
      mode = "richness"
    )
  })
  
  region_results <- if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(worker_args, function(a) {
      do.call(region_curve_worker, a)
    }, mc.cores = cores)
  } else {
    lapply(worker_args, function(a) {
      do.call(region_curve_worker, a)
    })
  }
  
  dplyr::bind_rows(region_results)
}

compute_region_global_mass_curves_ms <- function(
    fit_obj,
    f_grid,
    n_post_draws = 400,
    seed = 2026,
    cores = max(1L, parallel::detectCores(logical = FALSE) - 1L)) {
  
  set.seed(seed)
  
  prep <- fit_obj$prep
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  feature_names <- prep$features
  
  pi_draws <- extract_pi_draws(draws_df, feature_names)
  B <- nrow(pi_draws)
  draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
  
  pi_draws <- pi_draws[draw_idx, , drop = FALSE]
  
  has_shared_tau <- "tau" %in% names(draws_df)
  has_region_tau <- any(grepl("^tau_r\\[", names(draws_df)))
  
  tau_shared <- if (has_shared_tau) draws_df$tau[draw_idx] else NULL
  
  tau_r_mat <- NULL
  if (has_region_tau) {
    tau_r_cols <- grep("^tau_r\\[", names(draws_df), value = TRUE)
    tau_r_mat <- as.matrix(draws_df[draw_idx, tau_r_cols, drop = FALSE])
    colnames(tau_r_mat) <- prep$regions
  }
  
  if (!has_shared_tau && !has_region_tau) {
    stop("No tau or tau_r found in fit.")
  }
  
  region_names <- prep$regions
  
  worker_args <- lapply(seq_along(region_names), function(r) {
    region_name <- region_names[r]
    y_r <- as.numeric(prep$y[r, ])
    
    tau_vec <- if (has_shared_tau) {
      tau_shared
    } else {
      tau_r_mat[, region_name]
    }
    
    list(
      region_name = region_name,
      y_r = y_r,
      pi_draws = pi_draws,
      tau_vec = tau_vec,
      f_grid = f_grid,
      total_k = length(feature_names),
      mode = "mass"
    )
  })
  
  region_results <- if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(worker_args, function(a) {
      do.call(region_curve_worker, a)
    }, mc.cores = cores)
  } else {
    lapply(worker_args, function(a) {
      do.call(region_curve_worker, a)
    })
  }
  
  dplyr::bind_rows(region_results)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plotting functions
plot_global_freq_vs_pi_r <- function(
    fit_r,
    fit_r_label = "fit_r",
    n_post_draws = 1000,
    seed = 2026,
    conf_level = 0.95,
    pseudo_count = NULL,
    point_size = 2.2,
    point_alpha = 0.3,
    line_size = 0.8) {
  
  set.seed(seed)
  
  summarise_one_fit <- function(fit_obj, model_label) {
    prep <- fit_obj$prep
    features <- prep$features
    
    # safer novel filtering
    keep <- !features %in% prep$novel_label
    
    draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
    pi_draws <- as.matrix(extract_pi_draws(draws_df, features))
    
    if (!is.null(colnames(pi_draws)) && all(features %in% colnames(pi_draws))) {
      pi_draws <- pi_draws[, features, drop = FALSE]
    }
    
    # filter out novel category for plotting
    pi_draws <- pi_draws[, keep, drop = FALSE]
    features <- features[keep]
    
    B <- nrow(pi_draws)
    draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
    
    alpha <- (1 - conf_level) / 2
    pi_sub <- pi_draws[draw_idx, , drop = FALSE]
    
    pi_median <- matrixStats::colMedians(pi_sub, na.rm = TRUE)
    pi_lo <- matrixStats::colQuantiles(pi_sub, probs = alpha, na.rm = TRUE, drop = TRUE)
    pi_hi <- matrixStats::colQuantiles(pi_sub, probs = 1 - alpha, na.rm = TRUE, drop = TRUE)
    
    gf <- colSums(prep$y) / sum(prep$y)
    gf <- gf[keep]
    
    tibble::tibble(
      feature = features,
      model = model_label,
      actual_freq = as.numeric(gf),
      estimate_median = as.numeric(pi_median),
      estimate_lo = as.numeric(pi_lo),
      estimate_hi = as.numeric(pi_hi)
    )
  }
  
  if (is.null(pseudo_count)) {
    y <- fit_r$prep$y
    pseudo_count <- 0.5 / max(sum(y), 1)
  }
  
  df_r <- summarise_one_fit(fit_r, fit_r_label)
  
  plot_df <- df_r |>
    dplyr::mutate(
      actual_plot = pmax(actual_freq, pseudo_count),
      estimate_plot = pmax(estimate_median, pseudo_count),
      lo_plot = pmax(estimate_lo, pseudo_count),
      hi_plot = pmax(estimate_hi, pseudo_count)
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = actual_plot, y = estimate_plot)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", linewidth = line_size, color = "grey50"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo_plot, ymax = hi_plot),
      width = 0, alpha = point_alpha
    ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Actual global frequency",
      y = "Estimated pi (posterior median)",
      title = "Global frequency vs estimated prevalence"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


plot_regional_freq_vs_theta_r <- function(
    fit_r,
    fit_r_label = "fit_r",
    n_post_draws = 1000,
    seed = 2026,
    conf_level = 0.95,
    pseudo_count = NULL,
    point_size = 1.8,
    point_alpha = 0.3,
    line_size = 0.7,
    ncol = NULL) {
  
  summarise_one_fit <- function(fit_obj, model_label) {
    prep <- fit_obj$prep
    features <- prep$features
    regions <- prep$regions
    R <- length(regions)
    K <- length(features)
    
    draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
    pi_draws <- as.matrix(extract_pi_draws(draws_df, features))
    tau_info <- extract_tau_draws(draws_df, regions = regions)
    
    if (!is.null(colnames(pi_draws)) && all(features %in% colnames(pi_draws))) {
      pi_draws <- pi_draws[, features, drop = FALSE]
    }
    
    B <- nrow(pi_draws)
    draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
    
    alpha <- (1 - conf_level) / 2
    
    # actual regional frequencies: y / n_iso for each region
    actual_mat <- prep$y / rowSums(prep$y)
    colnames(actual_mat) <- features
    rownames(actual_mat) <- regions
    
    theta_arr <- array(NA_real_, dim = c(length(draw_idx), R, K))
    
    for (i in seq_along(draw_idx)) {
      b <- draw_idx[i]
      pi_b <- as.numeric(pi_draws[b, ])
      
      for (r in seq_along(regions)) {
        region_name <- regions[r]
        y_r <- as.numeric(prep$y[r, ])
        n_r <- as.integer(rowSums(prep$y)[region_name])
        
        tau_b <- if (!is.null(tau_info$shared)) {
          tau_info$shared[b]
        } else {
          tau_info$region[b, region_name]
        }
        
        theta_b <- rbeta(
          n = length(pi_b),
          shape1 = y_r + tau_b * pi_b,
          shape2 = (n_r - y_r) + tau_b * (1 - pi_b)
        )
        
        theta_arr[i, r, ] <- theta_b
      }
    }
    
    out <- vector("list", R)
    
    for (r in seq_along(regions)) {
      est_draws <- theta_arr[, r, , drop = FALSE][, 1, ]
      
      est_median <- matrixStats::colMedians(est_draws, na.rm = TRUE)
      est_lo <- matrixStats::colQuantiles(est_draws, probs = alpha, na.rm = TRUE, drop = TRUE)
      est_hi <- matrixStats::colQuantiles(est_draws, probs = 1 - alpha, na.rm = TRUE, drop = TRUE)
      
      out[[r]] <- tibble::tibble(
        region = regions[r],
        feature = features,
        model = model_label,
        actual_freq = as.numeric(actual_mat[r, ]),
        estimate_median = as.numeric(est_median),
        estimate_lo = as.numeric(est_lo),
        estimate_hi = as.numeric(est_hi)
      )
    }
    
    dplyr::bind_rows(out)
    
  }
  
  if (is.null(pseudo_count)) {
    pseudo_count <- 0.5 / max(rowSums(fit_r$prep$y))
  }
  
  df_r <- summarise_one_fit(fit_r, fit_r_label)
  
  plot_df <- df_r |>
    dplyr::mutate(
      actual_plot = pmax(actual_freq, pseudo_count),
      estimate_plot = pmax(estimate_median, pseudo_count),
      lo_plot = pmax(estimate_lo, pseudo_count),
      hi_plot = pmax(estimate_hi, pseudo_count)
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = actual_plot, y = estimate_plot)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", linewidth = line_size, color = "grey50"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo_plot, ymax = hi_plot),
      width = 0, alpha = point_alpha
    ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ region, ncol = ncol) +
    ggplot2::labs(
      x = "Actual regional frequency",
      y = "Estimated theta (posterior median)",
      title = "Regional frequencies vs estimated prevalence"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


# plot functions to visulaise resuls
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


# plot regional mass curve (- by global mass)sample plot functions for species richness and mass)
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


make_cov_plot_region <- function(df, ylab = "coverage", xlab = "Sample size", x_limit = c(10, 10000), pal = region_pal) {
  ggplot(df |> dplyr::filter(f > 0), aes(y = global_median, colour = region, fill = region)) +
    geom_ribbon(aes(x = min_sample_95, ymin = global_q2.5, ymax = global_q97.5), alpha = 0.1, colour = NA) +
    geom_line(aes(x = min_sample_95), linewidth = 1.2) +
    geom_hline(yintercept = 0.80, linetype = "dashed") +
    scale_fill_manual(values = region_pal) +
    scale_colour_manual(values = region_pal) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_log10(limits = x_limit,
                  breaks = c(10, 100, 1000, 10000),
                  labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    labs(x = xlab, y = ylab) +
    big_theme
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
length(unique(ecoli_bsi_count_region$escherichia__mlst_achtman__ST))
#View(ecoli_bsi_count_region)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.1ai Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~#
# regional tau
ecoli_mlst_regional_results <- run_sensitivity_grid_regional_resumable(
  df = ecoli_bsi_count_region,
  region_col = "region",
  feature_col = "escherichia__mlst_achtman__ST",
  count_col = "count",
  alpha_val_grid = c(0.1, 1, 10),
  alpha_novel_num_grid = c(0, 1, 10, 100, 500, 1000), 
  mu_tau_mean_grid =  c(0.1, 1, 10),
  mu_tau_sd_grid =  c(0.1, 1, 10),
  sigma_tau_rate_grid =  c(0.1, 1, 10),
  train_prop = 2/3,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  refresh = 100,#
  save_dir = "sens_results",
  checkpoint_file = "sens_results/ecoli_mlst_regional_checkpoint_with_diagnostics_novel.rds"
)

# Flag valid model runs
ecoli_mlst_regional_results <- flag_valid_fit(ecoli_mlst_regional_results)


# save results
write.csv(ecoli_mlst_regional_results, "sens_results/ecoli_mlst_regional_tau_sensitivity_results_fixed_holdout_lps.csv", row.names = FALSE)
#ecoli_mlst_regional_results <- read.csv("sens_results/ecoli_mlst_regional_tau_sensitivity_results_fixed_holdout_lps.csv")
#View(ecoli_mlst_regional_results)

# truncate extreme low-fit/high-fit values
upper <- quantile(ecoli_mlst_regional_results$holdout_log_score_per_isolate, 0.975, na.rm = TRUE)
upper # -1.595
ecoli_mlst_regional_results <- ecoli_mlst_regional_results |>
  mutate(fill_val = pmin(holdout_log_score_per_isolate, upper))

# # filter if needed
ecoli_mlst_regional_results_neat <- ecoli_mlst_regional_results |>
  filter(alpha_val %in% c(0.1, 1, 10),
         alpha_novel_num %in% c(1, 10, 100, 500, 1000), 
         mu_tau_mean %in%  c(0.1, 1, 10),
         mu_tau_sd %in%  c(0.1, 1, 10),
         sigma_tau_rate %in%  c(0.1, 1, 10))


#plot regional heatmap
ecoli_mlst_regional_tau_plot <- plot_regional_holdout_heatmap(ecoli_mlst_regional_results_neat, value_col = "fill_val")
ecoli_mlst_regional_tau_plot
ggsave("sens_results/ecoli_mlst_regional_tau_sensitivity_analysis_plot.png", ecoli_mlst_regional_tau_plot, units = "in", width = 10, height = 4.2, dpi = 300 )
# best fit is a1_mu0.01_sd0.01_sr10
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.1a.ii Run with single prior parameter combo ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set priors
df <- ecoli_bsi_mlst_df
counts <- df$count
N <- sum(df$count) # 1471
K <- length(df$count) # 263  # num unique features
f1 <- sum(df$count == 1)  # 161       # number of singletons
q_hat <- f1 / N  # 11% Good-Turing first-order - proportion of singletons

# Estimate novel/ unseen mass using Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(counts)
crp_fit$summary_df # mean: theta = 62.5 -> 4% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.0407

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$kleborate_mlst
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #12


# Run model single prior parameter combo
ecoli_mlst_prep <- prep_mlst_data(ecoli_bsi_count_region,
                        region_col = "region",
                        feature_col = "escherichia__mlst_achtman__ST",
                        count_col = "count",
                        alpha_novel_sum = alpha_novel_sum,
                        alpha_novel_num = alpha_novel_sum,
                        alpha_other = 1)
# fit model
fit_r <- fit_region_tau_model(ecoli_mlst_prep, 
                              mu_tau_mean = 1,
                              mu_tau_sd = 1 ,
                              sigma_tau_rate = 1,
                              iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# save
save_fit_bundle(fit_r, dir = "model_results", basename = "ecoli_bsi_mlst_region_tau")
#~~~~~~~~~~~~~~~~#
# read back in
#obj_r <- load_fit_bundle("model_results", "ecoli_bsi_mlst_region_tau")
#fit_r <- list(fit = obj_r$fit, prep = obj_r$prep, model = obj_r$meta$model,
#              stan_file = obj_r$meta$stan_file, exe_file = obj_r$meta$exe_file)
#~~~~~~~~~~~~~~~~#

# model summaries and diagnostics
print(fit_r$fit$summary(), n=300)
mean(fit_r$fit$summary()$rhat)
range(fit_r$fit$summary()$rhat)

mean(fit_r$fit$summary()$ess_bulk)
range(fit_r$fit$summary()$ess_bulk)

mean(fit_r$fit$summary()$ess_tail)
range(fit_r$fit$summary()$ess_tail)

fit_r$fit$cmdstan_summary()
fit_r$fit$cmdstan_diagnose()
mean(fit_r$fit$diagnostic_summary()$ebfmi)
fit_r$fit$diagnostic_summary()
fit_r$fit$sampler_diagnostics()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Cumulative curves with global and regional richness/ masses  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative species richness as proportion
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

curve_df_richness <- compute_region_global_richness_curves_ms(
  fit_obj = fit_r,    
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_df_richness)
#colnames(curve_df_richness)

# save
write.csv(curve_df_richness, "model_results/ecoli_bsi_mlst_richness_curves.csv", row.names = FALSE)
#curve_df_richness <- read.csv("model_results/ecoli_bsi_mlst_richness_curves.csv")

#~~~~~~~~~~~~~~~~~#
# cumulative mass
curve_df <- compute_region_global_mass_curves_ms(
  fit_obj = fit_r,     
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_df)
#colnames(curve_df)

# save
write.csv(curve_df, "model_results/ecoli_bsi_mlst_mass_curves.csv", row.names = FALSE)
#curve_df <- read.csv("model_results/ecoli_bsi_mlst_mass_curves.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * * 4.1.aiii Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# colourBrewerSet2 + Dark2 themed
region_pal <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
)
region_pal <- rev(region_pal)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#plots to assess fit
pi_plot <- plot_global_freq_vs_pi_r(fit_r,
                                    fit_r_label = "Regional tau")
pi_plot 
# save
ggsave("model_results/ecoli_bsi_mlst_overall_fit_hierarchical_model.png", 
       pi_plot, units = "in", width = 6, height = 4.5, dpi = 300)


theta_plot <- plot_regional_freq_vs_theta_r(fit_r,
                                            fit_r_label = "Regional tau")
theta_plot
# save 
ggsave("model_results/ecoli_bsi_mlst_regional_fit_hierarchical_model.png", 
       theta_plot, units = "in", width = 10, height = 7, dpi = 300)


#~~~~~~~~~~~~~~~~~~~#
# make individual plots:
# species richness
p1 <- make_mass_plot_region(
  curve_df_richness,
  xlab = "Regional MLST frequency (f)",
  ylab = "Proportion of unique MLSTs\n with regional frequency ≥f",
  xlim_min = min(curve_df_richness$f[curve_df_richness$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("model_results/ecoli_bsi_mlst_regional_cumulative_richness_curve_regionl_tau_bayes.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_df_richness, 
                           pal = region_pal,
                           x_limit = c(10, 10000),
                           ylab = "Proportion of unique\nMLSTs detected")
p2
# save
ggsave("model_results/ecoli_bsi_mlst_regional_richness_vs_ss_regionl_tau_bayes.png", p2, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~#
#coverage
# make individual plots:
p3 <- make_mass_plot_region(
  curve_df,
  xlab = "MLST frequency (f)",
  ylab = "Proportion of bacterial population\nwith MLST of frequency ≥f",
  xlim_min = min(curve_df$f[curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p3
ggsave("model_results/ecoli_bsi_mlst_regional_cumulative_mass_curve_regionl_tau_bayes.png", p3, units = "in", width = 6, height = 4, dpi = 300)


p4 <- make_cov_plot_region(curve_df, pal = region_pal)
p4
# save
ggsave("model_results/ecoli_bsi_mlst_regional_samplecoverage_vs_ss_regionl_tau_bayes.png", p4, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * * 4.1.aiv Summary tables ####
# for E. coli MLSTs by region
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

# add global cumulative richness and ss rows
national_df <- curve_df_richness |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_df_richness, national_df)

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
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_mlst_regional_bhm_summary_richness <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(ecoli_bsi_mlst_regional_bhm_summary_richness)

# tidy table
ecoli_bsi_mlst_regional_bhm_summary_richness <- ecoli_bsi_mlst_regional_bhm_summary_richness |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_mlst_regional_bhm_summary_richness)
#View(ecoli_bsi_mlst_regional_bhm_summary_richness)

# save
write.csv(ecoli_bsi_mlst_regional_bhm_summary_richness, "model_results/ecoli_bsi_mlst_regional_bhm_summary_richness.csv", row.names = FALSE)
#ecoli_bsi_mlst_regional_bhm_summary_richness <- read.csv("model_results/ecoli_bsi_mlst_regional_bhm_summary_richness.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
ecoli_bsi_mlst_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
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
#View(ecoli_bsi_mlst_regional_bhm_summary)

# save
write.csv(ecoli_bsi_mlst_regional_bhm_summary, "model_results/ecoli_bsi_mlst_regional_bhm_summary.csv", row.names = FALSE)
#ecoli_bsi_mlst_regional_bhm_summary <- read.csv("model_results/ecoli_bsi_mlst_regional_bhm_summary.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 4.1b Klebsiella MLSTs - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
kleb_bsi_count_region <- kleb_bsi_samples_metadata |>
  group_by(region, klebsiella_mlst_ST) |>
  summarise(count = n())
length(unique(kleb_bsi_count_region$klebsiella_mlst_ST))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * *  * 4.1bi Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Regional tau
kleb_mlst_regional_results <- run_sensitivity_grid_regional_resumable(
  df = kleb_bsi_count_region,
  region_col = "region",
  feature_col = "klebsiella_mlst_ST",
  count_col = "count",
  alpha_val_grid = c(0.1, 1, 10),
  alpha_novel_num_grid = c(0, 1, 10, 100, 500, 1000), 
  mu_tau_mean_grid =  c(0.1, 1, 10),
  mu_tau_sd_grid =  c(0.1, 1, 10),
  sigma_tau_rate_grid =  c(0.1, 1, 10),
  train_prop = 2/3,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  refresh = 100,
  save_dir = "sens_results",
  checkpoint_file = "sens_results/kleb_mlst_regional_checkpoint_with_diagnostics_novel.rds"
)

# flag valid fits
kleb_mlst_regional_results <- flag_valid_fit(kleb_mlst_regional_results)

# save results
write.csv(kleb_mlst_regional_results, "sens_results/kleb_mlst_regional_tau_sensitivity_results_fixed_holdout_lps.csv", row.names = FALSE)
#kleb_mlst_regional_results <- read.csv("sens_results/kleb_mlst_regional_tau_sensitivity_results_fixed_holdout_lps.csv")
#View(kleb_mlst_regional_results)

# truncate extreme low-fit/high-fit values
upper <- quantile(kleb_mlst_regional_results$holdout_log_score_per_isolate, 0.975, na.rm = TRUE)
upper # -3.068
kleb_mlst_regional_results <- kleb_mlst_regional_results |>
  mutate(fill_val = pmin(holdout_log_score_per_isolate, upper))

# filter if needed
kleb_mlst_regional_results_neat <- kleb_mlst_regional_results |>
  filter(alpha_val %in% c(0.1, 1, 10),
         alpha_novel_num %in% c(1, 10, 100, 500, 1000), 
         mu_tau_mean %in%  c(0.1, 1, 10),
         mu_tau_sd %in%  c(0.1, 1, 10),
         sigma_tau_rate %in%  c(0.1, 1, 10))

#plot regional heatmap
kleb_mlst_regional_tau_plot <- plot_regional_holdout_heatmap(kleb_mlst_regional_results_neat, value_col = "fill_val")
kleb_mlst_regional_tau_plot
ggsave("sens_results/kleb_mlst_regional_tau_sensitivity_analysis_plot.png", kleb_mlst_regional_tau_plot, units = "in", width = 10, height = 4.2, dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * *  * 4.1bii  Run with single prior parameter combo ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set priors
df <- kleb_bsi_mlst_df
counts <- df$count
N <- sum(df$count)
K <- nrow(df)
f1 <- sum(df$count == 1)  # 223       # number of singletons
q_hat <- f1 / N  # 47% Good-Turing first-order - proportion of singletons

# fit data to estimate Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(counts)
crp_fit$summary_df # mean: theta = 111.2 -> 19% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.192

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$kleborate_mlst
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #71

# run main models
kleb_mlst_prep <- prep_mlst_data(kleb_bsi_count_region,
                       region_col = "region",
                       feature_col = "klebsiella_mlst_ST",
                       count_col = "count",
                       alpha_novel_sum = alpha_novel_sum,
                       alpha_novel_num = alpha_novel_sum,
                       alpha_other = 1)
#View(prep)
#str(prep) # list of 6

# fit model
fit_r <- fit_region_tau_model(kleb_mlst_prep, 
                              mu_tau_mean = 1,
                              mu_tau_sd = 1,
                              sigma_tau_rate = 1,
                              iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# save
save_fit_bundle(fit_r, dir = "model_results", basename = "kleb_bsi_mlst_region_tau_1")
#~~~~~~~~~~~~~~~~#
# read back in
#obj_u <- load_fit_bundle("model_results", "kleb_bsi_mlst_shared_tau")
#obj_r <- load_fit_bundle("model_results", "kleb_bsi_mlst_region_tau")
#fit_u <- list(fit = obj_u$fit, prep = obj_u$prep, model = obj_u$meta$model,
#              stan_file = obj_u$meta$stan_file, exe_file = obj_u$meta$exe_file)
#fit_r <- list(fit = obj_r$fit, prep = obj_r$prep, model = obj_r$meta$model,
#              stan_file = obj_r$meta$stan_file, exe_file = obj_r$meta$exe_file)
#~~~~~~~~~~~~~~~~#
# model summaries and diagnostics
print(fit_r$fit$summary(), n=300)
mean(fit_r$fit$summary()$rhat)
range(fit_r$fit$summary()$rhat)

mean(fit_r$fit$summary()$ess_bulk)
range(fit_r$fit$summary()$ess_bulk)

mean(fit_r$fit$summary()$ess_tail)
range(fit_r$fit$summary()$ess_tail)

fit_r$fit$cmdstan_summary()
fit_r$fit$cmdstan_diagnose()
fit_r$fit$diagnostic_summary()
mean(fit_r$fit$diagnostic_summary()$ebfmi)
fit_r$fit$sampler_diagnostics()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Cumulative curves with global and regional richness/masses 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative species richness as proportion
# cumulative species richness as proportion
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

curve_df_richness <- compute_region_global_richness_curves_ms(
  fit_obj = fit_r,     
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_df_richness)

# save
write.csv(curve_df_richness, "model_results/kleb_bsi_mlst_richness_curves_1.csv", row.names = FALSE)
#curve_df_richness <-  read.csv("model_results/kleb_bsi_mlst_richness_curves_1.csv")

#~~~~~~~~~~~~~~~~~#
curve_df <- compute_region_global_mass_curves_ms(
  fit_obj = fit_r,     
  f_grid = f_grid_99,  
  n_post_draws = 400)
#View(curve_df)
#colnames(curve_df)
#unique(curve_df$region)
# save
write.csv(curve_df, "model_results/kleb_bsi_mlst_mass_curves_1.csv", row.names = FALSE)
#curve_df <- read.csv("model_results/kleb_bsi_mlst_mass_curves_1.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.1.bii Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#plots to assess fit
pi_plot <- plot_global_freq_vs_pi_r(fit_r,
                                    fit_r_label = "Regional tau")
pi_plot 
# save
ggsave("model_results/kleb_bsi_mlst_overall_fit_hierarchical_model_alpha_1.png", 
       pi_plot, units = "in", width = 6, height = 4.5, dpi = 300)


theta_plot <- plot_regional_freq_vs_theta_r(fit_r,
                                            fit_r_label = "Regional tau")
theta_plot
# save 
ggsave("model_results/kleb_bsi_mlst_regional_fit_hierarchical_model_alpha_1.png", 
       theta_plot, units = "in", width = 10, height = 7, dpi = 300)

#~~~~~~~~~~~~~~~~~~~#
# make individual plots:
# species richness
p5 <- make_mass_plot_region(
  curve_df_richness,
  xlab = "Regional MLST frequency (f)",
  ylab = "Proportion of unique MLSTs\n with regional frequency ≥f",
  xlim_min = min(curve_df_richness$f[curve_df_richness$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p5
ggsave("model_results/kleb_bsi_mlst_regional_cumulative_richness_curve_regionl_tau_bayes_1.png", p5, units = "in", width = 6, height = 4, dpi = 300)


p6 <- make_cov_plot_region(curve_df_richness, 
                           pal = region_pal,
                           x_limit = c(10, 10000),
                           ylab = "Proportion of unique\nMLSTs detected")
p6
# save
ggsave("model_results/kleb_bsi_mlst_regional_richness_vs_ss_regionl_tau_bayes_1.png", p6, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make individual plots:
p7 <- make_mass_plot_region(
  curve_df,
  xlab = "Regional MLST frequency (f)",
  ylab = "Proportion of bacterial population\nwith MLST of frequency ≥f",
  xlim_min = min(curve_df$f[curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p7
ggsave("model_results/kleb_bsi_mlst_regional_cumulative_mass_curve_regionl_tau_bayes_1.png", p7, units = "in", width = 6, height = 4, dpi = 300)


p8 <- make_cov_plot_region(curve_df, pal = region_pal)
p8
# save
ggsave("model_results/kleb_bsi_mlst_regional_samplecoverage_vs_ss_regionl_tau_bayes_1.png", p8, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.1.bii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative richness and ss rows
national_df <- curve_df_richness |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_df_richness, national_df)

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
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_mlst_regional_bhm_summary_richness <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(kleb_bsi_mlst_regional_bhm_summary_richness)

# tidy table
kleb_bsi_mlst_regional_bhm_summary_richness <- kleb_bsi_mlst_regional_bhm_summary_richness |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_mlst_regional_bhm_summary_richness)
#View(kleb_bsi_mlst_regional_bhm_summary_richness)

# save
write.csv(kleb_bsi_mlst_regional_bhm_summary_richness, "model_results/kleb_bsi_mlst_regional_bhm_summary_richness_1.csv", row.names = FALSE)
#kleb_bsi_mlst_regional_bhm_summary_richness <- read.csv("model_results/kleb_bsi_mlst_regional_bhm_summary_richness.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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


thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_mlst_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
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
#View(kleb_bsi_mlst_regional_bhm_summary)

write.csv(kleb_bsi_mlst_regional_bhm_summary, "model_results/kleb_bsi_mlst_regional_bhm_summary_1.csv", row.names = FALSE)
#kleb_bsi_mlst_regional_bhm_summary <- read.csv("model_results/kleb_bsi_mlst_regional_bhm_summary_1.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 4.2 fastbaps_L3 CLUSTERS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 4.2a E. coli fastbaps_L3 - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#colnames(ecoli_bsi_samples_metadata)
ecoli_bsi_count_region <- ecoli_bsi_samples_metadata |>
  group_by(region, Level.3) |>
  summarise(count = n())
length(unique(ecoli_bsi_count_region$Level.3))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.1ai Sensitivity analysis to find best parameter values by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# regional tau
ecoli_fastbaps_L3_regional_results <- run_sensitivity_grid_regional_resumable(
  df = ecoli_bsi_count_region,
  region_col = "region",
  feature_col = "Level.3",
  count_col = "count",
  alpha_val_grid = c(0.1, 1, 10),
  alpha_novel_num_grid = c(0, 1, 10, 100, 500, 1000), 
  mu_tau_mean_grid =  c(0.1, 1, 10),
  mu_tau_sd_grid =  c(0.1, 1, 10),
  sigma_tau_rate_grid =  c(0.1, 1, 10),
  train_prop = 2/3,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  refresh = 100,#
  save_dir = "sens_results",
  checkpoint_file = "sens_results/ecoli_fastbaps_L3_regional_checkpoint_with_diagnostics_novel.rds"
)

# flag valid fits
ecoli_fastbaps_L3_regional_results <- flag_valid_fit(ecoli_fastbaps_L3_regional_results)

# save results
write.csv(ecoli_fastbaps_L3_regional_results, "sens_results/ecoli_fastbaps_L3_regional_tau_sensitivity_results_fixed_holdout_lps.csv", row.names = FALSE)
#ecoli_fastbaps_L3_regional_results <- read.csv("sens_results/ecoli_fastbaps_L3_regional_tau_sensitivity_results_fixed_holdout_lps.csv")
#View(ecoli_fastbaps_L3_regional_results)

# truncate extreme low-fit/high-fit values
upper <- quantile(ecoli_fastbaps_L3_regional_results$holdout_log_score_per_isolate, 0.975, na.rm = TRUE)
upper # -1.19
ecoli_fastbaps_L3_regional_results <- ecoli_fastbaps_L3_regional_results |>
  mutate(fill_val = pmin(holdout_log_score_per_isolate, upper))

# filter if needed
ecoli_fastbaps_L3_regional_results_neat <- ecoli_fastbaps_L3_regional_results |>
  filter(alpha_val %in% c(0.1, 1, 10),
         alpha_novel_num %in% c(1, 10, 100, 500, 1000), 
         mu_tau_mean %in%  c(0.1, 1, 10),
         mu_tau_sd %in%  c(0.1, 1, 10),
         sigma_tau_rate %in%  c(0.1, 1, 10))


#plot regional heatmap
ecoli_fastbaps_L3_regional_tau_plot <- plot_regional_holdout_heatmap(ecoli_fastbaps_L3_regional_results_neat, value_col = "fill_val")
ecoli_fastbaps_L3_regional_tau_plot
ggsave("sens_results/ecoli_fastbaps_L3_regional_tau_sensitivity_analysis_plot.png", ecoli_fastbaps_L3_regional_tau_plot, units = "in", width = 10, height = 4.2, dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.1aii  Run with single prior parameter combo ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set priors:
df <- ecoli_bsi_fastbaps_L3
counts <- df$count
N <- sum(df$count)
K <- nrow(df)
f1 <- sum(df$count == 1)  # 41      # number of singletons
q_hat <- f1 / N  # 2.8% Good-Turing first-order - proportion of singletons

# fit data to estimate Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(counts)
crp_fit$summary_df # mean: theta = 33.5-> 2% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.022

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$Level.3
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #4

ecoli_fastbaps_L3_prep <- prep_mlst_data(ecoli_bsi_count_region,
                                           region_col = "region",
                                           feature_col = "Level.3",
                                           count_col = "count",
                                           alpha_novel_sum = alpha_novel_sum,
                                           alpha_novel_num = alpha_novel_sum,
                                           alpha_other = 1)

# fit models
fit_r <- fit_region_tau_model(ecoli_fastbaps_L3_prep,
                              mu_tau_mean = 1,
                              mu_tau_sd =  1,
                              sigma_tau_rate = 1,
                              iter_warmup = 1000, 
                              iter_sampling = 1000, 
                              chains = 4)

# save
save_fit_bundle(fit_r, dir = "model_results", basename = "ecoli_bsi_fastbaps_L3_region_tau")
#~~~~~~~~~~~~~~~~#
# read back in
#obj_r <- load_fit_bundle("model_results", "ecoli_bsi_fastbaps_L3_region_tau")
#fit_r <- list(fit = obj_r$fit, prep = obj_r$prep, model = obj_r$meta$model,
#stan_file = obj_r$meta$stan_file, exe_file = obj_r$meta$exe_file)
#~~~~~~~~~~~~~~~~#
# model summaries and diagnostics
print(fit_r$fit$summary(), n=300)
mean(fit_r$fit$summary()$rhat)
range(fit_r$fit$summary()$rhat)

mean(fit_r$fit$summary()$ess_bulk)
range(fit_r$fit$summary()$ess_bulk)

mean(fit_r$fit$summary()$ess_tail)
range(fit_r$fit$summary()$ess_tail)

fit_r$fit$cmdstan_summary()
fit_r$fit$cmdstan_diagnose()
fit_r$fit$diagnostic_summary()
mean(fit_r$fit$diagnostic_summary()$ebfmi)
fit_r$fit$sampler_diagnostics()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Cumulative curves with global and regional masses 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative species richness as proportion
curve_df_richness <- compute_region_global_richness_curves_ms(
  fit_obj = fit_r,    # use fit_u as this is slightly better
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_df_richness)

# save
write.csv(curve_df_richness, "model_results/ecoli_bsi_fastbaps_L3_richness_curves.csv", row.names = FALSE)
#curve_df_richness <- read.csv("model_results/ecoli_bsi_fastbaps_L3_richness_curves.csv")

#~~~~~~~~~~~~~~~~~#
#f_grid <- c(seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
curve_df <- compute_region_global_mass_curves_ms(
  fit_obj = fit_r,    
  f_grid = f_grid_99,  
  n_post_draws = 400)
View(curve_df)
colnames(curve_df)
unique(curve_df$region)
# save
write.csv(curve_df, "model_results/ecoli_bsi_fastbaps_L3_mass_curves.csv", row.names = FALSE)
#curve_df <- read.csv("model_results/ecoli_bsi_fastbaps_L3_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2.ai Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#plots to assess fit
pi_plot <- plot_global_freq_vs_pi_r(fit_r,
                                    fit_r_label = "Regional tau")
pi_plot 
# save
ggsave("model_results/ecoli_bsi_fastbaps_L3_overall_fit_hierarchical_model.png", 
       pi_plot, units = "in", width = 6, height = 4.5, dpi = 300)

theta_plot <- plot_regional_freq_vs_theta_r(fit_r,
                                            fit_r_label = "Regional tau")
theta_plot
# save 
ggsave("model_results/ecoli_bsi_fastbaps_L3_regional_fit_hierarchical_model.png", 
       theta_plot, units = "in", width = 10, height = 7, dpi = 300)
#~~~~~~~~~~~~~~~~~~~#
# make individual plots:
# species richness
p9 <- make_mass_plot_region(
  curve_df_richness,
  xlab = "Regional fastBAPS cluster frequency (f)",
  ylab = "Proportion of unique fastBAPS clusters\nwith regional frequency ≥f",
  xlim_min = min(curve_df_richness$f[curve_df_richness$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p9
ggsave("model_results/ecoli_bsi_fastbaps_L3_regional_cumulative_richness_curve_regionl_tau_bayes.png", p9, units = "in", width = 6, height = 4, dpi = 300)


p10 <- make_cov_plot_region(curve_df_richness, 
                           pal = region_pal,
                           x_limit = c(10, 10000),
                           ylab = "Proportion of unique\nfastBAPS clusters detected")
p10
# save
ggsave("model_results/ecoli_bsi_fastbaps_L3_regional_richness_vs_ss_regionl_tau_bayes.png", p10, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~#
# make individual plots:
p11 <- make_mass_plot_region(
  curve_df,
  xlab = "Regional fastBAPS cluster frequency (f)",
  ylab = "Proportion of bacterial population\nin fastBAPS cluster of frequency ≥f",
  xlim_min = min(curve_df$f[curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p11
ggsave("model_results/ecoli_bsi_fastbaps_L3_regional_cumulative_mass_curve_regionl_tau_bayes.png", p11, units = "in", width = 6, height = 4, dpi = 300)


p12 <- make_cov_plot_region(curve_df, pal = region_pal)
p12
# save
ggsave("model_results/ecoli_bsi_fastbaps_L3_regional_samplecoverage_vs_ss_regionl_tau_bayes.png", p12, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2.aii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative richness and ss rows
national_df <- curve_df_richness |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_df_richness, national_df)

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
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness)

# tidy table
ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness <- ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness)
#View(ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness)

# save
write.csv(ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness, "model_results/ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness.csv", row.names = FALSE)
#ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness <- read.csv("model_results/ecoli_bsi_fastbaps_L3_regional_bhm_summary_richness.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
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
#View(ecoli_bsi_fastbaps_L3_regional_bhm_summary)

write.csv(ecoli_bsi_fastbaps_L3_regional_bhm_summary, "model_results/ecoli_bsi_fastbaps_L3_regional_bhm_summary.csv", row.names = FALSE)
#ecoli_bsi_fastbaps_L3_regional_bhm_summary <- read.csv("model_results/ecoli_bsi_fastbaps_L3_regional_bhm_summary.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 4.2b Klebsiella fastbaps_L3 - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
kleb_bsi_count_region <- kleb_bsi_samples_metadata |>
  group_by(region, Level.3) |>
  summarise(count = n())
length(unique(kleb_bsi_count_region$Level.3))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2bi Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# regional tau
kleb_fastbaps_L3_regional_results <- run_sensitivity_grid_regional_resumable(
  df = kleb_bsi_count_region,
  region_col = "region",
  feature_col = "Level.3",
  count_col = "count",
  alpha_val_grid = c(0.1, 1, 10),
  alpha_novel_num_grid = c(0, 1, 10, 100, 500, 1000), 
  mu_tau_mean_grid =  c(0.1, 1, 10),
  mu_tau_sd_grid =  c(0.1, 1, 10),
  sigma_tau_rate_grid =  c(0.1, 1, 10),
  train_prop = 2/3,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  refresh = 100,
  save_dir = "sens_results",
  checkpoint_file = "sens_results/kleb_fastbaps_L3_regional_checkpoint_with_diagnostics_novel.rds"
)

# flag valid models
kleb_fastbaps_L3_regional_results <- flag_valid_fit(kleb_fastbaps_L3_regional_results)

# save results
write.csv(kleb_fastbaps_L3_regional_results, "sens_results/kleb_fastbaps_L3_regional_tau_sensitivity_results_fixed_holdout_lps.csv", row.names = FALSE)
#kleb_fastbaps_L3_regional_results <- read.csv("sens_results/kleb_fastbaps_L3_regional_tau_sensitivity_results_fixed_holdout_lps.csv")
#View(kleb_fastbaps_L3_regional_results)

# truncate extreme low-fit/high-fit values
upper <- quantile(kleb_fastbaps_L3_regional_results$holdout_log_score_per_isolate, 0.975, na.rm = TRUE)
upper # -0.6422
kleb_fastbaps_L3_regional_results <- kleb_fastbaps_L3_regional_results |>
  mutate(fill_val = pmin(holdout_log_score_per_isolate, upper))

# filter
kleb_fastbaps_L3_regional_results_neat <- kleb_fastbaps_L3_regional_results |>
  filter(alpha_val %in% c(0.1, 1, 10),
         alpha_novel_num %in% c(1, 10, 100, 500, 1000), 
         mu_tau_mean %in%  c(0.1, 1, 10),
         mu_tau_sd %in%  c(0.1, 1, 10),
         sigma_tau_rate %in%  c(0.1, 1, 10))

#plot regional heatmap
kleb_fastbaps_L3_regional_tau_plot <- plot_regional_holdout_heatmap(kleb_fastbaps_L3_regional_results_neat, value_col = "fill_val")
kleb_fastbaps_L3_regional_tau_plot
ggsave("sens_results/kleb_fastbaps_L3_regional_tau_sensitivity_analysis_plot.png", kleb_fastbaps_L3_regional_tau_plot, units = "in", width = 10, height = 4.2, dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2bii  Run with single prior parameter combo ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set priors
df <- kleb_bsi_fastbaps_L3_df
counts <- df$count
N <- sum(df$count)
K <- nrow(df)
f1 <- sum(df$count == 1)  # 11      # number of singletons
q_hat <- f1 / N  # 2.3% Good-Turing first-order - proportion of singletons

# fit data to estimate Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(counts)
crp_fit$summary_df # mean: theta = 5.2 -> 1% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.0108

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$Level.3
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat)) # round up
alpha_novel_sum #1


kleb_fastbaps_L3_prep <- prep_mlst_data(kleb_bsi_count_region,
                       region_col = "region",
                       feature_col = "Level.3",
                       count_col = "count",
                       alpha_novel_sum = alpha_novel_sum,
                       alpha_novel_num = alpha_novel_sum,
                       alpha_other = 1)

# fit model
fit_r <- fit_region_tau_model(kleb_fastbaps_L3_prep,
                              mu_tau_mean = 1,
                              mu_tau_sd = 1, 
                              sigma_tau_rate = 1,
                              iter_warmup = 1000, 
                              iter_sampling = 1000, 
                              chains = 4)

# save
save_fit_bundle(fit_r, dir = "model_results", basename = "kleb_bsi_fastbaps_L3_region_tau")
#~~~~~~~~~~~~~~~~#
# read back in
#obj_r <- load_fit_bundle("model_results", "kleb_bsi_fastbaps_L3_region_tau")
#fit_r <- list(fit = obj_r$fit, prep = obj_r$prep, model = obj_r$meta$model,
#stan_file = obj_r$meta$stan_file, exe_file = obj_r$meta$exe_file)
#~~~~~~~~~~~~~~~~#
# model summaries and diagnostics
print(fit_r$fit$summary(), n=300)
mean(fit_r$fit$summary()$rhat)
range(fit_r$fit$summary()$rhat)

mean(fit_r$fit$summary()$ess_bulk)
range(fit_r$fit$summary()$ess_bulk)

mean(fit_r$fit$summary()$ess_tail)
range(fit_r$fit$summary()$ess_tail)

fit_r$fit$cmdstan_summary()
fit_r$fit$cmdstan_diagnose()
fit_r$fit$diagnostic_summary()
mean(fit_r$fit$diagnostic_summary()$ebfmi)
fit_r$fit$sampler_diagnostics()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Mass curves with global and regional masses 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative species richness as proportion
n_grid_course <- c(seq(0, 10000, by = 10), seq(10010, 50000, by = 20))
f_grid_99 <- 1-((1-0.99)^(1/n_grid_course))

curve_df_richness_course <- compute_region_global_richness_curves_ms(
  fit_obj = fit_r,     
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_df_richness_course)

# save
write.csv(curve_df_richness_course, "model_results/kleb_bsi_fastbaps_L3_richness_curves_course.csv", row.names = FALSE)
#curve_df_richness_course <- read.csv("model_results/kleb_bsi_fastbaps_L3_richness_curves_course.csv")
#table(curve_df_richness_course$region)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# repeat for finer f_grid for summary tables
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

curve_df_richness <- compute_region_global_richness_curves_ms(
  fit_obj = fit_r,    
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_df_richness)

# save
write.csv(curve_df_richness, "model_results/kleb_bsi_fastbaps_L3_richness_curves_fine.csv", row.names = FALSE)
#curve_df_richness <- read.csv("model_results/kleb_bsi_fastbaps_L3_richness_curves_fine.csv")
#table(curve_df_richness$region)

#~~~~~~~~~~~~~~~~~#
curve_df <- compute_region_global_mass_curves_ms(
  fit_obj = fit_r,    
  f_grid = f_grid_99,  
  n_post_draws = 400)
#View(curve_df)

# save
write.csv(curve_df, "model_results/kleb_bsi_fastbaps_L3_mass_curves.csv", row.names = FALSE)
#curve_df <- read.csv("model_results/kleb_bsi_fastbaps_L3_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2.biii Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#plots to assess fit
pi_plot <- plot_global_freq_vs_pi_r(fit_r,
                                    fit_r_label = "Regional tau")
pi_plot 
# save
ggsave("model_results/kleb_bsi_fastbaps_L3_overall_fit_hierarchical_model.png", 
       pi_plot, units = "in", width = 6, height = 4.5, dpi = 300)

theta_plot <- plot_regional_freq_vs_theta_r(fit_r,
                                            fit_r_label = "Regional tau")
theta_plot
# save 
ggsave("model_results/kleb_bsi_fastbaps_L3_regional_fit_hierarchical_model.png", 
       theta_plot, units = "in", width = 10, height = 7, dpi = 300)


#~~~~~~~~~~~~~~~~~~~#
# make individual plots:
# species richness
p13 <- make_mass_plot_region(
  curve_df_richness_course,
  xlab = "Regional fastBAPS cluster frequency (f)",
  ylab = "Proportion of unique fastBAPS clusters\nwith regional frequency ≥f",
  xlim_min = min(curve_df_richness_course$f[curve_df_richness_course$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p13
ggsave("model_results/kleb_bsi_fastbaps_L3_regional_cumulative_richness_curve_regionl_tau_bayes.png", p13, units = "in", width = 6, height = 4, dpi = 300)


p14 <- make_cov_plot_region(curve_df_richness_course, 
                           pal = region_pal,
                           x_limit = c(10, 10000),
                           ylab = "Proportion of unique\nfastBAPS clusters detected")
p14
# save
ggsave("model_results/kleb_bsi_fastbaps_L3_regional_richness_vs_ss_regionl_tau_bayes.png", p14, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make individual plots:
p15 <- make_mass_plot_region(
  curve_df,
  xlab = "Regional fastBAPS cluster frequency (f)",
  ylab = "Proportion of bacterial population\nin fastBAPS cluster of frequency ≥f",
  xlim_min = min(curve_df$f[curve_df$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p15
ggsave("model_results/kleb_bsi_fastbaps_L3_regional_cumulative_mass_curve_regionl_tau_bayes.png", p15, units = "in", width = 6, height = 4, dpi = 300)


p16 <- make_cov_plot_region(curve_df, pal = region_pal)
p16
# save
ggsave("model_results/kleb_bsi_fastbaos_L3_regional_samplecoverage_vs_ss_regionl_tau_bayes.png", p16, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 4.2.biv Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative richness and ss rows
national_df <- curve_df_richness |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_df_richness, national_df)

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
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_fastbaps_L3_regional_bhm_summary_richness <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(kleb_bsi_fastbaps_L3_regional_bhm_summary_richness)

# tidy table
kleb_bsi_fastbaps_L3_regional_bhm_summary_richness <- kleb_bsi_fastbaps_L3_regional_bhm_summary_richness |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_fastbaps_L3_regional_bhm_summary_richness)
#View(kleb_bsi_fastbaps_L3_regional_bhm_summary_richness)

# save
write.csv(kleb_bsi_fastbaps_L3_regional_bhm_summary_richness, "model_results/kleb_bsi_fastbaps_L3_regional_bhm_summary_richness.csv", row.names = FALSE)
#kleb_bsi_fastbaps_L3_regional_bhm_summary_richness <- read.csv("model_results/kleb_bsi_fastbaps_L3_regional_bhm_summary_richness.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
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
#View(kleb_bsi_fastbaps_L3_regional_bhm_summary)

write.csv(kleb_bsi_fastbaps_L3_regional_bhm_summary, "model_results/kleb_bsi_fastbaps_L3_regional_bhm_summary.csv", row.names = FALSE)
#kleb_bsi_fastbaps_L3_regional_bhm_summary_cells_only <- read.csv("model_results/kleb_bsi_fastbaps_L3_regional_bhm_summary_cels_only.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. HIERARCHICAL BAYESIAN BETA-BINOMIAL MODEL FOR ARGS AND PLASMIDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# write feature-level bayesian bootstrapping functions, using Beta-Binomial hierarchical distributions
prep_subiso_data <- function(df,
                             region_col = "region",
                             isolate_id_col = "sample",
                             feature_cols,
                             presence_threshold = 0,
                             beta_novel_sum = 0,
                             beta_novel_num = 0,
                             novel_prefix = "Novel_") {
  stopifnot(is.data.frame(df))
  stopifnot(all(c(region_col, feature_cols) %in% names(df)))
  
  if (length(beta_novel_num) != 1L || is.na(beta_novel_num) || beta_novel_num < 0) {
    stop("beta_novel_num must be a single non-negative integer.")
  }
  beta_novel_num <- as.integer(beta_novel_num)
  
  if (beta_novel_num > 0 && (length(beta_novel_sum) != 1L || is.na(beta_novel_sum))) {
    stop("beta_novel_sum must be a single numeric value when beta_novel_num > 0.")
  }
  
  df2 <- df |>
    dplyr::mutate(
      dplyr::across(all_of(feature_cols), ~ as.integer(.x > presence_threshold))
    )
  
  regions <- sort(unique(df2[[region_col]]))
  observed_features <- feature_cols
  
  # number of isolates per region
  n_by_region <- df2 |>
    dplyr::group_by(.data[[region_col]]) |>
    dplyr::summarise(n_iso = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(match(.data[[region_col]], regions))
  
  # region x observed-feature counts of positive isolates
  y_obs <- matrix(
    0L,
    nrow = length(regions),
    ncol = length(observed_features),
    dimnames = list(regions, observed_features)
  )
  
  for (r in regions) {
    rows_r <- df2[[region_col]] == r
    y_obs[r, ] <- colSums(df2[rows_r, observed_features, drop = FALSE])
  }
  
  # append novel zero-count features
  novel_features <- character(0)
  y_mat <- y_obs
  
  if (beta_novel_num > 0L) {
    novel_features <- paste0(novel_prefix, seq_len(beta_novel_num))
    y_novel <- matrix(
      0L,
      nrow = length(regions),
      ncol = beta_novel_num,
      dimnames = list(regions, novel_features)
    )
    y_mat <- cbind(y_obs, y_novel)
  }
  
  features <- c(observed_features, novel_features)
  
  # global counts/frequencies over the expanded feature set
  global_pos <- colSums(y_mat)
  global_n <- sum(n_by_region$n_iso)
  global_freq <- global_pos / global_n
  
  list(
    y = y_mat,
    n_iso = setNames(n_by_region$n_iso, n_by_region[[region_col]]),
    regions = regions,
    features = features,
    observed_features = observed_features,
    novel_features = novel_features,
    novel_idx = if (length(novel_features) > 0L) match(novel_features, features) else integer(0),
    beta_novel_sum = beta_novel_sum,
    beta_novel_num = beta_novel_num,
    beta_novel_value = if (beta_novel_num > 0L) beta_novel_sum / beta_novel_num else numeric(0),
    global_pos = global_pos,
    global_n = global_n,
    global_freq = global_freq
  )
}

# region-spcific tau model
stan_region_tau_subiso <- '
data {
  int<lower=1> R;
  int<lower=1> K;
  array[R, K] int<lower=0> y;
  array[R] int<lower=1> n_iso;
  vector<lower=0>[K] alpha_prior;
  vector<lower=0>[K] beta_prior;
  real<lower=0> mu_tau_mean;
  real<lower=0> mu_tau_sd;
  real<lower=0> sigma_tau_rate;
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

  mu_log_tau ~ normal(log(mu_tau_mean), mu_tau_sd);
  sigma_log_tau ~ exponential(sigma_tau_rate);
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# wrappers to compile and save models more efficiently
# (Optional) make CmdStanR write to a persistent directory instead of temp
options(cmdstanr_output_dir = "cmdstan_output")

# Write Stan code only if needed
write_stan_file_if_needed <- function(code, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  code_lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  
  if (!file.exists(file) || !identical(readLines(file), code_lines)) {
    writeLines(code_lines, con = file)
  }
  normalizePath(file, mustWork = TRUE)
}

# Compile once and reuse in-session
.stan_model_cache <- new.env(parent = emptyenv())

get_compiled_model <- function(model_name,
                               code,
                               stan_file,
                               exe_dir = "stan_executables",
                               force_recompile = FALSE) {
  if (exists(model_name, envir = .stan_model_cache, inherits = FALSE) && !force_recompile) {
    cached <- get(model_name, envir = .stan_model_cache, inherits = FALSE)
    if (!is.null(cached$exe_file) && file.exists(cached$exe_file)) {
      return(cached)
    }
  }
  
  stan_file <- write_stan_file_if_needed(code, stan_file)
  dir.create(exe_dir, recursive = TRUE, showWarnings = FALSE)
  
  mod <- cmdstan_model(stan_file, compile = FALSE)
  mod$compile(dir = exe_dir, force_recompile = force_recompile)
  
  out <- list(
    model = mod,
    stan_file = stan_file,
    exe_file = mod$exe_file()
  )
  
  assign(model_name, out, envir = .stan_model_cache)
  out
}

# fit wrapper function
fit_region_tau_subiso <- function(prep,
                                  alpha_prior = NULL,
                                  beta_prior = NULL,
                                  beta_novel_sum = NULL,
                                  beta_novel_num = NULL,
                                  mu_tau_mean = 1,
                                  mu_tau_sd = 0.01,
                                  sigma_tau_rate = 0.1,
                                  iter_warmup = 1000,
                                  iter_sampling = 1000,
                                  chains = 4,
                                  seed = 2026,
                                  model_cache_dir = "stan_cache",
                                  exe_dir = "stan_executables",
                                  force_recompile = FALSE) {
  
  compiled <- get_compiled_model(
    model_name = "region_tau_subiso",
    code = stan_region_tau_subiso,
    stan_file = file.path(model_cache_dir, "region_tau_subiso.stan"),
    exe_dir = exe_dir,
    force_recompile = force_recompile
  )
  
  K <- ncol(prep$y)
  R <- nrow(prep$y)
  
  # infer novel counts from prep unless explicitly overridden
  if (is.null(beta_novel_num)) beta_novel_num <- length(prep$novel_idx)
  beta_novel_num <- as.integer(beta_novel_num)
  
  K_obs <- K - beta_novel_num
  if (K_obs < 0L) stop("beta_novel_num is larger than the number of features in prep$y.")
  
  # observed-feature priors
  if (is.null(alpha_prior)) {
    alpha_prior <- rep(1, K_obs)
  } else if (length(alpha_prior) == 1L) {
    alpha_prior <- rep(alpha_prior, K_obs)
  }
  if (length(alpha_prior) != K_obs) {
    stop("alpha_prior must have length 1 or match the number of observed features.")
  }
  
  if (is.null(beta_prior)) {
    beta_prior <- rep(50, K_obs)
  } else if (length(beta_prior) == 1L) {
    beta_prior <- rep(beta_prior, K_obs)
  }
  if (length(beta_prior) != K_obs) {
    stop("beta_prior must have length 1 or match the number of observed features.")
  }
  
  # append novel-feature beta priors
  if (beta_novel_num > 0L) {
    if (is.null(beta_novel_sum)) {
      beta_novel_sum <- prep$beta_novel_sum
    }
    if (is.null(beta_novel_sum)) {
      stop("beta_novel_sum must be provided when beta_novel_num > 0.")
    }
    beta_novel_value <- beta_novel_sum / beta_novel_num
    
    alpha_prior <- c(alpha_prior, rep(1, beta_novel_num))
    beta_prior  <- c(beta_prior, rep(beta_novel_value, beta_novel_num))
  }
  
  fit <- compiled$model$sample(
    data = list(
      R = R,
      K = K,
      y = prep$y,
      n_iso = as.integer(prep$n_iso[prep$regions]),
      alpha_prior = alpha_prior,
      beta_prior = beta_prior,
      mu_tau_mean = mu_tau_mean,
      mu_tau_sd = mu_tau_sd,
      sigma_tau_rate = sigma_tau_rate
    ),
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  
  list(
    fit = fit,
    prep = prep,
    model = "region_tau_subiso",
    stan_file = compiled$stan_file,
    exe_file = compiled$exe_file
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# wrapper to extract posterior draws
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# stratified isolate-levle split into training and test set
split_subiso_holdout <- function(df,
                                 region_col = "region",
                                 isolate_id_col = "sample",
                                 train_prop = 2/3,
                                 seed = 2026) {
  stopifnot(all(c(region_col, isolate_id_col) %in% names(df)))
  set.seed(seed)
  
  isolates <- df |>
    dplyr::distinct(.data[[isolate_id_col]], .data[[region_col]])
  
  split_list <- split(isolates, isolates[[region_col]], drop = TRUE)
  
  split_list <- lapply(split_list, function(.x) {
    n <- nrow(.x)
    
    if (n <= 1L) {
      .x$split <- "train"
      return(.x)
    }
    
    n_train <- floor(train_prop * n)
    n_train <- max(1L, min(n - 1L, n_train))
    
    idx <- sample.int(n)
    .x$split <- "test"
    .x$split[idx[seq_len(n_train)]] <- "train"
    .x
  })
  
  isolates_split <- bind_rows(split_list)
  
  df_split <- df |>
    dplyr::left_join(isolates_split, by = c(isolate_id_col, region_col))
  
  list(
    train = df_split |> dplyr::filter(split == "train"),
    test  = df_split |> dplyr::filter(split == "test"),
    split = df_split
  )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Held-out log predictive score for sensitivity analysis
# helper function to calculate PMF for beta-binomial
lbbinom_manual_vec <- function(x, size, alpha, beta) {
  n <- max(length(x), length(size), length(alpha), length(beta))
  x <- rep_len(x, n)
  size <- rep_len(size, n)
  alpha <- rep_len(alpha, n)
  beta <- rep_len(beta, n)
  
  out <- rep(-Inf, n)
  
  ok <- is.finite(alpha) & is.finite(beta) &
    alpha > 0 & beta > 0 &
    is.finite(x) & is.finite(size) &
    x >= 0 & x <= size
  
  if (!any(ok)) return(out)
  
  xo <- x[ok]
  so <- size[ok]
  ao <- alpha[ok]
  bo <- beta[ok]
  
  out[ok] <-
    lgamma(so + 1) - lgamma(xo + 1) - lgamma(so - xo + 1) +
    lgamma(xo + ao) + lgamma(so - xo + bo) - lgamma(so + ao + bo) -
    (lgamma(ao) + lgamma(bo) - lgamma(ao + bo))
  
  out
}


# held-out log predictive score for regional tau 
score_holdout_subiso_regional <- function(fit_obj, test_prep, exclude_novel = TRUE) {
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  
  pi_mat   <- extract_pi_draws_subiso(draws_df, test_prep$features)
  tau_info <- extract_tau_draws_subiso(draws_df, regions = test_prep$regions)
  tau_mat  <- tau_info$region
  
  keep <- if (exclude_novel && !is.null(test_prep$novel_idx) && length(test_prep$novel_idx) > 0) {
    rep(TRUE, ncol(test_prep$y))
  } else {
    rep(TRUE, ncol(test_prep$y))
  }
  
  # If novel columns are stored explicitly, drop them from scoring
  if (exclude_novel && !is.null(test_prep$novel_idx) && length(test_prep$novel_idx) > 0) {
    keep[test_prep$novel_idx] <- FALSE
  }
  
  if (!any(keep)) stop("All features removed when excluding novel categories.")
  
  pi_mat <- pi_mat[, keep, drop = FALSE]
  y_kept  <- test_prep$y[, keep, drop = FALSE]
  
  R <- nrow(y_kept)
  S <- nrow(pi_mat)
  
  n_iso_vec <- test_prep$n_iso[test_prep$regions]
  
  lp_draws <- numeric(S)
  
  for (s in seq_len(S)) {
    lp <- 0
    
    for (r in seq_len(R)) {
      y_rk <- y_kept[r, ]
      n_r  <- n_iso_vec[[r]]
      tau_r <- tau_mat[s, r]
      
      alpha <- tau_r * pi_mat[s, ]
      beta  <- tau_r * (1 - pi_mat[s, ])
      
      lp <- lp + sum(lbbinom_manual_vec(
        x = y_rk,
        size = n_r,
        alpha = alpha,
        beta = beta
      ))
    }
    
    lp_draws[s] <- lp
  }
  
  total_score <- log_mean_exp(lp_draws)
  n_test <- sum(y_kept)
  
  tibble::tibble(
    holdout_log_score_total = total_score,
    holdout_log_score_per_isolate = if (n_test > 0) total_score / n_test else NA_real_,
    n_test = n_test
  )
}


# diagnostics
make_combo_id_subiso <- function(alpha, beta, beta_novel_num,
                                 mu_tau_mean, mu_tau_sd, sigma_tau_rate) {
  paste0(
    "a", alpha,
    "_b", beta,
    "_bn", beta_novel_num,
    "_mu", mu_tau_mean,
    "_sd", mu_tau_sd,
    "_sr", sigma_tau_rate
  )
}


# function to extract model run diagnostics
extract_diagnostics <- function(fit_obj, iter_sampling = NULL, chains = NULL) {
  
  diag <- fit_obj$fit$diagnostic_summary()
  
  # aggregate chain-level diagnostics 
  num_divergent <- sum(diag$num_divergent)
  num_max_treedepth <- sum(diag$num_max_treedepth)
  ebfmi <- mean(diag$ebfmi, na.rm = TRUE)
  
  # compute fraction of treedepth hits
  if (!is.null(iter_sampling) && !is.null(chains)) {
    treedepth_frac <- num_max_treedepth / (iter_sampling * chains)
  } else {
    treedepth_frac <- NA_real_
  }
  
  # parameter summaries 
  summ <- fit_obj$fit$summary()
  
  rhat <- summ$rhat
  ess  <- summ$ess_bulk
  
  summarise_vec <- function(x, prefix) {
    tibble::tibble(
      !!paste0(prefix, "_min")   := min(x, na.rm = TRUE),
      !!paste0(prefix, "_q2.5")  := quantile(x, 0.025, na.rm = TRUE),
      !!paste0(prefix, "_q25")   := quantile(x, 0.25,  na.rm = TRUE),
      !!paste0(prefix, "_median"):= median(x, na.rm = TRUE),
      !!paste0(prefix, "_q75")   := quantile(x, 0.75,  na.rm = TRUE),
      !!paste0(prefix, "_q97.5") := quantile(x, 0.975, na.rm = TRUE),
      !!paste0(prefix, "_max")   := max(x, na.rm = TRUE)
    )
  }
  
  dplyr::bind_cols(
    tibble::tibble(
      num_divergent = num_divergent,
      num_max_treedepth = num_max_treedepth,
      treedepth_frac = treedepth_frac,
      ebfmi = ebfmi
    ),
    summarise_vec(rhat, "rhat"),
    summarise_vec(ess, "ess_bulk")
  )
}


empty_diagnostics <- function() {
  tibble::tibble(
    num_divergent = NA_real_,
    num_max_treedepth = NA_real_,
    treedepth_frac = NA_real_,
    ebfmi = NA_real_,
    
    rhat_min = NA_real_,
    rhat_q2.5 = NA_real_,
    rhat_q25 = NA_real_,
    rhat_median = NA_real_,
    rhat_q75 = NA_real_,
    rhat_q97.5 = NA_real_,
    rhat_max = NA_real_,
    
    ess_bulk_min = NA_real_,
    ess_bulk_q2.5 = NA_real_,
    ess_bulk_q25 = NA_real_,
    ess_bulk_median = NA_real_,
    ess_bulk_q75 = NA_real_,
    ess_bulk_q97.5 = NA_real_,
    ess_bulk_max = NA_real_
  )
}

run_one_regional_combo_holdout_subiso <- function(train_prep,
                                                  test_prep,
                                                  alpha,
                                                  beta,
                                                  beta_novel_num,
                                                  mu_tau_mean,
                                                  mu_tau_sd,
                                                  sigma_tau_rate,
                                                  iter_warmup,
                                                  iter_sampling,
                                                  chains,
                                                  seed) {
  tryCatch({
    
    K <- ncol(train_prep$y)
    K_obs <- K - as.integer(beta_novel_num)
    
    if (K_obs < 0L) stop("beta_novel_num is larger than the number of columns in prep$y.")
    
    fit <- fit_region_tau_subiso(
      prep = train_prep,
      alpha_prior = rep(alpha, K_obs),
      beta_prior  = rep(beta, K_obs),
      beta_novel_sum = beta_novel_num * beta,
      beta_novel_num = beta_novel_num,
      mu_tau_mean = mu_tau_mean,
      mu_tau_sd = mu_tau_sd,
      sigma_tau_rate = sigma_tau_rate,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = seed
    )
    
    score <- score_holdout_subiso_regional(fit, test_prep)
    
    diag <- extract_diagnostics(
      fit_obj = fit,
      iter_sampling = iter_sampling,
      chains = chains
    )
    
    dplyr::bind_cols(
      tibble(
        alpha = alpha,
        beta = beta,
        beta_novel_num = beta_novel_num,
        beta_novel_sum = beta_novel_num * beta,
        mu_tau_mean = mu_tau_mean,
        mu_tau_sd = mu_tau_sd,
        sigma_tau_rate = sigma_tau_rate,
        status = "ok",
        error = NA_character_
      ),
      score,
      diag
    )
    
  }, error = function(e) {
    dplyr::bind_cols(
      tibble(
        alpha = alpha,
        beta = beta,
        beta_novel_num = beta_novel_num,
        beta_novel_sum = beta_novel_num * beta,
        mu_tau_mean = mu_tau_mean,
        mu_tau_sd = mu_tau_sd,
        sigma_tau_rate = sigma_tau_rate,
        status = "failed",
        error = conditionMessage(e),
        holdout_log_score_total = NA_real_,
        holdout_log_score_per_isolate = NA_real_
      ),
      empty_diagnostics()
    )
  })
}

# run sensitivity grid -regional - fixed split, and resumable
run_sensitivity_regional_subiso_resumable <- function(df,
                                                      feature_cols,
                                                      alpha_grid,
                                                      beta_grid,
                                                      beta_novel_num_grid,
                                                      mu_tau_mean_grid,
                                                      mu_tau_sd_grid,
                                                      sigma_tau_rate_grid,
                                                      seed_split = 2026,
                                                      train_prop = 2/3,
                                                      iter_warmup = 500,
                                                      iter_sampling = 500,
                                                      chains = 1,
                                                      checkpoint_file = "sens_results/beta_binomial_regional_checkpoint.rds") {
  
  dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
  
  results <- if (file.exists(checkpoint_file)) readRDS(checkpoint_file) else tibble()
  
  split <- split_subiso_holdout(
    df = df,
    seed = seed_split,
    train_prop = train_prop
  )
  
  grid <- tidyr::crossing(
    alpha = alpha_grid,
    beta = beta_grid,
    beta_novel_num = beta_novel_num_grid,
    mu_tau_mean = mu_tau_mean_grid,
    mu_tau_sd = mu_tau_sd_grid,
    sigma_tau_rate = sigma_tau_rate_grid
  ) |>
    dplyr::mutate(
      combo_id = make_combo_id_subiso(
        alpha, beta, beta_novel_num, mu_tau_mean, mu_tau_sd, sigma_tau_rate
      )
    )
  
  if (nrow(results) > 0 && "combo_id" %in% names(results)) {
    grid <- grid |> dplyr::filter(!combo_id %in% results$combo_id)
  }
  
  if (nrow(grid) == 0L) return(results)
  
  # Build preps once per beta_novel_num
  prep_cache <- vector("list", length(unique(grid$beta_novel_num)))
  names(prep_cache) <- as.character(sort(unique(grid$beta_novel_num)))
  
  for (bn in sort(unique(grid$beta_novel_num))) {
    beta_novel_sum <- bn * unique(grid$beta[grid$beta_novel_num == bn])[1]
    
    prep_cache[[as.character(bn)]] <- list(
      train = prep_subiso_data(
        split$train,
        feature_cols = feature_cols,
        beta_novel_num = bn,
        beta_novel_sum = beta_novel_sum
      ),
      test = prep_subiso_data(
        split$test,
        feature_cols = feature_cols,
        beta_novel_num = bn,
        beta_novel_sum = beta_novel_sum
      )
    )
  }
  
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]
    
    preps <- prep_cache[[as.character(row$beta_novel_num)]]
    
    res <- run_one_regional_combo_holdout_subiso(
      train_prep = preps$train,
      test_prep = preps$test,
      alpha = row$alpha,
      beta = row$beta,
      beta_novel_num = row$beta_novel_num,
      mu_tau_mean = row$mu_tau_mean,
      mu_tau_sd = row$mu_tau_sd,
      sigma_tau_rate = row$sigma_tau_rate,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = seed_split + i
    ) |>
      dplyr::mutate(combo_id = row$combo_id)
    
    results <- dplyr::bind_rows(results, res)
    saveRDS(results, checkpoint_file)
  }
  
  results
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# same function as for isolate-level features
flag_valid_fit <- function(df, iterations = 250) {
  df |>
    mutate(num_divergent_frac = num_divergent / iterations) |>
    mutate(
      valid_fit = status == "ok" &
        is.finite(num_divergent) &
        is.finite(treedepth_frac) &
        is.finite(rhat_max) &
        is.finite(ess_bulk_min) &
        num_divergent_frac < 0.05 &
        treedepth_frac < 0.05 &
        rhat_q75 < 1.05 &
        ess_bulk_q25 > 100 
    )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
plot_regional_holdout_heatmap_subiso_one_combo <- function(regional_results,
                                                                beta_value,
                                                                sigma_value,
                                                                value_col = "holdout_log_score_per_isolate") {
  
  stopifnot(value_col %in% names(regional_results))
  stopifnot(all(c("alpha", "beta", "beta_novel_num",
                  "mu_tau_mean", "mu_tau_sd", "sigma_tau_rate",
                  "valid_fit") %in% names(regional_results)))
  
  fill_limits <- range(regional_results[[value_col]], na.rm = TRUE)
  fill_breaks <- pretty(fill_limits, n = 5)
  
  mu_mean_levels <- sort(unique(regional_results$mu_tau_mean))
  mu_sd_levels   <- sort(unique(regional_results$mu_tau_sd))
  
  mu_mean_labels <- paste0("mu[mu*log(tau)]==", mu_mean_levels)
  mu_sd_labels   <- paste0("sigma[mu*log(tau)]==", mu_sd_levels)
  
  df <- regional_results |>
    filter(beta == beta_value, sigma_tau_rate == sigma_value) |>
    group_by(alpha, beta, beta_novel_num, mu_tau_mean, mu_tau_sd, sigma_tau_rate) |>
    slice_tail(n = 1) |>
    ungroup() |>
    mutate(
      valid_fit = coalesce(valid_fit, FALSE),
      fill_val = .data[[value_col]],
      tile_alpha = if_else(valid_fit, 1, 0.35),
      alpha = factor(alpha, levels = sort(unique(regional_results$alpha))),
      beta_novel_num = factor(beta_novel_num, levels = sort(unique(regional_results$beta_novel_num))),
      mu_tau_mean = factor(
        paste0("mu[mu*log(tau)]==", mu_tau_mean),
        levels = mu_mean_labels
      ),
      mu_tau_sd = factor(
        paste0("sigma[mu*log(tau)]==", mu_tau_sd),
        levels = mu_sd_labels
      )
    )
  
  ggplot(df, aes(x = alpha, y = beta_novel_num, fill = fill_val, alpha = tile_alpha)) +
    geom_tile(width = 0.95, height = 0.95, colour = "white", linewidth = 0.25) +
    facet_grid(
      rows = vars(mu_tau_sd),
      cols = vars(mu_tau_mean),
      labeller = labeller(
        mu_tau_mean = label_parsed,
        mu_tau_sd = label_parsed
      )
    ) +
    scale_fill_viridis_c(
      name = "Hold-out score",
      limits = fill_limits,
      breaks = fill_breaks,
      oob = scales::squish,
      na.value = "grey90"
    ) +
    scale_alpha_identity(guide = "none") +
    labs(
      x = expression(beta),
      y = expression(m),
      title = bquote(beta * "’" == .(beta_value) ~ "," ~ lambda[sigma*log(tau)] == .(sigma_value))
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
}

plot_regional_holdout_heatmap_subiso <- function(regional_results,
                                                      value_col = "holdout_log_score_per_isolate") {
  
  beta_values <- sort(unique(regional_results$beta))
  sigma_values <- sort(unique(regional_results$sigma_tau_rate))
  
  plots <- list()
  
  for (b in beta_values) {
    for (s in sigma_values) {
      plots[[length(plots) + 1]] <- plot_regional_holdout_heatmap_subiso_one_combo(
        regional_results = regional_results,
        beta_value = b,
        sigma_value = s,
        value_col = value_col
      )
    }
  }
  
  wrap_plots(plots, ncol = length(sigma_values)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Helper functions for mass/ richness curve estimation
tail_count_desc <- function(x_sorted_desc, f_grid) {
  # number of elements >= each threshold in f_grid
  findInterval(-f_grid, -x_sorted_desc)
}

tail_sum_desc <- function(x_sorted_desc, f_grid) {
  # sum of elements >= each threshold in f_grid
  idx <- tail_count_desc(x_sorted_desc, f_grid)
  cs <- cumsum(x_sorted_desc)
  out <- numeric(length(f_grid))
  keep <- idx > 0L
  out[keep] <- cs[idx[keep]]
  out
}

summarise_draw_matrix <- function(mat, prefix) {
  tibble::tibble(
    !!paste0(prefix, "_median") := matrixStats::colMedians(mat, na.rm = TRUE),
    !!paste0(prefix, "_q2.5")   := matrixStats::colQuantiles(mat, probs = 0.025, na.rm = TRUE, drop = TRUE),
    !!paste0(prefix, "_q97.5")  := matrixStats::colQuantiles(mat, probs = 0.975, na.rm = TRUE, drop = TRUE)
  )
}

compute_global_all_matrix <- function(pi_draws, f_grid, mode = c("richness", "mass")) {
  mode <- match.arg(mode)
  
  nd <- nrow(pi_draws)
  nf <- length(f_grid)
  out <- matrix(NA_real_, nrow = nd, ncol = nf)
  
  for (i in seq_len(nd)) {
    pi_b <- as.numeric(pi_draws[i, ])
    
    if (mode == "richness") {
      pi_sorted <- sort(pi_b, decreasing = TRUE)
      out[i, ] <- tail_count_desc(pi_sorted, f_grid) / length(pi_b)
    } else {
      pi_sorted <- sort(pi_b, decreasing = TRUE)
      s <- sum(pi_sorted)
      if (s > 0) {
        out[i, ] <- tail_sum_desc(pi_sorted, f_grid) / s
      } else {
        out[i, ] <- 0
      }
    }
  }
  
  out
}

region_worker_richness <- function(
    region_name,
    y_r,
    n_r,
    pi_draws,
    tau_vec,
    f_grid,
    total_k,
    global_all_mat) {
  
  nd <- nrow(pi_draws)
  nf <- length(f_grid)
  
  y_r_non_zero <- sum(y_r > 0)
  if (y_r_non_zero == 0) {
    y_r_non_zero <- NA_real_
  }
  
  region_richness_mat <- matrix(NA_real_, nrow = nd, ncol = nf)
  regional_richness_of_global_mat <- matrix(NA_real_, nrow = nd, ncol = nf)
  
  for (i in seq_len(nd)) {
    pi_b <- as.numeric(pi_draws[i, ])
    tau_b <- tau_vec[i]
    
    theta_b <- rbeta(
      n = length(pi_b),
      shape1 = y_r + tau_b * pi_b,
      shape2 = (n_r - y_r) + tau_b * (1 - pi_b)
    )
    
    theta_sorted <- sort(theta_b, decreasing = TRUE)
    k_theta <- tail_count_desc(theta_sorted, f_grid)
    
    region_richness_mat[i, ] <- k_theta / y_r_non_zero
    regional_richness_of_global_mat[i, ] <- k_theta / total_k
  }
  
  tibble::tibble(
    region = region_name,
    f = f_grid
  ) |>
    dplyr::bind_cols(
      summarise_draw_matrix(region_richness_mat, "region"),
      summarise_draw_matrix(regional_richness_of_global_mat, "global"),
      summarise_draw_matrix(global_all_mat, "global_all")
    ) |>
    dplyr::mutate(
      min_sample_90 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.90) / log(1 - f)), NA_integer_),
      min_sample_95 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.95) / log(1 - f)), NA_integer_),
      min_sample_99 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.99) / log(1 - f)), NA_integer_)
    )
}

region_worker_mass <- function(
    region_name,
    y_r,
    n_r,
    pi_draws,
    tau_vec,
    f_grid,
    global_all_mat) {
  
  nd <- nrow(pi_draws)
  nf <- length(f_grid)
  
  region_mass_mat <- matrix(NA_real_, nrow = nd, ncol = nf)
  global_mass_mat <- matrix(NA_real_, nrow = nd, ncol = nf)
  
  for (i in seq_len(nd)) {
    pi_b <- as.numeric(pi_draws[i, ])
    tau_b <- tau_vec[i]
    
    theta_b <- rbeta(
      n = length(pi_b),
      shape1 = y_r + tau_b * pi_b,
      shape2 = (n_r - y_r) + tau_b * (1 - pi_b)
    )
    
    ord_theta <- order(theta_b, decreasing = TRUE)
    theta_sorted <- theta_b[ord_theta]
    pi_in_theta_order <- pi_b[ord_theta]
    
    k_theta <- tail_count_desc(theta_sorted, f_grid)
    
    sum_theta <- sum(theta_sorted)
    if (sum_theta > 0) {
      region_mass_mat[i, ] <- tail_sum_desc(theta_sorted, f_grid) / sum_theta
    } else {
      region_mass_mat[i, ] <- 0
    }
    
    sum_pi_sel <- numeric(nf)
    cs_pi_theta <- cumsum(pi_in_theta_order)
    keep <- k_theta > 0L
    sum_pi_sel[keep] <- cs_pi_theta[k_theta[keep]]
    
    sum_pi_all <- sum(pi_b)
    if (sum_pi_all > 0) {
      global_mass_mat[i, ] <- sum_pi_sel / sum_pi_all
    } else {
      global_mass_mat[i, ] <- 0
    }
  }
  
  tibble::tibble(
    region = region_name,
    f = f_grid
  ) |>
    dplyr::bind_cols(
      summarise_draw_matrix(region_mass_mat, "region"),
      summarise_draw_matrix(global_mass_mat, "global"),
      summarise_draw_matrix(global_all_mat, "global_all")
    ) |>
    dplyr::mutate(
      min_sample_90 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.90) / log(1 - f)), NA_integer_),
      min_sample_95 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.95) / log(1 - f)), NA_integer_),
      min_sample_99 = dplyr::if_else(f > 0 & f < 1, ceiling(log(1 - 0.99) / log(1 - f)), NA_integer_)
    )
}

# Fast richness curves 
compute_subiso_richness_curves_ms <- function(
    fit_obj,
    f_grid,
    n_post_draws = 1000,
    seed = 2026,
    cores = max(1L, parallel::detectCores(logical = FALSE) - 1L)) {
  
  set.seed(seed)
  
  prep <- fit_obj$prep
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  
  pi_draws <- extract_pi_draws_subiso(draws_df, prep$features)
  tau_info <- extract_tau_draws_subiso(draws_df, regions = prep$regions)
  
  B <- nrow(pi_draws)
  draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
  
  pi_draws <- pi_draws[draw_idx, , drop = FALSE]
  
  nd <- nrow(pi_draws)
  nf <- length(f_grid)
  
  # global_all is independent of region, so compute it once
  global_all_mat <- compute_global_all_matrix(pi_draws, f_grid, mode = "richness")
  
  region_specs <- lapply(seq_along(prep$regions), function(r) {
    region_name <- prep$regions[r]
    y_r <- as.numeric(prep$y[r, ])
    n_r <- as.integer(prep$n_iso[region_name])
    
    tau_vec <- if (!is.null(tau_info$shared)) {
      tau_info$shared[seq_len(nd)]
    } else {
      tau_info$region[seq_len(nd), region_name]
    }
    
    list(
      region_name = region_name,
      y_r = y_r,
      n_r = n_r,
      tau_vec = tau_vec
    )
  })
  
  worker_fun <- function(spec) {
    region_worker_richness(
      region_name = spec$region_name,
      y_r = spec$y_r,
      n_r = spec$n_r,
      pi_draws = pi_draws,
      tau_vec = spec$tau_vec,
      f_grid = f_grid,
      total_k = length(prep$features),
      global_all_mat = global_all_mat
    )
  }
  
  res <- if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(region_specs, worker_fun, mc.cores = cores)
  } else {
    lapply(region_specs, worker_fun)
  }
  
  dplyr::bind_rows(res)
}

# Fast mass curves 
compute_subiso_mass_curves_ms <- function(
    fit_obj,
    f_grid,
    n_post_draws = 400,
    seed = 2026,
    cores = max(1L, parallel::detectCores(logical = FALSE) - 1L)) {
  
  set.seed(seed)
  
  prep <- fit_obj$prep
  draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
  
  pi_draws <- extract_pi_draws_subiso(draws_df, prep$features)
  tau_info <- extract_tau_draws_subiso(draws_df, regions = prep$regions)
  
  B <- nrow(pi_draws)
  draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
  
  pi_draws <- pi_draws[draw_idx, , drop = FALSE]
  
  nd <- nrow(pi_draws)
  nf <- length(f_grid)
  
  # global_all is independent of region, so compute it once
  global_all_mat <- compute_global_all_matrix(pi_draws, f_grid, mode = "mass")
  
  region_specs <- lapply(seq_along(prep$regions), function(r) {
    region_name <- prep$regions[r]
    y_r <- as.numeric(prep$y[r, ])
    n_r <- as.integer(prep$n_iso[region_name])
    
    tau_vec <- if (!is.null(tau_info$shared)) {
      tau_info$shared[seq_len(nd)]
    } else {
      tau_info$region[seq_len(nd), region_name]
    }
    
    list(
      region_name = region_name,
      y_r = y_r,
      n_r = n_r,
      tau_vec = tau_vec
    )
  })
  
  worker_fun <- function(spec) {
    region_worker_mass(
      region_name = spec$region_name,
      y_r = spec$y_r,
      n_r = spec$n_r,
      pi_draws = pi_draws,
      tau_vec = spec$tau_vec,
      f_grid = f_grid,
      global_all_mat = global_all_mat
    )
  }
  
  res <- if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(region_specs, worker_fun, mc.cores = cores)
  } else {
    lapply(region_specs, worker_fun)
  }
  
  dplyr::bind_rows(res)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot code for model comparison to actual global frequencies
plot_global_freq_vs_pi_subiso_r <- function(
    fit_r,
    fit_r_label = "fit_r",
    n_post_draws = 400,
    seed = 2026,
    conf_level = 0.95,
    pseudo_count = NULL,
    point_size = 2.2,
    point_alpha = 0.3,
    line_size = 0.8) {
  
  set.seed(seed)
  
  summarise_one_fit <- function(fit_obj, model_label) {
    prep <- fit_obj$prep
    features <- prep$features
    draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
    pi_draws <- as.matrix(extract_pi_draws_subiso(draws_df, features))
    
    if (!is.null(colnames(pi_draws)) && all(features %in% colnames(pi_draws))) {
      pi_draws <- pi_draws[, features, drop = FALSE]
    }
    
    B <- nrow(pi_draws)
    draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
    
    alpha <- (1 - conf_level) / 2
    pi_sub <- pi_draws[draw_idx, , drop = FALSE]
    pi_median <- matrixStats::colMedians(pi_sub, na.rm = TRUE)
    pi_lo <- matrixStats::colQuantiles(pi_sub, probs = alpha, na.rm = TRUE, drop = TRUE)
    pi_hi <- matrixStats::colQuantiles(pi_sub, probs = 1 - alpha, na.rm = TRUE, drop = TRUE)
    
    gf <- prep$global_freq
    if (!is.null(names(gf)) && all(features %in% names(gf))) {
      gf <- gf[features]
    } else {
      gf <- as.numeric(gf)
      if (length(gf) != length(features)) {
        stop("prep$global_freq does not match the number of features.")
      }
    }
    
    tibble::tibble(
      feature = features,
      model = model_label,
      actual_freq = as.numeric(gf),
      estimate_median = as.numeric(pi_median),
      estimate_lo = as.numeric(pi_lo),
      estimate_hi = as.numeric(pi_hi)
    )
    
  }
  
  if (is.null(pseudo_count)) {
    pseudo_count <- 0.5 / max(sum(fit_r$prep$n_iso), 1)
  }
  
  df_r <- summarise_one_fit(fit_r, fit_r_label)
  
  plot_df <- df_r |>
    dplyr::mutate(
      actual_plot = pmax(actual_freq, pseudo_count),
      estimate_plot = pmax(estimate_median, pseudo_count),
      lo_plot = pmax(estimate_lo, pseudo_count),
      hi_plot = pmax(estimate_hi, pseudo_count)
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = actual_plot, y = estimate_plot)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", linewidth = line_size, color = "grey50"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo_plot, ymax = hi_plot),
      width = 0, alpha = point_alpha
    ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Actual global frequency",
      y = "Estimated pi (posterior median)",
      title = "Global frequency vs estimated prevalence"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

plot_regional_freq_vs_theta_subiso_r <- function(
    fit_r,
    fit_r_label = "fit_r",
    n_post_draws = 400,
    seed = 2026,
    conf_level = 0.95,
    pseudo_count = NULL,
    point_size = 1.8,
    point_alpha = 0.3,
    line_size = 0.7,
    ncol = NULL) {
  
  set.seed(seed)
  
  summarise_one_fit <- function(fit_obj, model_label) {
    prep <- fit_obj$prep
    features <- prep$features
    regions <- prep$regions
    R <- length(regions)
    K <- length(features)
    
    draws_df <- posterior::as_draws_df(fit_obj$fit$draws())
    pi_draws <- as.matrix(extract_pi_draws_subiso(draws_df, features))
    tau_info <- extract_tau_draws_subiso(draws_df, regions = regions)
    
    if (!is.null(colnames(pi_draws)) && all(features %in% colnames(pi_draws))) {
      pi_draws <- pi_draws[, features, drop = FALSE]
    }
    
    B <- nrow(pi_draws)
    draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
    
    alpha <- (1 - conf_level) / 2
    
    actual_mat <- sweep(prep$y, 1, prep$n_iso, "/")
    colnames(actual_mat) <- features
    rownames(actual_mat) <- regions
    
    theta_arr <- array(NA_real_, dim = c(length(draw_idx), R, K))
    
    for (i in seq_along(draw_idx)) {
      b <- draw_idx[i]
      pi_b <- as.numeric(pi_draws[b, ])
      
      for (r in seq_along(regions)) {
        region_name <- regions[r]
        y_r <- as.numeric(prep$y[r, ])
        n_r <- as.integer(prep$n_iso[region_name])
        
        tau_b <- if (!is.null(tau_info$shared)) {
          tau_info$shared[b]
        } else {
          tau_info$region[b, region_name]
        }
        
        theta_b <- rbeta(
          n = length(pi_b),
          shape1 = y_r + tau_b * pi_b,
          shape2 = (n_r - y_r) + tau_b * (1 - pi_b)
        )
        
        theta_arr[i, r, ] <- theta_b
      }
    }
    
    out <- vector("list", R)
    
    for (r in seq_along(regions)) {
      est_draws <- theta_arr[, r, , drop = FALSE][, 1, ]
      
      est_median <- matrixStats::colMedians(est_draws, na.rm = TRUE)
      est_lo <- matrixStats::colQuantiles(est_draws, probs = alpha, na.rm = TRUE, drop = TRUE)
      est_hi <- matrixStats::colQuantiles(est_draws, probs = 1 - alpha, na.rm = TRUE, drop = TRUE)
      
      out[[r]] <- tibble::tibble(
        region = regions[r],
        feature = features,
        model = model_label,
        actual_freq = as.numeric(actual_mat[r, ]),
        estimate_median = as.numeric(est_median),
        estimate_lo = as.numeric(est_lo),
        estimate_hi = as.numeric(est_hi)
      )
    }
    dplyr::bind_rows(out)
  }
  
  if (is.null(pseudo_count)) {
    pseudo_count <- 0.5 / max(max(rowSums(fit_r$prep$y)), 1)
  }
  
  df_r <- summarise_one_fit(fit_r, fit_r_label)
  
  plot_df <- df_r |>
    dplyr::mutate(
      actual_plot = pmax(actual_freq, pseudo_count),
      estimate_plot = pmax(estimate_median, pseudo_count),
      lo_plot = pmax(estimate_lo, pseudo_count),
      hi_plot = pmax(estimate_hi, pseudo_count)
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = actual_plot, y = estimate_plot)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", linewidth = line_size, color = "grey50"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo_plot, ymax = hi_plot),
      width = 0, alpha = point_alpha
    ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ region, ncol = ncol) +
    ggplot2::labs(
      x = "Actual regional frequency",
      y = "Estimated theta (posterior median)",
      title = "Regional frequencies vs estimated prevalence"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 5.1 ARGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 5.1a E. coli ARGs - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load and prepare data, if not already loaded
ecoli_bsi_amrfinder_metadata <- read.csv("neksus_ecoli_bsi_amrfinder_metadata.csv")
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
#View(ecoli_bsi_arg_presence_region)
#colnames(ecoli_bsi_arg_presence_region)
#length(unique(ecoli_bsi_arg_presence$sample)) # 1471
gene_cols <- setdiff(names(ecoli_bsi_arg_presence_region), c("sample", "region"))
prep_arg <- prep_subiso_data(
  df = ecoli_bsi_arg_presence_region,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.ai Sensitivity analysis of priors ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# regional tau
ecoli_arg_regional_tau_sensitivity <- run_sensitivity_regional_subiso_resumable(
  ecoli_bsi_arg_presence_region,
  feature_cols = gene_cols,
  alpha_grid = c(0.1, 1, 10),
  beta_grid = c(0.1, 1, 10, 100),
  beta_novel_num_grid = c(10, 100, 1000),
  mu_tau_mean_grid = c(0.1, 1, 10),
  mu_tau_sd_grid = c(0.1, 1, 10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/ecoli_arg_beta_binomial_regional_checkpoint_novel.rds")

# Top 5 and run with 10 different held-out splits
ecoli_arg_regional_tau_sensitivity <- flag_valid_fit(ecoli_arg_regional_tau_sensitivity)

# save results
write.csv(ecoli_arg_regional_tau_sensitivity, "sens_results/ecoli_arg_regional_tau_sensitivity_results_fixed_holdout_lps.csv", row.names = FALSE)
#ecoli_arg_regional_tau_sensitivity <- read.csv("sens_results/ecoli_arg_regional_tau_sensitivity_results_fixed_holdout_lps.csv")
#View(ecoli_arg_regional_tau_sensitivity)

# truncate extreme low-fit values
lower <- quantile(ecoli_arg_regional_tau_sensitivity$holdout_log_score_per_isolate, 0.025, na.rm = TRUE)
lower # -1.015
ecoli_arg_regional_tau_sensitivity <- ecoli_arg_regional_tau_sensitivity %>%
  mutate(fill_val = pmax(holdout_log_score_per_isolate, lower))


ecoli_arg_regional_tau_sensitivity_neat <- ecoli_arg_regional_tau_sensitivity |>
  filter(alpha %in% c(0.1, 1, 10),
         beta %in% c(0.1, 1, 10, 100),
         beta_novel_num %in% c(10, 100, 1000),
         mu_tau_mean %in% c(0.1, 1, 10),
         mu_tau_sd %in% c(0.1, 1, 10),
         sigma_tau_rate %in% c(0.1, 1, 10))
#View(ecoli_arg_regional_tau_sensitivity_neat)

#plot regional heatmap
ecoli_arg_regional_tau_plot <- plot_regional_holdout_heatmap_subiso(ecoli_arg_regional_tau_sensitivity_neat, value_col = "fill_val")
ecoli_arg_regional_tau_plot
ggsave("sens_results/ecoli_arg_regional_tau_sensitivity_analysis_plot.png", ecoli_arg_regional_tau_plot, units = "in", width = 10, height = 12, dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.aii Run models with single prior parameter combo ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set priors
df <- ecoli_bsi_arg_df |>
  group_by(Element.symbol) |>
  summarise(count = n()) |>
  arrange(count)
N_isolate <- length(unique(ecoli_bsi_arg_df$sample)) # number of genes
N_gene <- sum(df$count) # number of genes
K <- length(unique(ecoli_bsi_arg_df$Element.symbol)) # 189 # num unique features
f1 <- sum(df$count == 1) # 61 - no. singletons
u_hat <- f1 / N_isolate  # 0.041 - frequency of singletons
u_hat

beta_prior <- 1
alpha_prior <- 1
# set novel priors to preserve 'novel' mass, without and with including singletons 
beta_novel_num = ceiling((u_hat * K)/ (1 - u_hat)) # 1 - derived to make the prior novel frequency the same as in the data
beta_novel_sum = beta_novel_num * beta_prior
beta_novel_num
beta_novel_sum

# prep data
gene_cols <- setdiff(names(ecoli_bsi_arg_presence_region), c("sample", "region"))
prep_arg <- prep_subiso_data(
  df = ecoli_bsi_arg_presence_region,
  region_col = "region",
  isolate_id_col = "sample",
  beta_novel_sum = beta_novel_sum,
  beta_novel_num = beta_novel_num,
  feature_cols = gene_cols
)


#View(prep_arg$y)
fit_r_arg <- fit_region_tau_subiso(prep_arg, 
                                   alpha_prior =  c(rep(alpha_prior, length(gene_cols))),
                                   beta_prior = c(rep(beta_prior, length(gene_cols))), # use uniform beta prior for pi
                                   beta_novel_sum = beta_novel_sum,
                                   beta_novel_num = beta_novel_num,
                                   mu_tau_mean = 1,
                                   mu_tau_sd = 1,
                                   sigma_tau_rate = 1,
                                   iter_warmup = 1000, 
                                   iter_sampling = 1000, 
                                   chains = 4,
                                   seed = 2025)

save_fit_bundle(fit_r_arg, dir = "model_results", basename = "ecoli_bsi_arg_region_tau")
#~~~~~~~~~~~~~~~~#
# read back in
#obj_r <- load_fit_bundle("model_results", "ecoli_bsi_arg_region_tau")
#fit_r_arg <- list(fit = obj_r$fit, prep = obj_r$prep, model = obj_r$meta$model,
# stan_file = obj_r$meta$stan_file, exe_file = obj_r$meta$exe_file)
#~~~~~~~~~~~~~~~~#
# check model diagnostics
#print(fit_r_arg$fit$summary(), n=300)
mean(fit_r_arg$fit$summary()$rhat)
#range(fit_r_arg$fit$summary()$rhat)
mean(fit_r_arg$fit$summary()$ess_bulk)
#range(fit_r_arg$fit$summary()$ess_bulk)
mean(fit_r_arg$fit$summary()$ess_tail)
#range(fit_r_arg$fit$summary()$ess_tail)
#fit_r_arg$fit$cmdstan_summary()
fit_r_arg$fit$cmdstan_diagnose()
fit_r_arg$fit$diagnostic_summary()
mean(fit_r_arg$fit$diagnostic_summary()$ebfmi)
fit_r_arg$fit$sampler_diagnostics()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative richness and mass curves
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

# cumulative species richness curve
curve_arg_richness <- compute_subiso_richness_curves_ms(
  fit_obj = fit_r_arg,
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_arg_richness)
#colnames(curve_arg_richness)

#save
write.csv(curve_arg_richness, "model_results/ecoli_arg_subiso_richness_curves.csv", row.names = FALSE)
#curve_arg_richness <- read.csv("model_results/ecoli_arg_subiso_richness_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~#
# cumulative mass grid
curve_arg <- compute_subiso_mass_curves_ms(
  fit_obj = fit_r_arg,
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_arg)
#colnames(curve_arg)
#table(curve_arg$region)

write.csv(curve_arg, "model_results/ecoli_arg_subiso_mass_curves.csv", row.names = FALSE)
#curve_arg <- read.csv("model_results/ecoli_arg_subiso_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.aiii Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot actual vs modelled global frequencies
global_obs_vs_pred <- plot_global_freq_vs_pi_subiso_r(fit_r = fit_r_arg)
global_obs_vs_pred
#save
ggsave("model_results/ecoli_bsi_arg_global_obs_vs_est_f_logscale.png", 
       global_obs_vs_pred, units = "in", width = 6, height = 4.5, dpi = 300)

# plot actual vs modelled global frequencies
regional_obs_vs_pred <- plot_regional_freq_vs_theta_subiso_r(fit_r = fit_r_arg)
regional_obs_vs_pred
#save
ggsave("model_results/ecoli_bsi_arg_regional_obs_vs_est_f_logscale.png", 
       regional_obs_vs_pred, units = "in", width = 10, height = 7, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~#
# make individual plots
# species_richness:
p17 <- make_mass_plot_region(
  curve_arg_richness,
  xlab = "Regional ARG frequency (f)",
  ylab = "Proportion of unique ARGs\nwith regional frequency ≥f",
  xlim_min = min(curve_arg_richness$f[curve_arg_richness$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p17
ggsave("model_results/ecoli_bsi_ARG_regional_cumulative_richness_curve_regionl_tau_bayes.png", p17, units = "in", width = 6, height = 4, dpi = 300)


p18 <- make_cov_plot_region(curve_arg_richness, 
                           pal = region_pal,
                           ylab = "Proportion of unique\nARGs detected", 
                           xlab = "Sample size",
                           x_limit = c(10,10000))
p18
# save
ggsave("model_results/ecoli_bsi_ARG_regional_richness_vs_ss_regionl_tau_bayes.png", p18, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~#
# coverage:
p19 <- make_mass_plot_region(
  curve_arg,
  xlab = "Regional ARG frequency (f)",
  ylab = "Proportion of ARG population\nwith allele of frequency ≥f",
  xlim_min = min(curve_arg$f[curve_arg$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p19
ggsave("model_results/ecoli_bsi_ARG_regional_cumulative_mass_curve_regionl_tau_bayes.png", p19, units = "in", width = 6, height = 4, dpi = 300)


p20 <- make_cov_plot_region(curve_arg, pal = region_pal)
p20
# save
ggsave("model_results/ecoli_bsi_ARG_regional_samplecoverage_vs_ss_regionl_tau_bayes.png", p20, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.aiv Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_arg_richness |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_arg_richness, national_df)

# save
write.csv(curve_df2, "model_results/ecoli_arg_subiso_richness_curves_global.csv", row.names = FALSE)
#curve_df2 <- read.csv("model_results/ecoli_arg_subiso_richness_curves.csv")


curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", 
             "min_sample_95", 
             "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
                                 confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
                                 confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
                                 TRUE ~ NA_character_))
#View(curve_df_long)
#table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_arg_regional_bhm_summary_richness <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_arg_regional_bhm_summary_richness <- ecoli_bsi_arg_regional_bhm_summary_richness |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_arg_regional_bhm_summary_richness)
#View(ecoli_bsi_arg_regional_bhm_summary_richness)

# save
write.csv(ecoli_bsi_arg_regional_bhm_summary_richness, "model_results/ecoli_bsi_arg_regional_bhm_summary_richness.csv", row.names = FALSE)
#write.csv(ecoli_bsi_arg_regional_bhm_summary_richness, "model_results/ecoli_bsi_arg_regional_bhm_summary_richness.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
write.csv(curve_df2, "model_results/ecoli_arg_subiso_mass_curves_global.csv", row.names = FALSE)
#curve_df2 <- read.csv("model_results/ecoli_arg_subiso_mass_curves.csv")


curve_df_long <- curve_df2 |>
  pivot_longer(
    cols = c("min_sample_90", "min_sample_95", "min_sample_99"),
    names_to = "confidence_level",
    #names_sep = "_",
    values_to = "sample_size"
  ) |>
  mutate(region_conf = case_when(
    confidence_level == "min_sample_90" ~ paste0(region, " (90% conf.)"),
    confidence_level == "min_sample_95" ~ paste0(region, " (95% conf.)"),
    confidence_level == "min_sample_99" ~ paste0(region, " (99% conf)"),
    TRUE ~ NA_character_))
#View(curve_df_long)
#table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)


# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_arg_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_arg_regional_bhm_summary_cells_only <- ecoli_bsi_arg_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_arg_regional_bhm_summary)
#View(ecoli_bsi_arg_regional_bhm_summary)
# save
write.csv(ecoli_bsi_arg_regional_bhm_summary, "model_results/ecoli_bsi_arg_regional_bhm_summary.csv", row.names = FALSE)
#write.csv(ecoli_bsi_arg_regional_bhm_summary, "model_results/ecoli_bsi_arg_regional_bhm_summary.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 5.1b Klebsiella ARGs - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare data
kleb_bsi_amrfinder_metadata <- read.csv("neksus_kleb_bsi_amrfinder_metadata.csv")

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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1bi Sensitivity analysis of priors ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~#
# regional tau
kleb_arg_regional_tau_sensitivity <- run_sensitivity_regional_subiso_resumable(
  kleb_bsi_arg_presence_region,
  feature_cols = gene_cols,
  alpha_grid = c(0.1, 1, 10),
  beta_grid = c(0.1, 1, 10, 100),
  beta_novel_num_grid = c(10, 100, 1000),
  mu_tau_mean_grid = c(0.1, 1, 10),
  mu_tau_sd_grid = c(0.1, 1, 10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/kleb_arg_beta_binomial_regional_checkpoint_with_diagnostics_novel.rds")

# flag valid splits
kleb_arg_regional_tau_sensitivity <- flag_valid_fit(kleb_arg_regional_tau_sensitivity)

# save results
write.csv(kleb_arg_regional_tau_sensitivity, "sens_results/kleb_arg_regional_tau_sensitivity_results_fixed_holdout_lps.csv", row.names = FALSE)
#kleb_arg_regional_tau_sensitivity <- read.csv("sens_results/kleb_arg_regional_tau_sensitivity_results_fixed_holdout_lps.csv")
#View(kleb_arg_regional_tau_sensitivity)

# truncate extreme low-fit/high-fit values
lower <- quantile(kleb_arg_regional_tau_sensitivity$holdout_log_score_per_isolate, 0.025, na.rm = TRUE)
lower # -2.57
upper <- quantile(kleb_arg_regional_tau_sensitivity$holdout_log_score_per_isolate, 0.975, na.rm = TRUE)
upper # -0.95
kleb_arg_regional_tau_sensitivity <- kleb_arg_regional_tau_sensitivity |>
  mutate(fill_val = pmax(holdout_log_score_per_isolate, lower),
         fill_val = pmin(fill_val, upper))
         
# filter if needed
kleb_arg_regional_tau_sensitivity_neat <- kleb_arg_regional_tau_sensitivity |>
  filter(alpha %in% c(0.1, 1, 10),
         beta %in% c(0.1, 1, 10, 100),
         beta_novel_num %in% c(10, 100, 1000),
         mu_tau_mean %in% c(0.1, 1, 10),
         mu_tau_sd %in% c(0.1, 1, 10),
         sigma_tau_rate %in% c(0.1, 1, 10))

#plot regional heatmap
kleb_arg_regional_tau_plot <- plot_regional_holdout_heatmap_subiso(kleb_arg_regional_tau_sensitivity_neat, value_col = "fill_val")
kleb_arg_regional_tau_plot
ggsave("sens_results/kleb_arg_regional_tau_sensitivity_analysis_plot.png", kleb_arg_regional_tau_plot, units = "in", width = 10, height = 12, dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1bii Run with single prior parameter combo ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set priors
df <- kleb_bsi_arg_df |>
  group_by(Element.symbol) |>
  summarise(count = n()) |>
  arrange(count)
N_isolate <- length(unique(kleb_bsi_arg_df$sample)) # 468
N_gene <- sum(df$count) # number of genes 3258
K <- length(unique(kleb_bsi_arg_df$Element.symbol)) # num unique features 199
f1 <- sum(df$count == 1) # 64 - no. singletons
u_hat <- f1 / N_isolate # frequency of singletons
u_hat

beta_prior <- 1
alpha_prior <- 1
# set novel priors to preserve 'novel' mass, without and with including singletons 
beta_novel_num = ceiling((u_hat * K)/ (1 - u_hat)) # derived to make the prior novel frequency the same as in the data
beta_novel_sum = beta_novel_num * beta_prior # 4
beta_novel_sum
beta_novel_num

gene_cols <- setdiff(names(kleb_bsi_arg_presence_region), c("sample", "region"))
prep_arg <- prep_subiso_data(
  df = kleb_bsi_arg_presence_region,
  region_col = "region",
  isolate_id_col = "sample",
  beta_novel_sum = beta_novel_sum,
  beta_novel_num = beta_novel_num,
  feature_cols = gene_cols)



# fit model
fit_r_arg <- fit_region_tau_subiso(prep_arg, 
                                   alpha_prior =  c(rep(beta_prior, length(gene_cols))),
                                   beta_prior = c(rep(beta_prior, length(gene_cols))), # use uniform beta prior for pi
                                   beta_novel_sum = beta_novel_sum,
                                   beta_novel_num = beta_novel_num,
                                   mu_tau_mean = 1,
                                   mu_tau_sd = 1,
                                   sigma_tau_rate = 1,
                                   iter_warmup = 1000, 
                                   iter_sampling = 1000, 
                                   chains = 4,
                                   seed = 2025)
#save
save_fit_bundle(fit_r_arg, dir = "model_results", basename = "kleb_bsi_arg_region_tau")
#~~~~~~~~~~~~~~~~#
# read back in
#obj_r <- load_fit_bundle("model_results", "kleb_bsi_arg_region_tau")
#fit_r_arg <- list(fit = obj_r$fit, prep = obj_r$prep, model = obj_r$meta$model,
# stan_file = obj_r$meta$stan_file, exe_file = obj_r$meta$exe_file)
#~~~~~~~~~~~~~~~~#
# check diagnostics
#print(fit_r_arg$fit$summary(), n=300)
mean(fit_r_arg$fit$summary()$rhat)
#range(fit_r_arg$fit$summary()$rhat)
mean(fit_r_arg$fit$summary()$ess_bulk)
#range(fit_r_arg$fit$summary()$ess_bulk)
mean(fit_r_arg$fit$summary()$ess_tail)
#range(fit_r_arg$fit$summary()$ess_tail)
#fit_r_arg$fit$cmdstan_summary()
fit_r_arg$fit$cmdstan_diagnose()
fit_r_arg$fit$diagnostic_summary()
mean(fit_r_arg$fit$diagnostic_summary()$ebfmi)
fit_r_arg$fit$sampler_diagnostics()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative species richness curve
curve_arg_richness <- compute_subiso_richness_curves_ms(
  fit_obj = fit_r_arg,
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_arg_richness)
#colnames(curve_arg_richness)

#save
write.csv(curve_arg_richness, "model_results/kleb_arg_subiso_richness_curves.csv", row.names = FALSE)
#curve_arg_richness <- read.csv("model_results/kleb_arg_subiso_richness_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative mass grid
curve_arg <- compute_subiso_mass_curves_ms(
  fit_obj = fit_r_arg,
  f_grid = f_grid_99,
  n_post_draws = 400)

#View(curve_arg)
#colnames(curve_arg)
#table(curve_arg$region)
#save

write.csv(curve_arg, "model_results/kleb_arg_subiso_mass_curves.csv", row.names = FALSE)
#curve_arg <- read.csv("model_results/kleb_arg_subiso_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.bi Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot actual vs modelled global frequencies
global_obs_vs_pred <- plot_global_freq_vs_pi_subiso_r(fit_r = fit_r_arg)
global_obs_vs_pred
# save
ggsave("model_results/kleb_bsi_arg_global_obs_vs_est_f_logscale.png", 
       global_obs_vs_pred, units = "in", width = 6, height = 4.5, dpi = 300)

# plot actual vs modelled global frequencies
regional_obs_vs_pred <- plot_regional_freq_vs_theta_subiso_r(fit_r = fit_r_arg)
regional_obs_vs_pred
#save
ggsave("model_results/kleb_bsi_arg_regional_obs_vs_est_f_logscale.png", 
       regional_obs_vs_pred, units = "in", width = 10, height = 7, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~#
# make individual plots
# species_richness:
p21 <- make_mass_plot_region(
  curve_arg_richness,
  xlab = "Regional ARG frequency (f)",
  ylab = "Proportion of unique ARGs\nwith regional frequency ≥f",
  xlim_min = min(curve_arg_richness$f[curve_arg_richness$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p21
ggsave("model_results/kleb_bsi_ARG_regional_cumulative_richness_curve_regionl_tau_bayes.png", p21, units = "in", width = 6, height = 4, dpi = 300)


p22 <- make_cov_plot_region(curve_arg_richness, 
                           pal = region_pal,
                           ylab = "Proportion of unique\nARGs detected", 
                           xlab = "Sample size",
                           x_limit = c(10,10000))
p22
# save
ggsave("model_results/kleb_bsi_ARG_regional_richness_vs_ss_regionl_tau_bayes.png", p22, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~#
# coverage:
p23 <- make_mass_plot_region(
  curve_arg,
  xlab = "ARG frequency (f)",
  ylab = "Proportion of ARG population\nwith allele of frequency ≥f",
  xlim_min = min(curve_arg$f[curve_arg$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p23
ggsave("model_results/kleb_bsi_ARG_regional_cumulative_mass_curve_regionl_tau_bayes.png", p23, units = "in", width = 6, height = 4, dpi = 300)

p24 <- make_cov_plot_region(curve_arg, pal = region_pal)
p24
# save
ggsave("model_results/kleb_bsi_ARG_regional_samplecoverage_vs_ss_regionl_tau_bayes.png", p24, units = "in", width = 6, height = 4, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.1.bii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_arg_richness |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_arg_richness, national_df)

# save
write.csv(curve_df2, "model_results/kleb_arg_subiso_richness_curves_global.csv", row.names = FALSE)
#curve_df2 <- read.csv("model_results/kleb_arg_subiso_richness_curves.csv")


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
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_arg_regional_bhm_summary_richness <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_arg_regional_bhm_summary_richness <- kleb_bsi_arg_regional_bhm_summary_richness |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_arg_regional_bhm_summary_richness)
#View(kleb_bsi_arg_regional_bhm_summary_richness)

# save
write.csv(kleb_bsi_arg_regional_bhm_summary_richness, "model_results/kleb_bsi_arg_regional_bhm_summary_richness.csv", row.names = FALSE)
#write.csv(kleb_bsi_arg_regional_bhm_summary_richness, "model_results/kleb_bsi_arg_regional_bhm_summary_richness.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
write.csv(curve_df2, "kleb_arg_subiso_mass_curves_global.csv", row.names = FALSE)
#curve_df2 <- read.csv("kleb_arg_subiso_mass_curves.csv")


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
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_arg_regional_bhm_summary_cells_only <- kleb_bsi_arg_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_arg_regional_bhm_summary)
#View(kleb_bsi_arg_regional_bhm_summary)

# save
write.csv(kleb_bsi_arg_regional_bhm_summary, "model_results/kleb_bsi_arg_regional_bhm_summary.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 5.2a E. coli Plasmids - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare E. coli plasmids df
# master sample list (one row per sample)
all_samples <- ecoli_bsi_amrfinder_metadata |> distinct(sample)   

# detected pairs from PLING (1 if detected)
detected <- ecoli_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)
#View(detected)
nrow(detected)

# full grid of sample x feature using master sample list and the set of observed features
all_features <- detected |> pull(community_subcommunity) |> unique()
full_grid <- tidyr::expand_grid(sample = all_samples$sample,
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

counts <- presence_long |>
  group_by(community_subcommunity) |>
  summarise(count = sum(presence), .groups = "drop")
#table(counts$count)
sum(counts[counts$count ==1,]$count)

# add region data
ecoli_bsi_pling_df <- ecoli_bsi_pling_df |>
  left_join(ecoli_bsi_samples_metadata |> select(sequencing_id, region), by = c("sample" = "sequencing_id")) |>
  select(sample, region, everything())

#View(ecoli_bsi_pling_df)
length(unique(ecoli_bsi_pling_df$sample))# 1471
table(ecoli_bsi_pling_df$region)# 1471
sum(ecoli_bsi_pling_df[,3:827])# 3361

iso_df <- ecoli_bsi_pling_df
iso_df <- iso_df |> dplyr::ungroup()
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))

gene_cols <- setdiff(names(ecoli_bsi_pling_df), c("sample", "region"))
prep_pling <- prep_subiso_data(
  df = ecoli_bsi_pling_df,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2ai Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# regional tau #972
ecoli_pling_regional_tau_sensitivity <- run_sensitivity_regional_subiso_resumable(
  ecoli_bsi_pling_df,
  feature_cols = gene_cols,
  alpha_grid = c(0.1, 1, 10),
  beta_grid = c(0.1, 1, 10, 100),
  beta_novel_num_grid = c(10, 100, 1000),
  mu_tau_mean_grid = c(0.1, 1, 10),
  mu_tau_sd_grid = c(0.1, 1, 10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel.rds")


# a0.1_b100 = 81 combos
ecoli_pling_regional_tau_sensitivity_a0.1_b100 <- run_sensitivity_regional_subiso_resumable(
  ecoli_bsi_pling_df,
  feature_cols = gene_cols,
  alpha_grid = c(0.1),
  beta_grid = c(100),
  beta_novel_num_grid = c(10, 100, 1000),
  mu_tau_mean_grid = c(0.1, 1, 10),
  mu_tau_sd_grid = c(0.1, 1, 10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a0.1_b100.rds")



ecoli_pling_regional_tau_sensitivity_a0.1_b100_m1000 <- run_sensitivity_regional_subiso_resumable(
  ecoli_bsi_pling_df,
  feature_cols = gene_cols,
  alpha_grid = c(0.1),
  beta_grid = c(100),
  beta_novel_num_grid = c(1000),
  mu_tau_mean_grid = c(0.1, 1, 10),
  mu_tau_sd_grid = c(0.1, 1, 10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a0.1_b100_m1000.rds")

ecoli_pling_regional_tau_sensitivity_a0.1_b100_m1000_mu10 <- run_sensitivity_regional_subiso_resumable(
  ecoli_bsi_pling_df,
  feature_cols = gene_cols,
  alpha_grid = c(0.1),
  beta_grid = c(100),
  beta_novel_num_grid = c(1000),
  mu_tau_mean_grid = c(10),
  mu_tau_sd_grid = c(10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a0.1_b100_m1000_mu10.rds")


ecoli_pling_regional_tau_sensitivity_a1 <- run_sensitivity_regional_subiso_resumable(
  ecoli_bsi_pling_df,
  feature_cols = gene_cols,
  alpha_grid = c( 1),
  beta_grid = c(0.1, 1, 10, 100),
  beta_novel_num_grid = c(10, 100, 1000),
  mu_tau_mean_grid = c(0.1, 1, 10),
  mu_tau_sd_grid = c(0.1, 1, 10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a1.rds")


ecoli_pling_regional_tau_sensitivity_a10 <- run_sensitivity_regional_subiso_resumable(
  ecoli_bsi_pling_df,
  feature_cols = gene_cols,
  alpha_grid = c(10),
  beta_grid = c(0.1, 1, 10, 100),
  beta_novel_num_grid = c(10, 100, 1000),
  mu_tau_mean_grid = c(0.1, 1, 10),
  mu_tau_sd_grid = c(0.1, 1, 10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/ecoli_pling_beta_binomial_regional_checkpoint_novel_a10.rds")


# flag valid hits
ecoli_pling_regional_tau_sensitivity <- flag_valid_fit(ecoli_pling_regional_tau_sensitivity)

# save results
write.csv(ecoli_pling_regional_tau_sensitivity, "sens_results/ecoli_pling_regional_tau_sensitivity_results_fixed_holdout_lps.csv", row.names = FALSE)
#ecoli_pling_regional_tau_sensitivity <- read.csv("sens_results/ecoli_pling_regional_tau_sensitivity_results_fixed_holdout_lps.csv")
#View(ecoli_pling_regional_tau_sensitivity)

# skip truncation to truncate extreme high/low-fit values
#upper <- quantile(ecoli_pling_regional_tau_sensitivity$holdout_log_score_per_isolate, 0.975, na.rm = TRUE)
#upper # -2.15
#upper <- sort(unique(ecoli_pling_regional_tau_sensitivity$holdout_log_score_per_isolate), decreasing = TRUE)[2]
#upper # -2.09
#ecoli_pling_regional_tau_sensitivity <- ecoli_pling_regional_tau_sensitivity |>
 # mutate(fill_val = pmin(holdout_log_score_per_isolate, upper))


# filter
ecoli_pling_regional_tau_sensitivity_neat <- ecoli_pling_regional_tau_sensitivity |>
  filter(alpha %in% c(0.1, 1, 10),
         beta %in% c(0.1, 1, 10, 100),
         beta_novel_num %in% c(10, 100, 1000),
         mu_tau_mean %in% c(0.1, 1, 10),
         mu_tau_sd %in% c(0.1, 1, 10),
         sigma_tau_rate %in% c(0.1, 1, 10))

#plot regional heatmap
ecoli_pling_regional_tau_plot <- plot_regional_holdout_heatmap_subiso(ecoli_pling_regional_tau_sensitivity_neat)
ecoli_pling_regional_tau_plot
ggsave("sens_results/ecoli_pling_regional_tau_sensitivity_analysis_plot.png", ecoli_pling_regional_tau_plot, units = "in", width = 10, height = 12, dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2aii Run with single prior parameter combo ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set priors
N_gene <- sum(presence_long$presence) # 3361 number of plasmids
N_isolate <- length(unique(ecoli_bsi_pling_df$sample)) # number of genes
K <- length(unique(presence_long$community_subcommunity)) # 825 num unique features
f1 <- sum(counts$count == 1) #  622- no. singletons
u_hat <- f1 / N_isolate  #  mass of singletons
u_hat # 

# set priors
alpha_prior <- 1
beta_prior <- 50 # set slightly higher to have right-skewed prior (more low frequency) to improve model mixing as many plasmid subcommunities very rare
# set novel priors to preserve 'novel' mass, without and with including singletons 
beta_novel_num = ceiling((u_hat * K)/ (1 - u_hat)) # 188 - derived to make the prior novel frequency the same as in the data
beta_novel_sum = beta_novel_num * beta_prior 
beta_novel_num
beta_novel_sum
  
# prep data
gene_cols <- setdiff(names(ecoli_bsi_pling_df), c("sample", "region"))
prep_pling <- prep_subiso_data(
  df = ecoli_bsi_pling_df,
  region_col = "region",
  isolate_id_col = "sample",
  beta_novel_sum = beta_novel_sum,
  beta_novel_num = beta_novel_num,
  feature_cols = gene_cols
)

# fit model
fit_r_pling <- fit_region_tau_subiso(prep_pling, 
                                     alpha_prior = c(rep(alpha_prior, length(gene_cols))), # use uniform beta prior for pi
                                     beta_prior = c(rep(beta_prior, length(gene_cols))), # use uniform beta prior for pi
                                     mu_tau_mean = 1, 
                                     mu_tau_sd = 1, 
                                     sigma_tau_rate = 1, 
                                     beta_novel_sum = beta_novel_sum,
                                     beta_novel_num = beta_novel_num,
                                     iter_warmup = 1000, 
                                     iter_sampling = 1000, 
                                     chains = 4,
                                     seed = 2025)

save_fit_bundle(fit_r_pling, dir = "model_results", basename = "ecoli_bsi_pling_region_tau")
#~~~~~~~~~~~~~~~~#
# read back in
#obj_r <- load_fit_bundle("model_results", "ecoli_bsi_pling_region_tau")
#fit_r_pling <- list(fit = obj_r$fit, prep = obj_r$prep, model = obj_r$meta$model,
# stan_file = obj_r$meta$stan_file, exe_file = obj_r$meta$exe_file)
#~~~~~~~~~~~~~~~~#
# check diagnostics
# model summaries and diagnostics
#print(fit_r_pling$fit$summary(), n=300)
mean(fit_r_pling$fit$summary()$rhat)
mean(fit_r_pling$fit$summary()$ess_bulk)
mean(fit_r_pling$fit$summary()$ess_tail)
#fit_r_pling$fit$cmdstan_summary()
fit_r_pling$fit$cmdstan_diagnose()
fit_r_pling$fit$diagnostic_summary()
mean(fit_r_pling$fit$diagnostic_summary()$ebfmi)
fit_r_pling$fit$sampler_diagnostics()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative species richness curve
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

curve_pling_richness <- compute_subiso_richness_curves_ms(
  fit_obj = fit_r_pling,
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_pling_richness)
#colnames(curve_pling_richness)

#save
write.csv(curve_pling_richness, "model_results/ecoli_pling_subiso_richness_curves.csv", row.names = FALSE)
#curve_pling_richness <- read.csv("model_results/ecoli_pling_subiso_richness_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# calculate cumulative mass grid
curve_pling <- compute_subiso_mass_curves_ms(
  fit_obj = fit_r_pling,
  f_grid = f_grid_99,
  n_post_draws = 400)

#View(curve_pling)
#colnames(curve_pling)
#table(curve_pling$region)
#save
write.csv(curve_pling, "model_results/ecoli_pling_subiso_mass_curves.csv", row.names = FALSE)
#curve_pling <- read.csv("model_results/ecoli_pling_subiso_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2.ai Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot actual vs modelled global frequencies
global_obs_vs_pred <- plot_global_freq_vs_pi_subiso_r(fit_r = fit_r_pling)
global_obs_vs_pred
#save
ggsave("model_results/ecoli_bsi_pling_global_obs_vs_est_f_logscale.png", 
       global_obs_vs_pred, units = "in", width = 6, height = 4.5, dpi = 300)

# plot actual vs modelled global frequencies
regional_obs_vs_pred <- plot_regional_freq_vs_theta_subiso_r(fit_r = fit_r_pling)
regional_obs_vs_pred
#save
ggsave("model_results/ecoli_bsi_pling_regional_obs_vs_est_f_logscale.png", 
       regional_obs_vs_pred, units = "in", width = 10, height = 7, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~#
# make individual plots
# species_richness:
p25 <- make_mass_plot_region(
  curve_pling_richness,
  xlab = "Regional plasmid subcommunity frequency (f)",
  ylab = "Proportion of unique plasmid subcommunities\nwith regional frequency ≥f",
  xlim_min = min(curve_pling_richness$f[curve_pling_richness$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p25
ggsave("model_results/ecoli_bsi_pling_regional_cumulative_richness_curve_regionl_tau_bayes.png", p25, units = "in", width = 6, height = 4, dpi = 300)


p26 <- make_cov_plot_region(curve_pling_richness, 
                           pal = region_pal,
                           ylab = "Proportion of unique\nplasmid subcommunities detected", 
                           xlab = "Sample size",
                           x_limit = c(10,10000))
p26
# save
ggsave("model_results/ecoli_bsi_pling_regional_richness_vs_ss_regionl_tau_bayes.png", p26, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~#
# coverage:
p27 <- make_mass_plot_region(
  curve_pling,
  xlab = "Plasmid subcommunity frequency (f)",
  ylab = "Proportion of plasmid population\nin subcommunity of frequency ≥f",
  xlim_min = min(curve_pling$f[curve_pling$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p27
ggsave("model_results/ecoli_bsi_pling_regional_cumulative_mass_curve_regionl_tau_bayes.png", p27, units = "in", width = 6, height = 4, dpi = 300)

# make samaple size vs coverage plot
p28 <- make_cov_plot_region(curve_pling, pal = region_pal)
p28
# save
ggsave("model_results/ecoli_bsi_pling_regional_samplecoverage_vs_ss_regionl_tau_bayes.png", p28, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2.aii Summary tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_pling_richness |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_pling_richness, national_df)

# save
write.csv(curve_df2, "model_results/ecoli_pling_subiso_richness_curves_global.csv", row.names = FALSE)
#curve_df2 <- read.csv("model_results/ecoli_pling_subiso_richness_curves.csv")


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
#table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_pling_regional_bhm_summary_richness <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_pling_regional_bhm_summary_richness <- ecoli_bsi_pling_regional_bhm_summary_richness |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_pling_regional_bhm_summary_richness)
#View(ecoli_bsi_pling_regional_bhm_summary_richness)

# save
write.csv(ecoli_bsi_pling_regional_bhm_summary_richness, "model_results/ecoli_bsi_pling_regional_bhm_summary_richness.csv", row.names = FALSE)
#write.csv(ecoli_bsi_pling_regional_bhm_summary_richness, "model_results/ecoli_bsi_pling_regional_bhm_summary_richness.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

# save
write.csv(curve_df2, "model_results/ecoli_pling_subiso_mass_curves_global.csv", row.names = FALSE)
#curve_df2 <- read.csv("model_results/ecoli_pling_subiso_mass_curves.csv")


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
#table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
ecoli_bsi_pling_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_pling_regional_bhm_summary_cells_only <- ecoli_bsi_pling_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_pling_regional_bhm_summary)
#View(ecoli_bsi_pling_regional_bhm_summary)

# save
write.csv(ecoli_bsi_pling_regional_bhm_summary, "model_results/ecoli_bsi_pling_regional_bhm_summary.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 5.2b Klebsiella plasmids - by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare Kleb plasmids df
# master sample list (one row per sample)
all_samples <- kleb_bsi_amrfinder_metadata |> distinct(sample)   

# detected pairs from AMRFinder (1 if detected)
detected <- kleb_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)
nrow(detected)

# full grid of sample x feature using master sample list and the set of observed features
all_features <- detected |> pull(community_subcommunity) |> unique()
full_grid <- tidyr::expand_grid(sample = all_samples$sample,
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
counts <- presence_long |>
  group_by(community_subcommunity) |>
  summarise(count = sum(presence), .groups = "drop")
#table(counts$count)
sum(counts[counts$count ==1,]$count)

# add region data
kleb_bsi_pling_df <- kleb_bsi_pling_df |>
  left_join(kleb_bsi_samples_metadata |> select(sequencing_id, region), by = c("sample" = "sequencing_id")) |>
  select(sample, region, everything())

#View(kleb_bsi_pling_df)
length(unique(kleb_bsi_pling_df$sample))# 468
table(kleb_bsi_pling_df$region)# 

gene_cols <- setdiff(names(kleb_bsi_pling_df), c("sample", "region"))
prep_pling <- prep_subiso_data(
  df = kleb_bsi_pling_df,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2bi Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# regional tau
kleb_pling_regional_tau_sensitivity <- run_sensitivity_regional_subiso_resumable(
  kleb_bsi_pling_df,
  feature_cols = gene_cols,
  alpha_grid = c(0.1, 1, 10),
  beta_grid = c(0.1, 1, 10, 100),
  beta_novel_num_grid = c(10, 100, 1000),
  mu_tau_mean_grid = c(0.1, 1, 10),
  mu_tau_sd_grid = c(0.1, 1, 10),
  sigma_tau_rate_grid = c(0.1, 1, 10),
  seed_split = 2026,
  iter_warmup = 250,
  iter_sampling = 250,
  chains = 1,
  checkpoint_file = "sens_results/kleb_pling_beta_binomial_regional_checkpoint_with_diagnostics_novel.rds")

# Top 5 and run with 10 different held-out splits
kleb_pling_regional_tau_sensitivity <- flag_valid_fit(kleb_pling_regional_tau_sensitivity)

# save results
write.csv(kleb_pling_regional_tau_sensitivity, "sens_results/kleb_pling_regional_tau_sensitivity_results_fixed_holdout_lps.csv", row.names = FALSE)
#kleb_pling_regional_tau_sensitivity <- read.csv("sens_results/kleb_pling_regional_tau_sensitivity_results_fixed_holdout_lps.csv")
#View(kleb_pling_regional_tau_sensitivity)


# truncate extreme low/high-fit values if invalid
upper <- quantile(kleb_pling_regional_tau_sensitivity$holdout_log_score_per_isolate, 0.975, na.rm = TRUE)
upper # -3.5
kleb_pling_regional_tau_sensitivity <- kleb_pling_regional_tau_sensitivity |>
  mutate(fill_val = pmin(holdout_log_score_per_isolate, upper))

# filter if needed
kleb_pling_regional_tau_sensitivity_neat <- kleb_pling_regional_tau_sensitivity |>
  filter(alpha %in% c(0.1, 1, 10),
         beta %in% c(0.1, 1, 10, 100),
         beta_novel_num %in% c(10, 100, 1000),
         mu_tau_mean %in% c(0.1, 1, 10),
         mu_tau_sd %in% c(0.1, 1, 10),
         sigma_tau_rate %in% c(0.1, 1, 10))

#plot regional heatmap
kleb_pling_regional_tau_plot <- plot_regional_holdout_heatmap_subiso(kleb_pling_regional_tau_sensitivity_neat, value_col = "fill_val")
kleb_pling_regional_tau_plot
ggsave("sens_results/kleb_pling_regional_tau_sensitivity_analysis_plot.png", kleb_pling_regional_tau_plot, units = "in", width = 10, height = 12, dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2bii Run with single prior parameter combo ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set priors
N_isolate <- length(unique(kleb_bsi_pling_df$sample)) # number of genes
N_gene <- sum(presence_long$presence) #  number of plasmids
K <- length(unique(presence_long$community_subcommunity)) # num unique features
f1 <- sum(counts$count == 1) #   no. singletons
u_hat <- f1 / N_isolate  # 0.390 frequency of singletons
u_hat 

# set priors
alpha_prior <- 1
beta_prior <- 50 # set slightly higher to have right-skewed prior (more low frequency) to improve model mixing as many plasmid subcommunities very rare
# set novel priors to preserve 'novel' mass, without and with including singletons 
beta_novel_num = ceiling((u_hat * K)/ (1 - u_hat)) # 267 derived to make the prior novel frequency the same as in the data
beta_novel_sum = beta_novel_num * beta_prior # 
beta_novel_num
beta_novel_sum

# prep 
gene_cols <- setdiff(names(kleb_bsi_pling_df), c("sample", "region"))
prep_pling <- prep_subiso_data(
  df = kleb_bsi_pling_df,
  region_col = "region",
  isolate_id_col = "sample",
  beta_novel_sum = beta_novel_sum,
  beta_novel_num = beta_novel_num,
  feature_cols = gene_cols
)

# fit model
fit_r_pling <- fit_region_tau_subiso(prep_pling, 
                                     alpha_prior = c(rep(alpha_prior, length(gene_cols))),
                                     beta_prior = c(rep(beta_prior, length(gene_cols))),
                                     mu_tau_mean = 1,
                                     mu_tau_sd = 1,
                                     sigma_tau_rate = 1,
                                     beta_novel_sum = beta_novel_sum,
                                     beta_novel_num = beta_novel_num,
                                     iter_warmup = 1000, 
                                     iter_sampling = 1000, 
                                     chains = 4,
                                     seed = 2025)

save_fit_bundle(fit_r_pling, dir = "model_results", basename = "kleb_bsi_pling_region_tau")
#~~~~~~~~~~~~~~~~#
# read back in
#obj_r <- load_fit_bundle("model_results", "kleb_bsi_pling_region_tau")
#fit_r_pling <- list(fit = obj_r$fit, prep = obj_r$prep, model = obj_r$meta$model,
#                    stan_file = obj_r$meta$stan_file, exe_file = obj_r$meta$exe_file)
#~~~~~~~~~~~~~~~~#
# check diagnostics
#print(fit_r_pling$fit$summary(), n=300)
mean(fit_r_pling$fit$summary()$rhat)
#range(fit_r_pling$fit$summary()$rhat)

mean(fit_r_pling$fit$summary()$ess_bulk)
#range(fit_r_pling$fit$summary()$ess_bulk)

mean(fit_r_pling$fit$summary()$ess_tail)
#range(fit_r_pling$fit$summary()$ess_tail)

#fit_r_pling$fit$cmdstan_summary()
fit_r_pling$fit$cmdstan_diagnose()
fit_r_pling$fit$diagnostic_summary()
mean(fit_r_pling$fit$diagnostic_summary()$ebfmi)
fit_r_pling$fit$sampler_diagnostics()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# cumulative richness and mass curves
curve_pling_richness <- compute_subiso_richness_curves_ms(
  fit_obj = fit_r_pling,
  f_grid = f_grid_99,
  n_post_draws = 400)
#View(curve_pling_richness)
colnames(curve_pling_richness)

#save
write.csv(curve_pling_richness, "model_results/kleb_pling_subiso_richness_curves.csv", row.names = FALSE)
#curve_pling_richness <- read.csv("model_results/kleb_pling_subiso_richness_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# call cumulative mass grid
curve_pling <- compute_subiso_mass_curves_ms(
  fit_obj = fit_r_pling,
  f_grid = f_grid_99,
  n_post_draws = 400)

#View(curve_pling)
colnames(curve_pling)
#table(curve_pling$region)

write.csv(curve_pling, "model_results/kleb_pling_subiso_mass_curves.csv", row.names = FALSE)
#curve_pling <- read.csv("model_results/kleb_pling_subiso_mass_curves.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2.bi Plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot actual vs modelled global frequencies
global_obs_vs_pred <- plot_global_freq_vs_pi_subiso_r(fit_r = fit_r_pling)
global_obs_vs_pred
#save
ggsave("model_results/kleb_bsi_pling_global_obs_vs_est_f_logscale.png", 
       global_obs_vs_pred, units = "in", width = 6, height = 4.5, dpi = 300)

# plot actual vs modelled global frequencies
regional_obs_vs_pred <- plot_regional_freq_vs_theta_subiso_r(fit_r = fit_r_pling)
regional_obs_vs_pred
#save
ggsave("model_results/kleb_bsi_pling_regional_obs_vs_est_f_logscale.png", 
       regional_obs_vs_pred, units = "in", width = 10, height = 7, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~#
# make individual plots
# species_richness:
p29 <- make_mass_plot_region(
  curve_pling_richness,
  xlab = "Regional plasmid subcommunity frequency (f)",
  ylab = "Proportion of unique plasmid subcommunities\nwih regional frequency ≥ f",
  xlim_min = min(curve_pling_richness$f[curve_pling_richness$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p29
ggsave("model_results/kleb_bsi_pling_regional_cumulative_richness_curve_regionl_tau_bayes.png", p29, units = "in", width = 6, height = 4, dpi = 300)


p30 <- make_cov_plot_region(curve_pling_richness, 
                           pal = region_pal,
                           ylab = "Proportion of unique\nplasmid subcommunities detected", 
                           xlab = "Sample size",
                           x_limit = c(10,10000))
p30
# save
ggsave("model_results/kleb_bsi_pling_regional_richness_vs_ss_regionl_tau_bayes.png", p30, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~#
# coverage:
p31 <- make_mass_plot_region(
  curve_pling,
  xlab = "Plasmid subcommunity frequency (f)",
  ylab = "Proportion of plasmid population\nin subcommunity of frequency ≥ f",
  xlim_min = min(curve_pling$f[curve_pling$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p31

ggsave("model_results/kleb_bsi_pling_regional_cumulative_mass_curve_regionl_tau_bayes.png", p31, units = "in", width = 6, height = 4, dpi = 300)


p32 <- make_cov_plot_region(curve_pling, pal = region_pal)
p32
# save
ggsave("model_results/kleb_bsi_pling_regional_samplecoverage_vs_ss_regionl_tau_bayes.png", p32, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * 5.2.bii Summary table ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add global cumulative mass and ss rows
national_df <- curve_pling_richness |>
  dplyr::distinct(f, global_all_median, global_all_q2.5, global_all_q97.5, min_sample_90, min_sample_95, min_sample_99) |>
  dplyr::mutate(
    region = "UK National",
    f = f,
    global_median = global_all_median,
    global_q2.5 = global_all_q2.5,
    global_q97.5 = global_all_q97.5
  )

curve_df2 <- dplyr::bind_rows(curve_pling_richness, national_df)

# save
write.csv(curve_df2, "model_results/kleb_pling_subiso_richness_curves_global_r.csv", row.names = FALSE)
#curve_df2 <- read.csv("model_results/kleb_pling_subiso_richness_curves.csv")


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
#table(curve_df_long$region_conf)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_pling_regional_bhm_summary_richness <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_pling_regional_bhm_summary_richness <- kleb_bsi_pling_regional_bhm_summary_richness |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_pling_regional_bhm_summary_richness)
#View(kleb_bsi_pling_regional_bhm_summary_richness)

# save
write.csv(kleb_bsi_pling_regional_bhm_summary_richness, "model_results/kleb_bsi_pling_regional_bhm_summary_richness_r.csv", row.names = FALSE)
#write.csv(kleb_bsi_pling_regional_bhm_summary_richness, "model_results/kleb_bsi_pling_regional_bhm_summary_richness.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# save
write.csv(curve_df2, "model_results/kleb_pling_subiso_mass_curves_global_r.csv", row.names = FALSE)
#curve_df2 <- read.csv("model_results/kleb_pling_subiso_mass_curves.csv")


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
kleb_bsi_pling_regional_bhm_summary <- curve_df_long |>
  group_by(region_conf) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "global_median", t)
      lo_ss    <- min_sample_at_or_above(df, "global_q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "global_q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(region_conf, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_pling_regional_bhm_summary_cells_only <- kleb_bsi_pling_regional_bhm_summary |>
  select(region_conf, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_pling_regional_bhm_summary)
#View(kleb_bsi_pling_regional_bhm_summary)

# save
write.csv(kleb_bsi_pling_regional_bhm_summary, "model_results/kleb_bsi_pling_regional_bhm_summary_r.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. PANEL PLOT FOR ALL BY-REGION HIERARCHICAL MODELS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data - cumulative richness
ecoli_mlst_curve_df <- read.csv("model_results/ecoli_bsi_mlst_richness_curves.csv")
kleb_mlst_curve_df <- read.csv("model_results/kleb_bsi_mlst_richness_curves.csv")
ecoli_fastbaps_L3_curve_df <- read.csv("model_results/ecoli_bsi_fastbaps_L3_richness_curves.csv")
kleb_fastbaps_L3_curve_df <- read.csv("model_results/kleb_bsi_fastbaps_L3_richness_curves_course.csv")
ecoli_pling_curve_df <- read.csv("model_results/ecoli_pling_subiso_richness_curves.csv")
kleb_pling_curve_df <- read.csv("model_results/kleb_pling_subiso_richness_curves.csv")
ecoli_arg_curve_df <- read.csv("model_results/ecoli_arg_subiso_richness_curves.csv")
kleb_arg_curve_df <- read.csv("model_results/kleb_arg_subiso_richness_curves.csv")

# load data - cumulative mass
ecoli_mlst_curve_df <- read.csv("model_results/ecoli_bsi_mlst_mass_curves.csv")
kleb_mlst_curve_df <- read.csv("model_results/kleb_bsi_mlst_mass_curves.csv")
ecoli_fastbaps_L3_curve_df <- read.csv("model_results/ecoli_bsi_fastbaps_L3_mass_curves.csv")
kleb_fastbaps_L3_curve_df <- read.csv("model_results/kleb_bsi_fastbaps_L3_mass_curves.csv")
ecoli_pling_curve_df <- read.csv("model_results/ecoli_pling_subiso_mass_curves.csv")
kleb_pling_curve_df <- read.csv("model_results/kleb_pling_subiso_mass_curves.csv")
ecoli_arg_curve_df <- read.csv("model_results/ecoli_arg_subiso_mass_curves.csv")
kleb_arg_curve_df <- read.csv("model_results/kleb_arg_subiso_mass_curves.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Cumulative mass plots for ecoli and kleb by region ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# E. coli richness and mass -based sample sizes
combined_figure <-
  (p2 | p4) /
  (p10 | p12) /
  (p26 | p28) /
  (p18 | p20) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined_figure
# save
ggsave("model_results/ecoli_4x2_hierarchical_richness_and_mass_ss_curves_by_region.png",
       plot = combined_figure, width = 11, height = 16.5, units = "in", dpi = 300)


# Kleb richness and mass-based sample sizes.
combined_figure <-
  (p6 | p8) /
  (p14 | p16) /
  (p30 | p32) /
  (p22 | p24) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined_figure
# save
ggsave("model_results/kleb_4x2_hierarchical_richness_and_mass_ss_curves_by_region.png",
       plot = combined_figure, width = 11, height = 16.5, units = "in", dpi = 300)


# E.coli and Kleb mass-based sample sizes.
combined_figure <-
  (p4 | p8) /
  (p12 | p16) /
  (p28 | p32) /
  (p20 | p24) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined_figure
# save
ggsave("model_results/ecoli_kleb_4x2_hierarchical_sample_cov_vs_ss_by_region.png",
       plot = combined_figure, width = 11, height = 16.5, units = "in", dpi = 300)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# E.coli and kleb richness then E. coli and Kleb mass sample sizes
# add titles
p2 <- p2 + ggtitle("E. coli")
p6 <- p6 + ggtitle("Klebsiella")
p4 <- p4 + ggtitle("E. coli")
p8 <- p8 + ggtitle("Klebsiella")

combined_figure <-
  (p2 | p6 | p4 | p8 ) /
  (p10 | p14 | p12 | p16) /
  (p26 | p30 | p28 | p32) /
  (p18 | p22 | p20 | p24) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)",
                                      "j)", "k)", "l)", "m)", "n)", "o)", "p)"))) &
theme(
  legend.position = "bottom",
  plot.tag = element_text(family = "Arial", face = "bold", size = 16),
  plot.tag.position = c(0, 1),
  plot.title = element_text(family = "Arial", face = "bold", size = 18, hjust = 0.5)
)

combined_figure
# save
ggsave("model_results/combned_ecoli_kleb_4x4_hierarchical_richness_and_mass_ss_curves_by_region.png",
       plot = combined_figure, width = 16, height = 14, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# remove legends from plots, except p1, p3, p5, p7
# E. coli richness - both cumulative and ss plots
p2  <- p2  + guides(fill = "none", color = "none", alpha = "none")
p4  <- p4  + guides(fill = "none", color = "none", alpha = "none")
p6  <- p6  + guides(fill = "none", color = "none", alpha = "none")
p8  <- p8  + guides(fill = "none", color = "none", alpha = "none")
p9  <- p9  + guides(fill = "none", color = "none", alpha = "none")
p10 <- p10 + guides(fill = "none", color = "none", alpha = "none")
p11 <- p11 + guides(fill = "none", color = "none", alpha = "none")
p12 <- p12 + guides(fill = "none", color = "none", alpha = "none")
p13 <- p13 + guides(fill = "none", color = "none", alpha = "none")
p14 <- p14 + guides(fill = "none", color = "none", alpha = "none")
p15 <- p15 + guides(fill = "none", color = "none", alpha = "none")
p16 <- p16 + guides(fill = "none", color = "none", alpha = "none")
p17 <- p17 + guides(fill = "none", color = "none", alpha = "none")
p18 <- p18 + guides(fill = "none", color = "none", alpha = "none")
p19 <- p19 + guides(fill = "none", color = "none", alpha = "none")
p20 <- p20 + guides(fill = "none", color = "none", alpha = "none")
p21 <- p21 + guides(fill = "none", color = "none", alpha = "none")
p22 <- p22 + guides(fill = "none", color = "none", alpha = "none")
p23 <- p23 + guides(fill = "none", color = "none", alpha = "none")
p24 <- p24 + guides(fill = "none", color = "none", alpha = "none")
p25 <- p25 + guides(fill = "none", color = "none", alpha = "none")
p26 <- p26 + guides(fill = "none", color = "none", alpha = "none")
p27 <- p27 + guides(fill = "none", color = "none", alpha = "none")
p28 <- p28 + guides(fill = "none", color = "none", alpha = "none")
p29 <- p29 + guides(fill = "none", color = "none", alpha = "none")
p30 <- p30 + guides(fill = "none", color = "none", alpha = "none")
p31 <- p31 + guides(fill = "none", color = "none", alpha = "none")
p32 <- p32 + guides(fill = "none", color = "none", alpha = "none")


# E. coli species richness - cumulative and sample size
combined_figure <-
  (p1  | p2) /
  (p9  | p10) /
  (p25 | p26) /
  (p17 | p18) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_figure
# save
ggsave("model_results/ecoli_4x2_hierarchical_cumulative_richness_and_ss_curves_by_region.png",
       plot = combined_figure, width = 12, height = 16, units = "in", dpi = 300)


# E. coli population mass - both cumulative and ss plots
# combine plots:
combined_figure <-
  (p3 | p4) /
  (p11 | p12) /
  (p27 | p28) /
  (p19 | p20) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined_figure
# save
ggsave("model_results/ecoli_4x2_hierarchical_cumulative_mass_and_ss_curves_by_region.png",
       plot = combined_figure, width = 12, height = 16, units = "in", dpi = 300)



# Klebsiella richness - both cumulative and ss plots
# combine plots:
combined_figure <-
  (p5 | p6) /
  (p13 | p14) /
  (p29 | p30) /
  (p21 | p22) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined_figure
# save
ggsave("model_results/kleb_4x2_hierarchical_cumulative_richness_and_ss_curves_by_region.png",
       plot = combined_figure, width = 12, height = 16, units = "in", dpi = 300)


# Klebsiella mass - both cumulative and ss plots
# combine plots:
combined_figure <-
  (p7 | p8) /
  (p15 | p16) /
  (p31 | p32) /
  (p23 | p24) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
combined_figure
# save
ggsave("model_results/kleb_4x2_hierarchical_cumulative_mass_and_ss_curves_by_region.png",
       plot = combined_figure, width = 12, height = 16, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Bar plots for overall sample size - 80% coverage at 95% 95% confidence level for ecoli and kleb ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load overall data
mass_df_mlst <- read.csv("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
mass_df_fastbaps_L3 <- read.csv("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_frequency_mass_df.csv")
combined_bsi_pling_mass_df <- read.csv("model_results/combined_bsi_pling_bayes_mass_df.csv")
combined_bsi_arg_mass_df <- read.csv("model_results/combined_bsi_arg_bayes_mass_df.csv")

# load richness data
richness_df_mlst <- read.csv("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
richness_df_fastbaps_L3 <- read.csv("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_fine.csv")
combined_bsi_pling_richness_df <- read.csv("model_results/combined_bsi_pling_bayes_richness_df.csv")
combined_bsi_arg_richness_df <- read.csv("model_results/combined_bsi_arg_bayes_richness_df.csv")

# load hierarchical model data
ecoli_mlst_curve_df <- read.csv("model_results/ecoli_bsi_mlst_mass_curves.csv")
kleb_mlst_curve_df <- read.csv("model_results/kleb_bsi_mlst_mass_curves.csv")
ecoli_fastbaps_L3_curve_df <- read.csv("model_results/ecoli_bsi_fastbaps_L3_mass_curves.csv")
kleb_fastbaps_L3_curve_df <- read.csv("model_results/kleb_bsi_fastbaps_L3_mass_curves.csv")
ecoli_pling_curve_df <- read.csv("model_results/ecoli_pling_subiso_mass_curves.csv")
kleb_pling_curve_df <- read.csv("model_results/kleb_pling_subiso_mass_curves.csv")
ecoli_arg_curve_df <- read.csv("model_results/ecoli_arg_subiso_mass_curves.csv")
kleb_arg_curve_df <- read.csv("model_results/kleb_arg_subiso_mass_curves.csv")


# load richness data
ecoli_mlst_curve_df_richness <- read.csv("model_results/ecoli_bsi_mlst_richness_curves.csv")
kleb_mlst_curve_df_richness <- read.csv("model_results/kleb_bsi_mlst_richness_curves.csv")
ecoli_fastbaps_L3_curve_df_richness <- read.csv("model_results/ecoli_bsi_fastbaps_L3_richness_curves.csv")
kleb_fastbaps_L3_curve_df_richness <- read.csv("model_results/kleb_bsi_fastbaps_L3_richness_curves_fine.csv")
ecoli_pling_curve_df_richness <- read.csv("model_results/ecoli_pling_subiso_richness_curves.csv")
kleb_pling_curve_df_richness <- read.csv("model_results/kleb_pling_subiso_richness_curves.csv")
ecoli_arg_curve_df_richness <- read.csv("model_results/ecoli_arg_subiso_richness_curves.csv")
kleb_arg_curve_df_richness <- read.csv("model_results/kleb_arg_subiso_richness_curves.csv")


# format overall tables
overall_mlst <- mass_df_mlst |>
  mutate(feature = "MLST", model = "overall") |>
  left_join(richness_df_mlst, by = c("Genus" = "Genus", "f"="f"), suffix = c("_mass", "_richness")) |>
  rename(avg = median_mass) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, median_species_proportion, q2.5_species_proportion, q97.5_species_proportion, min_sample_95_richness, feature, model, Genus))

overall_fastbaps_L3 <- mass_df_fastbaps_L3 |>
  mutate(feature = "fastBAPS", model = "overall") |>
  left_join(richness_df_fastbaps_L3, by = c("Genus" = "Genus", "f"="f"), suffix = c("_mass", "_richness")) |>
  rename(avg = median_mass) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, median_species_proportion, q2.5_species_proportion, q97.5_species_proportion, min_sample_95_richness, feature, model, Genus))

overall_plasmid <- combined_bsi_pling_mass_df |>
  mutate(feature = "Plasmid", model = "overall") |>
  left_join(combined_bsi_pling_richness_df, by = c("Genus" = "Genus", "f"="f"), suffix = c("_mass", "_richness")) |>
  rename(avg = median_mass, q2.5_mass = q2.5, q97.5_mass = q97.5, q2.5_richness = q2.5_species_count, q97.5_richness = q97.5_species_count) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, median_species_proportion, q2.5_species_proportion, q97.5_species_proportion, min_sample_95_richness, feature, model, Genus))

overall_arg <- combined_bsi_arg_mass_df |>
  mutate(feature = "ARG", model = "overall") |> 
  left_join(combined_bsi_arg_richness_df, by = c("Genus" = "Genus", "f"="f"), suffix = c("_mass", "_richness")) |>
  rename(avg = median_mass, q2.5_mass = q2.5, q97.5_mass = q97.5, q2.5_richness = q2.5_species_count, q97.5_richness = q97.5_species_count) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, median_species_proportion, q2.5_species_proportion, q97.5_species_proportion, min_sample_95_richness, feature, model, Genus))



# reformat regionally-stratified data
ecoli_mlst_curve_national_df <- ecoli_mlst_curve_df |>
  left_join(ecoli_mlst_curve_df_richness, by = c("f" = "f", "region" = "region"), suffix = c("_mass", "_richness")) |>
  dplyr::distinct(f, global_all_median_mass, global_all_q2.5_mass, global_all_q97.5_mass, min_sample_95_mass, 
                  global_all_median_richness, global_all_q2.5_richness, global_all_q97.5_richness, min_sample_95_richness) |>
  dplyr::mutate(region = "UK National")

kleb_mlst_curve_national_df <- kleb_mlst_curve_df |>
  left_join(kleb_mlst_curve_df_richness, by = c("f" = "f", "region" = "region"), suffix = c("_mass", "_richness")) |>
  dplyr::distinct(f, global_all_median_mass, global_all_q2.5_mass, global_all_q97.5_mass, min_sample_95_mass, 
                  global_all_median_richness, global_all_q2.5_richness, global_all_q97.5_richness, min_sample_95_richness) |>
  dplyr::mutate(region = "UK National")

ecoli_fastbaps_L3_curve_national_df <- ecoli_fastbaps_L3_curve_df |>
  left_join(ecoli_fastbaps_L3_curve_df_richness, by = c("f" = "f", "region" = "region"), suffix = c("_mass", "_richness")) |>
  dplyr::distinct(f, global_all_median_mass, global_all_q2.5_mass, global_all_q97.5_mass, min_sample_95_mass, 
                  global_all_median_richness, global_all_q2.5_richness, global_all_q97.5_richness, min_sample_95_richness) |>
  dplyr::mutate(region = "UK National")

kleb_fastbaps_L3_curve_national_df <- kleb_fastbaps_L3_curve_df |>
  left_join(kleb_fastbaps_L3_curve_df_richness, by = c("f" = "f", "region" = "region"), suffix = c("_mass", "_richness")) |>
  dplyr::distinct(f, global_all_median_mass, global_all_q2.5_mass, global_all_q97.5_mass, min_sample_95_mass, 
                  global_all_median_richness, global_all_q2.5_richness, global_all_q97.5_richness, min_sample_95_richness) |>
  dplyr::mutate(region = "UK National")

ecoli_pling_curve_national_df <- ecoli_pling_curve_df |>
  left_join(ecoli_pling_curve_df_richness, by = c("f" = "f", "region" = "region"), suffix = c("_mass", "_richness")) |>
  dplyr::distinct(f, global_all_median_mass, global_all_q2.5_mass, global_all_q97.5_mass, min_sample_95_mass, 
                  global_all_median_richness, global_all_q2.5_richness, global_all_q97.5_richness, min_sample_95_richness) |>
  dplyr::mutate(region = "UK National")

kleb_pling_curve_national_df <- kleb_pling_curve_df |>
  left_join(kleb_pling_curve_df_richness, by = c("f" = "f", "region" = "region"), suffix = c("_mass", "_richness")) |>
  dplyr::distinct(f, global_all_median_mass, global_all_q2.5_mass, global_all_q97.5_mass, min_sample_95_mass, 
                  global_all_median_richness, global_all_q2.5_richness, global_all_q97.5_richness, min_sample_95_richness) |>
  dplyr::mutate(region = "UK National")

ecoli_arg_curve_national_df <- ecoli_arg_curve_df |>
  left_join(ecoli_arg_curve_df_richness, by = c("f" = "f", "region" = "region"), suffix = c("_mass", "_richness")) |>
  dplyr::distinct(f, global_all_median_mass, global_all_q2.5_mass, global_all_q97.5_mass, min_sample_95_mass, 
                  global_all_median_richness, global_all_q2.5_richness, global_all_q97.5_richness, min_sample_95_richness) |>
  dplyr::mutate(region = "UK National")

kleb_arg_curve_national_df <- kleb_arg_curve_df |>
  left_join(kleb_arg_curve_df_richness, by = c("f" = "f", "region" = "region"), suffix = c("_mass", "_richness")) |>
  dplyr::distinct(f, global_all_median_mass, global_all_q2.5_mass, global_all_q97.5_mass, min_sample_95_mass, 
                  global_all_median_richness, global_all_q2.5_richness, global_all_q97.5_richness, min_sample_95_richness) |>
  dplyr::mutate(region = "UK National")


# format regionally-stratified tables
ecoli_mlst_regional <- ecoli_mlst_curve_national_df |>
  filter(region == "UK National") |>
  mutate(feature = "MLST", model = "regional", Genus = "Escherichia") |>
  rename(avg = global_all_median_mass,
         q2.5_mass = global_all_q2.5_mass,
         q97.5_mass = global_all_q97.5_mass,
         median_species_proportion= global_all_median_richness,
         q2.5_species_proportion= global_all_q2.5_richness,
         q97.5_species_proportion = global_all_q97.5_richness) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, 
           median_species_proportion, q2.5_species_proportion, q97.5_species_proportion,
           min_sample_95_richness, 
           feature, model, Genus))

kleb_mlst_regional <- kleb_mlst_curve_national_df |>
  filter(region == "UK National") |>
  mutate(feature = "MLST", model = "regional", Genus = "Klebsiella") |>
  rename(avg = global_all_median_mass,
         q2.5_mass = global_all_q2.5_mass,
         q97.5_mass = global_all_q97.5_mass,
         median_species_proportion= global_all_median_richness,
         q2.5_species_proportion= global_all_q2.5_richness,
         q97.5_species_proportion = global_all_q97.5_richness) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, 
           median_species_proportion, q2.5_species_proportion, q97.5_species_proportion,
           min_sample_95_richness, 
           feature, model, Genus))

ecoli_fastbaps_L3_regional <- ecoli_fastbaps_L3_curve_national_df |>
  filter(region == "UK National") |>
  mutate(feature = "fastBAPS", model = "regional", Genus = "Escherichia") |>
  rename(avg = global_all_median_mass,
         q2.5_mass = global_all_q2.5_mass,
         q97.5_mass = global_all_q97.5_mass,
         median_species_proportion= global_all_median_richness,
         q2.5_species_proportion= global_all_q2.5_richness,
         q97.5_species_proportion = global_all_q97.5_richness) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, 
           median_species_proportion, q2.5_species_proportion, q97.5_species_proportion,
           min_sample_95_richness, 
           feature, model, Genus))

kleb_fastbaps_L3_regional <- kleb_fastbaps_L3_curve_national_df |>
  filter(region == "UK National") |>
  mutate(feature = "fastBAPS", model = "regional", Genus = "Klebsiella") |>
  rename(avg = global_all_median_mass,
         q2.5_mass = global_all_q2.5_mass,
         q97.5_mass = global_all_q97.5_mass,
         median_species_proportion= global_all_median_richness,
         q2.5_species_proportion= global_all_q2.5_richness,
         q97.5_species_proportion = global_all_q97.5_richness) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, 
           median_species_proportion, q2.5_species_proportion, q97.5_species_proportion,
           min_sample_95_richness, 
           feature, model, Genus))


ecoli_pling_regional <- ecoli_pling_curve_national_df |>
  filter(region == "UK National") |>
  mutate(feature = "Plasmid", model = "regional", Genus = "Escherichia") |>
  rename(avg = global_all_median_mass,
         q2.5_mass = global_all_q2.5_mass,
         q97.5_mass = global_all_q97.5_mass,
         median_species_proportion= global_all_median_richness,
         q2.5_species_proportion= global_all_q2.5_richness,
         q97.5_species_proportion = global_all_q97.5_richness) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, 
           median_species_proportion, q2.5_species_proportion, q97.5_species_proportion,
           min_sample_95_richness, 
           feature, model, Genus))

kleb_pling_regional <- kleb_pling_curve_national_df |>
  filter(region == "UK National") |>
  mutate(feature = "Plasmid", model = "regional", Genus = "Klebsiella") |>
  rename(avg = global_all_median_mass,
         q2.5_mass = global_all_q2.5_mass,
         q97.5_mass = global_all_q97.5_mass,
         median_species_proportion= global_all_median_richness,
         q2.5_species_proportion= global_all_q2.5_richness,
         q97.5_species_proportion = global_all_q97.5_richness) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, 
           median_species_proportion, q2.5_species_proportion, q97.5_species_proportion,
           min_sample_95_richness, 
           feature, model, Genus))

ecoli_arg_regional <- ecoli_arg_curve_national_df |>
  filter(region == "UK National") |>
  mutate(feature = "ARG", model = "regional", Genus = "Escherichia") |>
  rename(avg = global_all_median_mass,
         q2.5_mass = global_all_q2.5_mass,
         q97.5_mass = global_all_q97.5_mass,
         median_species_proportion= global_all_median_richness,
         q2.5_species_proportion= global_all_q2.5_richness,
         q97.5_species_proportion = global_all_q97.5_richness) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, 
           median_species_proportion, q2.5_species_proportion, q97.5_species_proportion,
           min_sample_95_richness, 
           feature, model, Genus))

kleb_arg_regional <- kleb_arg_curve_national_df |>
  filter(region == "UK National") |>
  mutate(feature = "ARG", model = "regional", Genus = "Klebsiella") |>
  rename(avg = global_all_median_mass,
         q2.5_mass = global_all_q2.5_mass,
         q97.5_mass = global_all_q97.5_mass,
         median_species_proportion= global_all_median_richness,
         q2.5_species_proportion= global_all_q2.5_richness,
         q97.5_species_proportion = global_all_q97.5_richness) |>
  select(c(f, avg, q2.5_mass, q97.5_mass, min_sample_95_mass, 
           median_species_proportion, q2.5_species_proportion, q97.5_species_proportion,
           min_sample_95_richness, 
           feature, model, Genus))

all_estimates <- rbind(overall_mlst, overall_fastbaps_L3, overall_plasmid, overall_arg,
      ecoli_mlst_regional, kleb_mlst_regional,
      ecoli_fastbaps_L3_regional, kleb_fastbaps_L3_regional,
      ecoli_pling_regional, kleb_pling_regional,
      ecoli_arg_regional, kleb_arg_regional) 

#View(overall_estimates)
#View(all_estimates)
# save
write.csv(all_estimates, "model_results/combined_bsi_overall_and_regional_sample_size_estimates.csv", row.names = FALSE)


first_crossing <- function(x, y, threshold = 0.8) {
  idx <- which(y >= threshold)[1]
  if (is.na(idx)) NA_real_ else x[idx]
}

get_threshold_proportions <- function(data, threshold) {
    data |>
    filter(f >= threshold) |>
    group_by(Genus) |>
    slice_min(
      order_by = f,
      n = 1,
      with_ties = FALSE
    ) |>
    ungroup() |>
    select(
      Genus,
      f,
      median_species_proportion,
      q2.5_species_proportion,
      q97.5_species_proportion
    )
}

get_threshold_proportions(richness_df_mlst, 0.0022)
get_threshold_proportions(richness_df_mlst, 0.0021)
get_threshold_proportions(richness_df_fastbaps_L3, 0.0022)
get_threshold_proportions(richness_df_fastbaps_L3, 0.0021)
get_threshold_proportions(combined_bsi_pling_richness_df, 0.0022)
get_threshold_proportions(combined_bsi_pling_richness_df, 0.0021)
get_threshold_proportions(combined_bsi_arg_richness_df, 0.0022)
get_threshold_proportions(combined_bsi_arg_richness_df, 0.0021)


sam_cov_summary <- all_estimates |>
  arrange(feature, model, Genus, min_sample_95_mass) |>
  group_by(feature, model, Genus) |>
  summarise(
    median = first_crossing(min_sample_95_mass, avg, threshold = 0.8),
    lower  = first_crossing(min_sample_95_mass, `q97.5_mass`, threshold = 0.8),
    upper  = first_crossing(min_sample_95_mass, `q2.5_mass`, threshold = 0.8),
    median_f = first_crossing(f, avg, threshold = 0.8),
    upper_f  = first_crossing(f, `q97.5_mass`, threshold = 0.8),
    lower_f  = first_crossing(f, `q2.5_mass`, threshold = 0.8),
    median_richness = first_crossing(median_species_proportion, avg, threshold = 0.8) ,
    lower_richness = first_crossing(q2.5_species_proportion, `q97.5_mass`, threshold = 0.8),
    upper_richness = first_crossing(q97.5_species_proportion, `q2.5_mass`, threshold = 0.8),
    .groups = "drop"
  )
# add pct of annual cases

sam_cov_summary <- sam_cov_summary |>
  mutate(total_pop = case_when(Genus == "Escherichia" ~ 42224,
                               Genus == "Klebsiella" ~ 13078),
         median_pct = median / total_pop * 100,
         lower_pct = lower / total_pop * 100,
         upper_pct = upper / total_pop * 100) |>
  rename(genus = Genus) |>
  mutate(sample_size = paste0(round(median, 0), " (", round(lower, 0), "-", round(upper, 0) , ")"),
         sample_size_pct = paste0(round(median_pct, 1), " (", round(lower_pct, 1), "-", round(upper_pct, 1) , ")"),
         richness = paste0(round(median_richness, 2)*100, " (", round(lower_richness, 2)*100, "-", round(upper_richness, 2)*100 , ")"),
         frequency = paste0(round(median_f, 4), " (", round(lower_f, 4), "-", round(upper_f, 4) , ")")
         )
#View(sam_cov_summary)
sam_cov_summary$feature <- factor(sam_cov_summary$feature, 
                                  levels = c("MLST", "fastBAPS","Plasmid","ARG"))

# save
write.csv(sam_cov_summary, "model_results/regional_and_global_sample_coverage_and_richness_summary.csv", row.names = FALSE)

# plot function
plot_sample_size_bars <- function(summary_df,
                                  feature_order = c(
                                    "MLST",
                                    "fastBAPS",
                                    "Plasmid",
                                    "ARG"
                                  ),
                                  genus_colors = c(
                                    Escherichia = "seagreen3",
                                    Klebsiella = "darkorange"
                                  ),
                                  model_alpha = c(
                                    overall = 1,
                                    regional = 0.4
                                  )) {
  summary_df <- summary_df |>
    mutate(
      feature = factor(feature, levels = feature_order),
      genus = factor(genus, levels = names(genus_colors)),
      model = factor(model, levels = names(model_alpha))
    )
  
  make_plot <- function(y, ymin, ymax, ylab) {
    ggplot(summary_df, aes(x = feature, y = .data[[y]], fill = genus, alpha = model)) +
      geom_col(
        aes(group = interaction(model, genus)),
        position = position_dodge(width = 0.75),
        width = 0.6,
        color = NA
      ) +
      geom_errorbar(
        aes(
          ymin = .data[[ymin]],
          ymax = .data[[ymax]],
          group = interaction(model, genus)
        ),
        position = position_dodge(width = 0.75),
        width = 0.2
      ) +
      scale_fill_manual(values = genus_colors, drop = FALSE) +
      scale_alpha_manual(values = model_alpha, drop = FALSE) +
      labs(x = NULL, y = ylab, fill = "Genus", alpha = "Model") +
      guides(
        fill = guide_legend(order = 1),
        alpha = guide_legend(order = 2)
      ) +
      theme_classic(base_size = 14)
  }
  
  list(
    count_plot = make_plot("median", "lower", "upper", "Count"),
    pct_plot = make_plot("median_pct", "lower_pct", "upper_pct",
                         "Percent (%) of annual BSIs to be sampled")
  )
}


plots <- plot_sample_size_bars(sam_cov_summary)

plots$count_plot
plots$pct_plot


ggsave("model_results/80_sam_cov_sample_size_barplot_count.png", plots$count_plot, width = 7, height = 4.5, dpi = 300)
ggsave("model_results/80_sam_cov_sample_size_barplot_pct.png", plots$pct_plot, width = 7, height = 4.5, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Bar plots for overall species richness (proportion of unique features detected) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Species richness summary
# load overall richness data
richness_df_mlst <- read.csv("model_results/combined_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
richness_df_fastbaps_L3 <- read.csv("model_results/combined_bsi_fastbaps_L3_bayes_cumulative_species_richness_df_fine.csv")
combined_bsi_pling_richness_df <- read.csv("model_results/combined_bsi_pling_bayes_richness_df.csv")
combined_bsi_arg_richness_df <- read.csv("model_results/combined_bsi_arg_bayes_richness_df.csv")

# format overall tables
overall_mlst <- richness_df_mlst |>
  mutate(feature = "MLST", model = "overall") |>
  rename(avg = median_species_proportion) |>
  select(c(f, avg, q2.5_species_proportion, q97.5_species_proportion, min_sample_95, feature, model, Genus)) |>
  rename(q2.5 = q2.5_species_proportion,
         q97.5 = q97.5_species_proportion)

overall_fastbaps_L3 <- richness_df_fastbaps_L3 |>
  mutate(feature = "fastbaps_L3", model = "overall") |>
  rename(avg = median_species_proportion) |>
  select(c(f, avg, q2.5_species_proportion, q97.5_species_proportion, min_sample_95, feature, model, Genus)) |>
  rename(q2.5 = q2.5_species_proportion,
       q97.5 = q97.5_species_proportion)

overall_pling <- combined_bsi_pling_richness_df |>
  mutate(feature = "Plasmid", model = "overall") |>
  rename(avg = median_species_proportion,
         q2.5 =  q2.5_species_proportion,
         q97.5 = q97.5_species_proportion) |>
  select(c(f, avg, q2.5, q97.5, min_sample_95, feature, model, Genus))

overall_arg <- combined_bsi_arg_richness_df |>
  mutate(feature = "ARG", model = "overall") |>
  rename(avg = median_species_proportion,
         q2.5 =  q2.5_species_proportion,
         q97.5 = q97.5_species_proportion) |>
  select(c(f, avg, q2.5, q97.5, min_sample_95, feature, model, Genus))

# rbind all
all_overall_species_richness_estimates <- rbind(overall_mlst, overall_fastbaps_L3, overall_pling, overall_arg)
#View(all_overall_species_richness_estimates)
# save
write.csv(all_overall_species_richness_estimates, "combined_bsi_overall_species_richness_and_ss_estimates.csv", row.names = FALSE)

# left join this onto sam cov summary
# round sample size estimates in both cases:
sam_cov_summary_overall <- sam_cov_summary |>
  filter(model == "overall") |>
  mutate(median = round(median),
         lower = round(lower),
         upper = round(upper))

# round sample size estimates in overall species richness table
all_overall_species_richness_estimates <- all_overall_species_richness_estimates |>
  mutate(min_sample_95 = round(min_sample_95)) |>
  group_by(feature, model, Genus, min_sample_95) |>
  slice_head() |>
  ungroup()

table(all_overall_species_richness_estimates$Genus, all_overall_species_richness_estimates$feature)


# left join median column
species_richness_for_sam_cov_estimates <- sam_cov_summary_overall |>
  left_join(all_overall_species_richness_estimates |> select(feature, model, Genus, min_sample_95, avg), 
            by = c("feature", "model", "genus" = "Genus", "median" = "min_sample_95")) |>
  rename(median_species_richness = avg)


# left join lower estimate column
species_richness_for_sam_cov_estimates <- species_richness_for_sam_cov_estimates |>
  left_join(all_overall_species_richness_estimates |> select(feature, model, Genus, min_sample_95, q2.5), 
            by = c("feature", "model", "genus" = "Genus", "lower" = "min_sample_95")) |>
  rename(lower_species_richness = q2.5)

# left join upper estimate column
species_richness_for_sam_cov_estimates <- species_richness_for_sam_cov_estimates |>
  left_join(all_overall_species_richness_estimates |> select(feature, model, Genus, min_sample_95, q97.5), 
            by = c("feature", "model", "genus" = "Genus", "upper" = "min_sample_95")) |>
  rename(upper_species_richness = q97.5)
#View(species_richness_for_sam_cov_estimates)

# save
write.csv(species_richness_for_sam_cov_estimates, "species_richness_for_80_pct_sam_cov.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. EXTERNAL VALIDATION ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 7.a NORM DATASET ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data 
NORM_supplementary <- read.csv("NORM_supplementary.csv")
#View(NORM_supplementary)
colnames(NORM_supplementary)
dim(NORM_supplementary) # 3254   = whole dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 7.ai Overall ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
NORM_ST_overall <- NORM_supplementary |>
  dplyr::group_by(ST) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(ST, count) |>
  dplyr::arrange(count, ST)
#View(NORM_ST_overall)
#table(NORM_ST_overall$count)
#nrow(NORM_ST_overall) # 263

# set prior parameters
df <- NORM_ST_overall
N <- sum(df$count) #3254
K <- length(df$count) # 403
f1 <- sum(df$count == 1)  # 259     # number of singletons
q_hat <- f1 / N  # 7.9% Good-Turing first-order - proportion of singletons

# Estimate novel/ unseen mass using Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(df$count)
crp_fit$summary_df # mean: theta = 86.4 -> 2.6% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.0258

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$ST
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #11

# apply simple bayes function
NORM_mlst_bayes <- isolate_bayes_dirichlet(NORM_ST_overall,
                                           feature_col = "ST",
                                           alpha_named = alpha_named,
                                           alpha_novel_sum = alpha_novel_sum,
                                           alpha_novel_num = alpha_novel_sum,
                                           B = 10000)

# save
saveRDS(NORM_mlst_bayes, "model_results/NORM_mlst_bayes.rds")

# mass and richness curves
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

# cumulative feature count (i.e. species richness) curve (number and proportion of MLSTs at least as frequent as f)
NORM_species_richness_df <- compute_species_richness_curve(
  prep_bootstrap_draws(NORM_mlst_bayes$draws), 
  f_grid = f_grid_99)
#View(NORM_species_richness_df)

# add sample sizes
NORM_species_richness_df <- NORM_species_richness_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  species_richness df
write.csv(NORM_species_richness_df, "model_results/NORM_species_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
NORM_mass_df <- compute_mass_curve(prep_bootstrap_draws(NORM_mlst_bayes$draws), f_grid = f_grid_99)
#View(NORM_mass_df)

# add sample sizes
NORM_mass_df <- NORM_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(NORM_mass_df, "model_results/NORM_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#NORM_mass_df <- read.csv("model_results/NORM_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# summary tables - population species_richness
species_richness_df_mlst_long <- NORM_species_richness_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))
#View(species_richness_df_mlst)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
norm_mlst_species_richness_summary_table <- species_richness_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_species_proportion", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5_species_proportion", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5_species_proportion", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(norm_mlst_species_richness_summary_table)

# tidy table - extract columns named like "75%_cell", "80%_cell", ...
norm_mlst_species_richness_summary_table <- norm_mlst_species_richness_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(norm_mlst_species_richness_summary_table)
#View(norm_mlst_species_richness_summary_table)
# save
write.csv(norm_mlst_species_richness_summary_table, "model_results/norm_mlst_species_richness_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~#
# summary tables 
# transform to long
mass_df_mlst_long <- NORM_mass_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
NORM_mlst_sample_coverage_summary_table <- mass_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
NORM_mlst_sample_coverage_summary_table <- NORM_mlst_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(NORM_mlst_sample_coverage_summary_table)
View(NORM_mlst_sample_coverage_summary_table)
# save
write.csv(NORM_mlst_sample_coverage_summary_table, "model_results/NORM_mlst_sample_coverage_summary_table.csv", row.names = FALSE)
# ~1,8000 for 1,600-2.100 over all years

#
#View(NORM_mlst_bayes$summary_df)
NORM_ST_frequency <- NORM_ST_overall |>
  mutate(actual = count / sum(count)) |>
  mutate(ST = as.character(ST))
#View(NORM_ST_frequency)
NORM_plot_df <- NORM_mlst_bayes$summary_df |>
  left_join(NORM_ST_frequency, by = c("feature" = "ST")) |>
  mutate(alpha_named = 1,
         alpha_novel_num = 1,
         alpha_novel_sum = 1,
         estimate = median,
         lo = q2.5,
         hi = q97.5) 

# plot observed vs expected
NORM_pred_vs_obs <- plot_obs_vs_pred(NORM_plot_df, alpha_novel_num_vals = c(1),  alpha_novel_sum_vals = c(1), 1)
NORM_pred_vs_obs
ggsave("model_results/NORM_pred_vs_obs.png", NORM_pred_vs_obs, width = 6, height =5, dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 7.b BSAC DATASET ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
BSAC_supplementary <- read_csv("BSAC_Supplemental_Data_S1.csv")
#View(BSAC_supplementary)
colnames(BSAC_supplementary)
#start just with MLST
dim(BSAC_supplementary) # 1509 x 407
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * 7.bi Overall ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
BSAC_ST_overall <- BSAC_supplementary |>
  dplyr::group_by(MLST) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(MLST, count) |>
  dplyr::arrange(count, MLST)
#View(BSAC_ST_overall)
#table(BSAC_ST_overall$count)
#nrow(BSAC_ST_overall) # 227

# set prior parameters
df <- BSAC_ST_overall
N <- sum(df$count)
K <- length(df$count) # 227
f1 <- sum(df$count == 1)  # 145    # number of singletons
q_hat <- f1 / N  # Good-Turing first-order - proportion of singletons
q_hat # 0.096

# Estimate novel/ unseen mass using Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(df$count)
crp_fit$summary_df # mean: theta = 51 -> 3.3% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.3286

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$MLST
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #8

# apply simple bayes function
BSAC_mlst_bayes <- isolate_bayes_dirichlet(BSAC_ST_overall,
                                           feature_col = "MLST",
                                           alpha_named = alpha_named,
                                           alpha_novel_sum = alpha_novel_sum,
                                           alpha_novel_num = alpha_novel_sum,
                                           B = 10000)

# save
saveRDS(BSAC_mlst_bayes, "model_results/BSAC_mlst_bayes.rds")

# mass and richness curves
n_grid <- c(seq(0, 10000, by = 1), seq(10005, 50000, by = 5), seq(50010, 100000, by = 10))
f_grid_99 <- 1-((1-0.99)^(1/n_grid))

# cumulative feature count (i.e. species richness) curve (number and proportion of MLSTs at least as frequent as f)
BSAC_species_richness_df <- compute_species_richness_curve(prep_bootstrap_draws(BSAC_mlst_bayes$draws), f_grid = f_grid_99)
#View(BSAC_species_richness_df)

# add sample sizes
BSAC_species_richness_df <- BSAC_species_richness_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  species_richness df
write.csv(BSAC_species_richness_df, "model_results/BSAC_species_richness_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
BSAC_mass_df <- compute_mass_curve(prep_bootstrap_draws(BSAC_mlst_bayes$draws), f_grid = f_grid_99)
#View(BSAC_mass_df)

# add sample sizes
BSAC_mass_df <- BSAC_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(BSAC_mass_df, "model_results/BSAC_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#BSAC_mass_df <- read.csv("model_results/BSAC_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# summary tables - population species_richness
species_richness_df_mlst_long <- BSAC_species_richness_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))
#View(species_richness_df_mlst)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
BSAC_mlst_species_richness_summary_table <- species_richness_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_species_proportion", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5_species_proportion", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5_species_proportion", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(BSAC_mlst_species_richness_summary_table)

# tidy table - extract columns named like "75%_cell", "80%_cell", ...
BSAC_mlst_species_richness_summary_table <- BSAC_mlst_species_richness_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(BSAC_mlst_species_richness_summary_table)
#View(BSAC_mlst_species_richness_summary_table)
# save
write.csv(BSAC_mlst_species_richness_summary_table, "model_results/BSAC_mlst_species_richness_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~#
# summary tables 
# transform to long
mass_df_mlst_long <- BSAC_mass_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
BSAC_mlst_sample_coverage_summary_table <- mass_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
BSAC_mlst_sample_coverage_summary_table <- BSAC_mlst_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(BSAC_mlst_sample_coverage_summary_table)
#View(BSAC_mlst_sample_coverage_summary_table)
# save
write.csv(BSAC_mlst_sample_coverage_summary_table, "model_results/BSAC_mlst_sample_coverage_summary_table.csv", row.names = FALSE)
# ~1,8000 for 1,600-2.100 over all years

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#View(BSAC_mlst_bayes$summary_df)
BSAC_ST_frequency <- BSAC_ST_overall |>
  mutate(actual = count / sum(count)) |>
  mutate(MLST = as.character(MLST))
#View(BSAC_ST_frequency)
BSAC_plot_df <- BSAC_mlst_bayes$summary_df |>
  left_join(BSAC_ST_frequency, by = c("feature" = "MLST")) |>
  mutate(alpha_named = 1,
         alpha_novel_sum = 1,
         alpha_novel_num = 1,
         estimate = median,
         lo = q2.5,
         hi = q97.5) 

# plot observed vs expected
BSAC_pred_vs_obs <- plot_obs_vs_pred(BSAC_plot_df, alpha_novel_num_vals = c(1),  alpha_novel_sum_vals = c(1), 1)
BSAC_pred_vs_obs
ggsave("model_results/BSAC_pred_vs_obs.png", BSAC_pred_vs_obs, width = 6, height =5, dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 7.c Overall comparative plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# combine NEKSUS, NORM and NEKSUS plots
# load data
# richness dfs
neksus_species_richness_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
NORM_species_richness_df <- read.csv("model_results/NORM_species_richness_df.csv")
BSAC_species_richness_df <- read.csv("model_results/BSAC_species_richness_df.csv")

neksus_species_richness_df <- neksus_species_richness_df |>   mutate(dataset = "NEKSUS")
NORM_species_richness_df <- NORM_species_richness_df |>   mutate(dataset = "NORM")
BSAC_species_richness_df <- BSAC_species_richness_df |>   mutate(dataset = "BSAC")

# combine
combined_richness_dfs <- rbind(neksus_species_richness_df, NORM_species_richness_df, BSAC_species_richness_df)
#View(combined_richness_dfs)

# mss dfs
neksus_mass_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
NORM_mass_df <- read.csv("model_results/NORM_mass_df.csv")
BSAC_mass_df <- read.csv("model_results/BSAC_mass_df.csv")

neksus_mass_df <- neksus_mass_df |> mutate(dataset = "NEKSUS")
NORM_mass_df <- NORM_mass_df |> mutate(dataset = "NORM")
BSAC_mass_df <- BSAC_mass_df |> mutate(dataset = "BSAC")

# combine
combined_mass_dfs <- rbind(neksus_mass_df, NORM_mass_df, BSAC_mass_df)
#View(combined_mass_dfs)

#~~~~~~~~~~~~~~#
dataset_colours <- c(
  NEKSUS = "#009E73",  # bluish green
  NORM   = "#CC79A7",  # reddish purple
  BSAC   = "#D55E00"   # vermillion
)

# Plot (cumulative species_richness of population at frequency >= f)
NEKSUS_NORM_BSAC_cumulative_species_richness_plot <- ggplot(combined_richness_dfs, aes(x = f, y = median_species_proportion, colour =  dataset, fill = dataset)) +
  geom_line(data = combined_richness_dfs, aes(x = f, y = median_species_proportion, colour =  dataset)) +
  geom_ribbon(data = combined_richness_dfs , aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Dataset", values = dataset_colours) +
  scale_colour_manual(name = "Dataset", values = dataset_colours) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(combined_richness_dfs$f), 1)) +
  labs(x = "MLST frequency (f) (log scale)",
       y = "Proportion of unique MLSTs of frequency ≥ f"
  ) +
  
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
NEKSUS_NORM_BSAC_cumulative_species_richness_plot
# save
ggsave("model_results/NEKSUS_NORM_BSAC_cumulative_species_richness_proportion_plot.png", plot = NEKSUS_NORM_BSAC_cumulative_species_richness_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Cumulative fraction with sample size plot 
bayes_sample_coverage_plot <- ggplot(combined_richness_dfs, aes(y = median_species_proportion, colour = dataset, fill = dataset)) +
  #geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  #geom_line(aes(x = min_sample_99)) +
  #geom_ribbon(aes(x = min_sample_90, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.35, colour = NA) +
  #geom_ribbon(aes(x = min_sample_99, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = dataset_colours) +
  scale_colour_manual(values = dataset_colours) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_x_log10(
    limits = c(10, 10000),
    breaks = c(1, 10, 100, 1000, 10000)
  ) +
  scale_y_continuous() +
  labs(x = "Sample size",
       y = "coverage") +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/NEKSUS_NORM_BSAC_mlst_bayes_species_richness_vs_sample_size_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#

# Plot (cumulative mass of population at frequency >= f)
NEKSUS_NORM_BSAC_cumulative_mass_plot <- ggplot(combined_mass_dfs, aes(x = f, y = median_mass, colour =  dataset, fill = dataset)) +
  geom_line(data = combined_mass_dfs, aes(x = f, y = median_mass, colour =  dataset)) +
  geom_ribbon(data = combined_mass_dfs , aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Dataset", values = dataset_colours) +
  scale_colour_manual(name = "Dataset", values = dataset_colours) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(combined_mass_dfs$f), 1)) +
  labs(x = "MLST frequency (f) (log scale)",
       y = "Proportion of population belonging to MLSTs of frequency ≥ f"
  ) +
  
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
NEKSUS_NORM_BSAC_cumulative_mass_plot
# save
ggsave("model_results/NEKSUS_NORM_BSAC_cumulative_mass_plot.png", plot = NEKSUS_NORM_BSAC_cumulative_mass_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Cumulative fraction with sample size plot 
bayes_sample_coverage_plot <- ggplot(combined_mass_dfs, aes(y = median_mass, colour = dataset, fill = dataset)) +
  #geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  #geom_line(aes(x = min_sample_99)) +
  #geom_ribbon(aes(x = min_sample_90, ymin = q2.5, ymax = q97.5), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5, ymax = q97.5), alpha = 0.35, colour = NA) +
  #geom_ribbon(aes(x = min_sample_99, ymin = q2.5, ymax = q97.5), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = dataset_colours) +
  scale_colour_manual(values = dataset_colours) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_x_log10(
    limits = c(10, 10000),
    breaks = c(1, 10, 100, 1000, 10000)
  ) +
  scale_y_continuous() +
  labs(x = "Sample size",
       y = "coverage") +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/NEKSUS_NORM_BSAC_mlst_bayes_sample_coverage_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 7.d NORM Date Filtered ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
NORM_ST_overall_post2011 <- NORM_supplementary |>
  dplyr::filter(year >= 2011) |>
  dplyr::group_by(ST) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(ST, count) |>
  dplyr::arrange(count, ST)
#View(NORM_ST_overall_post2011)
sum(NORM_ST_overall_post2011$count) # 1977
nrow(NORM_ST_overall_post2011) #  297
write.csv(NORM_ST_overall_post2011, "NORM_ST_overall_post2011_mlst_counts.csv", row.names = FALSE, quote = FALSE)

# set prior parameters
df <- NORM_ST_overall_post2011
N <- sum(df$count) #1977
K <- length(df$count) # 297
f1 <- sum(df$count == 1)  # 191     # number of singletons
q_hat <- f1 / N  # Good-Turing first-order - proportion of singletons
q_hat # 9.7%

# Estimate novel/ unseen mass using Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(df$count)
crp_fit$summary_df # mean: theta = 67 -> 3.3 prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.033

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$ST
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #11

# apply simple bayes function
NORM_mlst_bayes_post2011 <- isolate_bayes_dirichlet(NORM_ST_overall_post2011,
                                           feature_col = "ST",
                                           alpha_named = alpha_named,
                                           alpha_novel_sum = alpha_novel_sum,
                                           alpha_novel_num = alpha_novel_sum,
                                           B = 10000)

# save
saveRDS(NORM_mlst_bayes_post2011, "model_results/NORM_mlst_bayes_post2011.rds")

# mass and richness curves
NORM_species_richness_df <- compute_species_richness_curve(prep_bootstrap_draws(NORM_mlst_bayes_post2011$draws), f_grid = f_grid_99)
#View(NORM_species_richness_df)

# add sample sizes
NORM_species_richness_df <- NORM_species_richness_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  species_richness df
write.csv(NORM_species_richness_df, "model_results/NORM_species_richness_df_post2011.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
ORM_mass_df <- compute_mass_curve(prep_bootstrap_draws(NORM_mlst_bayes_post2011$draws), f_grid = f_grid_99)
#View(NORM_mass_df)

# add sample sizes
NORM_mass_df <- NORM_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(NORM_mass_df, "model_results/NORM_mass_df_post2011.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#NORM_mass_df <- read.csv("model_results/NORM_mass_df_post2011.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# summary tables - population species_richness
species_richness_df_mlst_long <- NORM_species_richness_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))
#View(species_richness_df_mlst)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
norm_mlst_species_richness_summary_table <- species_richness_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_species_proportion", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5_species_proportion", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5_species_proportion", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(norm_mlst_species_richness_summary_table)

# tidy table - extract columns named like "75%_cell", "80%_cell", ...
norm_mlst_species_richness_summary_table <- norm_mlst_species_richness_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(norm_mlst_species_richness_summary_table)
#View(norm_mlst_species_richness_summary_table)
# save
write.csv(norm_mlst_species_richness_summary_table, "model_results/norm_mlst_species_richness_summary_table_post2011.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~#
# summary tables 
# transform to long
mass_df_mlst_long <- NORM_mass_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
NORM_mlst_sample_coverage_summary_table <- mass_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
NORM_mlst_sample_coverage_summary_table <- NORM_mlst_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(NORM_mlst_sample_coverage_summary_table)
#View(NORM_mlst_sample_coverage_summary_table)
# save
write.csv(NORM_mlst_sample_coverage_summary_table, "model_results/NORM_mlst_sample_coverage_summary_table_post2011.csv", row.names = FALSE)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#View(NORM_mlst_bayes_post2011$summary_df)
NORM_ST_frequency <- NORM_ST_overall_post2011 |>
  mutate(actual = count / sum(count)) |>
  mutate(ST = as.character(ST))
#View(NORM_ST_frequency)
NORM_plot_df <- NORM_mlst_bayes_post2011$summary_df |>
  left_join(NORM_ST_frequency, by = c("feature" = "ST")) |>
  mutate(alpha_named = 1,
         alpha_novel_sum = 1,
         alpha_novel_num = 1,
         estimate = median,
         lo = q2.5,
         hi = q97.5) 

# plot observed vs expected
NORM_pred_vs_obs <- plot_obs_vs_pred(NORM_plot_df, alpha_novel_sum_vals = c(1), alpha_novel_num_vals = c(1), 1)
NORM_pred_vs_obs
ggsave("model_results/NORM_pred_vs_obs_post2011.png", NORM_pred_vs_obs, width = 6, height =5, dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 7.e BSAC date filtered ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
BSAC_ST_overall_post2003 <- BSAC_supplementary |>
  dplyr::filter(Year_of_isolation >= 2003) |>
  dplyr::group_by(MLST) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(MLST, count) |>
  dplyr::arrange(count, MLST)
#View(BSAC_ST_overall_post2003)
sum(BSAC_ST_overall_post2003$count)
nrow(BSAC_ST_overall_post2003) # 171
write.csv(BSAC_ST_overall_post2003, "BSAC_ST_overall_post2003_mlst_counts.csv", row.names = FALSE, quote = FALSE)

# set prior parameters
df <- BSAC_ST_overall_post2003
N <- sum(df$count)
K <- length(df$count) # 203
f1 <- sum(df$count == 1)  # 132    # number of singletons
q_hat <- f1 / N  # Good-Turing first-order - proportion of singletons
q_hat # 0.099

# Estimate novel/ unseen mass using Chinese restaurant process (CRP) theta parameter
crp_fit <- fit_crp_theta_bayes(df$count)
crp_fit$summary_df # mean: theta = 46 -> 3.3% prob that next isolate is novel
# estimate for total unseen/ novel probability mass:
novelty_prob_hat <- median(crp_fit$novelty_prob_draws)
novelty_prob_hat # # 0.334

# set priors
alpha_named <- rep(1, K)  # set uninformative priors
names(alpha_named) <- df$MLST
alpha_novel_sum <- ceiling(set_alpha_novel_sum_from_crp(alpha_named, novelty_prob_hat))
alpha_novel_sum #8

# apply simple bayes function
BSAC_mlst_bayes_post2003 <- isolate_bayes_dirichlet(BSAC_ST_overall_post2003,
                                           feature_col = "MLST",
                                           alpha_named = alpha_named,
                                           alpha_novel_sum = alpha_novel_sum,
                                           alpha_novel_num = alpha_novel_sum,
                                           B = 10000)

# save
saveRDS(BSAC_mlst_bayes_post2003, "model_results/BSAC_mlst_bayes_post2003.rds")

# mass and richness curves
BSAC_species_richness_df <- compute_species_richness_curve(prep_bootstrap_draws(BSAC_mlst_bayes_post2003$draws), f_grid = f_grid_99)
#View(BSAC_species_richness_df)

# add sample sizes
BSAC_species_richness_df <- BSAC_species_richness_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  species_richness df
write.csv(BSAC_species_richness_df, "model_results/BSAC_species_richness_df_post2003.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
BSAC_mass_df <- compute_mass_curve(prep_bootstrap_draws(BSAC_mlst_bayes_post2003$draws), f_grid = f_grid_99)
#View(BSAC_mass_df)

# add sample sizes
BSAC_mass_df <- BSAC_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(BSAC_mass_df, "model_results/BSAC_mass_df_post2003.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# Reload saved data 
#BSAC_mass_df <- read.csv("model_results/BSAC_mass_df_post2003.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# summary tables - population species_richness
species_richness_df_mlst_long <- BSAC_species_richness_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))
#View(species_richness_df_mlst)

# thresholds to evaluate
thresh <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# main pipeline: compute per-Estimator x threshold cells
BSAC_mlst_species_richness_summary_table <- species_richness_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_species_proportion", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5_species_proportion", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5_species_proportion", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )
#View(BSAC_mlst_species_richness_summary_table)

# tidy table - extract columns named like "75%_cell", "80%_cell", ...
BSAC_mlst_species_richness_summary_table <- BSAC_mlst_species_richness_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(BSAC_mlst_species_richness_summary_table)
#View(BSAC_mlst_species_richness_summary_table)
# save
write.csv(BSAC_mlst_species_richness_summary_table, "model_results/BSAC_mlst_species_richness_summary_table_post2003.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~#
# summary tables 
# transform to long
mass_df_mlst_long <- BSAC_mass_df |>
  pivot_longer(cols = c("min_sample_90", "min_sample_95", "min_sample_99") , names_to = "Estimator" , values_to = "sample_size") |>
  mutate(Genus_Estimator = paste0(Genus, "|", Estimator))

# thresholds to evaluate
thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
BSAC_mlst_sample_coverage_summary_table <- mass_df_mlst_long |>
  group_by(Genus_Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      median_ss  <- min_sample_at_or_above(df, "median_mass", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t) # switch 97.5 and 2.5ht percentiles as lowest sample size 
      hi_ss    <- min_sample_at_or_above(df, "q2.5", t)
      
      # Format: "median (lower - upper)". If none available, return NA string.
      formatted <- if (is.na(median_ss) && is.na(lo_ss) && is.na(hi_ss)) {
        NA_character_
      } else {
        # Replace NA components with "NA" in the string or use >max indicator if preferred
        median_txt <- if (is.na(median_ss)) "NA" else formatC(median_ss, format = "d", big.mark = ",")
        lo_txt   <- if (is.na(lo_ss))   "NA" else formatC(lo_ss, format = "d", big.mark = ",")
        hi_txt   <- if (is.na(hi_ss))   "NA" else formatC(hi_ss, format = "d", big.mark = ",")
        str_c(median_txt, " (", lo_txt, " - ", hi_txt, ")")
      }
      
      tibble(threshold = t, cell = formatted,
             median_ss = median_ss, lo_ss = lo_ss, hi_ss = hi_ss)
    })
    out
  }, .keep = TRUE) |>
  ungroup() |>
  # pivot thresholds into columns named "75%", "80%", ...
  mutate(threshold_label = paste0(as.integer(threshold * 100), "%")) |>
  select(Genus_Estimator, threshold_label, cell, median_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, median_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
BSAC_mlst_sample_coverage_summary_table <- BSAC_mlst_sample_coverage_summary_table |>
  select(Genus_Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(BSAC_mlst_sample_coverage_summary_table)
#View(BSAC_mlst_sample_coverage_summary_table)
# save
write.csv(BSAC_mlst_sample_coverage_summary_table, "model_results/BSAC_mlst_sample_coverage_summary_table_post2003.csv", row.names = FALSE)
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#View(BSAC_mlst_bayes$summary_df)
BSAC_ST_frequency <- BSAC_ST_overall |>
  mutate(actual = count / sum(count)) |>
  mutate(MLST = as.character(MLST))
#View(BSAC_ST_frequency)
BSAC_plot_df <- BSAC_mlst_bayes$summary_df |>
  left_join(BSAC_ST_frequency, by = c("feature" = "MLST")) |>
  mutate(alpha_named = 1,
         alpha_novel_sum = 1,
         alpha_novel_num = 1,
         estimate = median,
         lo = q2.5,
         hi = q97.5) 

# plot observed vs expected
BSAC_pred_vs_obs <- plot_obs_vs_pred(BSAC_plot_df, alpha_novel_num_vals = c(1), alpha_novel_sum_vals = c(1), 1)
BSAC_pred_vs_obs
ggsave("model_results/BSAC_pred_vs_obs.png", BSAC_pred_vs_obs, width = 6, height =5, dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * 7.f Overall comparative plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# combine NEKSUS, NORM and NEKSUS plots
# load data
# richness dfs
neksus_species_richness_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_species_richness_df.csv")
NORM_species_richness_df <- read.csv("model_results/NORM_species_richness_df_post2011.csv")
BSAC_species_richness_df <- read.csv("model_results/BSAC_species_richness_df_post2003.csv")

neksus_species_richness_df <- neksus_species_richness_df |> mutate(dataset = "NEKSUS")
NORM_species_richness_df <- NORM_species_richness_df |> mutate(dataset = "NORM")
BSAC_species_richness_df <- BSAC_species_richness_df |> mutate(dataset = "BSAC")

# combine
combined_richness_dfs <- rbind(neksus_species_richness_df, NORM_species_richness_df, BSAC_species_richness_df)
#View(combined_richness_dfs)

# mss dfs
neksus_mass_df <- read.csv("model_results/ecoli_bsi_kleborate_mlst_bayes_cumulative_frequency_mass_df.csv")
NORM_mass_df <- read.csv("model_results/NORM_mass_df_post2011.csv")
BSAC_mass_df <- read.csv("model_results/BSAC_mass_df_post2003.csv")

neksus_mass_df <- neksus_mass_df |>  mutate(dataset = "NEKSUS")
NORM_mass_df <- NORM_mass_df |>  mutate(dataset = "NORM")
BSAC_mass_df <- BSAC_mass_df |>  mutate(dataset = "BSAC")

# combine
combined_mass_dfs <- rbind(neksus_mass_df, NORM_mass_df, BSAC_mass_df)
#View(combined_mass_dfs)

#~~~~~~~~~~~~~~#
dataset_colours <- c(
  NEKSUS = "#009E73",  # bluish green
  NORM   = "#CC79A7",  # reddish purple
  BSAC   = "#D55E00"   # vermillion
)

# Plot (cumulative species_richness of population at frequency >= f)
NEKSUS_NORM_BSAC_cumulative_species_richness_plot <- ggplot(combined_richness_dfs, aes(x = f, y = median_species_proportion, colour =  dataset, fill = dataset)) +
  geom_line(data = combined_richness_dfs, aes(x = f, y = median_species_proportion, colour =  dataset)) +
  geom_ribbon(data = combined_richness_dfs , aes(ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Dataset", values = dataset_colours) +
  scale_colour_manual(name = "Dataset", values = dataset_colours) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(combined_richness_dfs$f), 1)) +
  labs(x = "MLST frequency (f) (log scale)",
       y = "Proportion of unique MLSTs of frequency ≥ f"
  ) +
  
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
NEKSUS_NORM_BSAC_cumulative_species_richness_plot
# save
ggsave("model_results/NEKSUS_NORM_BSAC_cumulative_species_richness_proportion_plot.png", plot = NEKSUS_NORM_BSAC_cumulative_species_richness_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Cumulative fraction with sample size plot 
bayes_sample_coverage_plot <- ggplot(combined_richness_dfs, aes(y = median_species_proportion, colour = dataset, fill = dataset)) +
  #geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  #geom_line(aes(x = min_sample_99)) +
  #geom_ribbon(aes(x = min_sample_90, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.35, colour = NA) +
  #geom_ribbon(aes(x = min_sample_99, ymin = q2.5_species_proportion, ymax = q97.5_species_proportion), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = dataset_colours) +
  scale_colour_manual(values = dataset_colours) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_x_log10(
    limits = c(10, 10000),
    breaks = c(1, 10, 100, 1000, 10000)
  ) +
  scale_y_continuous() +
  labs(x = "Sample size",
       y = "coverage") +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/NEKSUS_NORM_BSAC_mlst_bayes_species_richness_vs_sample_size_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#

# Plot (cumulative mass of population at frequency >= f)
NEKSUS_NORM_BSAC_cumulative_mass_plot <- ggplot(combined_mass_dfs, aes(x = f, y = median_mass, colour =  dataset, fill = dataset)) +
  geom_line(data = combined_mass_dfs, aes(x = f, y = median_mass, colour =  dataset)) +
  geom_ribbon(data = combined_mass_dfs , aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Dataset", values = dataset_colours) +
  scale_colour_manual(name = "Dataset", values = dataset_colours) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(combined_mass_dfs$f), 1)) +
  labs(x = "MLST frequency (f) (log scale)",
       y = "Proportion of population belonging to MLSTs of frequency ≥ f"
  ) +
  
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
NEKSUS_NORM_BSAC_cumulative_mass_plot
# save
ggsave("model_results/NEKSUS_NORM_BSAC_cumulative_mass_plot.png", plot = NEKSUS_NORM_BSAC_cumulative_mass_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# Cumulative fraction with sample size plot 
bayes_sample_coverage_plot <- ggplot(combined_mass_dfs, aes(y = median_mass, colour = dataset, fill = dataset)) +
  #geom_line(aes(x = min_sample_90)) +
  geom_line(aes(x = min_sample_95) ) +
  #geom_line(aes(x = min_sample_99)) +
  #geom_ribbon(aes(x = min_sample_90, ymin = q2.5, ymax = q97.5), alpha = 0.1, colour = NA) +
  geom_ribbon(aes(x = min_sample_95, ymin = q2.5, ymax = q97.5), alpha = 0.35, colour = NA) +
  #geom_ribbon(aes(x = min_sample_99, ymin = q2.5, ymax = q97.5), alpha = 0.6, colour = NA) +
  scale_fill_manual(values = dataset_colours) +
  scale_colour_manual(values = dataset_colours) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_x_log10(
    limits = c(10, 10000),
    breaks = c(1, 10, 100, 1000, 10000)
  ) +
  scale_y_continuous() +
  labs(x = "Sample size",
       y = "coverage") +
  theme_minimal()
bayes_sample_coverage_plot
ggsave("model_results/NEKSUS_NORM_BSAC_mlst_bayes_sample_coverage_plot_80.png", plot = bayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#

