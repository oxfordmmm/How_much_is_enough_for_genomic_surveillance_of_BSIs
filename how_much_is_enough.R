#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# RAREFACTION ANALYSIS ####
# for NEKSUS stuy overall 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * E. coli Overall ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages
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

# set workign directory
setwd("~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/")

# load E.coli and Kleb MLST data
ecoli_bsi_samples_metadata <- read.csv("rarefaction/ecoli_bsi_samples_metadata.csv")
kleb_bsi_samples_metadata <- read.csv("rarefaction/kleb_bsi_samples_metadata.csv")
#colnames(ecoli_bsi_samples_metadata)
#colnames(kleb_bsi_samples_metadata)

per_region_ecoli <- ecoli_bsi_samples_metadata |>
  group_by(region) |>
  summarise(n = n(),
    n_mlsts = n_distinct(escherichia__mlst_achtman__ST),
    n_clonal_complex = n_distinct(escherichia__mlst_achtman__clonal_complex),
    n_L1 = n_distinct(Level.1),
    n_L2 = n_distinct(Level.2),
    n_L3 = n_distinct(Level.3)) |>
  # mutate cols to give proportion of variation relative to regional sample size
  mutate(mlst_ration = n_mlsts / n,
         prop_clonal_complex = n_clonal_complex / n,
         prop_L1 = n_L1 / n,
         prop_L2 = n_L2 / n,
         prop_L3 = n_L3 / n
         )


View(per_region_ecoli)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load amrfinder data
# this df already had PLING plasmid "community_subcommunity" and "rep_types_whole_plasmid"
amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
#colnames(amrfinder_metadata_updated)
table(amrfinder_metadata_updated$kraken2_species)

# filter only ecoli BSI from amrfinder df
ecoli_bsi_amrfinder_metadata <- amrfinder_metadata_updated |>
  filter(grepl("Escherichia", kraken2_species) & (sampletype == "blood" | sampletype == "unknown" | is.na(sampletype))) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
length(unique(ecoli_bsi_amrfinder_metadata$sample)) #1471 non-duplicates
# filter only Kleb BSI from amrfinder df
kleb_bsi_amrfinder_metadata <- amrfinder_metadata_updated |>
  filter(grepl("Klebsiella", kraken2_species) & (sampletype == "blood" | sampletype == "unknown" | is.na(sampletype))) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
length(unique(kleb_bsi_amrfinder_metadata$sample)) #468 non-duplicates


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load fastbaps clusters (from fasatbaps_clustering.R script) and merge with dfs
all_ecoli_partition <- read.csv("ecoli_fastbaps_clusters.csv")
all_kleb_partition <- read.csv("kleb_by_sc_fastbaps_clusters.csv")

#View(all_ecoli_partition)
#nrow(all_ecoli_partition)
#nrow(ecoli_bsi_samples_metadata)
#View(ecoli_bsi_samples_metadata)
#table(all_ecoli_partition$run)
colnames(all_ecoli_partition)
colnames(ecoli_bsi_samples_metadata)
intersect(colnames(ecoli_bsi_samples_metadata), colnames(all_ecoli_partition))
#"run"             "hospital"        "kraken2_species" "sampletype_cat"      

# merge onto metadata
ecoli_bsi_samples_metadata <- ecoli_bsi_samples_metadata |>
  left_join(all_ecoli_partition |> select(c(sample, run, Level.1, Level.2, Level.3)), by = c("sequencing_id" = "sample", "run"= "run"))
# merge onto amrfinder
nrow(ecoli_bsi_samples_metadata)
sum(is.na(ecoli_bsi_samples_metadata$Level.1)) # 0
sum(is.na(ecoli_bsi_samples_metadata$hospital)) # 0
#save
write.csv(ecoli_bsi_samples_metadata, "rarefaction/ecoli_bsi_samples_metadata.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in ####
#ecoli_bsi_samples_metadata <- read.csv("rarefaction/ecoli_bsi_samples_metadata.csv")
#~~~~~~~~~~~~~~~#
#kleb_bsi_samples_metadata <- kleb_bsi_samples_metadata[,C(1:56)]
#nrow(kleb_bsi_samples_metadata)
#colnames(kleb_bsi_samples_metadata)

# merge kleb metadata
#View(all_kleb_partition)
#nrow(all_kleb_partition)
#nrow(kleb_bsi_samples_metadata)
#View(kleb_bsi_samples_metadata)
#table(all_kleb_partition$run)
#nrow(all_kleb_partition)
colnames(all_kleb_partition)
colnames(kleb_bsi_samples_metadata)
intersect(colnames(kleb_bsi_samples_metadata), colnames(all_kleb_partition))
#"run"             "hospital"        "kraken2_species" "sampletype_cat"      


kleb_bsi_samples_metadata <- kleb_bsi_samples_metadata |>
  left_join(all_kleb_partition |> select(c(sample, run, `Level 1`, `Level 2`, `Level 3`)), by = c("sequencing_id" = "sample", "run"= "run"))
nrow(kleb_bsi_samples_metadata) # 468
sum(is.na(kleb_bsi_samples_metadata$`Level 1`)) #0
sum(is.na(kleb_bsi_samples_metadata$hospital)) #0
#save
write.csv(kleb_bsi_samples_metadata, "rarefaction/kleb_bsi_samples_metadata.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in ####
#kleb_bsi_samples_metadata <- read.csv("rarefaction/kleb_bsi_samples_metadata.csv")
#~~~~~~~~~~~~~~~#

# merge amrfinder with fastbaps clusters info
# ecoli
intersect(colnames(amrfinder_metadata_updated), colnames(all_ecoli_partition))
nrow(ecoli_bsi_amrfinder_metadata)

ecoli_bsi_amrfinder_metadata <- ecoli_bsi_amrfinder_metadata |>
  left_join(all_ecoli_partition |> select(c(sample, run, Level.1, Level.2, Level.3)), by = c("sample"="sample", "run"= "run"))
nrow(ecoli_bsi_amrfinder_metadata) # 74106
length(unique(ecoli_bsi_amrfinder_metadata$sample)) # 1471
sum(is.na(ecoli_bsi_amrfinder_metadata$`Level.1`))
sum(is.na(ecoli_bsi_amrfinder_metadata$mlst))
colnames(ecoli_bsi_amrfinder_metadata)

#save
write.csv(ecoli_bsi_amrfinder_metadata, "rarefaction/ecoli_bsi_amrfinder_metadata.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in ####
#ecoli_bsi_amrfinder_metadata <- read.csv("rarefaction/ecoli_bsi_amrfinder_metadata.csv")
#~~~~~~~~~~~~~~~#

#kleb
kleb_bsi_amrfinder_metadata <- kleb_bsi_amrfinder_metadata |>
  left_join(all_kleb_partition |> select(c(sample, run, `Level 1`, `Level 2`, `Level 3`)), by = c("sample"="sample", "run"= "run"))
nrow(kleb_bsi_amrfinder_metadata) # 14825
length(unique(kleb_bsi_amrfinder_metadata$sample)) # 468
sum(is.na(kleb_bsi_amrfinder_metadata$`Level.1`))
sum(is.na(kleb_bsi_amrfinder_metadata$mlst))
colnames(kleb_bsi_amrfinder_metadata)

#save
write.csv(kleb_bsi_amrfinder_metadata, "rarefaction/kleb_bsi_amrfinder_metadata.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in ####
#kleb_bsi_amrfinder_metadata <- read.csv("rarefaction/kleb_bsi_amrfinder_metadata.csv")
#~~~~~~~~~~~~~~~#

# read in kleborate data
#~~~~~~~~~~~~#
# * Entry point to read in cleaned df ####
kleborate_escherichia_cleaned <- read.csv("kleborate_escherichia_cleaned.csv")
View(kleborate_escherichia_cleaned)
length(unique(kleborate_escherichia_cleaned$sample)) #1559 so all non-duplicate ecoli. 
colnames(kleborate_escherichia_cleaned)
table(kleborate_escherichia_cleaned$escherichia__mlst_achtman__ST)
table(kleborate_escherichia_cleaned$escherichia__mlst_achtman__clonal_complex)
#~~~~~~~~~~~~#

# merge ecoli kleborate data
all_ecoli_partition <- all_ecoli_partition |>
  mutate(run = gsub("microbesng_failed", "microbesng_failed_assembly", run))
# merge onto metadata
kleborate_escherichia_cleaned <- kleborate_escherichia_cleaned |>
  left_join(all_ecoli_partition |> select(c(sample, run, Level.1, Level.2, Level.3)), by = c("sample" = "sample", "run"= "run"))
# merge onto amrfinder
nrow(kleborate_escherichia_cleaned)
#View(kleborate_escherichia_cleaned)
sum(is.na(kleborate_escherichia_cleaned$Level.1)) # 0
sum(is.na(kleborate_escherichia_cleaned$escherichia__mlst_achtman__ST)) # 0
table(kleborate_escherichia_cleaned$Level.1, kleborate_escherichia_cleaned$escherichia__mlst_achtman__clonal_complex) # ST131 is split into 2 partitions, but otherwise perfect correspondence
table(kleborate_escherichia_cleaned$Level.1, kleborate_escherichia_cleaned$escherichia__mlst_achtman__ST) # ST131 is split into 2 partitions, but otherwise perfect correspondence
table(kleborate_escherichia_cleaned$Level.2, kleborate_escherichia_cleaned$escherichia__mlst_achtman__clonal_complex) # this already splits clonal complexes
table(kleborate_escherichia_cleaned$Level.2, kleborate_escherichia_cleaned$escherichia__mlst_achtman__ST) # ST131 is split into 2 partitions, but otherwise perfect correspondence

kleborate_escherichia_cleaned <- kleborate_escherichia_cleaned |>
  mutate(run = gsub("microbesng_failed_assembly", "microbesng_failed", run))

#save
write.csv(kleborate_escherichia_cleaned, "kleborate_escherichia_cleaned_with_fastbaps_clusters.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in ####
#kleborate_escherichia_cleaned <- read.csv("kleborate_escherichia_cleaned_with_fastbaps_clusters.csv")

#~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#
# * Entry point to read in cleaned df ####
kleborate_kpsc_cleaned <- read.csv("kleborate_kpsc_cleaned.csv")
kleborate_kosc_cleaned <- read.csv("kleborate_kosc_cleaned.csv")

#View(kleborate_kpsc_cleaned)
#View(kleborate_kosc_cleaned)
length(unique(kleborate_kpsc_cleaned$sample)) #406
length(unique(kleborate_kosc_cleaned$sample)) #91
colnames(kleborate_kpsc_cleaned)
colnames(kleborate_kosc_cleaned)
table(kleborate_kpsc_cleaned$klebsiella_pneumo_complex__mlst__ST) # no clonal complex column
intersect(colnames(kleborate_kpsc_cleaned), colnames(kleborate_kosc_cleaned))
setdiff(colnames(kleborate_kpsc_cleaned), colnames(kleborate_kosc_cleaned))
setdiff(colnames(kleborate_kosc_cleaned), colnames(kleborate_kpsc_cleaned))

intersect(kleborate_kosc_cleaned$klebsiella_oxytoca_complex__mlst__ST, kleborate_kpsc_cleaned$klebsiella_pneumo_complex__mlst__ST)

# get kpsc df into correct format
kleborate_kpsc <- kleborate_kpsc_cleaned |>
   mutate(klebsiella_mlst_ST = paste0("kpsc ", klebsiella_pneumo_complex__mlst__ST),
         klebsiella_mlst_clonal_complex = paste0("kpsc ", NA),
         klebsiella_mlst_gene1 = paste0("kpsc gapA(", klebsiella_pneumo_complex__mlst__gapA, ")"),
         klebsiella_mlst_gene2 = paste0("kpsc infB(", klebsiella_pneumo_complex__mlst__infB, ")"),
         klebsiella_mlst_gene3 = paste0("kpsc mdh(", klebsiella_pneumo_complex__mlst__mdh, ")"),
         klebsiella_mlst_gene4 = paste0("kpsc pgi(", klebsiella_pneumo_complex__mlst__pgi, ")"),
         klebsiella_mlst_gene5_1 = paste0("kpsc phoE(", klebsiella_pneumo_complex__mlst__phoE, ")"),
         klebsiella_mlst_gene5_2 = paste0("kpsc rpoB(", klebsiella_pneumo_complex__mlst__rpoB, ")"),
         klebsiella_mlst_gene6 = paste0("kpsc tonB(", klebsiella_pneumo_complex__mlst__tonB, ")")
         ) |>
  select(c(sample, run, sampletype, enterobacterales__species__species, enterobacterales__species__species_match, 
           general__contig_stats__contig_count, general__contig_stats__N50, general__contig_stats__largest_contig, 
           general__contig_stats__total_size, general__contig_stats__ambiguous_bases, general__contig_stats__QC_warnings,
           klebsiella_mlst_ST, klebsiella_mlst_clonal_complex, klebsiella_mlst_gene1, klebsiella_mlst_gene2, klebsiella_mlst_gene3,
           klebsiella_mlst_gene4, klebsiella_mlst_gene5_1, klebsiella_mlst_gene5_2, klebsiella_mlst_gene6,
           everything())) |>
  select(-c(klebsiella_pneumo_complex__mlst__ST, klebsiella_pneumo_complex__mlst__gapA, klebsiella_pneumo_complex__mlst__infB,
            klebsiella_pneumo_complex__mlst__mdh, klebsiella_pneumo_complex__mlst__pgi, klebsiella_pneumo_complex__mlst__phoE, 
            klebsiella_pneumo_complex__mlst__rpoB, klebsiella_pneumo_complex__mlst__tonB))

colnames(kleborate_kpsc)

intersect(colnames(kleborate_kpsc_cleaned), colnames(kleborate_kosc_cleaned))
setdiff(colnames(kleborate_kpsc), colnames(kleborate_kosc_cleaned))
setdiff(colnames(kleborate_kosc_cleaned), colnames(kleborate_kpsc))

# get kosc df into correct format 
kleborate_kosc <- kleborate_kosc_cleaned |>
  # make standardised mlst columns
  mutate(klebsiella_mlst_ST = paste0("kosc ", klebsiella_oxytoca_complex__mlst__ST),
         klebsiella_mlst_clonal_complex = paste0("kosc ", klebsiella_oxytoca_complex__mlst__clonal_complex),
         klebsiella_mlst_gene1 = paste0("kosc gapA(", klebsiella_oxytoca_complex__mlst__gapA, ")"),
         klebsiella_mlst_gene2 = paste0("kosc infB(", klebsiella_oxytoca_complex__mlst__infB, ")"),
         klebsiella_mlst_gene3 = paste0("kosc mdh(", klebsiella_oxytoca_complex__mlst__mdh, ")"),
         klebsiella_mlst_gene4 = paste0("kosc pgi(", klebsiella_oxytoca_complex__mlst__pgi, ")"),
         klebsiella_mlst_gene5_1 = paste0("kosc phoE(", klebsiella_oxytoca_complex__mlst__phoE, ")"),
         klebsiella_mlst_gene5_2 = paste0("kosc rpoB(", klebsiella_oxytoca_complex__mlst__rpoB, ")"),
         klebsiella_mlst_gene6 = paste0("kosc tonB(", klebsiella_oxytoca_complex__mlst__tonB, ")")
         ) |>
  # remove other columns
  select(-c(klebsiella_oxytoca_complex__mlst__ST, klebsiella_oxytoca_complex__mlst__clonal_complex, klebsiella_oxytoca_complex__mlst__gapA, klebsiella_oxytoca_complex__mlst__infB,
            klebsiella_oxytoca_complex__mlst__mdh, klebsiella_oxytoca_complex__mlst__pgi, klebsiella_oxytoca_complex__mlst__phoE, 
            klebsiella_oxytoca_complex__mlst__rpoB, klebsiella_oxytoca_complex__mlst__tonB)) |>
  select(c(sample, run, sampletype, enterobacterales__species__species, enterobacterales__species__species_match, 
           general__contig_stats__contig_count, general__contig_stats__N50, general__contig_stats__largest_contig, 
           general__contig_stats__total_size, general__contig_stats__ambiguous_bases, general__contig_stats__QC_warnings,
           klebsiella_mlst_ST, klebsiella_mlst_clonal_complex, klebsiella_mlst_gene1, klebsiella_mlst_gene2, klebsiella_mlst_gene3,
           klebsiella_mlst_gene4, klebsiella_mlst_gene5_1, klebsiella_mlst_gene5_2, klebsiella_mlst_gene6,
           everything())) 

# check which columns still missing
kpsc_only_cols <- setdiff(colnames(kleborate_kpsc), colnames(kleborate_kosc))
setdiff(colnames(kleborate_kosc), colnames(kleborate_kpsc)) # 0

# add them to kosc and fill with NA
kleborate_kosc[, kpsc_only_cols] <- NA

# bind together
kleborate_kpsc_kosc <- rbind(kleborate_kpsc, kleborate_kosc)
#View(kleborate_kpsc_kosc)
nrow(kleborate_kpsc_kosc)
length(unique(kleborate_kpsc_kosc$sample)) # 497
kpsc_kosc_samples <- kleborate_kpsc_kosc$sample
  
table(amrfinder_metadata_updated$typing_scheme, amrfinder_metadata_updated$kraken2_species, useNA = "ifany")
table(amrfinder_metadata_updated$typing_scheme_cleaned, amrfinder_metadata_updated$kraken2_species, useNA = "ifany")
table(amrfinder_metadata_updated$typing_scheme)
table(amrfinder_metadata_updated$typing_scheme)
  

raoultella <- amrfinder_metadata_updated |>
  filter(grepl("Raoultella", kraken2_species)) |>
  group_by(sample, run) |>
  slice_head() |> 
  pull(sample, sampletype_cat)
raoultella # both missing form the kleborate kpsc and kosc outputs :()


# filter non-kleb oxy and non-kleb pneumo from the amrfinder metadataframe
all_kleb__metadata <- amrfinder_metadata_updated |>
  group_by(sample, run) |>
  arrange(desc(size)) |>
  slice_head() |>
  ungroup() |>
  filter(typing_scheme %in% c("kaerogenes", "klebsiella", "koxytoca")) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
nrow(all_kleb__metadata) # 520
length(unique(all_kleb__metadata$sample)) # 520
table(all_kleb__metadata$typing_scheme)
#kaerogenes klebsiella   koxytoca 
#21        406         93 

# which are the 2 that are in this df but not in the clustering df
setdiff(all_kleb__metadata$sample, all_kleb_partition$sample) #  "AJ269" "AJ278" # these are the 2 Raoultella (should be under koxytoca)

# get kleb aerogenes into correct format
kleb_aerogenes_metadata <- all_kleb__metadata |>
  filter(typing_scheme == "kaerogenes")
#View(kleb_aerogenes_metadata) # 21
nrow(kleb_aerogenes_metadata) # 21
colnames(kleborate_kpsc_kosc)
colnames(kleb_aerogenes_metadata)
table(kleb_aerogenes_metadata$kraken2_species, useNA = "ifany")
table(kleb_aerogenes_metadata$sampletype, useNA = "ifany")

# get kleb aerogenes into correct format
kleborate_kaerogenes <- kleb_aerogenes_metadata |>
  mutate(enterobacterales__species__species = kraken2_species,
         enterobacterales__species__species_match = NA , 
         general__contig_stats__contig_count = Total_Contigs ,
         general__contig_stats__N50 = Contig_N50, 
         general__contig_stats__largest_contig = Max_Contig_Length,
         general__contig_stats__total_size = Genome_Size ,
         general__contig_stats__ambiguous_bases = "no",
         general__contig_stats__QC_warnings =  "-",
         klebsiella_mlst_ST = case_when(mlst != "-" ~ paste0("ST", mlst),
                                        TRUE ~ mlst),
         klebsiella_mlst_clonal_complex = "kaerogenes",
         klebsiella_mlst_gene1 = gene1 ,
         klebsiella_mlst_gene2 = gene2,
         klebsiella_mlst_gene3 = gene3,
         klebsiella_mlst_gene4 = gene4,
         klebsiella_mlst_gene5_1 = gene5_1,
         klebsiella_mlst_gene5_2 = gene5_2,
         klebsiella_mlst_gene6 = gene6
         ) |>
  select(c(sample, run, sampletype, enterobacterales__species__species, enterobacterales__species__species_match, 
           general__contig_stats__contig_count, general__contig_stats__N50, general__contig_stats__largest_contig, 
           general__contig_stats__total_size, general__contig_stats__ambiguous_bases, general__contig_stats__QC_warnings,
           klebsiella_mlst_ST, klebsiella_mlst_clonal_complex, klebsiella_mlst_gene1, klebsiella_mlst_gene2, klebsiella_mlst_gene3,
           klebsiella_mlst_gene4, klebsiella_mlst_gene5_1, klebsiella_mlst_gene5_2, klebsiella_mlst_gene6))
  
# get missing kaero cols
kaero_missing_cols <- setdiff(colnames(kleborate_kpsc_kosc), colnames(kleborate_kaerogenes))
# add other columns
kleborate_kaerogenes[, kaero_missing_cols] <- NA
length(kleborate_kaerogenes)  #116
length(kleborate_kpsc_kosc)  #116
setdiff(colnames(kleborate_kaerogenes), colnames(kleborate_kpsc_kosc)) #0
setdiff( colnames(kleborate_kpsc_kosc), colnames(kleborate_kaerogenes)) #0
  
# bind rows together
kleborate_all_kleb <- rbind(kleborate_kpsc_kosc, kleborate_kaerogenes)
nrow(kleborate_all_kleb) # 518

# merge wit fastbasp cluster info
all_kleb_partition <- all_kleb_partition |>
  mutate(run = gsub("microbesng_failed", "microbesng_failed_assembly", run))
#View(all_kleb_partition)

kleborate_all_kleb_with_fastbaps <- kleborate_all_kleb |>
  left_join(all_kleb_partition |> select(c(sample, run, `Level 1`, `Level 2`, `Level 3`)), by = c("sample" = "sample", "run"= "run"))
nrow(kleborate_all_kleb_with_fastbaps) # 518
#View(kleborate_all_kleb_with_fastbaps)
sum(is.na(kleborate_all_kleb_with_fastbaps$`Level 1`)) # 0 NA
sum(is.na(kleborate_all_kleb_with_fastbaps$general__contig_stats__contig_count)) # 0 NA

table(kleborate_all_kleb_with_fastbaps$`Level 1`, kleborate_all_kleb_with_fastbaps$klebsiella_mlst_clonal_complex) # kosc clonal complexes probably mean very little- 2 Kosc,  3 kpsc, 1 kaerogenes
table(kleborate_all_kleb_with_fastbaps$`Level 1`, kleborate_all_kleb_with_fastbaps$klebsiella_mlst_ST) # 
table(kleborate_all_kleb_with_fastbaps$`Level 2`, kleborate_all_kleb_with_fastbaps$klebsiella_mlst_clonal_complex) # kosc in 6, kaerogenes n 2 (with 1 singleton), and kpsc across 7, also incl 1 singleton
table(kleborate_all_kleb_with_fastbaps$`Level 2`, kleborate_all_kleb_with_fastbaps$klebsiella_mlst_ST) # 

# modify mlst column so no "-" characters
kleborate_all_kleb_with_fastbaps <- kleborate_all_kleb_with_fastbaps |>
  mutate(klebsiella_mlst_ST = case_when(klebsiella_mlst_ST  == "-" ~ paste0(klebsiella_mlst_clonal_complex, "|", 
                                                                            klebsiella_mlst_gene1, "|",
                                                                            klebsiella_mlst_gene2, "|",
                                                                            klebsiella_mlst_gene3, "|",
                                                                            klebsiella_mlst_gene4, "|",
                                                                            klebsiella_mlst_gene5_1, "|",
                                                                            klebsiella_mlst_gene5_2, "|",
                                                                            klebsiella_mlst_gene6, "|"),
                                        TRUE ~ klebsiella_mlst_ST))|>
  mutate(run = gsub("microbesng_failed_assembly", "microbesng_failed", run))

#save
write.csv(kleborate_all_kleb_with_fastbaps, "kleborate_all_kleb_with_fastbaps.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in ####
#kleborate_all_kleb_with_fastbaps <- read.csv("kleborate_all_kleb_with_fastbaps.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# merge ecoli bsi data with kleborate mlsts, and fastbaps partitions 
colnames(ecoli_bsi_samples_metadata)
colnames(kleborate_escherichia_cleaned)

setdiff(ecoli_bsi_samples_metadata$sequencing_id, kleborate_escherichia_cleaned$sample)
setdiff( kleborate_escherichia_cleaned$sample, ecoli_bsi_samples_metadata$sequencing_id)


ecoli_bsi_samples_metadata <- ecoli_bsi_samples_metadata |>
  left_join(kleborate_escherichia_cleaned |> select(c(sample, run,
                                                      escherichia__mlst_achtman__ST, 
                                                      escherichia__mlst_achtman__clonal_complex,
                                                      escherichia__mlst_achtman__adk,
                                                      escherichia__mlst_achtman__fumC,
                                                      escherichia__mlst_achtman__gyrB,
                                                      escherichia__mlst_achtman__icd,
                                                      escherichia__mlst_achtman__mdh,
                                                      escherichia__mlst_achtman__purA,
                                                      escherichia__mlst_achtman__recA)), by = c("sequencing_id" = "sample", "run" = "run"))

colnames(ecoli_bsi_samples_metadata)
nrow(ecoli_bsi_samples_metadata)
#View(ecoli_bsi_samples_metadata)
sum(is.na(ecoli_bsi_samples_metadata$escherichia__mlst_achtman__ST))
sum(is.na(ecoli_bsi_samples_metadata$escherichia__mlst_achtman__clonal_complex)) # 367
# merge kleb bsi data with kleborate mlsts, and fastbaps partitions 
table(ecoli_bsi_samples_metadata$mlst, ecoli_bsi_samples_metadata$escherichia__mlst_achtman__ST , useNA = "ifany")

#save
write.csv(ecoli_bsi_samples_metadata, "rarefaction/ecoli_bsi_samples_metadata.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in ####
#ecoli_bsi_samples_metadata <- read.csv("rarefaction/ecoli_bsi_samples_metadata.csv")
#~~~~~~~~~~~~~~~#

# also merge kleb kleborate mlst data
colnames(kleb_bsi_samples_metadata)
colnames(kleborate_all_kleb_with_fastbaps)
#setdiff(kleb_bsi_samples_metadata$sequencing_id, kleborate_all_kleb_with_fastbaps$sample) # 0
setdiff( kleborate_all_kleb_with_fastbaps$sample, kleb_bsi_samples_metadata$sequencing_id) # non-BSI samples missing here

kleb_bsi_samples_metadata <- kleb_bsi_samples_metadata |>
  left_join(kleborate_all_kleb_with_fastbaps |> select(c(sample, run,
                                                         klebsiella_mlst_ST, 
                                                         klebsiella_mlst_clonal_complex,
                                                         klebsiella_mlst_gene1,
                                                         klebsiella_mlst_gene2,
                                                         klebsiella_mlst_gene3,
                                                         klebsiella_mlst_gene4,
                                                         klebsiella_mlst_gene5_1,
                                                         klebsiella_mlst_gene5_2,
                                                         klebsiella_mlst_gene6)), 
            by = c("sequencing_id" = "sample", "run" = "run"))

colnames(kleb_bsi_samples_metadata)
nrow(kleb_bsi_samples_metadata) # 468
#View(kleb_bsi_samples_metadata)
sum(is.na(kleb_bsi_samples_metadata$klebsiella_mlst_ST)) # 0
sum(is.na(kleb_bsi_samples_metadata$klebsiella_mlst_clonal_complex)) #0
#table(kleb_bsi_samples_metadata$mlst, kleb_bsi_samples_metadata$klebsiella_mlst_ST , useNA = "ifany")
table(kleb_bsi_samples_metadata$klebsiella_mlst_clonal_complex, kleb_bsi_samples_metadata$kraken2_species, useNA = "ifany")
table(kleb_bsi_samples_metadata$klebsiella_mlst_clonal_complex, kleb_bsi_samples_metadata$typing_scheme, useNA = "ifany")

#save
write.csv(kleb_bsi_samples_metadata, "rarefaction/kleb_bsi_samples_metadata.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# read in ####
#kleb_bsi_samples_metadata <- read.csv("rarefaction/kleb_bsi_samples_metadata.csv")
#~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Frequency distribution histograms of MLSTs, fastbaps clusters (level 1,2,3) ARGs, and plasmids (inc; pling plasmid subcommunities) ####
 # Set df and column names
# Your Escherichia objects (you already have these names)
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

# 1) Counts histogram (counts on x, number of clusters/STs on y)
counts_plot <- ggplot(all_counts, aes(x = value, fill = genus, colour = NULL)) +
  geom_histogram(bins = n_bins, position = "identity", alpha = 0.6, closed = "right") +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  scale_fill_manual(values = genus_cols, name = "Genus") +
  labs(x = "Number of isolates per unit", y = "Number of unique units") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")
counts_plot
ggsave("panel_counts_histograms_six_panel.png", counts_plot, width = 6, height = 6, dpi = 300)

# 1.b) Counts histogram (counts on x, number of clusters/STs on y) - only fastbaps L2
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
freq_plot_by_genus <- ggplot(all_counts, aes(x = value, fill = genus, group = genus)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins = n_bins,
                 position = "identity", alpha = 0.6, colour = NA) +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(Escherichia = "seagreen3", Klebsiella = "darkorange")) +
  labs(x = "Size (number of isolates per unit)", y = "Proportion of units (within genus)", fill = "Genus") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))
freq_plot_by_genus
# not sure the frequencies are correct here
ggsave("panel_freq_histograms_six_panel.png", freq_plot_by_genus, width = 6, height = 6, dpi = 300)

# 2.b) Frequency histogram (counts on x, number of clusters/STs on y)
freq_plot_by_genus <- ggplot(all_counts[!all_counts$metric %in% c("FastBAPS clusters L1", "FastBAPS clusters L2"), ], aes(x = value, fill = genus, group = genus)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins = n_bins,
                 position = "identity", alpha = 0.6, colour = NA) +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(Escherichia = "seagreen3", Klebsiella = "darkorange")) +
  labs(x = "Size (number of isolates per unit)", y = "Proportion of units (within genus)", fill = "Genus") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))
freq_plot_by_genus
# not sure the frequencies are correct here
ggsave("panel_freq_histograms_four_panel.png", freq_plot_by_genus, width = 6, height = 4.5, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# sankey diagram to visualse concordance of MLSTs (mlst and kleborate) and fastBAPS clusters
# install if needed
# install.packages(c("ggplot2","ggalluvial","dplyr"))
library(dplyr)
library(ggplot2)
install.packages("ggalluvial")
library(ggalluvial)


# for ggalluvial, provide a data.frame where each row is an observation and columns are axes
df2 <- ecoli_bsi_samples_metadata
colnames(df2)
df2$Level.1 <- as.factor(df2$Level.1)
df2$Level.2 <- as.factor(df2$Level.2)
df2$Level.3 <- as.factor(df2$Level.3)

# Optionally aggregate identical trajectories and get counts
df_counts <- df2 %>%
  count(escherichia__mlst_achtman__clonal_complex, 
        escherichia__mlst_achtman__ST,
        mlst,
        Level.1, Level.2, Level.3, name = "freq") |>
  mutate(alluvium_id = interaction(
    escherichia__mlst_achtman__clonal_complex,
    escherichia__mlst_achtman__ST,
    mlst,
    Level.1, Level.2, Level.3,
    drop = TRUE
  ))
#View(df_counts)

# Draw alluvial
ecoli_mlst_sankey <- ggplot(df_counts,
       aes(axis1 = escherichia__mlst_achtman__clonal_complex, axis2 = escherichia__mlst_achtman__ST, axis3 = Level.1, axis4 = Level.2, axis5 = Level.3,
           y = freq)) +
  scale_x_discrete(labels = c("Clonal\nComplex","MLST (Kleborate)", "MLST","Level 1","Level 2","Level 3"),
                   expand = c(.05, .05)) +
  geom_alluvium(aes(fill = Level.1), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
 # geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_viridis_d(option = "viridis") +  # or "viridis", "plasma", "magma"  
  theme_minimal() +
  labs(y = "Count", x = "")
ecoli_mlst_sankey
ggsave("ecoli_mlst_fastbaps_sankey.png", ecoli_mlst_sankey, width = 15, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# repeat for Klebsiella
# for ggalluvial, provide a data.frame where each row is an observation and columns are axes
df2 <- kleb_bsi_samples_metadata
colnames(df2)
df2$Level.1 <- as.factor(df2$Level.1)
df2$Level.2 <- as.factor(df2$Level.2)
df2$Level.3 <- as.factor(df2$Level.3)

# Optionally aggregate identical trajectories and get counts
df_counts <- df2 |>
  count(klebsiella_mlst_clonal_complex, 
        klebsiella_mlst_ST,
        mlst,
        Level.1, Level.2, Level.3, name = "freq")
df_counts <- df_counts |>
  arrange(klebsiella_mlst_ST, Level.3)
#View(df_counts)

# Draw alluvial
kleb_mlst_sankey <- ggplot(df_counts,
                            aes(axis1 = klebsiella_mlst_clonal_complex, axis2 = klebsiella_mlst_ST, axis3 = Level.1, axis4 = Level.2, axis5 = Level.3,
                                y = freq)) +
  scale_x_discrete(labels = c("Clonal\nComplex","MLST","Level 1","Level 2","Level 3"),
                   expand = c(.05, .05)) +
  geom_alluvium(aes(fill = Level.3), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  # geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_viridis_d(option = "plasma") +  # or "viridis", "plasma", "magma"  
  
  theme_minimal() +
  labs(y = "Count", x = "")
kleb_mlst_sankey

ggsave("kleb_mlst_fastbaps_sankey.png", kleb_mlst_sankey, width = 16, height = 10, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * MLST rarefaction E.coli ####

# Function that performs rarefaction from NEKSUS sample and calculates species richness, shannon and simpson diversity
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

# Set Parameters
n_reps <- 100  # number of bootstrap replicates
max_samples <- length(unique(ecoli_bsi_samples_metadata$isolateid)) # 1471

# Run rarefaction with and without replacement
# with replacement
ecoli_bsi_mlst_rarefaction_replacement <- rarefy_species(df = ecoli_bsi_samples_metadata,
                                feature = "mlst_profile", max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = 12, 
                                #seed = 42
                                )

# without replacement
ecoli_bsi_mlst_rarefaction_no_replacement <- rarefy_species(df = ecoli_bsi_samples_metadata,
                               feature = "mlst_profile", max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = 12, 
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
total_mlst_diversity <- length(unique(ecoli_bsi_samples_metadata$mlst_profile)) #278
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
# save output table
write.csv(ecoli_bsi_mlst_rarefaction_combined, file = "rarefaction/ecoli_bsi_mlst_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#ecoli_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_mlst_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~#

missing<- amrfinder_metadata_updated |>
  filter(sample %in% c("AG111", "AG130", "AG139"))
View(missing)
View(amrfinder_metadata_updated)
table(missing$duplicate_assembly_qc_pass)

#~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~#

#reorder metircs:
ecoli_bsi_mlst_rarefaction_combined$Metric <- factor(ecoli_bsi_mlst_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
ecoli_bsi_mlst_rarefaction_combined$Replacement <- factor(ecoli_bsi_mlst_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# define custom colours
pal <- c("grey90", "seagreen3" )

# plot 2 graphs
p <- ggplot(data = ecoli_bsi_mlst_rarefaction_combined, aes(x = sample_size, y = est, fill = Replacement, colour = Replacement)) +
  geom_line( linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3, colour = NA) +
  # scale viridis colours
  scale_colour_manual(name = "Replacement", values = pal) +
  scale_fill_manual(name = "Replacement", values = pal) +
  facet_wrap(~Metric, ncol = 4, scales = "free_y") +
   theme_minimal() +
  labs(title = "E.coli BSI MLST rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)

ggsave("rarefaction/ecoli_bsi_mlst_rarefaction_combined.png", plot = p, width = 14, height = 5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Key Gene Rarefaction ####
# function to perform rarefaction using isolate featues like ARGs or plasmids
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
table(ecoli_bsi_arg$Type, useNA = "ifany")
table(ecoli_bsi_amrfinder_metadata$Type, useNA = "ifany")
# check all isolates present in bith df
length(unique(ecoli_bsi_amrfinder_metadata$sample))
length(unique(ecoli_bsi_arg$sample))

# define run params
max_samples <- length(unique(ecoli_bsi_amrfinder_metadata$sample)) # 1471
n_reps <- 1000
n_cores <- 12

# Run with replacement
ecoli_bsi_arg_rarefaction_replacement <- rarefy_genes(ecoli_bsi_arg, sample_col = "sample", gene_col = "Element.symbol",
                          max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = n_cores #, seed = 1
                          )
# run without replacement
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
# save output table
write.csv(ecoli_bsi_arg_rarefaction_combined, file = "rarefaction/ecoli_bsi_arg_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#ecoli_bsi_arg_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_arg_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#

#reorder metircs:
ecoli_bsi_arg_rarefaction_combined$Metric <- factor(ecoli_bsi_arg_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
ecoli_bsi_arg_rarefaction_combined$Replacement <- factor(ecoli_bsi_arg_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# define custom colours
pal <- c("grey90", "seagreen3" )

# plot 2 graphs
p <- ggplot(data = ecoli_bsi_arg_rarefaction_combined, aes(x = sample_size, y = est, fill = Replacement, colour = Replacement)) +
  geom_line( linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, colour = NA) +
  # scale viridis colours
  scale_colour_manual(name = "Replacement", values = pal) +
  scale_fill_manual(name = "Replacement", values = pal) +
  facet_wrap(~Metric, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "E.coli BSI ARG rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
#save
ggsave("rarefaction/ecoli_bsi_arg_rarefaction_combined.png", plot = p, width = 14, height = 5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Plasmid Inc/rep types Rarefaction ####
# fix rep types - change non-plasmids to NAs for better handling by feature rarefaction function
ecoli_bsi_amrfinder_metadata <- ecoli_bsi_amrfinder_metadata |>
  mutate(rep_types_whole_plasmid = gsub("-,", "", rep_types_whole_plasmid)) |>
  mutate(rep_types_whole_plasmid = case_when(rep_types_whole_plasmid == "-" ~ NA,
                                             TRUE ~ rep_types_whole_plasmid))

ecoli_bsi_amrfinder_metadata$rep_types_whole_plasmid
  
# define run params
max_samples <- length(unique(ecoli_bsi_amrfinder_metadata$sample)) # 1471
n_reps <- 1000
n_cores <- 12

# Run with replacement
ecoli_bsi_inc_rarefaction_replacement <- rarefy_genes(ecoli_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "rep_types_whole_plasmid",
                                                      max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = n_cores #, seed = 1
)
# run without replacement
ecoli_bsi_inc_rarefaction_no_replacement <- rarefy_genes(ecoli_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "rep_types_whole_plasmid",
                                                         max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = n_cores #, seed = 1
)

# Clean, rescale and combine tables
# rbind 2 dfs after adding column to identify replacement vs no replacement
ecoli_bsi_inc_rarefaction_replacement <- ecoli_bsi_inc_rarefaction_replacement |>
  mutate(Replacement = "With replacement")
ecoli_bsi_inc_rarefaction_no_replacement <- ecoli_bsi_inc_rarefaction_no_replacement |>
  mutate(Replacement = "No replacement")
ecoli_bsi_inc_rarefaction_combined <- rbind(ecoli_bsi_inc_rarefaction_replacement,
                                            ecoli_bsi_inc_rarefaction_no_replacement)


# combine dfs
# add proportion of diversity captured to output df
total_inc_diversity <- length(unique(ecoli_bsi_amrfinder_metadata[!is.na(ecoli_bsi_amrfinder_metadata$rep_types_whole_plasmid), ]$rep_types_whole_plasmid)) #331
max_shannon <- round(max(ecoli_bsi_inc_rarefaction_replacement[ecoli_bsi_inc_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

ecoli_bsi_inc_rarefaction_combined <- ecoli_bsi_inc_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_inc_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor)|>
  mutate(Genus = "Escherichia")
#View(ecoli_bsi_inc_rarefaction_combined)
# save output table
write.csv(ecoli_bsi_inc_rarefaction_combined, file = "rarefaction/ecoli_bsi_inc_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#ecoli_bsi_inc_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_inc_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#
# Plot
#reorder metircs:
ecoli_bsi_inc_rarefaction_combined$Metric <- factor(ecoli_bsi_inc_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
ecoli_bsi_inc_rarefaction_combined$Replacement <- factor(ecoli_bsi_inc_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# define custom colours
pal <- c("grey90", "seagreen3" )

# plot 2 graphs
p <- ggplot(data = ecoli_bsi_inc_rarefaction_combined, aes(x = sample_size, y = est, fill = Replacement, colour = Replacement)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, colour = NA) +
  # scale viridis colours
  scale_colour_manual(name = "Replacement", values = pal) +
  scale_fill_manual(name = "Replacement", values = pal) +
  facet_wrap(~Metric, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "E.coli BSI inc rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
#save
ggsave("rarefaction/ecoli_bsi_inc_rarefaction_combined.png", plot = p, width = 14, height = 5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Plasmid subcommunities Rarefaction ####
# check column
ecoli_bsi_amrfinder_metadata$community_subcommunity

# define run params
max_samples <- length(unique(ecoli_bsi_amrfinder_metadata$sample)) # 1471
n_reps <- 1000
n_cores <- 12

# Run with replacement
ecoli_bsi_pling_rarefaction_replacement <- rarefy_genes(ecoli_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "community_subcommunity",
                                                      max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = n_cores #, seed = 1
)
# run without replacement
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
# save output table
write.csv(ecoli_bsi_pling_rarefaction_combined, file = "rarefaction/ecoli_bsi_pling_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#ecoli_bsi_pling_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_pling_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#
# Plot
#reorder metircs:
ecoli_bsi_pling_rarefaction_combined$Metric <- factor(ecoli_bsi_pling_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
ecoli_bsi_pling_rarefaction_combined$Replacement <- factor(ecoli_bsi_pling_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# define custom colours
pal <- c("grey90", "seagreen3" )

# plot 2 graphs
p <- ggplot(data = ecoli_bsi_pling_rarefaction_combined, aes(x = sample_size, y = est, fill = Replacement, colour = Replacement)) +
  geom_line( linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, colour = NA) +
  # scale viridis colours
  scale_colour_manual(name = "Replacement", values = pal) +
  scale_fill_manual(name = "Replacement", values = pal) +
  facet_wrap(~Metric, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "E.coli BSI pling rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
#save
ggsave("rarefaction/ecoli_bsi_pling_rarefaction_combined.png", plot = p, width = 14, height = 5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Klebsiella BSIs Overall Rarefaction ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Quick look at frequency distributions for MLSTs, ARGs, and plasmids
# plot histogram of species frequency
kleb_bsi_hist_data <- kleb_bsi_samples_metadata |>
  group_by(mlst_profile) |>
  summarise(count = n())

ggplot(data = kleb_bsi_hist_data, aes(x = count)) +
  geom_histogram()

# get amr gene frequency histogram
kleb_bsi_arg_hist_data <- kleb_bsi_amrfinder_metadata |>
  filter(!is.na(Element.symbol)) |>
  filter(Type == "AMR") |>
  group_by(Element.symbol) |>
  summarise(count = n())
#View(kleb_bsi_arg_hist_data)
ggplot(data = kleb_bsi_arg_hist_data, aes(x = count)) +
  geom_histogram()

# plasmid Inc type frequency plot
kleb_bsi_inc_hist_plot <- kleb_bsi_amrfinder_metadata |>
  filter(!is.na(rep_types_whole_plasmid) & rep_types_whole_plasmid != "-") |>
  group_by(sample, run, Contig.id) |> 
  slice_head() |>
  ungroup() |>
  mutate(rep_types_whole_plasmid = gsub("-,", "", rep_types_whole_plasmid)) |>
  group_by(rep_types_whole_plasmid) |>
  summarise(count = n())
#View(kleb_bsi_inc_hist_plot)
ggplot(data = kleb_bsi_inc_hist_plot, aes(x = count)) +
  geom_histogram()
# also decays very quickly

# plasmid subcommunity frequency/ size plot
kleb_bsi_subcommunity_hist_plot <- kleb_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  group_by(sample, run, Contig.id) |> 
  slice_head() |>
  ungroup() |>
  group_by(community_subcommunity) |>
  summarise(count = n())
View(kleb_bsi_subcommunity_hist_plot)
ggplot(data = kleb_bsi_subcommunity_hist_plot, aes(x = count)) +
  geom_histogram()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Klebsiella Rarefaction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Parameters
n_reps <- 1000  # number of bootstrap replicates
max_samples <- length(unique(kleb_bsi_samples_metadata$isolateid)) #468
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * MLST rarefaction for Kleb BSI ####
kleb_bsi_mlst_rarefaction_replacement <- rarefy_species(df = kleb_bsi_samples_metadata,
                                                         feature = "mlst_profile", max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = 12, 
                                                         #seed = 42
                                                        )

# run rarefactionwithout replacement
kleb_bsi_mlst_rarefaction_no_replacement <- rarefy_species(df = kleb_bsi_samples_metadata,
                                                            feature = "mlst_profile", max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = 12, 
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
total_mlst_diversity <- length(unique(kleb_bsi_samples_metadata$mlst_profile)) #278
max_shannon <- round(max(kleb_bsi_mlst_rarefaction_replacement[kleb_bsi_mlst_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

kleb_bsi_mlst_rarefaction_combined <- kleb_bsi_mlst_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_mlst_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor) |>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_mlst_rarefaction_combined)
# save output table
write.csv(kleb_bsi_mlst_rarefaction_combined, file = "rarefaction/kleb_bsi_mlst_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#kleb_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_mlst_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#

#reorder metircs:
kleb_bsi_mlst_rarefaction_combined$Metric <- factor(kleb_bsi_mlst_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
kleb_bsi_mlst_rarefaction_combined$Replacement <- factor(kleb_bsi_mlst_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# define custom colours
pal <- c("grey90", "darkorange")


# plot 2 graphs
p <- ggplot(data = kleb_bsi_mlst_rarefaction_combined, aes(x = sample_size, y = est, fill = Replacement, colour = Replacement)) +
  geom_line( linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, colour = NA) +
  # scale viridis colours
  scale_colour_manual(name = "Replacement", values = pal) +
  scale_fill_manual(name = "Replacement", values = pal) +
  facet_wrap(~Metric, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "Klebsiella BSI MLST rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)

ggsave("rarefaction/kleb_bsi_mlst_rarefaction_combined.png", plot = p, width = 14, height = 5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Combined Kleb and E.coli plot ##### 
# MLSTs
# load data
#ecoli_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_mlst_rarefaction_combined.csv")
#kleb_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_mlst_rarefaction_combined.csv")

combined_bsi_mlst_rarefaction <- rbind(ecoli_bsi_mlst_rarefaction_combined, kleb_bsi_mlst_rarefaction_combined)
combined_bsi_mlst_rarefaction$Replacement <- factor(combined_bsi_mlst_rarefaction$Replacement, levels = c("With replacement", "No replacement"))
combined_bsi_mlst_rarefaction <- combined_bsi_mlst_rarefaction |>
  mutate(ReplacementGenus = interaction(Replacement, Genus, sep = "_"))

pal <- c(
  "No replacement_Escherichia"   = "seagreen3",
  "No replacement_Klebsiella"   = "darkorange",
  "With replacement_Escherichia"= "grey90",
  "With replacement_Klebsiella" = "grey90"
)

# plot combined data
p <- ggplot(data = combined_bsi_mlst_rarefaction, aes(x = sample_size, y = est, fill = ReplacementGenus, colour = ReplacementGenus)) +
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
    )
  ) +
  scale_fill_manual(
    name = "Sampling scheme",
    values = pal,
    labels = c(
      "No replacement_Escherichia"    = "No replacement – Escherichia",
      "No replacement_Klebsiella"    = "No replacement – Klebsiella",
      "With replacement_Escherichia" = "With replacement – Escherichia",
      "With replacement_Klebsiella"  = "With replacement – Klebsiella"
    )
  ) +
  facet_wrap(~ Metric , ncol = 2,  scales = "free") +
  theme_minimal(base_size = 16) +
  labs(title = "E. coli and Klebsiella BSI MLST rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
# save
#ggsave("rarefaction/combined_bsi_mlst_rarefaction.png", plot = p, width = 12, height = 7, units = "in", dpi = 300)
ggsave("rarefaction/combined_bsi_mlst_rarefaction_same_axes.png", plot = p, width = 12, height = 9, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * ARG diversity Klebsiella ####
# filter amrfinder df for AMR genes only
kleb_bsi_arg <- kleb_bsi_amrfinder_metadata |>
  filter(Type == "AMR"| is.na(Type))
table(kleb_bsi_arg$Type, useNA = "ifany")
table(kleb_bsi_amrfinder_metadata$Type, useNA = "ifany")
# check all isolates present in bith df
length(unique(kleb_bsi_amrfinder_metadata$sample))
length(unique(kleb_bsi_arg$sample))


# define run params
max_samples <- length(unique(kleb_bsi_amrfinder_metadata$sample)) # 1471
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
total_ARG_diversity <- length(unique(kleb_bsi_amrfinder_metadata[!is.na(kleb_bsi_amrfinder_metadata$Element.symbol), ]$Element.symbol)) #277
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
# save output table
write.csv(kleb_bsi_arg_rarefaction_combined, file = "rarefaction/kleb_bsi_arg_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#kleb_bsi_arg_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_arg_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#

#reorder metircs:
kleb_bsi_arg_rarefaction_combined$Metric <- factor(kleb_bsi_arg_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
kleb_bsi_arg_rarefaction_combined$Replacement <- factor(kleb_bsi_arg_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# define custom colours
#pal <- viridis(n = 6, option = "D")
"#440154FF" "#414487FF" "#2A788EFF" "#22A884FF" "#7AD151FF" "#FDE725FF"
#"#440154FF" "#2A788EFF" "#7AD151FF" 
#pal <- c( pal[4], pal[1], pal[1], pal[6]) # select 4 colours
"#7AD151FF" "#2A788EFF" "#440154FF" "#FDE725FF"
pal <- c("grey90", "darkorange" )


# plot 2 graphs
p <- ggplot(data = kleb_bsi_arg_rarefaction_combined, aes(x = sample_size, y = est, fill = Replacement, colour = Replacement)) +
  geom_line( linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3, colour = NA) +
  # scale viridis colours
  scale_colour_manual(name = "Replacement", values = pal) +
  scale_fill_manual(name = "Replacement", values = pal) +
  facet_wrap(~Metric, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "Klebsiella BSI ARG rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
#save
ggsave("rarefaction/kleb_bsi_arg_rarefaction_combined.png", plot = p, width = 14, height = 5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Combined Kleb and E.coli plot ####
# ARGs
# load data
#ecoli_bsi_arg_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_arg_rarefaction_combined.csv")
#kleb_bsi_arg_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_arg_rarefaction_combined.csv")

combined_bsi_arg_rarefaction <- rbind(ecoli_bsi_arg_rarefaction_combined, kleb_bsi_arg_rarefaction_combined)
combined_bsi_arg_rarefaction$Metric <- factor(combined_bsi_arg_rarefaction$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
combined_bsi_arg_rarefaction$Replacement <- factor(combined_bsi_arg_rarefaction$Replacement, levels = c("With replacement", "No replacement"))
combined_bsi_arg_rarefaction <- combined_bsi_arg_rarefaction |>
  mutate(ReplacementGenus = interaction(Replacement, Genus, sep = "_"))

pal <- c(
  "No replacement_Escherichia"   = "seagreen3",
  "No replacement_Klebsiella"   = "darkorange",
  "With replacement_Escherichia"= "grey90",
  "With replacement_Klebsiella" = "grey90"
)

# plot combined data
p <- ggplot(data = combined_bsi_arg_rarefaction, aes(x = sample_size, y = est, fill = ReplacementGenus, colour = ReplacementGenus)) +
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
    )
  ) +
  scale_fill_manual(
    name = "Sampling scheme",
    values = pal,
    labels = c(
      "No replacement_Escherichia"    = "No replacement – Escherichia",
      "No replacement_Klebsiella"    = "No replacement – Klebsiella",
      "With replacement_Escherichia" = "With replacement – Escherichia",
      "With replacement_Klebsiella"  = "With replacement – Klebsiella"
    )
  ) +
  facet_wrap(~ Metric , ncol = 2,  scales = "free") +
  theme_minimal(base_size = 16) +
  labs(title = "E. coli and Klebsiella BSI ARG rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
# save
#ggsave("rarefaction/combined_bsi_arg_rarefaction.png", plot = p, width = 12, height = 7, units = "in", dpi = 300)
ggsave("rarefaction/combined_bsi_arg_rarefaction_same_axes.png", plot = p, width = 12, height = 9, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Plasmid Inc types ####
# fix rep types
kleb_bsi_amrfinder_metadata <- kleb_bsi_amrfinder_metadata |>
  mutate(rep_types_whole_plasmid = gsub("-,", "", rep_types_whole_plasmid)) |>
  mutate(rep_types_whole_plasmid = case_when(rep_types_whole_plasmid == "-" ~ NA,
                                             TRUE ~ rep_types_whole_plasmid))

kleb_bsi_amrfinder_metadata$rep_types_whole_plasmid

# define run params
max_samples <- length(unique(kleb_bsi_amrfinder_metadata$sample)) # 468
n_reps <- 1000
n_cores <- 12

# Run
kleb_bsi_inc_rarefaction_replacement <- rarefy_genes(kleb_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "rep_types_whole_plasmid",
                                                      max_samples = max_samples, reps = n_reps, replacement = TRUE, n_cores = n_cores #, seed = 1
)
kleb_bsi_inc_rarefaction_no_replacement <- rarefy_genes(kleb_bsi_amrfinder_metadata, sample_col = "sample", gene_col = "rep_types_whole_plasmid",
                                                         max_samples = max_samples, reps = n_reps, replacement = FALSE, n_cores = n_cores #, seed = 1
)

# Clean, rescale and combine tables
# rbind 2 dfs after adding column to identify replacement vs no replacement
kleb_bsi_inc_rarefaction_replacement <- kleb_bsi_inc_rarefaction_replacement |>
  mutate(Replacement = "With replacement")
kleb_bsi_inc_rarefaction_no_replacement <- kleb_bsi_inc_rarefaction_no_replacement |>
  mutate(Replacement = "No replacement")
kleb_bsi_inc_rarefaction_combined <- rbind(kleb_bsi_inc_rarefaction_replacement,
                                            kleb_bsi_inc_rarefaction_no_replacement)


# combine dfs
# add proportion of diversity captured to output df
total_inc_diversity <- length(unique(kleb_bsi_amrfinder_metadata[!is.na(kleb_bsi_amrfinder_metadata$rep_types_whole_plasmid), ]$rep_types_whole_plasmid)) #331
max_shannon <- round(max(kleb_bsi_inc_rarefaction_replacement[kleb_bsi_inc_rarefaction_replacement$Metric == "Shannon diversity", ]$est))
max_simpson <- 1 # by definition
max_coverage <- 1 # by definition

kleb_bsi_inc_rarefaction_combined <- kleb_bsi_inc_rarefaction_combined |>
  mutate(scale_factor = case_when(Metric == "Species richness" ~ total_inc_diversity,
                                  Metric == "Shannon diversity" ~ max_shannon,
                                  TRUE ~ 1),
         est_scaled = est / scale_factor,
         lcl_scaled = lcl / scale_factor,
         ucl_scaled = ucl / scale_factor)|>
  mutate(Genus = "Klebsiella")
#View(kleb_bsi_inc_rarefaction_combined)
# save output table
write.csv(kleb_bsi_inc_rarefaction_combined, file = "rarefaction/kleb_bsi_inc_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#kleb_bsi_inc_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_inc_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#
# Plot
#reorder metircs:
kleb_bsi_inc_rarefaction_combined$Metric <- factor(kleb_bsi_inc_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
kleb_bsi_inc_rarefaction_combined$Replacement <- factor(kleb_bsi_inc_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# define custom colours
pal <- c("grey90", "darkorange" )

# plot 2 graphs
p <- ggplot(data = kleb_bsi_inc_rarefaction_combined, aes(x = sample_size, y = est, fill = Replacement, colour = Replacement)) +
  geom_line( linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, colour = NA) +
  # scale viridis colours
  scale_colour_manual(name = "Replacement", values = pal) +
  scale_fill_manual(name = "Replacement", values = pal) +
  facet_wrap(~Metric, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "E.coli BSI inc rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
#save
ggsave("rarefaction/kleb_bsi_inc_rarefaction_combined.png", plot = p, width = 14, height = 5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Combined Kleb and E.coli plot ##### 
# Plasmid Inc types
# load data
#ecoli_bsi_inc_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_inc_rarefaction_combined.csv")
#kleb_bsi_inc_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_inc_rarefaction_combined.csv")

combined_bsi_inc_rarefaction <- rbind(ecoli_bsi_inc_rarefaction_combined, kleb_bsi_inc_rarefaction_combined)
combined_bsi_inc_rarefaction$Replacement <- factor(combined_bsi_inc_rarefaction$Replacement, levels = c("With replacement", "No replacement"))
combined_bsi_inc_rarefaction <- combined_bsi_inc_rarefaction |>
  mutate(ReplacementGenus = interaction(Replacement, Genus, sep = "_"))

pal <- c(
  "No replacement_Escherichia"   = "seagreen3",
  "No replacement_Klebsiella"   = "darkorange",
  "With replacement_Escherichia"= "grey90",
  "With replacement_Klebsiella" = "grey90"
)

# plot combined data
p <- ggplot(data = combined_bsi_inc_rarefaction, aes(x = sample_size, y = est, fill = ReplacementGenus, colour = ReplacementGenus)) +
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
    )
  ) +
  scale_fill_manual(
    name = "Sampling scheme",
    values = pal,
    labels = c(
      "No replacement_Escherichia"    = "No replacement – Escherichia",
      "No replacement_Klebsiella"    = "No replacement – Klebsiella",
      "With replacement_Escherichia" = "With replacement – Escherichia",
      "With replacement_Klebsiella"  = "With replacement – Klebsiella"
    )
  ) +
  facet_wrap(~ Metric , ncol = 2,  scales = "free") +
  theme_minimal(base_size = 16) +
  labs(title = "E. coli and Klebsiella BSI Plasmid Inc rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
# save
#ggsave("rarefaction/combined_bsi_inc_rarefaction.png", plot = p, width = 12, height = 7, units = "in", dpi = 300)
ggsave("rarefaction/combined_bsi_inc_rarefaction_same_axes.png", plot = p, width = 12, height = 9, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Plasmid subcommunities Rarefaction ####
# check column
kleb_bsi_amrfinder_metadata$community_subcommunity

# define run params
max_samples <- length(unique(kleb_bsi_amrfinder_metadata$sample)) # 468
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
# save output table
write.csv(kleb_bsi_pling_rarefaction_combined, file = "rarefaction/kleb_bsi_pling_rarefaction_combined.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~#
# Entry point to read in data
#kleb_bsi_pling_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_pling_rarefaction_combined.csv")
#~~~~~~~~~~~~~~~~~~~#
# Plot
#reorder metircs:
kleb_bsi_pling_rarefaction_combined$Metric <- factor(kleb_bsi_pling_rarefaction_combined$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))
kleb_bsi_pling_rarefaction_combined$Replacement <- factor(kleb_bsi_pling_rarefaction_combined$Replacement, levels = c("With replacement", "No replacement"))

# define custom colours
pal <- c("grey90", "darkorange" )

# plot 2 graphs
p <- ggplot(data = kleb_bsi_pling_rarefaction_combined, aes(x = sample_size, y = est, fill = Replacement, colour = Replacement)) +
  geom_line( linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5, colour = NA) +
  # scale viridis colours
  scale_colour_manual(name = "Replacement", values = pal) +
  scale_fill_manual(name = "Replacement", values = pal) +
  facet_wrap(~Metric, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "E.coli BSI pling rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
#save
ggsave("rarefaction/kleb_bsi_pling_rarefaction_combined.png", plot = p, width = 14, height = 5, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Combined Kleb and E.coli plot ##### 
# PLING community - subcommunities
# load data
#ecoli_bsi_pling_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_pling_rarefaction_combined.csv")
#kleb_bsi_pling_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_pling_rarefaction_combined.csv")

combined_bsi_pling_rarefaction <- rbind(ecoli_bsi_pling_rarefaction_combined, kleb_bsi_pling_rarefaction_combined)
combined_bsi_pling_rarefaction$Replacement <- factor(combined_bsi_pling_rarefaction$Replacement, levels = c("With replacement", "No replacement"))
combined_bsi_pling_rarefaction <- combined_bsi_pling_rarefaction |>
  mutate(ReplacementGenus = interaction(Replacement, Genus, sep = "_"))

pal <- c(
  "No replacement_Escherichia"   = "seagreen3",
  "No replacement_Klebsiella"   = "darkorange",
  "With replacement_Escherichia"= "grey90",
  "With replacement_Klebsiella" = "grey90"
)

# plot combined data
p <- ggplot(data = combined_bsi_pling_rarefaction, aes(x = sample_size, y = est, fill = ReplacementGenus, colour = ReplacementGenus)) +
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
    )
  ) +
  scale_fill_manual(
    name = "Sampling scheme",
    values = pal,
    labels = c(
      "No replacement_Escherichia"    = "No replacement – Escherichia",
      "No replacement_Klebsiella"    = "No replacement – Klebsiella",
      "With replacement_Escherichia" = "With replacement – Escherichia",
      "With replacement_Klebsiella"  = "With replacement – Klebsiella"
    )
  ) +
  facet_wrap(~ Metric , ncol = 2,  scales = "free") +
  theme_minimal(base_size = 16) +
  labs(title = "E. coli and Klebsiella BSI Plasmid PLING subcommunity rarefaction",
       x = "Sample Size",
       y = "Diversity Metric") +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14))
print(p)
# save
#ggsave("rarefaction/combined_bsi_pling_rarefaction.png", plot = p, width = 12, height = 7, units = "in", dpi = 300)
ggsave("rarefaction/combined_bsi_pling_rarefaction_same_axes.png", plot = p, width = 12, height = 9, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Combined plot with E.coli and Klebsiella ####
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
combined_bsi_rarefaction$Feature <- factor(combined_bsi_rarefaction$Feature, levels = c("MLST", "AMR genes", "PLING plasmid subcommunity"))
combined_bsi_rarefaction$Replacement <- factor(combined_bsi_rarefaction$Replacement, levels = c("With replacement", "No replacement"))
combined_bsi_rarefaction$Metric <- factor(combined_bsi_rarefaction$Metric, levels = c("Simpson diversity", "Shannon diversity", "Species richness", "Sample coverage"))

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
install.packages("ggh4x")

library(ggh4x)

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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Simulate species incursion ####
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
                                       feature = "mlst_profile", sample_col = "isolateid",
                                       sample_sizes = NULL, 
                                       max_samples = max_samples,
                                       k_interval = 10,
                                       reps = n_reps, new_species_label = "NEW_SPECIES",
                                       n_cores = n_cores, seed = NULL, verbose = TRUE)


# run for Kleb
max_samples <- length(unique(kleb_bsi_samples_metadata$isolateid))
kleb_bsi_mlst_incursion <- simulate_species_incursion(kleb_bsi_samples_metadata,
                                       feature = "mlst_profile", sample_col = "isolateid",
                                       sample_sizes = NULL, 
                                       max_samples = max_samples,
                                       k_interval = 10,
                                       reps = n_reps, new_species_label = "NEW_SPECIES",
                                       n_cores = n_cores, seed = NULL, verbose = TRUE)


ecoli_bsi_mlst_incursion <- ecoli_bsi_mlst_incursion |>
  mutate(Genus = "Escherichia")
kleb_bsi_mlst_incursion <- kleb_bsi_mlst_incursion |>
  mutate(Genus = "Klebsiella")

cobmined_bsi_mlst_incursion <- rbind(ecoli_bsi_mlst_incursion, kleb_bsi_mlst_incursion)

#set stratifying features as factors:
cobmined_bsi_mlst_incursion$Metric <- factor(cobmined_bsi_mlst_incursion$Metric, levels = c("Simpson_diversity", "Shannon_diversity", "Species_richness", "Sample_coverage"))
cobmined_bsi_mlst_incursion$sample_size <- as.factor(cobmined_bsi_mlst_incursion$sample_size)

# --- PREPROCESS: create ordered sample_size factor and per-genus colour palettes ---
# assume df is your data and genus_pal_list is defined as before
df <- cobmined_bsi_mlst_incursion
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

# --- PLOT ---
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

str(cobmined_bsi_mlst_incursion$sample_size)
unique(cobmined_bsi_mlst_incursion$sample_size)
sample_levels
cobmined_bsi_mlst_incursion$sample_size <- factor(cobmined_bsi_mlst_incursion$sample_size, levels  = rev(sample_levels))

# plot
p <- ggplot(data = cobmined_bsi_mlst_incursion, aes(x = new_species_frequency, y = est, colour = sample_size, fill = sample_size)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3, colour = NA) +
  scale_fill_viridis_d(option = "mako") +
  scale_colour_viridis_d(option = "mako") +
  scale_x_log10() +
  facet_wrap(Genus ~ Metric, scales= "free", ncol = 4) +
  theme_minimal(base_size = 16)

print(p)
# save
ggsave("rarefaction/combined_bsi_mlst_incursion_logscale.png", plot = p, width = 15, height = 10, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Combined plots with both ecoli and kleb together ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# add species of origin and metric column
ecoli_bsi_mlst_rarefaction_no_replacement <- ecoli_bsi_mlst_rarefaction_no_replacement |>
  mutate(species_group = "Escherichia",
         metric = "MLST")
ecoli_gene_rarefaction_AMR_only_no_replacement <- ecoli_gene_rarefaction_AMR_only_no_replacement |>
  mutate(species_group = "Escherichia",
         metric = "ARG")
ecoli_bsi_plasmid_inc_rarefaction_no_replacement <- ecoli_bsi_plasmid_inc_rarefaction_no_replacement |>
  mutate(species_group = "Escherichia",
         metric = "Plasmid Inc types")
ecoli_bsi_plasmid_subcommunity_rarefaction_no_replacement <- ecoli_bsi_plasmid_subcommunity_rarefaction_no_replacement |>
  mutate(species_group = "Escherichia",
         metric = "Plasmid subcommunities") |>
  select(-c(X))

colnames(ecoli_bsi_mlst_rarefaction_no_replacement)
colnames(ecoli_gene_rarefaction_AMR_only_no_replacement)
colnames(ecoli_bsi_plasmid_inc_rarefaction_no_replacement)
colnames(ecoli_bsi_plasmid_subcommunity_rarefaction_no_replacement)


# add species of origin and metric column
kleb_mlst_rarefaction_no_replacement <- kleb_mlst_rarefaction_no_replacement |>
  mutate(species_group = "Klebsiella",
         metric = "MLST")
kleb_gene_rarefaction_AMR_only_no_replacement <- kleb_gene_rarefaction_AMR_only_no_replacement |>
  mutate(species_group = "Klebsiella",
         metric = "ARG")
kleb_bsi_plasmid_inc_rarefaction_no_replacement <- kleb_bsi_plasmid_inc_rarefaction_no_replacement |>
  mutate(species_group = "Klebsiella",
         metric = "Plasmid Inc types")
kleb_bsi_plasmid_subcommunity_rarefaction_no_replacement <- kleb_bsi_plasmid_subcommunity_rarefaction_no_replacement |>
  mutate(species_group = "Klebsiella",
         metric = "Plasmid subcommunities")

# rbind all together
combined_overall_bsi_rarefaction_no_replacement <- rbind(ecoli_bsi_mlst_rarefaction_no_replacement,
                                                         ecoli_gene_rarefaction_AMR_only_no_replacement,
                                                         ecoli_bsi_plasmid_inc_rarefaction_no_replacement, 
                                                         ecoli_bsi_plasmid_subcommunity_rarefaction_no_replacement, 
                                                        
                                                         kleb_mlst_rarefaction_no_replacement, 
                                                         kleb_gene_rarefaction_AMR_only_no_replacement, 
                                                         kleb_bsi_plasmid_inc_rarefaction_no_replacement, 
                                                         kleb_bsi_plasmid_subcommunity_rarefaction_no_replacement
                                                         )
View(combined_overall_bsi_rarefaction_no_replacement)
# set order of metrics
combined_overall_bsi_rarefaction_no_replacement$metric <- factor(combined_overall_bsi_rarefaction_no_replacement$metric,
                                                                 levels = c("MLST",
                                                                            "ARG", 
                                                                            "Plasmid Inc types", 
                                                                            "Plasmid subcommunities"))

#save
write.csv(combined_overall_bsi_rarefaction_no_replacement, "rarefaction/combined_overall_bsi_rarefaction_no_replacement.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~#
# * * read in combined cleaned dataset ####
#combined_overall_bsi_rarefaction_no_replacement <- read.csv("rarefaction/combined_overall_bsi_rarefaction_no_replacement.csv")
#~~~~~~~~~~~~~~~#


# plot overall combined plot
combined_overall_bsi_rarefaction_no_replacement_plot <- ggplot(combined_overall_bsi_rarefaction_no_replacement, aes(x = sample_size, y = mean_unique, color = species_group)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = species_group), alpha = 0.2, color = NA) +
  scale_color_manual(values = species_group_colours) +
  scale_fill_manual(values = species_group_colours) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  labs(x = "Sample size",
       y = "Number of unique features",
       fill = "Genus",
       colour = "Genus") +
  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
combined_overall_bsi_rarefaction_no_replacement_plot
ggsave(filename = "rarefaction/combined_overall_bsi_rarefaction_no_replacement_plot_count.png", plot = combined_overall_bsi_rarefaction_no_replacement_plot, width = 8, height = 6, units = "in", dpi = 300 )

# plot proportion diversity
combined_overall_bsi_rarefaction_no_replacement <- combined_overall_bsi_rarefaction_no_replacement |>
  group_by(species_group) |>
  mutate(max_species_sample_size = max(sample_size)) |>
  ungroup() |>
  mutate(sample_size_pct = sample_size / max_species_sample_size * 100) |>
  rename(prop_diversity_lower = pro_diversity_lower)

colnames(combined_overall_bsi_rarefaction_no_replacement)

# plot by pct
combined_overall_bsi_rarefaction_no_replacement_plot_pct <- ggplot(combined_overall_bsi_rarefaction_no_replacement, 
                                                                   aes(x = sample_size_pct, y = prop_diversity, color = species_group)) +
  geom_line() +
  geom_ribbon(aes(ymin = prop_diversity_lower, ymax = prop_diversity_upper, fill = species_group), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "black") +
  scale_color_manual(values = species_group_colours) +
  scale_fill_manual(values = species_group_colours) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  labs(x = "Sample size",
       y = "Number of unique features",
       fill = "Genus",
       colour = "Genus") +
  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
combined_overall_bsi_rarefaction_no_replacement_plot_pct
ggsave(filename = "rarefaction/combined_overall_bsi_rarefaction_no_replacement_plot_pct.png", plot = combined_overall_bsi_rarefaction_no_replacement_plot_pct, width = 8, height = 6, units = "in", dpi = 300 )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# UPSAMPLING ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ecoli BSI upsampling with preseq
library(dplyr)
library(ggplot2)
#install.packages("PreseqR") # this doesn't work as removed from CRAN :(
install.packages("fs")
find.package("cli")
packageVersion("cli")
remove.packages("cli")
install.packages("cli")

install.packages("devtools")
devtools::install_github("smithlabcode/PreseqR")
library(preseqR)


# attempt 2
# 1. Make sure you have a clean R session
.rs.restartR()   # RStudio; or restart R manually

# 2. Install the missing package(s)
install.packages("polynom")

# 3. Reinstall preseqR from GitHub (ask to install dependencies too)
# Use remotes (lighter) or devtools; remotes is recommended for single installs
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("smithlabcode/PreseqR", dependencies = TRUE)

library(preseqR)


# Count profile frequencies
profile_counts <- ecoli_bsi_samples_metadata |>
  mutate(mlst_profile = as.character(mlst_profile)) |>
  group_by(mlst_profile) |>
  summarise(n = n()) |>
  group_by(n) |>
  summarise(no_mlsts = n_distinct(mlst_profile))
#View(profile_counts)
profile_counts <- as.matrix(profile_counts)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Species accumulation curve
# Question: how many unique MLSTs occur at least r times? (as a function of sample size)

# define current and estiamted upsampling sample sizes for E.coli BSIs:
current_n <- nrow(ecoli_bsi_samples_metadata) #1471
estimated_total_n <- current_n / 0.65 #2263
estimated_national_n <- 18732

#unique mlsts observed currently
actual_unique_mlsts <- length(unique(ecoli_bsi_samples_metadata$mlst_profile)) # 287

# define sample sizes
relative_sample_sizes <- c(seq(0,15, 0.1), (estimated_total_n/current_n), (estimated_national_n / current_n))
relative_sample_sizes <- sort(relative_sample_sizes)
#length(relative_sample_sizes)
actual_sample_sizes <- relative_sample_sizes * current_n
actual_sample_sizes <- sort(actual_sample_sizes)
# define parameters, including range of r (in the range 1 - 10)minimum number of times it is important to detect an MLST).
r_values <- c(seq(1,10,1))
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
View(all_results_df)
# save
write.csv(all_results_df, "rarefaction/ecoli_bsi_upsampling_by_r.csv", row.names = FALSE)
#~~~~~~~~~~~#
# Entry point for plotting
#all_results_df <- read.csv("rarefaction/ecoli_bsi_upsampling_by_r.csv")
#~~~~~~~~~~~#

# * * plot ####
# plot all r values on 1 plot
#pal_colors <- viridis::viridis(n = length(r_values), option = "plasma")
#"#0D0887FF" "#9C179EFF" "#ED7953FF" "#F0F921FF"

#plot all r curves on same plot
p <- ggplot(data = all_results_df[all_results_df$r %in% c(1,2,3,4,5,10),], aes(x = actual_sample_size, y = est, colour = r, fill = r)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  geom_vline(xintercept = current_n, linetype = "dashed", colour = "black") + # actual sample size 
  geom_vline(xintercept = estimated_total_n, linetype = "dashed", colour = "black") + # estimates total sample size if assuming 65% sampling rate
  geom_vline(xintercept = estimated_national_n, linetype = "dashed", colour = "black") + # estimated total number of E.coli BSIs in England in this time frame (6 months)
  scale_colour_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma") +
  scale_fill_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma", guide = "none") +
  labs(
    x = "Actual Sample size",
    y = "Estimated unique MLSTs",
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

p
ggsave(filename = "rarefaction/ecoli_bsi_upsampling_by_r.png", plot = p, width = 7, height = 6, units = "in", dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
 # * *  make summary table of proportion of unique MLSTs captured at each r value and sample size ####

#unique mlsts observed currently
actual_unique_mlsts_r1 <- length(unique(ecoli_bsi_samples_metadata$mlst_profile)) # 287
actual_unique_mlsts_r2 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 2, ]$mlst_profile)) # 96
actual_unique_mlsts_r3 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 3, ]$mlst_profile)) # 54
actual_unique_mlsts_r4 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 4, ]$mlst_profile)) # 44
actual_unique_mlsts_r5 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 5, ]$mlst_profile)) # 49
actual_unique_mlsts_r6 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 6, ]$mlst_profile)) # 33
actual_unique_mlsts_r7 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 7, ]$mlst_profile)) # 28
actual_unique_mlsts_r8 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 8, ]$mlst_profile)) # 26
actual_unique_mlsts_r9 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 9, ]$mlst_profile)) # 25
actual_unique_mlsts_r10 <- length(unique(ecoli_bsi_samples_metadata[ecoli_bsi_samples_metadata$ST_freq >= 10, ]$mlst_profile)) # 22

unique_mlsts_by_r <- c(actual_unique_mlsts_r1,
                         actual_unique_mlsts_r2,
                         actual_unique_mlsts_r3,
                         actual_unique_mlsts_r4,
                         actual_unique_mlsts_r5,
                         actual_unique_mlsts_r6,
                         actual_unique_mlsts_r7,
                         actual_unique_mlsts_r8,
                         actual_unique_mlsts_r9,
                         actual_unique_mlsts_r10)

# make df of r values and unique mlsts observed at each r (minimum number of occurences of an MLST in dataset)
# NOTE that estimates of unique MLSTs at the actual sample size are acceptably close to the actual only in the range r = 0-7
# for r=8-10, the estiamtes of unique MLSTs occurring >8 times is smaller than the observed. 
# this is probably due to the frequency distribution not being estimated correctly for higher frequency MLSTs??
r_unique_mlsts_df <- as.data.frame(cbind(r_values, unique_mlsts_by_r))
r_unique_mlsts_df <- r_unique_mlsts_df |>
  mutate(r_values = factor(r_values, levels = sort(unique(r_values))))


ecoli_bsi_sac_summary <- all_results_df |>
  filter(actual_sample_size %in% c(current_n, estimated_total_n, estimated_national_n)) |>
  dplyr::select(-c(relative_sample_size))|>
  pivot_wider(names_from = actual_sample_size, values_from = c(est, se, lower, upper))|>
  left_join(r_unique_mlsts_df, by = c("r"="r_values")) |>
  # divide each estimate value by actual number of STs (occuring at least r times) oberved in sample 
  mutate(prop_regional_est = round(unique_mlsts_by_r/est_2263.07692307692,3) ,
         prop_regional_lower = round(unique_mlsts_by_r/lower_2263.07692307692,3),
         prop_regional_upper = round(unique_mlsts_by_r/upper_2263.07692307692,3) ,
         prop_national_est = round(unique_mlsts_by_r/est_18732,3) ,
         prop_national_lower = round(unique_mlsts_by_r/lower_18732,3) ,
         prop_national_upper = round(unique_mlsts_by_r/upper_18732,3),
         est_regional = round(est_2263.07692307692, 0),
         lower_regional = round(lower_2263.07692307692, 0),
         upper_regional = round(upper_2263.07692307692, 0),
         est_national = round(est_18732, 0),
         lower_national = round(lower_18732, 0),
         upper_national = round(upper_18732, 0)
         ) |>
  #make neat summaries in the form est (lower CI - upper CI)
  mutate(regional_estimate = paste0(est_regional, " (", lower_regional, " - ", upper_regional, ")"),
         national_estimate = paste0(est_national, " (", lower_national, "-", upper_national, ")"),
         regional_proportion = paste0(prop_regional_est, " (", prop_regional_upper, " - ", prop_regional_lower, ")"),
         national_proportion = paste0(prop_national_est, " (", prop_national_upper, "-", prop_national_lower, ")")) |>
  dplyr::select(c(r, unique_mlsts_by_r,
                  regional_estimate, national_estimate, regional_proportion, national_proportion,
                  est_regional, lower_regional, upper_regional,
                  est_national, lower_national, upper_national,
                  prop_regional_est, prop_regional_lower, prop_regional_upper, 
                  prop_national_est, prop_national_lower, prop_national_upper)
         )
#View(ecoli_bsi_sac_summary)

# save csv
write.csv(ecoli_bsi_sac_summary, "rarefaction/ecoli_bsi_upsampling_summary_table_by_r.csv", row.names = FALSE)


#~~~~~~~~~~~#
# * * plot ####
# now plot with y axis as a proportion of actual observed MLSTs
all_results_df_proportions <- all_results_df |>
  left_join(r_unique_mlsts_df, by = c("r"="r_values")) |>
  # add propotiyon column by dividing estimates by actual observed MLST number at each r
  mutate(prop_est = est/unique_mlsts_by_r,
         prop_lower = lower/unique_mlsts_by_r,
         prop_upper = upper/unique_mlsts_by_r,
         est_observed_prop = unique_mlsts_by_r / est ,
         lower_observed_prop = unique_mlsts_by_r / lower,
         upper_observed_prop = unique_mlsts_by_r/ upper
         ) 
#View(all_results_df_proportions)


p <- ggplot(all_results_df_proportions[all_results_df_proportions$r %in% c(1,2,3,4,5,10),], aes(x = actual_sample_size, y = prop_est, colour = r, fill = r)) +
  geom_ribbon(aes(ymin = prop_lower, ymax = prop_upper), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  geom_vline(xintercept = current_n, linetype = "dashed", colour = "black") + # actual sample size 
  geom_vline(xintercept = estimated_total_n, linetype = "dashed", colour = "black") + # estimates total sample size if assuming 65% sampling rate
  geom_vline(xintercept = estimated_national_n, linetype = "dashed", colour = "black") + # estimated total number of E.coli BSIs in England in this time frame (6 months)
  scale_colour_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma") +
  scale_fill_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma", guide = "none") +
  labs(
    x = "Total population size of E.coli BSIs",
    y = "Unique MLSTs as proportion of observed in NEKSUS sample",
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

p
ggsave(filename = "rarefaction/ecoli_bsi_upsampling_by_r_mlst_proportion.png", plot = p, width = 7, height = 6, units = "in", dpi = 300 )


# * * plot ####
# plot the proportion we have observed assuming various sampling fractions
# filter df to remove too small sample sizes
all_results_df_proportions_filtered <- all_results_df_proportions |>
  dplyr::filter(relative_sample_size >= 1)

p <- ggplot(all_results_df_proportions_filtered[all_results_df_proportions_filtered$r %in% c(1,2,3,4,5,10),], aes(x = actual_sample_size, y = est_observed_prop, colour = r, fill = r)) +
  geom_ribbon(aes(ymin = lower_observed_prop, ymax = upper_observed_prop), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  geom_vline(xintercept = current_n, linetype = "dashed", colour = "black") + # actual sample size 
  geom_vline(xintercept = estimated_total_n, linetype = "dashed", colour = "black") + # estimates total sample size if assuming 65% sampling rate
  geom_vline(xintercept = estimated_national_n, linetype = "dashed", colour = "black") + # estimated total number of E.coli BSIs in England in this time frame (6 months)
  scale_colour_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma") +
  scale_fill_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma", guide = "none") +
  scale_y_log10(breaks = c(seq(0.1, 0.9, by = 0.1), 1:10),
                labels = c(seq(0.1, 0.9, by = 0.1), 1:10)) +
  labs(
    x = "Total population size of E.coli BSIs",
    y = "Proportion of MLSTs observed in NEKSUS sample",
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

p
ggsave(filename = "rarefaction/ecoli_bsi_upsampling_by_r_mlst_proportion_observed.png", plot = p, width = 7, height = 6, units = "in", dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#SAMPLE COVERAGE
# Question: what is the probability of observing an MLST at least r times? (as a function of sample size)
# aka: what proportion of all indiviudals in the community belong to a species (mlst) we have already observed 
# repeat for a range of r values
r_values <- c(seq(1,10,1))

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
# factor r for plotting
all_results_sample_cov_df <- all_results_sample_cov_df |>
  mutate(r = factor(r, levels = sort(unique(r))))
#View(all_results_sample_cov_df)

# save
write.csv(all_results_sample_cov_df, "rarefaction/ecoli_bsi_upsampling_sample_cov_by_r.csv", row.names = FALSE)
#~~~~~~~~~~~#
# Entry point for plotting
#all_results_sample_cov_df <- read.csv("rarefaction/ecoli_bsi_upsampling_sample_cov_by_r.csv")
#~~~~~~~~~~~#

# * * plot overall ####
p <- ggplot(all_results_sample_cov_df[all_results_sample_cov_df$r %in% c(1,2,3,4,5,10),], aes(x = actual_sample_size, y = est, colour = r, fill = r)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_vline(xintercept = 1471, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 2263, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 18732, linetype = "dashed", colour = "black") +
  scale_colour_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma") +
  scale_fill_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma", guide = "none") +
  labs(
    x = "Actual Sample size",
    y = "Probability of observing an MLST at least r times in a sample ",
    # title = "Probability of observing MLSTs at least r times"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p
ggsave(filename = "rarefaction/ecoli_bsi_probability_upsampling_by_r.png", plot = p, width = 7, height = 5, units = "in", dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Summary tables ####

ecoli_bsi_sample_coverage_summary <- all_results_sample_cov_df |>
  filter(actual_sample_size %in% c(current_n, estimated_total_n, estimated_national_n)) |>
  dplyr::select(-c(relative_sample_size))|>
  pivot_wider(names_from = actual_sample_size, values_from = c(est, se, lower, upper))|>
  left_join(r_unique_mlsts_df, by = c("r"="r_values")) |>
  # divide each estimate value by actual number of STs (occuring at least r times) oberved in sample 
  mutate(
         est_actual = round(est_1471, 3),
         lower_actual = round(lower_1471, 3),
         upper_actual = round(upper_1471, 3),
         est_regional = round(est_2263.07692307692, 3),
         lower_regional = round(lower_2263.07692307692, 3),
         upper_regional = round(upper_2263.07692307692, 3),
         est_national = round(est_18732, 3),
         lower_national = round(lower_18732, 3),
         upper_national = round(upper_18732, 3)
           ) |>
  #make neat summaries in the form est (lower CI - upper CI)
  mutate(
    actual_estimate = paste0(est_actual, " (", lower_actual, " - ", upper_actual, ")"),
    regional_estimate = paste0(est_regional, " (", lower_regional, " - ", upper_regional, ")"),
    national_estimate = paste0(est_national, " (", lower_national, "-", upper_national, ")")
        ) |>
  dplyr::select(c(r, unique_mlsts_by_r,
                  actual_estimate, regional_estimate, national_estimate,
                  est_regional, lower_regional, upper_regional,
                  est_national, lower_national, upper_national)
  )
#View(ecoli_bsi_sample_coverage_summary)

# save csv
write.csv(ecoli_bsi_sample_coverage_summary, "rarefaction/ecoli_bsi_upsampling_summary_table_by_r.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# try optimal sampling function
# Question: what is the optimal sampling fraction?? (developped for sc WGS from a shallow seq experiment)
#ecoli_bsi_optimal_seq <- preseqR.optimal.sequencing(n = profile_counts, efficiency=0.1, bin=1, r=1, mt=20, times=10, conf=0.95)
# not that useful as cost-benefit trade-off (efficiency) is arbitrary.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# KLEBSIELLAE BSI UP-SAMPLING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# count mlst profile frequencies
profile_counts <- kleb_bsi_samples_metadata |>
  mutate(mlst_profile = as.character(mlst_profile)) |>
  group_by(mlst_profile) |>
  summarise(n = n()) |>
  group_by(n) |>
  dplyr::summarise(no_mlsts = n_distinct(mlst_profile))
#View(profile_counts)
profile_counts <- as.matrix(profile_counts)


# Run species accumulation curve
# define current and estiamted upsampling sample sizes for Kleb BSIs:
current_n <- nrow(kleb_bsi_samples_metadata) #468
estimated_total_n <- current_n / 0.65 #720
estimated_national_n <- 4933

# Species accumulation curve
# Question: how many unique MLSTs occur at least r times? (as a function of sample size)

#unique mlsts observed currently
actual_unique_mlsts <- length(unique(kleb_bsi_samples_metadata$mlst_profile)) # 297

# define sample sizes
relative_sample_sizes <- c(seq(0,30, 0.1), (estimated_total_n/current_n), (estimated_national_n / current_n))
relative_sample_sizes <- sort(relative_sample_sizes)
length(relative_sample_sizes)
actual_sample_sizes <- relative_sample_sizes * current_n
actual_sample_sizes <- sort(actual_sample_sizes)
# define parameters, including range of r (in the range 1 - 10)minimum number of times it is important to detect an MLST).
r_values <- c(seq(1,10,1))
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
View(all_results_df)
# save
write.csv(all_results_df, "rarefaction/kleb_bsi_upsampling_by_r.csv", row.names = FALSE)
#~~~~~~~~~~~#
# Entry point for plotting
#all_results_df <- read.csv("rarefaction/kleb_bsi_upsampling_by_r.csv")
#~~~~~~~~~~~#

# * * plot ####
#plot all r curves on same plot
p <- ggplot(all_results_df[all_results_df$r %in% c(1,2,3,4,5,10),], aes(x = actual_sample_size, y = est, colour = r, fill = r)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  geom_vline(xintercept = current_n, linetype = "dashed", colour = "black") + # actual sample size 
  geom_vline(xintercept = estimated_total_n, linetype = "dashed", colour = "black") + # estimates total sample size if assuming 65% sampling rate
  geom_vline(xintercept = estimated_national_n, linetype = "dashed", colour = "black") + # estimated total number of E.coli BSIs in England in this time frame (6 months)
  scale_colour_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma") +
  scale_fill_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma", guide = "none") +
  labs(
    x = "Actual Sample size",
    y = "Estimated unique MLSTs",
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

p
ggsave(filename = "rarefaction/kleb_bsi_upsampling_by_r.png", plot = p, width = 7, height = 6, units = "in", dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * *  make summary table of proportion of unique MLSTs captured at each r value and sample size ####

#unique mlsts observed currently
actual_unique_mlsts_r1 <- length(unique(kleb_bsi_samples_metadata$mlst_profile)) # 297
actual_unique_mlsts_r2 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 2, ]$mlst_profile)) # 70
actual_unique_mlsts_r3 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 3, ]$mlst_profile)) # 29
actual_unique_mlsts_r4 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 4, ]$mlst_profile)) # 18
actual_unique_mlsts_r5 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 5, ]$mlst_profile)) # 14
actual_unique_mlsts_r6 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 6, ]$mlst_profile)) # 11
actual_unique_mlsts_r7 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 7, ]$mlst_profile)) # 9
actual_unique_mlsts_r8 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 8, ]$mlst_profile)) # 7
actual_unique_mlsts_r9 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 9, ]$mlst_profile)) # 3
actual_unique_mlsts_r10 <- length(unique(kleb_bsi_samples_metadata[kleb_bsi_samples_metadata$ST_freq >= 10, ]$mlst_profile)) #3

unique_mlsts_by_r <- c(actual_unique_mlsts_r1,
                       actual_unique_mlsts_r2,
                       actual_unique_mlsts_r3,
                       actual_unique_mlsts_r4,
                       actual_unique_mlsts_r5,
                       actual_unique_mlsts_r6,
                       actual_unique_mlsts_r7,
                       actual_unique_mlsts_r8,
                       actual_unique_mlsts_r9,
                       actual_unique_mlsts_r10)

# make df of r values and unique mlsts observed at each r (minimum number of occurences of an MLST in dataset)
# NOTE that estimates of unique MLSTs at the actual sample size are acceptably close to the actual only in the range r = 0-7
# for r=8-10, the estiamtes of unique MLSTs occurring >8 times is smaller than the observed. 
# this is probably due to the frequency distribution not being estimated correctly for higher frequency MLSTs??
r_unique_mlsts_df <- as.data.frame(cbind(r_values, unique_mlsts_by_r))
r_unique_mlsts_df <- r_unique_mlsts_df |>
  mutate(r_values = factor(r_values, levels = sort(unique(r_values))))


kleb_bsi_sac_summary <- all_results_df |>
  filter(actual_sample_size %in% c(current_n, estimated_total_n, estimated_national_n)) |>
  dplyr::select(-c(relative_sample_size))|>
  pivot_wider(names_from = actual_sample_size, values_from = c(est, se, lower, upper))|>
  left_join(r_unique_mlsts_df, by = c("r"="r_values"))|>
  # divide each estimate value by actual number of STs (occuring at least r times) oberved in sample 
  mutate(prop_regional_est = round(unique_mlsts_by_r/est_720,3) ,
         prop_regional_lower = round(unique_mlsts_by_r/lower_720,3),
         prop_regional_upper = round(unique_mlsts_by_r/upper_720,3) ,
         prop_national_est = round(unique_mlsts_by_r/est_4933,3) ,
         prop_national_lower = round(unique_mlsts_by_r/lower_4933,3) ,
         prop_national_upper = round(unique_mlsts_by_r/upper_4933,3),
         est_regional = round(est_720, 0),
         lower_regional = round(lower_720, 0),
         upper_regional = round(upper_720, 0),
         est_national = round(est_4933, 0),
         lower_national = round(lower_4933, 0),
         upper_national = round(upper_4933, 0)
  ) |>
  #make neat summaries in the form est (lower CI - upper CI)
  mutate(regional_estimate = paste0(est_regional, " (", lower_regional, " - ", upper_regional, ")"),
         national_estimate = paste0(est_national, " (", lower_national, "-", upper_national, ")"),
         regional_proportion = paste0(prop_regional_est, " (", prop_regional_upper, " - ", prop_regional_lower, ")"),
         national_proportion = paste0(prop_national_est, " (", prop_national_upper, "-", prop_national_lower, ")")) |>
  dplyr::select(c(r, unique_mlsts_by_r,
                  regional_estimate, national_estimate, regional_proportion, national_proportion,
                  est_regional, lower_regional, upper_regional,
                  est_national, lower_national, upper_national,
                  prop_regional_est, prop_regional_lower, prop_regional_upper, 
                  prop_national_est, prop_national_lower, prop_national_upper)
  )
#View(kleb_bsi_sac_summary)

# save csv
write.csv(kleb_bsi_sac_summary, "rarefaction/kleb_bsi_upsampling_summary_table_by_r.csv", row.names = FALSE)

#~~~~~~~~~~~#
# * * plot ####
# now plot with y axis as a proportion of actual observed MLSTs
all_results_df_proportions <- all_results_df |>
  left_join(r_unique_mlsts_df, by = c("r"="r_values")) |>
  # add propotiyon column by dividing estimates by actual observed MLST number at each r
  mutate(prop_est = est/unique_mlsts_by_r,
         prop_lower = lower/unique_mlsts_by_r,
         prop_upper = upper/unique_mlsts_by_r,
         est_observed_prop = unique_mlsts_by_r / est ,
         lower_observed_prop = unique_mlsts_by_r / lower,
         upper_observed_prop = unique_mlsts_by_r/ upper
  ) 
#View(all_results_df_proportions)


p <- ggplot(all_results_df_proportions[all_results_df_proportions$r %in% c(1,2,3,4,5,10),], aes(x = actual_sample_size, y = prop_est, colour = r, fill = r)) +
  geom_ribbon(aes(ymin = prop_lower, ymax = prop_upper), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  geom_vline(xintercept = current_n, linetype = "dashed", colour = "black") + # actual sample size 
  geom_vline(xintercept = estimated_total_n, linetype = "dashed", colour = "black") + # estimates total sample size if assuming 65% sampling rate
  geom_vline(xintercept = estimated_national_n, linetype = "dashed", colour = "black") + # estimated total number of E.coli BSIs in England in this time frame (6 months)
  scale_colour_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma") +
  scale_fill_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma", guide = "none") +
  labs(
    x = "Total population size of E.coli BSIs",
    y = "Unique MLSTs as proportion of observed in NEKSUS sample",
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

p
ggsave(filename = "rarefaction/kleb_bsi_upsampling_by_r_mlst_proportion.png", plot = p, width = 7, height = 6, units = "in", dpi = 300 )


# * * plot ####
# plot the proportion we have observed assuming various sampling fractions
# filter df to remove too small sample sizes
all_results_df_proportions_filtered <- all_results_df_proportions |>
  dplyr::filter(relative_sample_size >= 1)

p <- ggplot(all_results_df_proportions_filtered[all_results_df_proportions_filtered$r %in% c(1,2,3,4,5,10),], aes(x = actual_sample_size, y = est_observed_prop, colour = r, fill = r)) +
  geom_ribbon(aes(ymin = lower_observed_prop, ymax = upper_observed_prop), alpha = 0.3, colour = NA) +
  geom_line(size = 1) +
  geom_vline(xintercept = current_n, linetype = "dashed", colour = "black") + # actual sample size 
  geom_vline(xintercept = estimated_total_n, linetype = "dashed", colour = "black") + # estimates total sample size if assuming 65% sampling rate
  geom_vline(xintercept = estimated_national_n, linetype = "dashed", colour = "black") + # estimated total number of E.coli BSIs in England in this time frame (6 months)
  scale_colour_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma") +
  scale_fill_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma", guide = "none") +
  scale_y_log10(breaks = c(seq(0.1, 0.9, by = 0.1), 1:10),
                labels = c(seq(0.1, 0.9, by = 0.1), 1:10)) +
  labs(
    x = "Total population size of E.coli BSIs",
    y = "Proportion of MLSTs observed in NEKSUS sample",
    title = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

p
ggsave(filename = "rarefaction/kleb_bsi_upsampling_by_r_mlst_proportion_observed.png", plot = p, width = 7, height = 6, units = "in", dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#SAMPLE COVERAGE
# Question: what is the probability of observing an MLST at least r times? (as a function of sample size)
# repeat for a range of r values
r_values <- c(1,2,3,4,5,6,7,8,9)

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

# do manually for r = 10 as getting error message 'system is computationally singular: reciprocal condition number = 1.32704e-18'
sac <- preseqR.sample.cov.bootstrap(n = profile_counts, r = 10, mt = 20,
                                    times = bootstrap_times, conf = conf_level)

results_r10 <- tibble(r = 10, 
       relative_sample_size = relative_sample_sizes ,
       actual_sample_size = actual_sample_sizes ,
       est = sac$f(relative_sample_sizes ),
       se = sac$se(relative_sample_sizes ),
       lower = sac$lb(relative_sample_sizes ),
       upper = sac$ub(relative_sample_sizes )
)

results_list_sample_cov[[10]] <- results_r10

all_results_sample_cov_df <- bind_rows(results_list_sample_cov)
# factor r for plotting
all_results_sample_cov_df <- all_results_sample_cov_df |>
  mutate(r = factor(r, levels = sort(unique(r))))
#View(all_results_sample_cov_df)


 save
write.csv(all_results_sample_cov_df, "rarefaction/kleb_bsi_upsampling_sample_cov_by_r.csv", row.names = FALSE)
#~~~~~~~~~~~#
# Entry point for plotting
#all_results_sample_cov_df <- read.csv("rarefaction/kleb_bsi_upsampling_sample_cov_by_r.csv")
#~~~~~~~~~~~#
# * * plot overall ####
p <- ggplot(all_results_sample_cov_df[all_results_sample_cov_df$r %in% c(1,2,3,4,5,10),], aes(x = actual_sample_size, y = est, colour = r, fill = r)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_vline(xintercept = current_n, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = estimated_total_n, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = estimated_national_n, linetype = "dashed", colour = "black") +
  scale_colour_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma") +
  scale_fill_viridis_d(name = "MLSTs occur\nat least r times", option = "plasma", guide = "none") +
  labs(
    x = "Actual Sample size",
    y = "Probability of observing an MLST at least r times in a sample ",
    # title = "Probability of observing MLSTs at least r times"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p
ggsave(filename = "rarefaction/kleb_bsi_probability_upsampling_by_r.png", plot = p, width = 7, height = 5, units = "in", dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Summary tables ####

kleb_bsi_sample_coverage_summary <- all_results_sample_cov_df |>
  filter(actual_sample_size %in% c(current_n, estimated_total_n, estimated_national_n)) |>
  dplyr::select(-c(relative_sample_size))|>
  pivot_wider(names_from = actual_sample_size, values_from = c(est, se, lower, upper))|>
  left_join(r_unique_mlsts_df, by = c("r"="r_values")) |>
  # divide each estimate value by actual number of STs (occuring at least r times) oberved in sample 
  mutate(
    est_actual = round(est_468, 3),
    lower_actual = round(lower_468, 3),
    upper_actual = round(upper_468, 3),
    est_regional = round(est_720, 3),
    lower_regional = round(lower_720, 3),
    upper_regional = round(upper_720, 3),
    est_national = round(est_4933, 3),
    lower_national = round(lower_4933, 3),
    upper_national = round(upper_4933, 3)
  ) |>
  #make neat summaries in the form est (lower CI - upper CI)
  mutate(
    actual_estimate = paste0(est_actual, " (", lower_actual, " - ", upper_actual, ")"),
    regional_estimate = paste0(est_regional, " (", lower_regional, " - ", upper_regional, ")"),
    national_estimate = paste0(est_national, " (", lower_national, "-", upper_national, ")")
  ) |>
  dplyr::select(c(r, unique_mlsts_by_r,
                  actual_estimate, regional_estimate, national_estimate,
                  est_regional, lower_regional, upper_regional,
                  est_national, lower_national, upper_national)
  )
#View(kleb_bsi_sample_coverage_summary)

# save csv
write.csv(kleb_bsi_sample_coverage_summary, "rarefaction/kleb_bsi_upsampling_summary_table_by_r_sample_cov.csv", row.names = FALSE)



profile_counts[,1] %*% profile_counts[,2] # 468 = total number of Klebsiella isolates 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# iNEXT ####
install.packages("iNEXT")
library(iNEXT)

# * Abundance data (when counting how may of MLSTs present) ####
# * * sample size based sampling curve ####
length(unique(ecoli_bsi_samples_metadata$mlst_profile)) # 287 unique mlsts
# prepare input data as matrix with rownames as mlst_profiles and column names as sites, with first column being overall. 
ecoli_bsi_mlst_abundance_df <- ecoli_bsi_samples_metadata |>
  group_by(mlst_profile, region) |>
  summarise(n = n()) |>
  pivot_wider(names_from = region, values_from = n, values_fill = 0) |>
  mutate(Total = sum(Midlands, `North East A`, `North East B`, `North West`, `South West`, 
                                `East`, `London`,`South East A`, `South East B`, `South East C`)) |>
  dplyr::select(mlst_profile, Total, everything())
ecoli_bsi_mlst_abundance_df <- data.frame(ecoli_bsi_mlst_abundance_df)
rownames(ecoli_bsi_mlst_abundance_df) <- ecoli_bsi_mlst_abundance_df[,1]
ecoli_bsi_mlst_abundance_df <- ecoli_bsi_mlst_abundance_df[,-1]
#View(ecoli_bsi_mlst_abundance_df)
sum(ecoli_bsi_mlst_abundance_df$Total) # check

# get vector of sample sizes for which to estimate diversity
# reuse the actual_sample_sizes integer vector from above
# define current and estiamted upsampling sample sizes for E.coli BSIs:
current_n <- nrow(ecoli_bsi_samples_metadata) #1471
estimated_total_n <- current_n / 0.65 #2263
estimated_national_n <- 18732

#unique mlsts observed currently
actual_unique_mlsts <- length(unique(ecoli_bsi_samples_metadata$mlst_profile)) # 287

# define sample sizes
relative_sample_sizes <- c(seq(0,15, 0.1), (estimated_total_n/current_n), (estimated_national_n / current_n))
relative_sample_sizes <- sort(relative_sample_sizes)
#length(relative_sample_sizes)
actual_sample_sizes <- relative_sample_sizes * current_n
actual_sample_sizes <- sort(actual_sample_sizes)
str(actual_sample_sizes)
# make integer
actual_sample_sizes <- as.integer(round(actual_sample_sizes))

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
saveRDS(ecoli_bsi_iNEXT, "rarefaction/ecoli_bsi_iNEXT.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in ####
#ecoli_bsi_iNEXT <- readRDS("rarefaction/ecoli_bsi_iNEXT.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# covert p to factor for plotting
ecoli_bsi_iNEXT_size_based <- as.data.frame(ecoli_bsi_iNEXT$iNextEst$size_based)
ecoli_bsi_iNEXT_size_based$Order.q <- factor(ecoli_bsi_iNEXT_size_based$Order.q, levels = sort(unique(ecoli_bsi_iNEXT_size_based$Order.q)))
ecoli_bsi_iNEXT_coverage_based <- as.data.frame(ecoli_bsi_iNEXT$iNextEst$coverage_based)
ecoli_bsi_iNEXT_coverage_based$Order.q <- factor(ecoli_bsi_iNEXT_coverage_based$Order.q)
ecoli_bsi_asy <- as.data.frame(ecoli_bsi_iNEXT$AsyEst)
ecoli_bsi_asy$Diversity <- factor(ecoli_bsi_asy$Diversity)
View(ecoli_bsi_iNEXT_size_based)


# split into regional and total data
ecoli_bsi_iNEXT_size_based_total <- ecoli_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total" & m != 0 & m != 1)
ecoli_bsi_iNEXT_size_based_regional <- ecoli_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage != "Total" & m != 0 & m != 1)

ecoli_bsi_iNEXT_coverage_based_total <- ecoli_bsi_iNEXT_coverage_based |>
  dplyr::filter(Assemblage == "Total" & m != 0 & m != 1 )
ecoli_bsi_iNEXT_coverage_based_regional <- ecoli_bsi_iNEXT_coverage_based |>
  dplyr::filter(Assemblage != "Total" & m != 0& m != 1)

ecoli_bsi_asy <- ecoli_bsi_asy |>
  mutate(Order.q = case_when(Diversity == "Species richness" ~ 0,
                             Diversity == "Shannon diversity" ~ 1,
                             Diversity == "Simpson diversity" ~ 2))
ecoli_bsi_asy$Order.q <- factor(ecoli_bsi_asy$Order.q, levels = sort(unique(ecoli_bsi_asy$Order.q)))
ecoli_bsi_asy_total <- ecoli_bsi_asy |>
  dplyr::filter(Assemblage == "Total")
ecoli_bsi_asy_regional <- ecoli_bsi_asy |>
  dplyr::filter(Assemblage != "Total")


# ChatGPT summary or coverage based vs size based CIs:
# “Coverage-based confidence intervals for MLST richness were substantially wider
# than size-based intervals, reflecting high uncertainty in the number of rare,
# unobserved MLSTs. This effect was pronounced for species richness (q = 0) but 
# minimal for Shannon (q = 1) and Simpson (q = 2) diversity, which are dominated 
# by common MLSTs and are therefore less sensitive to incomplete sampling.”

# ChatGPT definition of sampling coverage:
#Sampling coverage answers this question:
#  “What fraction of the total probability mass (i.e. individuals) in the community
#  belongs to species that we have already observed?”
# Equivalently:
#   “If I take one more individual at random from the population, what is the 
#   probability it belongs to a species I’ve already seen?”

# ensure Order.q is a factor with nice labels
order_labels <- c("0" = "Species richness (q=0)",
                  "1" = "Shannon diversity (q=1)",
                  "2" = "Simpson diversity (q=2)")

# compute a scale factor to map coverage (0..1) onto the qD scale for the RHS axis
# we take the maximum upper CI of qD from the size-based df as the top of primary axis
max_qD <- max(ecoli_bsi_iNEXT_size_based_total$qD.UCL, ecoli_bsi_iNEXT_size_based_total$qD, na.rm = TRUE)
# ensure max_qD > 0
if (max_qD <= 0 || is.infinite(max_qD) || is.na(max_qD)) max_qD <- max(size_df$qD, na.rm = TRUE)

# the coverage values typically are in 0..1; multiply them by this factor
cov_scale <- max_qD / 1     # equivalently = max_qD

# colour palette (discrete, colorblind-friendly)
#n_q <- length(unique(ecoli_bsi_iNEXT_size_based$Order.q))
pal <- viridis(n = 6, option = "D")
"#440154FF" "#414487FF" "#2A788EFF" "#22A884FF" "#7AD151FF" "#FDE725FF"


# plot total with size based and coverage based CIs
p <- ggplot(data = ecoli_bsi_iNEXT_size_based_total) +
  # size-based ribbons and lines (primary axis)
  geom_line( aes(x = m, y = qD, colour = Order.q, group = Order.q), size = 1) +
  geom_ribbon( aes(x = m, ymin = qD.LCL, ymax = qD.UCL, fill = Order.q, group = Order.q),
    alpha = 0.35, colour = NA) +

  # coverage-based (overlay) lighter ribbons and dashed lines
  geom_ribbon(aes(x = ecoli_bsi_iNEXT_coverage_based_total$m, 
                  ymin = ecoli_bsi_iNEXT_coverage_based_total$qD.LCL, 
                  ymax = ecoli_bsi_iNEXT_coverage_based_total$qD.UCL, fill = Order.q, 
                  group = Order.q),
    alpha = 0.18, colour = NA) +
  
  # add horizontal dashed asymptotes from AsyEst (for Total assemblage)
  geom_hline(data = ecoli_bsi_asy_total, aes(yintercept = Estimator, colour = Order.q),
    linetype = "dashed", size = 0.6) +
  geom_text(data = ecoli_bsi_asy_total,aes(x = 0, y = Estimator, label = round(Estimator), colour = Order.q),
    hjust = 1, vjust = -0.3, size = 3.5, fontface = "bold", show.legend = FALSE )+
  
  # add the scaled coverage line (RHS axis mapping) from the size-based df
  geom_line(aes(x = m, y = SC * cov_scale), colour = "black", size = 0.8, inherit.aes = FALSE) +
  geom_ribbon(aes(x = m, ymin = SC.LCL * cov_scale, ymax = SC.UCL * cov_scale),
   fill = "grey60", alpha = 0.3, inherit.aes = FALSE) +
  
  # vertical reference lines for sample sizes
  geom_vline(xintercept = 1471, linetype = "longdash", colour = "grey40", size = 1) +
  geom_vline(xintercept = 2263, linetype = "longdash", colour = "grey40" , size = 1) +
  geom_vline(xintercept = 18732, linetype = "longdash", colour = "grey40", size = 1) +
  
  # scales: primary y for qD, secondary y for coverage
  scale_colour_manual(name = "Diversity order", values = pal, labels = order_labels) +
  scale_fill_manual(name = "Diversity order", values = pal, labels = order_labels, guide = guide_legend(override.aes = list(alpha = 0.4))) +
  
  scale_y_log10( name = "Diversity (qD)",
   # limits = c(0, max_qD * 1.05),
    sec.axis = sec_axis(~ . / cov_scale, name = "Sample coverage (SC)")) +
  scale_x_continuous(name = "Sample size (m)") +
  labs(
    title = "iNEXT size- and coverage-based diversity curves",
    subtitle = "Dark ribbon= size-based CI; light= coverage-based; black/grey = sample coverage (RHS axis)",
    caption = "Horizontal dashed asymptote estimates from iNEXT$AsyEst"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5)
  )
  
p

#save
ggsave(filename = "rarefaction/ecoli_bsi_iNEXT_mlst_total_logscale.png", plot = p, width = 10, height = 7, units = "in", dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Compare Preseq and iNEXT ####
# plot just total species richness (q=0) on same axes as preseqR estimate
ecoli_bsi_iNEXT_subset <- ecoli_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total", Order.q == 0) |>
  mutate(Method = "iNEXT", r = "iNEXT") |>
  rename(actual_sample_size = m ,
         est = qD,
         lower = qD.LCL, 
         upper = qD.UCL) |>
  dplyr::select(c(Method, r, actual_sample_size, est, lower, upper))

# load preseq data
ecoli_bsi_mlst_preseq <- read.csv("rarefaction/ecoli_bsi_upsampling_by_r.csv")
#View(ecoli_bsi_mlst_preseq)
ecoli_bsi_mlst_preseq_filtered <- ecoli_bsi_mlst_preseq |>
  dplyr::filter(r == 1) |>
  mutate(Method = "preseqR") |>
  dplyr::select(Method, r, actual_sample_size, est, lower, upper)
#rbind
preseq_iNEXT_combined <- rbind(ecoli_bsi_iNEXT_subset, ecoli_bsi_mlst_preseq_filtered)
#View(preseq_iNEXT_combined)

# plot
p <- ggplot(data  = preseq_iNEXT_combined, aes(x = actual_sample_size, y = est, colour = Method, fill = Method)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_ribbon(data = ecoli_bsi_iNEXT_coverage_based_total |> dplyr::filter(Order.q == 0) , aes(x = m, y = qD, ymin = qD.LCL, ymax = qD.UCL), alpha = 0.3, colour = NA, fill = "salmon") +
  # vertical reference lines for sample sizes
  geom_vline(xintercept = 1471, linetype = "longdash", colour = "grey20", size = 1) +
  geom_vline(xintercept = 2263, linetype = "longdash", colour = "grey20" , size = 1) +
  geom_vline(xintercept = 18732, linetype = "longdash", colour = "grey20", size = 1) +
  #scale_colour_viridis_d(name = "Method", option = "B") +
  #scale_fill_viridis_d(name = "Method", option = "B", guide = "none") +
  labs(
    x = "Sample Size",
    y = "Species Richness"
  ) +
  theme_minimal()
p
#save
ggsave(filename = "rarefaction/ecoli_bsi_mlst_iNEXT_vs_preseq_count.png", plot = p, width = 7, height = 5, units = "in", dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# summary table for species richness, Shannon and Simpson diversity, at important sample sizes. 
ecoli_bsi_iNEXT_summary <- ecoli_bsi_iNEXT_size_based |>
  dplyr::filter(m %in% c(current_n, round(estimated_total_n), estimated_national_n)) |>
  mutate(Method = "iNEXT") |>
  dplyr::select(c(Method, Assemblage, Order.q, m, qD, qD.LCL, qD.UCL)) |>
  rename(actual_sample_size = m ,
         est = qD,
         lower = qD.LCL, 
         upper = qD.UCL) |>
  mutate(Diversity = case_when(Order.q == 0 ~ "Species richness (q=0)",
                           Order.q == 1 ~ "Shannon diversity (q=1)",
                           Order.q == 2 ~ "Simpson diversity (q=2)")) |>
  mutate(Estimated_MLST_count = paste0(round(est), " (", round(lower), " - ", round(upper), ")")) |>
  dplyr::select(c(Method, Assemblage, Order.q, Diversity, actual_sample_size, Estimated_MLST_count, est, lower, upper))
View(ecoli_bsi_iNEXT_summary)

#save
write.csv(ecoli_bsi_iNEXT_summary, "rarefaction/ecoli_bsi_iNEXT_mlst_count_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# repeat for sampling coverage
ecoli_bsi_iNEXT_subset_samcov <- ecoli_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total", Order.q == 0) |>
  mutate(Method = "iNEXT") |>
  rename(actual_sample_size = m ,
         est = SC,
         lower = SC.LCL, 
         upper = SC.UCL) |>
  dplyr::select(c(Method, actual_sample_size, est, lower, upper))
View(ecoli_bsi_iNEXT_subset_samcov)

# load preseq data
all_results_sample_cov_df <- read.csv("rarefaction/ecoli_bsi_upsampling_sample_cov_by_r.csv")
#View(all_results_sample_cov_df)
ecoli_bsi_mlst_preseq_samcov_filtered <- all_results_sample_cov_df |>
  dplyr::filter(r == 1) |>
  mutate(Method = "preseqR") |>
  dplyr::select(Method, actual_sample_size, est, lower, upper)
#rbind
preseq_iNEXT_combined_samcov <- rbind(ecoli_bsi_iNEXT_subset_samcov, ecoli_bsi_mlst_preseq_samcov_filtered)
#View(preseq_iNEXT_combined_samcov)

# plot
p <- ggplot(data  = preseq_iNEXT_combined_samcov, aes(x = actual_sample_size, y = est, colour = Method, fill = Method)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  # vertical reference lines for sample sizes
  geom_vline(xintercept = 1471, linetype = "longdash", colour = "grey20", size = 1) +
  geom_vline(xintercept = 2263, linetype = "longdash", colour = "grey20" , size = 1) +
  geom_vline(xintercept = 18732, linetype = "longdash", colour = "grey20", size = 1) +
  #scale_colour_viridis_d(name = "Method", option = "B") +
  #scale_fill_viridis_d(name = "Method", option = "B", guide = "none") +
  #scale_y_log10() +
  labs(
    x = "Sample Size",
    y = "Sample Coverage\nProportion of all individuals belonging to observed MLSTs"
  ) +
  theme_minimal()
p
#save
ggsave(filename = "rarefaction/ecoli_bsi_mlst_iNEXT_vs_preseq_sample_coverage.png", plot = p, width = 7, height = 5, units = "in", dpi = 300 )

# * * Summary table for sampling coverage joint for preseQ (r=1) and iNEXT ####
# At various sampling depths, what proportion of the total community will belong to species we have already sampled?

preseq_iNEXT_combined_samcov$actual_sample_size <- round(preseq_iNEXT_combined_samcov$actual_sample_size)
ecoli_bsi_preseq_iNEXT_samcov_summary <- preseq_iNEXT_combined_samcov |>
  dplyr::filter(actual_sample_size %in% c(current_n, round(estimated_total_n), estimated_national_n))|>
  mutate(Estimated_sampling_coverage = paste0(round(est*100, 1), " (", round(lower*100, 1), " - ", round(upper*100,1), ")")) |>
  dplyr::select(c(Method, actual_sample_size, Estimated_sampling_coverage)) |>
  pivot_wider(id_cols = actual_sample_size, names_from = Method, values_from = Estimated_sampling_coverage)
View(ecoli_bsi_preseq_iNEXT_samcov_summary)

#save
write.csv(ecoli_bsi_preseq_iNEXT_samcov_summary, "rarefaction/ecoli_bsi_preseq_iNEXT_samcov_summary.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~#
# * Repeat for incidence data ####
# (using this as ARG presence/ absence, and bacterium is the 'sampling unit')
# load amrfinder metadata
amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
View(amrfinder_metadata_updated)
# Filter ecoli BSI specific data
ecoli_bsi_amrfinder_metadata <- amrfinder_metadata_updated |>
  filter(grepl("Escherichia", kraken2_species) & (sampletype == "blood" | sampletype == "unknown" | is.na(sampletype))) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
length(unique(ecoli_bsi_amrfinder_metadata$sample)) #1471 non-duplicates
length(unique(ecoli_bsi_amrfinder_metadata$Element.symbol)) #332 = 331 unique ARGs and NA
# 189 different AMR genes
# consider filtering for identity and coverage, but in practice, most >95% for both


ecoli_bsi_arg_presence_absence_matrix <- ecoli_bsi_amrfinder_metadata |>
  dplyr::filter(Type == "AMR") |>
  group_by(sample, Element.symbol) |>
  summarise(n = n()) |>
  mutate(presence = ifelse(n >= 1, 1, 0)) |>
  ungroup() |>
  pivot_wider(id_cols = Element.symbol, names_from = sample, values_from = presence, values_fill = 0) |>
  dplyr::filter(!is.na(Element.symbol)) # remove NAs
#View(ecoli_bsi_arg_presence_absence_matrix) 
ecoli_bsi_arg_presence_absence_matrix <- as.data.frame(ecoli_bsi_arg_presence_absence_matrix)
rownames(ecoli_bsi_arg_presence_absence_matrix) <- ecoli_bsi_arg_presence_absence_matrix$Element.symbol
ecoli_bsi_arg_presence_absence_matrix <- ecoli_bsi_arg_presence_absence_matrix[,-1]

#transform to incidence frequency input vectror needed for iNEXT
ecoli_bsi_arg_incfreq <- as.incfreq(ecoli_bsi_arg_presence_absence_matrix)
#View(ecoli_bsi_arg_incfreq)

ecoli_bsi_arg_iNEXT <- iNEXT(ecoli_bsi_arg_incfreq,
  q = c(0,1,2),  datatype = "incidence_freq", size = actual_sample_sizes, 
 # endpoint = NULL, knots = 40,
  se = TRUE, conf = 0.95, nboot = 100)

# inspect results
ecoli_bsi_arg_iNEXT$DataInfo # 99.46% sample coverage of observed 189 ARGs
ecoli_bsi_arg_iNEXT$iNextEst$size_based
ecoli_bsi_arg_iNEXT$iNextEst$coverage_based
ecoli_bsi_arg_iNEXT$AsyEst

#save
saveRDS(ecoli_bsi_arg_iNEXT, "rarefaction/ecoli_bsi_arg_iNEXT.rds")
#~~~~~~~~~~~~~~~#
# read in saved
#ecoli_bsi_arg_iNEXT <- readRDS("ecoli_bsi_arg_iNEXT.rds")
#~~~~~~~~~~~~~~~#

# get outputs into manageable data frame form
ecoli_bsi_arg_size_based <- as.data.frame(ecoli_bsi_arg_iNEXT$iNextEst$size_based)
ecoli_bsi_arg_size_based$Order.q <- factor(ecoli_bsi_arg_size_based$Order.q, levels = sort(unique(ecoli_bsi_arg_size_based$Order.q)))
ecoli_bsi_arg_size_based <- ecoli_bsi_arg_size_based |>
  dplyr::filter(t != 0 & t != 1) # remove zero and one sample size rows

ecoli_bsi_arg_coverage_based <- as.data.frame(ecoli_bsi_arg_iNEXT$iNextEst$coverage_based)
ecoli_bsi_arg_coverage_based$Order.q <- factor(ecoli_bsi_arg_coverage_based$Order.q, levels = sort(unique(ecoli_bsi_arg_coverage_based$Order.q)))
ecoli_bsi_arg_coverage_based <- ecoli_bsi_arg_coverage_based |>
  dplyr::filter(t != 0 & t != 1) 

ecoli_bsi_arg_asy <- as.data.frame(ecoli_bsi_arg_iNEXT$AsyEst)
ecoli_bsi_arg_asy$Diversity <- rownames(ecoli_bsi_arg_asy)
ecoli_bsi_arg_asy <- ecoli_bsi_arg_asy |>
  mutate(Order.q = case_when(Diversity == "Species Richness" ~ 0,
                             Diversity == "Shannon diversity" ~ 1,
                             Diversity == "Simpson diversity" ~ 2))
ecoli_bsi_arg_asy$Order.q <- factor(ecoli_bsi_arg_asy$Order.q, levels = sort(unique(ecoli_bsi_arg_asy$Order.q)))

# plot counts of ARG by sample size, both size-based and coverage based CIs
# colour palette (discrete, colorblind-friendly)
#n_q <- length(unique(ecoli_bsi_iNEXT_size_based$Order.q))
pal <- viridis(n = 6, option = "D")
"#440154FF" "#414487FF" "#2A788EFF" "#22A884FF" "#7AD151FF" "#FDE725FF"

# refine max on y-axis as a scale factor to plot sample coverage on same scale
max_qD <- max(ecoli_bsi_arg_size_based$qD.UCL, na.rm = TRUE)
cov_scale <- max_qD / 1     # equivalently = max_qD

# plot total with size based and coverage based CIs
p <- ggplot(data = ecoli_bsi_arg_size_based) +
  # size-based ribbons and lines (primary axis)
  geom_line( aes(x = t, y = qD, colour = Order.q, group = Order.q), size = 1) +
  geom_ribbon( aes(x = t, ymin = qD.LCL, ymax = qD.UCL, fill = Order.q, group = Order.q),
               alpha = 0.35, colour = NA) +
  # coverage-based (overlay) lighter ribbons and dashed lines
  geom_ribbon(aes(x = ecoli_bsi_arg_coverage_based$t, 
                  ymin = ecoli_bsi_arg_coverage_based$qD.LCL, 
                  ymax = ecoli_bsi_arg_coverage_based$qD.UCL, fill = Order.q, 
                  group = Order.q),
              alpha = 0.18, colour = NA) +
  # add horizontal dashed asymptotes from AsyEst (for Total assemblage)
  geom_hline(data = ecoli_bsi_arg_asy, aes(yintercept = Estimator, colour = Order.q),
             linetype = "dashed", size = 0.6) +
  geom_text(data = ecoli_bsi_arg_asy, aes(x = 0, y = Estimator, label = round(Estimator), colour = Order.q),
            hjust = 1, vjust = -0.3, size = 3.5, fontface = "bold", show.legend = FALSE )+
  # add the scaled coverage line (RHS axis mapping) from the size-based df
  geom_line(aes(x = t, y = SC * cov_scale), colour = "black", size = 0.8, inherit.aes = FALSE) +
  geom_ribbon(aes(x = t, ymin = SC.LCL * cov_scale, ymax = SC.UCL * cov_scale),
              fill = "grey60", alpha = 0.3, inherit.aes = FALSE) +
  # vertical reference lines for sample sizes
  geom_vline(xintercept = 1471, linetype = "longdash", colour = "grey40", size = 1) +
  geom_vline(xintercept = 2263, linetype = "longdash", colour = "grey40" , size = 1) +
  geom_vline(xintercept = 18732, linetype = "longdash", colour = "grey40", size = 1) +
  # scales: primary y for qD, secondary y for coverage
  scale_colour_manual(name = "Diversity order", values = pal, labels = order_labels) +
  scale_fill_manual(name = "Diversity order", values = pal, labels = order_labels, guide = guide_legend(override.aes = list(alpha = 0.4))) +
  scale_y_log10( name = "ARG Diversity (qD)",
                 # limits = c(0, max_qD * 1.05),
                 sec.axis = sec_axis(~ . / cov_scale, name = "ARG Sample coverage (SC)")) +
  scale_x_continuous(name = "Sample size (m)") +
  labs(title = "iNEXT ARG count accummulation curves",
       subtitle = "Dark ribbon= size-based CI; light= coverage-based; black/grey = sample coverage (RHS axis)",
       caption = "Horizontal dashed asymptote estimates from iNEXT$AsyEst"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5)
  )
p

#save
ggsave(filename = "rarefaction/ecoli_bsi_iNEXT_arg_logscale.png", plot = p, width = 10, height = 7, units = "in", dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Klebsiella iNEXT repeat ####
# * Abundance data (when counting how may of MLSTs present) ####
# * * sample size based sampling curve ####
length(unique(kleb_bsi_samples_metadata$mlst_profile)) # 297 unique mlsts
# prepare input data as matrix with rownames as mlst_profiles and column names as sites, with first column being overall. 
kleb_bsi_mlst_abundance_df <- kleb_bsi_samples_metadata |>
  group_by(mlst_profile, region) |>
  summarise(n = n()) |>
  pivot_wider(names_from = region, values_from = n, values_fill = 0) |>
  mutate(Total = sum(Midlands, `North East A`, `North East B`, `North West`, `South West`, 
                                `East`, `London`,`South East A`, `South East B`, `South East C`)) |>
  dplyr::select(mlst_profile, Total, everything())
kleb_bsi_mlst_abundance_df <- data.frame(kleb_bsi_mlst_abundance_df)
rownames(kleb_bsi_mlst_abundance_df) <- kleb_bsi_mlst_abundance_df[,1]
kleb_bsi_mlst_abundance_df <- kleb_bsi_mlst_abundance_df[,-1]
#View(kleb_bsi_mlst_abundance_df)
sum(kleb_bsi_mlst_abundance_df$Total) # 468 check

# get vector of sample sizes for which to estimate diversity
# reuse the actual_sample_sizes integer vector from above
# define current and estiamted upsampling sample sizes for E.coli BSIs:
current_n <- nrow(kleb_bsi_samples_metadata) #468
estimated_total_n <- current_n / 0.65 #720
estimated_national_n <- 4933

#unique mlsts observed currently
actual_unique_mlsts <- length(unique(kleb_bsi_samples_metadata$mlst_profile)) # 297

# define sample sizes
relative_sample_sizes <- c(seq(0,15, 0.1), (estimated_total_n/current_n), (estimated_national_n / current_n))
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
saveRDS(kleb_bsi_iNEXT, "rarefaction/kleb_bsi_iNEXT.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in ####
#kleb_bsi_iNEXT <- readRDS("rarefaction/kleb_bsi_iNEXT.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# covert p to factor for plotting
kleb_bsi_iNEXT_size_based <- as.data.frame(kleb_bsi_iNEXT$iNextEst$size_based)
kleb_bsi_iNEXT_size_based$Order.q <- factor(kleb_bsi_iNEXT_size_based$Order.q, levels = sort(unique(kleb_bsi_iNEXT_size_based$Order.q)))
kleb_bsi_iNEXT_coverage_based <- as.data.frame(kleb_bsi_iNEXT$iNextEst$coverage_based)
kleb_bsi_iNEXT_coverage_based$Order.q <- factor(kleb_bsi_iNEXT_coverage_based$Order.q)
kleb_bsi_asy <- as.data.frame(kleb_bsi_iNEXT$AsyEst)
kleb_bsi_asy$Diversity <- factor(kleb_bsi_asy$Diversity)
View(kleb_bsi_iNEXT_size_based)


# split into regional and total data
kleb_bsi_iNEXT_size_based_total <- kleb_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total" & m != 0 & m != 1)
kleb_bsi_iNEXT_size_based_regional <- kleb_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage != "Total" & m != 0 & m != 1)

kleb_bsi_iNEXT_coverage_based_total <- kleb_bsi_iNEXT_coverage_based |>
  dplyr::filter(Assemblage == "Total" & m != 0 & m != 1 )
kleb_bsi_iNEXT_coverage_based_regional <- kleb_bsi_iNEXT_coverage_based |>
  dplyr::filter(Assemblage != "Total" & m != 0& m != 1)

kleb_bsi_asy <- kleb_bsi_asy |>
  mutate(Order.q = case_when(Diversity == "Species richness" ~ 0,
                             Diversity == "Shannon diversity" ~ 1,
                             Diversity == "Simpson diversity" ~ 2))
kleb_bsi_asy$Order.q <- factor(kleb_bsi_asy$Order.q, levels = sort(unique(kleb_bsi_asy$Order.q)))
kleb_bsi_asy_total <- kleb_bsi_asy |>
  dplyr::filter(Assemblage == "Total")
kleb_bsi_asy_regional <- kleb_bsi_asy |>
  dplyr::filter(Assemblage != "Total")


# ChatGPT summary or coverage based vs size based CIs:
# “Coverage-based confidence intervals for MLST richness were substantially wider
# than size-based intervals, reflecting high uncertainty in the number of rare,
# unobserved MLSTs. This effect was pronounced for species richness (q = 0) but 
# minimal for Shannon (q = 1) and Simpson (q = 2) diversity, which are dominated 
# by common MLSTs and are therefore less sensitive to incomplete sampling.”

# ChatGPT definition of sampling coverage:
#Sampling coverage answers this question:
#  “What fraction of the total probability mass (i.e. individuals) in the community
#  belongs to species that we have already observed?”
# Equivalently:
#   “If I take one more individual at random from the population, what is the 
#   probability it belongs to a species I’ve already seen?”

# ensure Order.q is a factor with nice labels
order_labels <- c("0" = "Species richness (q=0)",
                  "1" = "Shannon diversity (q=1)",
                  "2" = "Simpson diversity (q=2)")

# compute a scale factor to map coverage (0..1) onto the qD scale for the RHS axis
# we take the maximum upper CI of qD from the size-based df as the top of primary axis
max_qD <- max(kleb_bsi_iNEXT_size_based_total$qD.UCL, kleb_bsi_iNEXT_size_based_total$qD, na.rm = TRUE)
# ensure max_qD > 0
if (max_qD <= 0 || is.infinite(max_qD) || is.na(max_qD)) max_qD <- max(size_df$qD, na.rm = TRUE)

# the coverage values typically are in 0..1; multiply them by this factor
cov_scale <- max_qD / 1     # equivalently = max_qD

# colour palette (discrete, colorblind-friendly)
#n_q <- length(unique(kleb_bsi_iNEXT_size_based$Order.q))
pal <- viridis(n = 6, option = "D")
"#440154FF" "#414487FF" "#2A788EFF" "#22A884FF" "#7AD151FF" "#FDE725FF"


# plot total with size based and coverage based CIs
p <- ggplot(data = kleb_bsi_iNEXT_size_based_total) +
  # size-based ribbons and lines (primary axis)
  geom_line( aes(x = m, y = qD, colour = Order.q, group = Order.q), size = 1) +
  geom_ribbon( aes(x = m, ymin = qD.LCL, ymax = qD.UCL, fill = Order.q, group = Order.q),
    alpha = 0.35, colour = NA) +

  # coverage-based (overlay) lighter ribbons and dashed lines
  geom_ribbon(aes(x = kleb_bsi_iNEXT_coverage_based_total$m, 
                  ymin = kleb_bsi_iNEXT_coverage_based_total$qD.LCL, 
                  ymax = kleb_bsi_iNEXT_coverage_based_total$qD.UCL, fill = Order.q, 
                  group = Order.q),
    alpha = 0.18, colour = NA) +
  
  # add horizontal dashed asymptotes from AsyEst (for Total assemblage)
  geom_hline(data = kleb_bsi_asy_total, aes(yintercept = Estimator, colour = Order.q),
    linetype = "dashed", size = 0.6) +
  geom_text(data = kleb_bsi_asy_total,aes(x = 0, y = Estimator, label = round(Estimator), colour = Order.q),
    hjust = 1, vjust = -0.3, size = 3.5, fontface = "bold", show.legend = FALSE )+
  
  # add the scaled coverage line (RHS axis mapping) from the size-based df
  geom_line(aes(x = m, y = SC * cov_scale), colour = "black", size = 0.8, inherit.aes = FALSE) +
  geom_ribbon(aes(x = m, ymin = SC.LCL * cov_scale, ymax = SC.UCL * cov_scale),
   fill = "grey60", alpha = 0.3, inherit.aes = FALSE) +
  
  # vertical reference lines for sample sizes
  geom_vline(xintercept = 1471, linetype = "longdash", colour = "grey40", size = 1) +
  geom_vline(xintercept = 2263, linetype = "longdash", colour = "grey40" , size = 1) +
  geom_vline(xintercept = 18732, linetype = "longdash", colour = "grey40", size = 1) +
  
  # scales: primary y for qD, secondary y for coverage
  scale_colour_manual(name = "Diversity order", values = pal, labels = order_labels) +
  scale_fill_manual(name = "Diversity order", values = pal, labels = order_labels, guide = guide_legend(override.aes = list(alpha = 0.4))) +
  
  scale_y_log10( name = "Diversity (qD)",
   # limits = c(0, max_qD * 1.05),
    sec.axis = sec_axis(~ . / cov_scale, name = "Sample coverage (SC)")) +
  scale_x_continuous(name = "Sample size (m)") +
  labs(
    title = "iNEXT size- and coverage-based diversity curves",
    subtitle = "Dark ribbon= size-based CI; light= coverage-based; black/grey = sample coverage (RHS axis)",
    caption = "Horizontal dashed asymptote estimates from iNEXT$AsyEst"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5)
  )
  
p

#save
ggsave(filename = "rarefaction/kleb_bsi_iNEXT_mlst_total_logscale.png", plot = p, width = 10, height = 7, units = "in", dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Compare Preseq and iNEXT ####
# plot just total species richness (q=0) on same axes as preseqR estimate
kleb_bsi_iNEXT_subset <- kleb_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total", Order.q == 0) |>
  mutate(Method = "iNEXT", r = "iNEXT") |>
  rename(actual_sample_size = m ,
         est = qD,
         lower = qD.LCL, 
         upper = qD.UCL) |>
  dplyr::select(c(Method, r, actual_sample_size, est, lower, upper))

# load preseq data
kleb_bsi_mlst_preseq <- read.csv("rarefaction/kleb_bsi_upsampling_by_r.csv")
#View(kleb_bsi_mlst_preseq)
kleb_bsi_mlst_preseq_filtered <- kleb_bsi_mlst_preseq |>
  dplyr::filter(r == 1) |>
  mutate(Method = "preseqR") |>
  dplyr::select(Method, r, actual_sample_size, est, lower, upper)
#rbind
preseq_iNEXT_combined <- rbind(kleb_bsi_iNEXT_subset, kleb_bsi_mlst_preseq_filtered)
#View(preseq_iNEXT_combined)

# plot
p <- ggplot(data  = preseq_iNEXT_combined, aes(x = actual_sample_size, y = est, colour = Method, fill = Method)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_ribbon(data = kleb_bsi_iNEXT_coverage_based_total |> dplyr::filter(Order.q == 0) , aes(x = m, y = qD, ymin = qD.LCL, ymax = qD.UCL), alpha = 0.3, colour = NA, fill = "salmon") +
  # vertical reference lines for sample sizes
  geom_vline(xintercept = 1471, linetype = "longdash", colour = "grey20", size = 1) +
  geom_vline(xintercept = 2263, linetype = "longdash", colour = "grey20" , size = 1) +
  geom_vline(xintercept = 18732, linetype = "longdash", colour = "grey20", size = 1) +
  #scale_colour_viridis_d(name = "Method", option = "B") +
  #scale_fill_viridis_d(name = "Method", option = "B", guide = "none") +
  labs(
    x = "Sample Size",
    y = "Species Richness"
  ) +
  theme_minimal()
p
#save
ggsave(filename = "rarefaction/kleb_bsi_mlst_iNEXT_vs_preseq_count.png", plot = p, width = 7, height = 5, units = "in", dpi = 300 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# summary table for species richness, Shannon and Simpson diversity, at important sample sizes. 
kleb_bsi_iNEXT_summary <- kleb_bsi_iNEXT_size_based |>
  dplyr::filter(m %in% c(current_n, round(estimated_total_n), estimated_national_n)) |>
  mutate(Method = "iNEXT") |>
  dplyr::select(c(Method, Assemblage, Order.q, m, qD, qD.LCL, qD.UCL)) |>
  rename(actual_sample_size = m ,
         est = qD,
         lower = qD.LCL, 
         upper = qD.UCL) |>
  mutate(Diversity = case_when(Order.q == 0 ~ "Species richness (q=0)",
                           Order.q == 1 ~ "Shannon diversity (q=1)",
                           Order.q == 2 ~ "Simpson diversity (q=2)")) |>
  mutate(Estimated_MLST_count = paste0(round(est), " (", round(lower), " - ", round(upper), ")")) |>
  dplyr::select(c(Method, Assemblage, Order.q, Diversity, actual_sample_size, Estimated_MLST_count, est, lower, upper))
View(kleb_bsi_iNEXT_summary)

#save
write.csv(kleb_bsi_iNEXT_summary, "rarefaction/kleb_bsi_iNEXT_mlst_count_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# repeat for sampling coverage
kleb_bsi_iNEXT_subset_samcov <- kleb_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total", Order.q == 0) |>
  mutate(Method = "iNEXT") |>
  rename(actual_sample_size = m ,
         est = SC,
         lower = SC.LCL, 
         upper = SC.UCL) |>
  dplyr::select(c(Method, actual_sample_size, est, lower, upper))
View(kleb_bsi_iNEXT_subset_samcov)

# load preseq data
all_results_sample_cov_df <- read.csv("rarefaction/kleb_bsi_upsampling_sample_cov_by_r.csv")
#View(all_results_sample_cov_df)
kleb_bsi_mlst_preseq_samcov_filtered <- all_results_sample_cov_df |>
  dplyr::filter(r == 1) |>
  mutate(Method = "preseqR") |>
  dplyr::select(Method, actual_sample_size, est, lower, upper)
#rbind
preseq_iNEXT_combined_samcov <- rbind(kleb_bsi_iNEXT_subset_samcov, kleb_bsi_mlst_preseq_samcov_filtered)
#View(preseq_iNEXT_combined_samcov)

# plot
p <- ggplot(data  = preseq_iNEXT_combined_samcov, aes(x = actual_sample_size, y = est, colour = Method, fill = Method)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  # vertical reference lines for sample sizes
  geom_vline(xintercept = 1471, linetype = "longdash", colour = "grey20", size = 1) +
  geom_vline(xintercept = 2263, linetype = "longdash", colour = "grey20" , size = 1) +
  geom_vline(xintercept = 18732, linetype = "longdash", colour = "grey20", size = 1) +
  #scale_colour_viridis_d(name = "Method", option = "B") +
  #scale_fill_viridis_d(name = "Method", option = "B", guide = "none") +
  #scale_y_log10() +
  labs(
    x = "Sample Size",
    y = "Sample Coverage\nProportion of all individuals belonging to observed MLSTs"
  ) +
  theme_minimal()
p
#save
ggsave(filename = "rarefaction/kleb_bsi_mlst_iNEXT_vs_preseq_sample_coverage.png", plot = p, width = 7, height = 5, units = "in", dpi = 300 )

# * * Summary table for sampling coverage joint for preseQ (r=1) and iNEXT ####
# At various sampling depths, what proportion of the total community will belong to species we have already sampled?

preseq_iNEXT_combined_samcov$actual_sample_size <- round(preseq_iNEXT_combined_samcov$actual_sample_size)
kleb_bsi_preseq_iNEXT_samcov_summary <- preseq_iNEXT_combined_samcov |>
  dplyr::filter(actual_sample_size %in% c(current_n, round(estimated_total_n), estimated_national_n))|>
  mutate(Estimated_sampling_coverage = paste0(round(est*100, 1), " (", round(lower*100, 1), " - ", round(upper*100,1), ")")) |>
  dplyr::select(c(Method, actual_sample_size, Estimated_sampling_coverage)) |>
  pivot_wider(id_cols = actual_sample_size, names_from = Method, values_from = Estimated_sampling_coverage)
View(kleb_bsi_preseq_iNEXT_samcov_summary)

#save
write.csv(kleb_bsi_preseq_iNEXT_samcov_summary, "rarefaction/kleb_bsi_preseq_iNEXT_samcov_summary.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~#
# * Repeat for incidence data ####
# (using this as ARG presence/ absence, and bacterium is the 'sampling unit')
# load amrfinder metadata
amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
View(amrfinder_metadata_updated)
# Filter kleb BSI specific data
kleb_bsi_amrfinder_metadata <- amrfinder_metadata_updated |>
  filter(grepl("Escherichia", kraken2_species) & (sampletype == "blood" | sampletype == "unknown" | is.na(sampletype))) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
length(unique(kleb_bsi_amrfinder_metadata$sample)) #1471 non-duplicates
length(unique(kleb_bsi_amrfinder_metadata$Element.symbol)) #332 = 331 unique ARGs and NA
# 189 different AMR genes
# consider filtering for identity and coverage, but in practice, most >95% for both


kleb_bsi_arg_presence_absence_matrix <- kleb_bsi_amrfinder_metadata |>
  dplyr::filter(Type == "AMR") |>
  group_by(sample, Element.symbol) |>
  summarise(n = n()) |>
  mutate(presence = ifelse(n >= 1, 1, 0)) |>
  ungroup() |>
  pivot_wider(id_cols = Element.symbol, names_from = sample, values_from = presence, values_fill = 0) |>
  dplyr::filter(!is.na(Element.symbol)) # remove NAs
#View(kleb_bsi_arg_presence_absence_matrix) 
kleb_bsi_arg_presence_absence_matrix <- as.data.frame(kleb_bsi_arg_presence_absence_matrix)
rownames(kleb_bsi_arg_presence_absence_matrix) <- kleb_bsi_arg_presence_absence_matrix$Element.symbol
kleb_bsi_arg_presence_absence_matrix <- kleb_bsi_arg_presence_absence_matrix[,-1]

#transform to incidence frequency input vectror needed for iNEXT
kleb_bsi_arg_incfreq <- as.incfreq(kleb_bsi_arg_presence_absence_matrix)
#View(kleb_bsi_arg_incfreq)

kleb_bsi_arg_iNEXT <- iNEXT(kleb_bsi_arg_incfreq,
  q = c(0,1,2),  datatype = "incidence_freq", size = actual_sample_sizes, 
 # endpoint = NULL, knots = 40,
  se = TRUE, conf = 0.95, nboot = 100)

# inspect results
kleb_bsi_arg_iNEXT$DataInfo # 99.46% sample coverage of observed 189 ARGs
kleb_bsi_arg_iNEXT$iNextEst$size_based
kleb_bsi_arg_iNEXT$iNextEst$coverage_based
kleb_bsi_arg_iNEXT$AsyEst

#save
saveRDS(kleb_bsi_arg_iNEXT, "rarefaction/kleb_bsi_arg_iNEXT.rds")
#~~~~~~~~~~~~~~~#
# read in saved
#kleb_bsi_arg_iNEXT <- readRDS("kleb_bsi_arg_iNEXT.rds")
#~~~~~~~~~~~~~~~#

# get outputs into manageable data frame form
kleb_bsi_arg_size_based <- as.data.frame(kleb_bsi_arg_iNEXT$iNextEst$size_based)
kleb_bsi_arg_size_based$Order.q <- factor(kleb_bsi_arg_size_based$Order.q, levels = sort(unique(kleb_bsi_arg_size_based$Order.q)))
kleb_bsi_arg_size_based <- kleb_bsi_arg_size_based |>
  dplyr::filter(t != 0 & t != 1) # remove zero and one sample size rows

kleb_bsi_arg_coverage_based <- as.data.frame(kleb_bsi_arg_iNEXT$iNextEst$coverage_based)
kleb_bsi_arg_coverage_based$Order.q <- factor(kleb_bsi_arg_coverage_based$Order.q, levels = sort(unique(kleb_bsi_arg_coverage_based$Order.q)))
kleb_bsi_arg_coverage_based <- kleb_bsi_arg_coverage_based |>
  dplyr::filter(t != 0 & t != 1) 

kleb_bsi_arg_asy <- as.data.frame(kleb_bsi_arg_iNEXT$AsyEst)
kleb_bsi_arg_asy$Diversity <- rownames(kleb_bsi_arg_asy)
kleb_bsi_arg_asy <- kleb_bsi_arg_asy |>
  mutate(Order.q = case_when(Diversity == "Species Richness" ~ 0,
                             Diversity == "Shannon diversity" ~ 1,
                             Diversity == "Simpson diversity" ~ 2))
kleb_bsi_arg_asy$Order.q <- factor(kleb_bsi_arg_asy$Order.q, levels = sort(unique(kleb_bsi_arg_asy$Order.q)))

# plot counts of ARG by sample size, both size-based and coverage based CIs
# colour palette (discrete, colorblind-friendly)
#n_q <- length(unique(kleb_bsi_iNEXT_size_based$Order.q))
pal <- viridis(n = 6, option = "D")
"#440154FF" "#414487FF" "#2A788EFF" "#22A884FF" "#7AD151FF" "#FDE725FF"

# refine max on y-axis as a scale factor to plot sample coverage on same scale
max_qD <- max(kleb_bsi_arg_size_based$qD.UCL, na.rm = TRUE)
cov_scale <- max_qD / 1     # equivalently = max_qD

# plot total with size based and coverage based CIs
p <- ggplot(data = kleb_bsi_arg_size_based) +
  # size-based ribbons and lines (primary axis)
  geom_line( aes(x = t, y = qD, colour = Order.q, group = Order.q), size = 1) +
  geom_ribbon( aes(x = t, ymin = qD.LCL, ymax = qD.UCL, fill = Order.q, group = Order.q),
               alpha = 0.35, colour = NA) +
  # coverage-based (overlay) lighter ribbons and dashed lines
  geom_ribbon(aes(x = kleb_bsi_arg_coverage_based$t, 
                  ymin = kleb_bsi_arg_coverage_based$qD.LCL, 
                  ymax = kleb_bsi_arg_coverage_based$qD.UCL, fill = Order.q, 
                  group = Order.q),
              alpha = 0.18, colour = NA) +
  # add horizontal dashed asymptotes from AsyEst (for Total assemblage)
  geom_hline(data = kleb_bsi_arg_asy, aes(yintercept = Estimator, colour = Order.q),
             linetype = "dashed", size = 0.6) +
  geom_text(data = kleb_bsi_arg_asy, aes(x = 0, y = Estimator, label = round(Estimator), colour = Order.q),
            hjust = 1, vjust = -0.3, size = 3.5, fontface = "bold", show.legend = FALSE )+
  # add the scaled coverage line (RHS axis mapping) from the size-based df
  geom_line(aes(x = t, y = SC * cov_scale), colour = "black", size = 0.8, inherit.aes = FALSE) +
  geom_ribbon(aes(x = t, ymin = SC.LCL * cov_scale, ymax = SC.UCL * cov_scale),
              fill = "grey60", alpha = 0.3, inherit.aes = FALSE) +
  # vertical reference lines for sample sizes
  geom_vline(xintercept = 1471, linetype = "longdash", colour = "grey40", size = 1) +
  geom_vline(xintercept = 2263, linetype = "longdash", colour = "grey40" , size = 1) +
  geom_vline(xintercept = 18732, linetype = "longdash", colour = "grey40", size = 1) +
  # scales: primary y for qD, secondary y for coverage
  scale_colour_manual(name = "Diversity order", values = pal, labels = order_labels) +
  scale_fill_manual(name = "Diversity order", values = pal, labels = order_labels, guide = guide_legend(override.aes = list(alpha = 0.4))) +
  scale_y_log10( name = "ARG Diversity (qD)",
                 # limits = c(0, max_qD * 1.05),
                 sec.axis = sec_axis(~ . / cov_scale, name = "ARG Sample coverage (SC)")) +
  scale_x_continuous(name = "Sample size (m)") +
  labs(title = "iNEXT ARG count accummulation curves",
       subtitle = "Dark ribbon= size-based CI; light= coverage-based; black/grey = sample coverage (RHS axis)",
       caption = "Horizontal dashed asymptote estimates from iNEXT$AsyEst"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5)
  )
p

#save
ggsave(filename = "rarefaction/kleb_bsi_iNEXT_arg_logscale.png", plot = p, width = 10, height = 7, units = "in", dpi = 300 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * RichnEst ####

# prepare occupancy vector as input for python-baced CLI implementation
occupancy_vector <- ecoli_bsi_samples_metadata |>
  group_by(mlst_profile) |>
  summarise(n = n()) |>
  arrange(mlst_profile)
# fill in the 0s
occupancy_vector$n <- factor(occupancy_vector$n, levels = c(seq(1, max(occupancy_vector$n), 1)))

# tabulate frequency of each occurence frequency
tab <- table(occupancy_vector$n)
# make sure it's numeric counts with numeric names
k_vals <- as.integer(names(tab))    # observed k values (1,2,3,...,251 in your example)
counts <- as.integer(tab)          # corresponding o_k values

# build occupancy vector o where o[k] = number of objects with count k
max_k <- max(k_vals)
o <- integer(max_k)                # initializes zeros for k = 1..max_k
o[k_vals] <- counts

# check results
o           # occupancy vector: position k -> o[k]
length(o)   # equals max_k # 251
# total number of distinct objects observed (species)
S <- sum(o)
# total sample size (sum of counts n_i) = sum(k * o_k)
T <- sum(seq_along(o) * o)
# check total unique MLST and total sample size
S; T

# write file as one-line vector, space separated
cat(o, file = "rarefaction/occupancy_oneline.txt", sep = " ")

#~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~#
# SAMPLE COVERAGE ESTIMATORS ####
# Quick look at different Estimators for sample Coverage 
# Good-Turing estimator
1-(182/1471) # close enough to 88%

# and Chao1 estimators (when f2 != 0)
278 + (182 * 182)/(2 * 42) #= 475 # 672
# and so empirircal coverage is 

# and chao1 bias corrected estimator (when f2 ==0)
278 + (182 * (182-1))/(2*(42+1)) # 661

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# BAYESIAN BOOTSTRAPPING ####
# ESTIMATES POSTERIOR FREQUENCY DISTRIBUTION OF MLSTS
# use prior of a dirichlet distribution
# install package
#install.packages("bayesboot")
# load package
library(bayesboot)
library(dplyr)
library(tidyr)
library(ggplot2)
#plotting
install.packages("forcats")
install.packages("patchwork")
install.packages("scales")
library(forcats)
library(patchwork)
library(scales)

# load data
ecoli_bsi_samples_metadata <- read.csv("rarefaction/ecoli_bsi_samples_metadata.csv")
kleb_bsi_samples_metadata <- read.csv("rarefaction/kleb_bsi_samples_metadata.csv")


# prepare data - E. coli
ecoli_bsi_mlst_df <- ecoli_bsi_samples_metadata |>
  dplyr::rename(kleborate_mlst = escherichia__mlst_achtman__ST) |>
  dplyr::group_by(kleborate_mlst) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(kleborate_mlst, count) |>
  dplyr::arrange(count, kleborate_mlst)
#View(ecoli_bsi_mlst_df)

# prepare data - Kleb
kleb_bsi_mlst_df <- kleb_bsi_samples_metadata |>
  dplyr::rename(kleborate_mlst = klebsiella_mlst_ST) |>
  dplyr::group_by(kleborate_mlst) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(kleborate_mlst, count) |>
  dplyr::arrange(count, kleborate_mlst)
#View(kleb_bsi_mlst_df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# look at non-bootstrapped empirical frequency distributions of E. coli (before bootstrapping) combined with the power calculation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot power calculation table
N_isolates <- n_distinct(ecoli_bsi_samples_metadata$isolateid)

freq_df <- ecoli_bsi_samples_metadata |>
  group_by(kleborate_mlst) |>
  dplyr::summarise(count = n(), .groups = "drop" ) |>
  dplyr::mutate(frequency = count/N_isolates) |>
  dplyr::arrange(count) |>
  mutate(
    cumulative_frequency = sapply(count, function(c) sum(frequency[count >= c])),
    cumulative_mlst_no = sapply(count, function(c) sum(count >=c)),
    cumulative_isolate_count = sapply(count, function(c) sum(count[count >= c])),
    empirical_sample_coverage = cumulative_isolate_count / N_isolates,
    # min sample sizes to observe this MLST at least once with given prob
    min_sample_90 = log(1-0.90)/log(1-frequency),
    min_sample_95 = log(1-0.95)/log(1-frequency),
    min_sample_99 = log(1-0.99)/log(1-frequency)
  )  
#View(freq_df)

# plot frequency vs cumulative frequency
freq_vs_cum_freq <- ggplot(data = freq_df, aes(x = frequency, y = cumulative_frequency)) +
  geom_line() +
  geom_point() +
  theme_minimal()
freq_vs_cum_freq

ss_vs_freq <- ggplot(data = freq_df, aes(x = frequency, y = min_sample_95)) +
  geom_line(aes(y = min_sample_90), colour = "black") +
  geom_line(aes(y = min_sample_95), colour = "blue") +
  geom_line(aes(y = min_sample_99), colour = "lightblue") +
  geom_point(aes(y = min_sample_90), colour = "black") +
  geom_point(aes(y = min_sample_95), colour = "blue") +
  geom_point(aes(y = min_sample_99), colour = "lightblue") +
  theme_minimal()

ss_vs_freq

ss_vs_sample_cov <- ggplot(data = freq_df, aes(x = empirical_sample_coverage, y = min_sample_95)) +
  geom_line(aes(y = min_sample_90), colour = "black") +
  geom_line(aes(y = min_sample_95), colour = "blue") +
  geom_line(aes(y = min_sample_99), colour = "lightblue") +
  geom_point(aes(y = min_sample_90), colour = "black") +
  geom_point(aes(y = min_sample_95), colour = "blue") +
  geom_point(aes(y = min_sample_99), colour = "lightblue") +
  theme_minimal()
ss_vs_sample_cov # wobbly line as noise from limited NEKSUS sample



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
#  * * check MLST count distribution  ####
# prep data
ecoli_bsi_mlst_df <- ecoli_bsi_samples_metadata |>
  dplyr::rename(kleborate_mlst = escherichia__mlst_achtman__ST) |>
  dplyr::group_by(kleborate_mlst) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(kleborate_mlst, count) |>
  dplyr::arrange(count, kleborate_mlst)
#View(ecoli_bsi_mlst_df)
table(ecoli_bsi_mlst_df$count)
nrow(ecoli_bsi_mlst_df) # 263
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
saveRDS(ecoli_bsi_mlst_bayesboot_non_exact, "rarefaction/ecoli_bsi_mlst_bayesboot_non_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in ####
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
# * * Reload saved data ####
#ecoli_mass_df <- read.csv("ecoli_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv")
#ecoli_mass_df_non_exact <- read.csv("ecoli_bsi_kelborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * repeat for level 1 fastbaps partition  ####
ecoli_bsi_fastbaps_L1 <- ecoli_bsi_samples_metadata |>
  dplyr::group_by(Level.1) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(Level.1, count) |>
  dplyr::arrange(count, Level.1)
#View(ecoli_bsi_fastbaps_L1)

# check count distribution
table(ecoli_bsi_fastbaps_L1$count)
nrow(ecoli_bsi_fastbaps_L1) # 21
# estimates for alpha novel
# prior for alpha_novel
# df has columns mlst_profile and count
df <- ecoli_bsi_fastbaps_L1
N <- sum(df$count)
f1 <- sum(df$count == 1)  #       # number of singletons
q_hat <- f1 / N  # Good-Turing first-order - proportion of singletons
alpha_1 <- 1

# Preferred / exact approach (supports non-integer alphas):
ecoli_bsi_fastbaps_L1_bayesboot_exact <- bayesboot_mlst(ecoli_bsi_fastbaps_L1,
                                                 feature_col = "Level.1",
                                                 alpha_named = NULL, alpha_novel = alpha_1, 
                                                 B = 10000, use_exact_dirichlet = TRUE)
print(ecoli_bsi_fastbaps_L1_bayesboot_exact$summary_df)

# save
saveRDS(ecoli_bsi_fastbaps_L1_bayesboot_exact, "rarefaction/ecoli_bsi_fastbaps_L1_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in ####
#ecoli_bsi_fastbaps_L1_bayesboot_exact <- readRDS("rarefaction/ecoli_bsi_fastbaps_L1_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#

# compute mass >= f curve
ecoli_mass_df <- compute_mass_curve(ecoli_bsi_fastbaps_L1_bayesboot_exact$draws, f_grid = f_grid)

# add sample sizes
ecoli_mass_df <- ecoli_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(ecoli_mass_df, "ecoli_bsi_fastbaps_L1_bayesboot_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Reload saved data ####
#ecoli_mass_df <- read.csv("ecoli_bsi_fastbaps_L1_bayesboot_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * repeat for level 2 fastbaps partition  ####
ecoli_bsi_fastbaps_L2 <- ecoli_bsi_samples_metadata |>
  dplyr::group_by(Level.2) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(Level.2, count) |>
  dplyr::arrange(count, Level.2)
#View(ecoli_bsi_fastbaps_L2)

# check count distribution
table(ecoli_bsi_fastbaps_L2$count)
nrow(ecoli_bsi_fastbaps_L2) # 75
# estimates for alpha novel
# prior for alpha_novel
# df has columns mlst_profile and count
df <- ecoli_bsi_fastbaps_L2
N <- sum(df$count)
f1 <- sum(df$count == 1)  #       # number of singletons
q_hat <- f1 / N  # Good-Turing first-order - proportion of singletons
alpha_1 <- 1

# Preferred / exact approach (supports non-integer alphas):
ecoli_bsi_fastbaps_L2_bayesboot_exact <- bayesboot_mlst(ecoli_bsi_fastbaps_L2,
                                                 feature_col = "Level.2",
                                                 alpha_named = NULL, alpha_novel = alpha_1, 
                                                 B = 10000, use_exact_dirichlet = TRUE)
print(ecoli_bsi_fastbaps_L2_bayesboot_exact$summary_df)

# save
saveRDS(ecoli_bsi_fastbaps_L2_bayesboot_exact, "rarefaction/ecoli_bsi_fastbaps_L2_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in ####
#ecoli_bsi_fastbaps_L2_bayesboot_exact <- readRDS("rarefaction/ecoli_bsi_fastbaps_L2_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#

# compute mass >= f curve
ecoli_mass_df <- compute_mass_curve(ecoli_bsi_fastbaps_L2_bayesboot_exact$draws, f_grid = f_grid)

# add sample sizes
ecoli_mass_df <- ecoli_mass_df |>
  dplyr::mutate(min_sample_90 = log(1-0.90)/log(1-f),
                min_sample_95 = log(1-0.95)/log(1-f),
                min_sample_99 = log(1-0.99)/log(1-f),
                Genus = "Escherichia",
                Exact = "Exact") |>
  filter(f != 0)

# save  mass df
write.csv(ecoli_mass_df, "ecoli_bsi_fastbaps_L2_bayesboot_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Reload saved data ####
#ecoli_mass_df <- read.csv("ecoli_bsi_fastbaps_L2_bayesboot_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * repeat for level 3 fastbaps partition  ####
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
# df has columns mlst_profile and count
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
# Entry point top read prediction back in ####
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
# * * Reload saved data ####
#ecoli_mass_df <- read.csv("ecoli_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Klebsiella BSI MLSTs repeat  ####
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
# Entry point top read prediction back in ####
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
# * * Reload saved data ####
#kleb_mass_df <- read.csv("kleb_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv")
#kleb_mass_df_non_exact <- read.csv("kleb_bsi_kelborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * repeat for fastbaps level 1 clusters ####
# prepare data - Kleb
kleb_bsi_fastbaps_L1_df <- kleb_bsi_samples_metadata |>
  dplyr::group_by(Level.1) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(Level.1, count) |>
  dplyr::arrange(count, Level.1)
#View(kleb_bsi_fastbaps_L1_df)

# check count distribution
table(kleb_bsi_fastbaps_L1_df$count)
nrow(kleb_bsi_fastbaps_L1_df) # 6
# estimates for alpha novel
# prior for alpha_novel
# df has columns mlst_profile and count
df <- kleb_bsi_fastbaps_L1_df
alpha_1 <- 1
alpha_0.1 <- 0.1

# Preferred / exact approach (supports non-integer alphas):
kleb_bsi_fastbaps_L1_bayesboot_exact <- bayesboot_mlst(kleb_bsi_fastbaps_L1_df,
                                             feature_col = "Level.1",
                                             alpha_named = NULL, alpha_novel = 1, 
                                             B = 10000, use_exact_dirichlet = TRUE)

print(kleb_bsi_fastbaps_L1_bayesboot_exact$summary_df)

# save
saveRDS(kleb_bsi_fastbaps_L1_bayesboot_exact, "rarefaction/kleb_bsi_fastbaps_L1_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in ####
#kleb_bsi_fastbaps_L1_bayesboot_exact <- readRDS("rarefaction/kleb_bsi_fastbaps_L1_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
kleb_mass_df <- compute_mass_curve(kleb_bsi_fastbaps_L1_bayesboot_exact$draws, f_grid = f_grid)
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
View(kleb_mass_df)
write.csv(kleb_mass_df, "kleb_bsi_fastbaps_L1_bayesboot_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Reload saved data ####
#kleb_mass_df <- read.csv("kleb_bsi_fastbaps_L1_bayesboot_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * repeat for fastbaps level 2 clusters ####
# prepare data - Kleb
kleb_bsi_fastbaps_L2_df <- kleb_bsi_samples_metadata |>
  dplyr::group_by(Level.2) |>
  dplyr::summarise(count = n()) |>
  dplyr::select(Level.2, count) |>
  dplyr::arrange(count, Level.2)
#View(kleb_bsi_fastbaps_L2_df)

# check count distribution
table(kleb_bsi_fastbaps_L2_df$count)
nrow(kleb_bsi_fastbaps_L2_df) # 15
# estimates for alpha novel
# prior for alpha_novel
# df has columns mlst_profile and count
df <- kleb_bsi_fastbaps_L2_df
N <- sum(df$count) # 468
alpha_1 <- 1

# Preferred / exact approach (supports non-integer alphas):
kleb_bsi_fastbaps_L2_bayesboot_exact <- bayesboot_mlst(kleb_bsi_fastbaps_L2_df,
                                             feature_col = "Level.2",
                                             alpha_named = NULL, alpha_novel = alpha_1, 
                                             B = 10000, use_exact_dirichlet = TRUE)

print(kleb_bsi_fastbaps_L2_bayesboot_exact$summary_df)

# save
saveRDS(kleb_bsi_fastbaps_L2_bayesboot_exact, "rarefaction/kleb_bsi_fastbaps_L2_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# Entry point top read prediction back in ####
#kleb_bsi_fastbaps_L2_bayesboot_exact <- readRDS("rarefaction/kleb_bsi_fastbaps_L2_bayesboot_exact.rds")
#~~~~~~~~~~~~~~~~~~~~~#
# compute mass >= f curve
kleb_mass_df <- compute_mass_curve(kleb_bsi_fastbaps_L2_bayesboot_exact$draws, f_grid = f_grid)
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
write.csv(kleb_mass_df, "kleb_bsi_fastbaps_L2_bayesboot_exact_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Reload saved data ####
#kleb_mass_df <- read.csv("kleb_bsi_fastbaps_L2_bayesboot_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * repeat for fastbaps level 3 clusters ####
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
# Entry point top read prediction back in ####
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
# * * Reload saved data ####
#kleb_mass_df <- read.csv("kleb_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Merge and Save dfs for E.coli and Klebsiella ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * MLST ####
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
# * * Reload saved data ####
#mass_df_mlst <- read.csv("combined_bsi_kleborate_mlst_bayesboot_cumulative_frequency_mass_df.csv")
#mass_df_mlst_non_exact <- read.csv("combined_bsi_kleborate_mlst_bayesboot_non_exact_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * FastBAPS level 1 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#
ecoli_mass_df <- read.csv("ecoli_bsi_fastbaps_L1_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df <- read.csv("kleb_bsi_fastbaps_L1_bayesboot_exact_cumulative_frequency_mass_df.csv")
mass_df_fastbaps_L1 <- rbind(ecoli_mass_df, kleb_mass_df)
#View(mass_df_fastbaps_L1)
# save merged mass df
write.csv(mass_df_fastbaps_L1, "combined_bsi_fastbaps_L1_bayesboot_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Reload saved data ####
#mass_df_fastbaps_L1 <- read.csv("combined_bsi_fastbaps_L1_bayesboot_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * FastBAPS level 2 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#
ecoli_mass_df <- read.csv("ecoli_bsi_fastbaps_L2_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df <- read.csv("kleb_bsi_fastbaps_L2_bayesboot_exact_cumulative_frequency_mass_df.csv")
mass_df_fastbaps_L2 <- rbind(ecoli_mass_df, kleb_mass_df)
#View(mass_df_fastbaps_L2)
# save merged mass df
write.csv(mass_df_fastbaps_L2, "combined_bsi_fastbaps_L2_bayesboot_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Reload saved data ####
#mass_df_fastbaps_L2 <- read.csv("combined_bsi_fastbaps_L2_bayesboot_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * FastBAPS level 3 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~#
ecoli_mass_df <- read.csv("ecoli_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df <- read.csv("kleb_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
mass_df_fastbaps_L3 <- rbind(ecoli_mass_df, kleb_mass_df)
#View(mass_df_fastbaps_L3)
# save merged mass df
write.csv(mass_df_fastbaps_L3, "combined_bsi_fastbaps_L3_bayesboot_cumulative_frequency_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Reload saved data ####
#mass_df_fastbaps_L3 <- read.csv("combined_bsi_fastbaps_L3_bayesboot_cumulative_frequency_mass_df.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# plot
# cumulative mass_df
# frequency and sample_size
# combined plot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * histogram of MLST frequencies ####
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
# * * Cumulative fraction with sample size plot ####
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * fastbaps L1 clustering ####
# Plot (cumulative mass of population at frequency >= f)
cumulative_mass_plot <- ggplot(mass_df_fastbaps_L1, aes(x = f, y = mean_mass, colour =  Genus, fill = Genus)) +
  geom_line(aes(x = f, y = mean_mass, colour =  Genus)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.4, colour = NA) +
  scale_fill_manual(name = "Genus", values = genus_colours) +
  scale_colour_manual(name = "Genus", values = genus_colours) +
  scale_x_log10( breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(xlim = c(min(mass_df_fastbaps_L1$f), 1)) +
  labs(x = "FastBAPS cluster frequency (f) (log scale)",
       y = "Proportion of population belonging to cluster of frequency ≥ f"
      ) +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = -0.5, size = 14)) +
  theme_minimal(base_size = 14)
cumulative_mass_plot
# save
ggsave("rarefaction/combined_bsi_fastbaps_L1_bayesboot_cumulative_mass_plot.png", plot = cumulative_mass_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# * * Cumulative arg fraction with sample size plot ####
bootbayes_sample_coverage_plot <- ggplot(mass_df_fastbaps_L1, aes(x = mean_mass, colour = Genus, fill = Genus)) +
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
ggsave("rarefaction/combined_bsi_fastbaps_L1_bootbayes_sample_coverage_plot_95.png", plot = bootbayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * fastbaps L2 clustering ####
# Plot (cumulative mass of population at frequency >= f)
cumulative_mass_plot <- ggplot(mass_df_fastbaps_L2, aes(x = f, y = mean_mass, colour =  Genus, fill = Genus)) +
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
ggsave("rarefaction/combined_bsi_fastbaps_L2_bayesboot_cumulative_mass_plot.png", plot = cumulative_mass_plot, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#
# * * Cumulative arg fraction with sample size plot ####
bootbayes_sample_coverage_plot <- ggplot(mass_df_fastbaps_L2, aes(x = mean_mass, colour = Genus, fill = Genus)) +
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
ggsave("rarefaction/combined_bsi_fastbaps_L2_bootbayes_sample_coverage_plot_95.png", plot = bootbayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * fastbaps L3 clustering ####
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
# * * Cumulative arg fraction with sample size plot ####
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# frequency vs min sample size to detect that frequency of isolate with 90, 95 or 99% probability
# this is generic plot, independent of the data
label_df <- mass_df |>
  filter(f == min(f), Genus == "Escherichia") |>
  select(f, min_sample_90, min_sample_95, min_sample_99) |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "level",
    values_to = "min_sample"
  ) |>
  mutate(
    label = case_when(
      level == "min_sample_90" ~ "90%",
      level == "min_sample_95" ~ "95%",
      level == "min_sample_99" ~ "99%"
    )
  )

freq_vs_min_sss <- ggplot(mass_df, aes(x = f, y = min_sample_99)) +
  geom_line(aes(y = min_sample_90), colour = "grey80") +
  geom_line(aes(y = min_sample_95), colour = "grey40") +
  geom_line(aes(y = min_sample_99), colour = "grey10") +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.2) +
  # line labels
  geom_text(data = label_df, aes(x = f, y = min_sample, label = label),
    hjust = 1, vjust = 0, size = 4, colour = "black" ) +
  scale_x_log10( breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
  coord_cartesian(clip = "off") +  # allow labels to extend slightly past plot area
  
  labs(x = "Frequency (f) (log scale)",
       y = "Sample size",
       title = "Minimum sample size to detect MLSTs at frequency f",
       ) + 
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 14)) +
  theme_minimal(base_size = 14)
freq_vs_min_sss
# save
ggsave("rarefaction/freq_vs_ss_power_calc.png", plot = freq_vs_min_sss, width = 8, height = 6, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Combined sample coverage estimators 
# Bayesian bootstrapped 
# rarefaction for Good's coverage 
# iNEXT rarefaction/ extrapolation for sample civerage
# preseqR esstimates for sample coverage

# load data
# bayesboot data
#write.csv(mass_df, "combined_bsi_mlst_bayesboot_cumulative_frequency_mass_df.csv", row.names = FALSE)
# load ecoli_mass_df
ecoli_mass_df <- read.csv("ecoli_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv")
ecoli_mass_df_long <- ecoli_mass_df |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(mean_sc = mean_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "MLST Bayesian bootstrap. 90%",
                               Estimator == "min_sample_95" ~ "MLST Bayesian bootstrap. 95%",
                               Estimator == "min_sample_99" ~ "MLST Bayesian bootstrap. 99%", 
                               TRUE ~ Estimator))
# add fastbasp cluster 1
ecoli_mass_df_fastbaps_L1 <- read.csv("ecoli_bsi_fastbaps_L1_bayesboot_exact_cumulative_frequency_mass_df.csv")
ecoli_mass_df_fastaps_L1_long <- ecoli_mass_df_fastbaps_L1 |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(mean_sc = mean_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "fastBAPS L1 Bayesian bootstrap. 90%",
                               Estimator == "min_sample_95" ~ "fastBAPS L1 Bayesian bootstrap. 95%",
                               Estimator == "min_sample_99" ~ "fastBAPS L1 Bayesian bootstrap. 99%", 
                               TRUE ~ Estimator))
# add fastbasp cluster 2
ecoli_mass_df_fastbaps_L2 <- read.csv("ecoli_bsi_fastbaps_L2_bayesboot_exact_cumulative_frequency_mass_df.csv")
ecoli_mass_df_fastaps_L2_long <- ecoli_mass_df_fastbaps_L2 |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(mean_sc = mean_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "fastBAPS L2 Bayesian bootstrap. 90%",
                               Estimator == "min_sample_95" ~ "fastBAPS L2 Bayesian bootstrap. 95%",
                               Estimator == "min_sample_99" ~ "fastBAPS L2 Bayesian bootstrap. 99%", 
                               TRUE ~ Estimator))

# add fastbasp cluster 3
ecoli_mass_df_fastbaps_L3 <- read.csv("ecoli_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
ecoli_mass_df_fastaps_L3_long <- ecoli_mass_df_fastbaps_L3 |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(mean_sc = mean_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "fastBAPS L3 Bayesian bootstrap. 90%",
                               Estimator == "min_sample_95" ~ "fastBAPS L3 Bayesian bootstrap. 95%",
                               Estimator == "min_sample_99" ~ "fastBAPS L3 Bayesian bootstrap. 99%", 
                               TRUE ~ Estimator))

# rarefaction data 
ecoli_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/ecoli_bsi_mlst_rarefaction_combined.csv")
ecoli_bsi_mlst_rarefaction_sc <- ecoli_bsi_mlst_rarefaction_combined |>
  filter(Metric == "Sample coverage" & Replacement == "No replacement") |>
  rename(mean_sc = est, q2.5 = lcl, q97.5 = ucl) |>
  mutate(Estimator = "rarefaction") |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5)

# preseqR data
all_results_sample_cov_df <- read.csv("rarefaction/ecoli_bsi_upsampling_sample_cov_by_r.csv")
all_results_sample_cov_df_r1 <- all_results_sample_cov_df |>  
  filter(r == 1) |>
  rename(sample_size = actual_sample_size, mean_sc = est, q2.5 = lower, q97.5 = upper) |>
  mutate(Estimator = "preseqR") |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5)
#View(all_results_sample_cov_df_r1)

#iNEXT data
ecoli_bsi_iNEXT <- readRDS("rarefaction/ecoli_bsi_iNEXT.rds")
ecoli_bsi_iNEXT_size_based <- as.data.frame(ecoli_bsi_iNEXT$iNextEst$size_based)
ecoli_bsi_iNEXT_size_based$Order.q <- factor(ecoli_bsi_iNEXT_size_based$Order.q, levels = sort(unique(ecoli_bsi_iNEXT_size_based$Order.q)))
ecoli_bsi_iNEXT_size_based_total <- ecoli_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total" & Order.q == 0 ) |>
  rename( sample_size = m, mean_sc = SC, q2.5 = SC.LCL, q97.5 = SC.UCL) |>
  mutate(Estimator = "iNEXT") |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5)
#View(ecoli_bsi_iNEXT_size_based_total)

# bind all rows together
ecoli_bsi_combined_sample_cov_estimators <- rbind(ecoli_mass_df_long,
                                                  ecoli_mass_df_fastaps_L1_long,
                                                  ecoli_mass_df_fastaps_L2_long,
                                                  ecoli_mass_df_fastaps_L3_long,
                                                  ecoli_bsi_mlst_rarefaction_sc,
                                                  all_results_sample_cov_df_r1,
                                                  ecoli_bsi_iNEXT_size_based_total)

# set order of Estimators:
ecoli_bsi_combined_sample_cov_estimators$Estimator <- factor(ecoli_bsi_combined_sample_cov_estimators$Estimator, 
                                                             levels = c("MLST Bayesian bootstrap. 99%",
                                                                        "MLST Bayesian bootstrap. 95%",
                                                                        "MLST Bayesian bootstrap. 90%",
                                                                        "fastBAPS L1 Bayesian bootstrap. 99%",
                                                                        "fastBAPS L1 Bayesian bootstrap. 95%",
                                                                        "fastBAPS L1 Bayesian bootstrap. 90%",
                                                                        "fastBAPS L2 Bayesian bootstrap. 99%",
                                                                        "fastBAPS L2 Bayesian bootstrap. 95%",
                                                                        "fastBAPS L2 Bayesian bootstrap. 90%",
                                                                        "fastBAPS L3 Bayesian bootstrap. 99%",
                                                                        "fastBAPS L3 Bayesian bootstrap. 95%",
                                                                        "fastBAPS L3 Bayesian bootstrap. 90%",
                                                                        "rarefaction",
                                                                        "iNEXT",
                                                                        "preseqR"))
#View(ecoli_bsi_combined_sample_cov_estimators)
unique(ecoli_bsi_combined_sample_cov_estimators$Estimator)
# fix colours
#add key
# make summarybtable 
sam_cov_colours <- c(
  `MLST Bayesian bootstrap. 99%` = "seagreen4",
  `MLST Bayesian bootstrap. 95%` = "seagreen3",
  `MLST Bayesian bootstrap. 90%` = "seagreen1",
  rarefaction = "black", 
  iNEXT = "grey40",
  preseqR = "grey70"
)
names(sam_cov_colours)

# define max sample size
max_actual_sample_coverage <- ecoli_bsi_combined_sample_cov_estimators |>
  filter(Estimator == "rarefaction") |>
  filter(sample_size == max(sample_size)) |>
  pull(mean_sc)

# filter out non-MLST metrics
plot_df <- ecoli_bsi_combined_sample_cov_estimators |>
  filter(Estimator %in% names(sam_cov_colours))

# plot
sample_coverage_comparison_plot <- ggplot(plot_df, aes(x = mean_sc, y = sample_size, colour = Estimator, fill = Estimator)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(xmin = q2.5, xmax = q97.5), alpha = 0.4, colour = NA) +
  # add vertical line for 90%
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  # add solid line for actual collected sample
  geom_vline(xintercept = max_actual_sample_coverage, linetype = "solid") +
  annotate("label", x = 0.72, y = 10000, label = "80% sample coverage", angle = 0, vjust = -0.1, hjust = 0.1, size = 4, colour = "black", fill = "white") +
  annotate("label", x = 0.92, y = 10000, label = "NEKSUS sample collection", angle = 0, vjust = -0.1, hjust = 0.5,  size = 4, colour = "black", fill = "white") +
  scale_fill_manual(name = "Estimator", values = sam_cov_colours) +
  scale_colour_manual(name = "Estimator", values = sam_cov_colours) +
  scale_y_continuous(limits = c(0,10000)) +
  scale_x_continuous(limits = c(0.5, 1.0)) +
  labs(x = "Sample coverage\n(proportion of population belonging to an MLST that has been observed in the sample)",
       y = "Minimum sample size",
       title = "Comparison of sample coverage estimators for E. coli BSI MLSTs",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(5.5, 20, 5.5, 5.5)  # space for vertical labels
  )
sample_coverage_comparison_plot
# save
ggsave("rarefaction/ecoli_bsi_mlst_sample_coverage_comparison_plot.png", plot = sample_coverage_comparison_plot, width = 10, height = 6, units = "in", dpi = 300)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# get summary table of median, and 95% CIS fo sample sizes needed to detect samples at minimum frequency x, and how this translates to population covered.
# 70, 75, 80, 85, 90, 95, 99% of populaiton for both E.coli ans Klebsiella
#View(ecoli_bsi_combined_sample_cov_estimators)
ecoli_bsi_sample_coverage_summary_table <- ecoli_bsi_combined_sample_cov_estimators |>
  group_by(Estimator)

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
ecoli_bsi_sample_coverage_summary_table <- ecoli_bsi_combined_sample_cov_estimators |>
  group_by(Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "mean_sc", t)
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
  select(Estimator, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
ecoli_bsi_sample_coverage_summary_table_cells_only <- ecoli_bsi_sample_coverage_summary_table |>
  select(Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(ecoli_bsi_sample_coverage_summary_table)
print(ecoli_bsi_sample_coverage_summary_table_cells_only)
View(ecoli_bsi_sample_coverage_summary_table)
View(ecoli_bsi_sample_coverage_summary_table_cells_only)

write.csv(ecoli_bsi_sample_coverage_summary_table_cells_only, "rarefaction/ecoli_bsi_sample_coverage_summary_table_cells_only.csv", row.names = FALSE)
write.csv(ecoli_bsi_sample_coverage_summary_table, "rarefaction/ecoli_bsi_sample_coverage_summary_table.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * repeat for klebsiella ####

# bayesboot data
# load kleb_mass_df
kleb_mass_df <- read.csv("kleb_bsi_kleborate_mlst_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df_long <- kleb_mass_df |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(mean_sc = mean_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "MLST Bayesian bootstrap. 90%",
                               Estimator == "min_sample_95" ~ "MLST Bayesian bootstrap. 95%",
                               Estimator == "min_sample_99" ~ "MLST Bayesian bootstrap. 99%", 
                               TRUE ~ Estimator))
# add fastbasp cluster 1
kleb_mass_df_fastbaps_L1 <- read.csv("kleb_bsi_fastbaps_L1_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df_fastaps_L1_long <- kleb_mass_df_fastbaps_L1 |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(mean_sc = mean_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "fastBAPS L1 Bayesian bootstrap. 90%",
                               Estimator == "min_sample_95" ~ "fastBAPS L1 Bayesian bootstrap. 95%",
                               Estimator == "min_sample_99" ~ "fastBAPS L1 Bayesian bootstrap. 99%", 
                               TRUE ~ Estimator))
# add fastbasp cluster 2
kleb_mass_df_fastbaps_L2 <- read.csv("kleb_bsi_fastbaps_L2_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df_fastaps_L2_long <- kleb_mass_df_fastbaps_L2 |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(mean_sc = mean_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "fastBAPS L2 Bayesian bootstrap. 90%",
                               Estimator == "min_sample_95" ~ "fastBAPS L2 Bayesian bootstrap. 95%",
                               Estimator == "min_sample_99" ~ "fastBAPS L2 Bayesian bootstrap. 99%", 
                               TRUE ~ Estimator))

# add fastbasp cluster 3
kleb_mass_df_fastbaps_L3 <- read.csv("kleb_bsi_fastbaps_L3_bayesboot_exact_cumulative_frequency_mass_df.csv")
kleb_mass_df_fastaps_L3_long <- kleb_mass_df_fastbaps_L3 |>
  tidyr::pivot_longer(
    cols = starts_with("min_sample"),
    names_to = "Estimator",
    values_to = "sample_size"
  ) |>
  rename(mean_sc = mean_mass, q2.5 = q2.5, q97.5 = q97.5) |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5) |>
  mutate(Estimator = case_when(Estimator == "min_sample_90" ~ "fastBAPS L3 Bayesian bootstrap. 90%",
                               Estimator == "min_sample_95" ~ "fastBAPS L3 Bayesian bootstrap. 95%",
                               Estimator == "min_sample_99" ~ "fastBAPS L3 Bayesian bootstrap. 99%", 
                               TRUE ~ Estimator))

# rarefaction data 
kleb_bsi_mlst_rarefaction_combined <- read.csv("rarefaction/kleb_bsi_mlst_rarefaction_combined.csv")
kleb_bsi_mlst_rarefaction_sc <- kleb_bsi_mlst_rarefaction_combined |>
  filter(Metric == "Sample coverage" & Replacement == "No replacement") |>
  rename(mean_sc = est, q2.5 = lcl, q97.5 = ucl) |>
  mutate(Estimator = "rarefaction") |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5)

# preseqR data
all_results_sample_cov_df <- read.csv("rarefaction/kleb_bsi_upsampling_sample_cov_by_r.csv")
all_results_sample_cov_df_r1 <- all_results_sample_cov_df |>  
  filter(r == 1) |>
  rename(sample_size = actual_sample_size, mean_sc = est, q2.5 = lower, q97.5 = upper) |>
  mutate(Estimator = "preseqR") |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5)
#View(all_results_sample_cov_df_r1)

#iNEXT data
kleb_bsi_iNEXT <- readRDS("rarefaction/kleb_bsi_iNEXT.rds")
kleb_bsi_iNEXT_size_based <- as.data.frame(kleb_bsi_iNEXT$iNextEst$size_based)
kleb_bsi_iNEXT_size_based$Order.q <- factor(kleb_bsi_iNEXT_size_based$Order.q, levels = sort(unique(kleb_bsi_iNEXT_size_based$Order.q)))
kleb_bsi_iNEXT_size_based_total <- kleb_bsi_iNEXT_size_based |>
  dplyr::filter(Assemblage == "Total" & Order.q == 0 ) |>
  rename( sample_size = m, mean_sc = SC, q2.5 = SC.LCL, q97.5 = SC.UCL) |>
  mutate(Estimator = "iNEXT") |>
  dplyr::select(Estimator, sample_size, mean_sc, q2.5, q97.5)
#View(kleb_bsi_iNEXT_size_based_total)

# bind all rows together
kleb_bsi_combined_sample_cov_estimators <- rbind(kleb_mass_df_long,
                                                  kleb_mass_df_fastaps_L1_long,
                                                  kleb_mass_df_fastaps_L2_long,
                                                  kleb_mass_df_fastaps_L3_long,
                                                  kleb_bsi_mlst_rarefaction_sc,
                                                  all_results_sample_cov_df_r1,
                                                  kleb_bsi_iNEXT_size_based_total)

# set order of Estimators:
kleb_bsi_combined_sample_cov_estimators$Estimator <- factor(kleb_bsi_combined_sample_cov_estimators$Estimator, 
                                                             levels = c("MLST Bayesian bootstrap. 99%",
                                                                        "MLST Bayesian bootstrap. 95%",
                                                                        "MLST Bayesian bootstrap. 90%",
                                                                        "fastBAPS L1 Bayesian bootstrap. 99%",
                                                                        "fastBAPS L1 Bayesian bootstrap. 95%",
                                                                        "fastBAPS L1 Bayesian bootstrap. 90%",
                                                                        "fastBAPS L2 Bayesian bootstrap. 99%",
                                                                        "fastBAPS L2 Bayesian bootstrap. 95%",
                                                                        "fastBAPS L2 Bayesian bootstrap. 90%",
                                                                        "fastBAPS L3 Bayesian bootstrap. 99%",
                                                                        "fastBAPS L3 Bayesian bootstrap. 95%",
                                                                        "fastBAPS L3 Bayesian bootstrap. 90%",
                                                                        "rarefaction",
                                                                        "iNEXT",
                                                                        "preseqR"))
#View(kleb_bsi_combined_sample_cov_estimators)
unique(kleb_bsi_combined_sample_cov_estimators$Estimator)
# fix colours
#add key
# make summary table 
sam_cov_colours <- c(
  `MLST Bayesian bootstrap. 99%` = "darkorange4",
  `MLST Bayesian bootstrap. 95%` = "darkorange3",
  `MLST Bayesian bootstrap. 90%` = "darkorange1",
  rarefaction = "black", 
  iNEXT = "grey40",
  preseqR = "grey70"
)


max_actual_sample_coverage <- kleb_bsi_combined_sample_cov_estimators |>
  filter(Estimator == "rarefaction") |>
  filter(sample_size == max(sample_size)) |>
  pull(mean_sc)


# filter out non-MLST metrics
plot_df <- kleb_bsi_combined_sample_cov_estimators |>
  filter(Estimator %in% names(sam_cov_colours))


sample_coverage_comparison_plot <- ggplot(plot_df, aes(x = mean_sc, y = sample_size, colour = Estimator, fill = Estimator)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(xmin = q2.5, xmax = q97.5), alpha = 0.4, colour = NA) +
  # add vertical line and label at top of line for 90%
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  # add solid line for actual collected sample
  geom_vline(xintercept = max_actual_sample_coverage, linetype = "solid") +
  # labels for vertical lines 
  annotate("label", x = 0.815, y = 10000, label = "80% sample coverage", angle = 0, vjust = -0.1, size = 4, colour = "black", fill = "white") +
  annotate("label", x = max_actual_sample_coverage, y = 10000, label = "NEKSUS sample collection", angle = 0, vjust = -0.1, size = 4, colour = "black", fill = "white") +
  scale_fill_manual(name = "Estimator", values = sam_cov_colours) +
  scale_colour_manual(name = "Estimator", values = sam_cov_colours) +
  scale_y_continuous(limits = c(0,10000)) +
  scale_x_continuous() +
  labs(x = "Sample coverage\n(proportion of population belonging to an MLST that has been observed in the sample)",
       y = "Minimum sample size",
       title = "Comparison of sample coverage estimators for Klebsiella BSI MLSTs",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(5.5, 20, 5.5, 5.5)  # space for vertical labels
  )
sample_coverage_comparison_plot
# save
ggsave("rarefaction/kleb_bsi_mlst_sample_coverage_comparison_plot.png", plot = sample_coverage_comparison_plot, width = 10, height = 6, units = "in", dpi = 300)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# sumary tables for Klebsiella
thresh <- c(0.5, 0.6, 0.7, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

# main pipeline: compute per-Estimator x threshold cells
kleb_bsi_sample_coverage_summary_table <- kleb_bsi_combined_sample_cov_estimators |>
  group_by(Estimator) |>
  group_modify(~ {
    df <- .x
    # for each threshold produce a row: threshold and formatted string
    out <- map_dfr(thresh, function(t) {
      mean_ss  <- min_sample_at_or_above(df, "mean_sc", t)
      lo_ss    <- min_sample_at_or_above(df, "q97.5", t)
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
  select(Estimator, threshold_label, cell, mean_ss, lo_ss, hi_ss) |>
  pivot_wider(
    names_from = threshold_label,
    values_from = c(cell, mean_ss, lo_ss, hi_ss),
    names_glue = "{threshold_label}_{.value}"
  )

# If you prefer a tidy wide table with only the formatted cells (no numeric subcolumns),
# extract columns named like "75%_cell", "80%_cell", ...
kleb_bsi_sample_coverage_summary_table_cells_only <- kleb_bsi_sample_coverage_summary_table |>
  select(Estimator, ends_with("_cell")) |>
  rename_with(~ str_remove(., "_cell"), ends_with("_cell"))

# View result
print(kleb_bsi_sample_coverage_summary_table)
View(kleb_bsi_sample_coverage_summary_table)
print(kleb_bsi_sample_coverage_summary_table_cells_only)
View(kleb_bsi_sample_coverage_summary_table_cells_only)

write.csv(kleb_bsi_sample_coverage_summary_table, "rarefaction/kleb_bsi_sample_coverage_summary_table.csv", row.names = FALSE)
write.csv(kleb_bsi_sample_coverage_summary_table_cells_only, "rarefaction/kleb_bsi_sample_coverage_summary_table_cells_only.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SENSITIVITY ANALYSIS ####
#for varying alpha_initial (counts between 0.001), and alpha_novel
# initial starting values are uninformative uniform prions for alpha_inital (all 1s)
# and varying starting values of alpha_novel
# exact and non-exact dirichlet
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Lightweight sensitivity analysis: only keep mass_df, not draws ----------------
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(future.apply)   # optional; will be used only if use_parallel = TRUE
#install.packages("future.apply")
# ---------- User settings -----------------------------------------------------
B <- 1000                        # posterior draws per run (reduce for speed)
f_grid <- c(0, seq(0.0001, 0.001, by = 0.00001),
            seq(0.0011, 0.01, by = 0.0001),
            seq(0.011, 1, by = 0.001))

alpha_named_values <- c(0.01, 0.1, 1)
use_parallel <- FALSE            # set TRUE to use parallel workers (requires future plan)
n_workers <- 12                  # only used when use_parallel = TRUE

# ---------- helper that runs one scenario but only returns mass_df ----------------
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

# ---------- build scenario table (same logic as before) ------------------------
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

message("Will run ", nrow(scenarios), " scenarios.")

# ---------- run scenarios (serial or parallel) --------------------------------
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

# ---------- optional plotting -------------------------------------------------
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
# * Non-exact ####
# repeat sensitivtiy analysis for non-exact. set alpha to 1. Vary alpha_novel still (skip 0.5 value).

# ---------- User settings -----------------------------------------------------
B <- 1000                        # posterior draws per run (reduce for speed)
f_grid <- c(0, seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))

alpha_named_values <- c(1)
use_parallel <- FALSE            # set TRUE to use parallel workers (requires future plan)
n_workers <- 12                  # only used when use_parallel = TRUE


# ---------- build scenario table (same logic as before) ------------------------
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

# ---------- run scenarios (serial or parallel) --------------------------------
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

# ---------- optional plotting -------------------------------------------------
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
ggsave("rarefaction/cumulative_mass_sensitivity_mass_only_non_exact.png", cumulative_mass_plot_sensitivity_non_exact, width = 16, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * GENE and PLASMID BAYESBOOT ####
# account for within-isolate co-occurrence and distribution of features
# perform isolate-level Bayesian bootstrapping
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load Required packages
#library(dplyr)
#library(purrr)
#library(tidyr)

# load data
amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
#colnames(amrfinder_metadata_updated)
# filter only ecoli BSI from amrfinder df
ecoli_bsi_amrfinder_metadata <- amrfinder_metadata_updated |>
  filter(grepl("Escherichia", kraken2_species) & (sampletype == "blood" | sampletype == "unknown" | is.na(sampletype))) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
#length(unique(ecoli_bsi_amrfinder_metadata$sample)) #1471 non-duplicates
#View(ecoli_bsi_amrfinder_metadata) 

# filter only Kleb BSI from amrfinder df
kleb_bsi_amrfinder_metadata <- amrfinder_metadata_updated |>
  filter(grepl("Klebsiella", kraken2_species) & (sampletype == "blood" | sampletype == "unknown" | is.na(sampletype))) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
#length(unique(kleb_bsi_amrfinder_metadata$sample)) #468 non-duplicates

# prepare E.coli and Klebsiella data into format where 1 isolate per row, and genes are columns
ecoli_bsi_arg_df <- ecoli_bsi_amrfinder_metadata |>
  filter(Type =="AMR" | is.na(Type)) |>
  group_by(sample, Element.symbol) |> 
  summarise(count = n(),
            presence = case_when(count >0 ~ 1,
                                 count <=0 ~ 0,
                                 TRUE ~ 0)) 
ecoli_bsi_arg_count <- ecoli_bsi_arg_df|>
  pivot_wider(id_cols = sample, names_from = Element.symbol, values_from = count, values_fill =0) |>
  select(-c(`NA`)) |>
  ungroup()
ecoli_bsi_arg_presence <- ecoli_bsi_arg_df|>
  pivot_wider(id_cols = sample, names_from = Element.symbol, values_from = presence, values_fill =0) |>
  select(-c(`NA`)) |>
  ungroup()
#View(ecoli_bsi_arg_count)
#View(ecoli_bsi_arg_presence)
length(unique(ecoli_bsi_arg_presence$sample))

# prepare Klebsiella data into format where 1 isolate per row, and genes are columns
kleb_bsi_arg_df <- kleb_bsi_amrfinder_metadata |>
  filter(Type =="AMR" | is.na(Type)) |>
  group_by(sample, Element.symbol) |> 
  summarise(count = n(),
            presence = case_when(count >0 ~ 1,
                                 count <=0 ~ 0,
                                 TRUE ~ 0)) 
kleb_bsi_arg_count <- kleb_bsi_arg_df|>
  pivot_wider(id_cols = sample, names_from = Element.symbol, values_from = count, values_fill =0) |>
  select(-c(`NA`)) |>
  ungroup()
kleb_bsi_arg_presence <- kleb_bsi_arg_df|>
  pivot_wider(id_cols = sample, names_from = Element.symbol, values_from = presence, values_fill =0) |>
  select(-c(`NA`)) |>
  ungroup()
#View(kleb_bsi_arg_count)
#View(kleb_bsi_arg_presence)

# Define isolate bayesboot function
# Bayesian bootstrap on isolates: function that returns posterior draws of prevalences and other stats
isolate_bayesboot <- function(iso_df, gene_cols, B = 2000, alpha_isolate = 1, include_counts = FALSE, copy_cols = NULL) {
  # iso_df: rows = isolates, columns: isolate_id, species (optional), gene_X (0/1), gene_Y (0/1), ...
  # You can also have integer copy counts columns or use binary presnece/absence.
  N <- nrow(iso_df)
  # Precompute indicator matrix
  Z <- as.matrix(iso_df[, gene_cols, drop = FALSE])
  r_vec <- iso_df$r
  out <- vector("list", B)
  for (b in seq_len(B)) {
    # draw gamma weights (posterior Dirichlet with shapes = alpha_isolate)
    Gdraws <- rgamma(N, shape = alpha_isolate, rate = 1)
    w <- Gdraws / sum(Gdraws)
    # marginal prevalence per gene:
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

# Set params
iso_df <- ecoli_bsi_arg_presence
iso_df <- iso_df |> dplyr::ungroup()
gene_cols <- setdiff(colnames(iso_df), "sample")
iso_df <- iso_df |>
  dplyr::mutate(r = rowSums(dplyr::across(all_of(gene_cols)) > 0))
#View(iso_df)
#length(unique(iso_df$sample)) # 1471, so all isolates represented
B <- 1000
# run
ecoli_bsi_bbs <- isolate_bayesboot(iso_df, gene_cols, B = B, alpha_isolate = 1)

# save
saveRDS(ecoli_bsi_bbs, "rarefaction/ecoli_bsi_ARG_bayesboot.rds")
#~~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_bbs <- read_rds("rarefaction/ecoli_bsi_ARG_bayesboot.rds")
#~~~~~~~~~~~~~~~#

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

# * * Prevalence historgam ####
ecoli_bsi_prevalence_summary <- extract_prevalence_df(ecoli_bsi_bbs, gene_cols)
print(ecoli_bsi_prevalence_summary)

# plot histogram of detection probabilities:
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Cumulative population mass curves with increasing frequency ####
# * * * 1. get posterior prevalence matrix (rows = draws, cols = genes/features)
# assumes ecoli_bsi_pling_bbs is a list of length B, each element has $p_g (named vector)
get_p_matrix <- function(bbs) {
  p_list <- lapply(bbs, function(x) x$p_g)
  p_mat <- do.call(rbind, p_list)  # B x G
  colnames(p_mat) <- names(p_list[[1]])
  return(as.matrix(p_mat))
}

p_mat <- get_p_matrix(ecoli_bsi_bbs)   # B x G
B <- nrow(p_mat)
G <- ncol(p_mat)

# * * 2. define frequency grid
# avoid 0 because log scale; include a tiny positive near-zero if you like
f_grid <- c(0, seq(0.00001, 0.0001, by = 0.000001), seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.0001))
f_grid <- unique(sort(f_grid))

# compute total mass per draw once 
total_mass_vec <- rowSums(p_mat)   # length B
# guard: if any total_mass == 0, set to NA to avoid division by zero later
total_mass_vec[total_mass_vec == 0] <- NA_real_

# corrected: compute normalised mass proportion for each f ----
# returns matrix: rows = draws (B), cols = length(f_grid)
prop_mat <- sapply(f_grid, function(f) {
  mask <- p_mat >= f                 # B x G logical
  masked_mass <- rowSums(p_mat * mask) # B-length numeric
  prop <- masked_mass / total_mass_vec
  prop
}) # result is B x length(f_grid) matrix (as columns)
#View(prop_mat)

# summarise across draws for each f ----
mean_mass <- colMeans(prop_mat, na.rm = TRUE)
q2.5 <- apply(prop_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
q97.5 <- apply(prop_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

ecoli_bsi_mass_df <- tibble(
  f = f_grid,
  mean_mass = mean_mass,
  q2.5 = q2.5,
  q97.5 = q97.5
) |> filter(!is.na(mean_mass))

# ---- optional: min sample sizes using f as feature frequency ----
# (these are unchanged; they depend on f not on normalization)
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

# check min and max mass
signif(mass_df$mean_mass[1], 4)  # should be ≈ 1
signif(mass_df$mean_mass[nrow(mass_df)], 4)  # should be ≈ 0

# plot (mean + 95% CI) 
cumulative_mass_plot <- ggplot(ecoli_bsi_mass_df, aes(x = f)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, fill = "seagreen3", colour = NA) +
  geom_line(aes(y = mean_mass), size = 1.05, colour = "seagreen3") +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(min(mass_df$f), 1)) +
  labs(
    x = "Feature frequency f (log scale)",
    y = "Proportion of population with features of frequency ≥ f",
    title = "Cumulative mass curve: proportion of population at frequency ≥ f",
    subtitle = "Mean (line) and 95% posterior interval (ribbon)"
  ) +
  theme_minimal(base_size = 13)
cumulative_mass_plot
#save
ggsave("rarefaction/ecoli_bsi_arg_cumulative_mass.png", cumulative_mass_plot, width =16, height = 8, dpi = 300)


# * * Cumulative sampling fraction with sample size plot ####
bootbayes_sample_coverage_plot <- ggplot(ecoli_bsi_mass_df, aes(x = mean_mass, colour = "seagreen3", fill = "seagreen3")) +
  geom_line(aes(y = min_sample_90), colour = "seagreen1") +
  geom_line(aes(y = min_sample_95), colour = "seagreen3") +
  geom_line(aes(y = min_sample_99), colour = "seagreen4") +
  geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "seagreen1", colour = NA) +
  geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "seagreen3", colour = NA) +
  geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "seagreen4", colour = NA) +
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(0,5000)) +
  labs(x = "Sample coverage",
       y = "Minimum sample size",
       title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bootbayes_sample_coverage_plot
ggsave("rarefaction/ecoli_bsi_ARG_bootbayes_sample_coverage_plot.png", plot = bootbayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Detection probability curves for a given gene ####
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
det_curve <- compute_detection_curve(bbs, "blaTEM1", m_values = 1:200)
head(det_curve)

# plot detection probability
detection_prob_plot <- ggplot(data = det_curve) +
  geom_line(aes(x = m, y = detect_prob_mean), colour = "seagreen3") +
  geom_ribbon(aes(x = m, ymin = detect_prob_lo, ymax = detect_prob_hi), alpha = 0.3, fill = "seagreen3") +
  theme_minimal() +
  labs(title = "Detection probability for blaTEM with increasing sample size",
    x = "sample size",
       y = "detection probability")
detection_prob_plot

#~~~~~~~~~~~~~~~~~~~#
# Get posterior prevalence matrix (B x G) from bbs with helper function
# Returns a matrix with columns = genes, rows = draws
get_p_matrix <- function(bbs) {
  # each bbs[[b]]$p_g is a named numeric vector (genes)
  p_list <- lapply(bbs, function(x) x$p_g)
  # convert to matrix: rows = draws, cols = genes
  p_mat <- do.call(rbind, p_list)   # B x G
  return(as.matrix(p_mat))
}

# Top-n genes by posterior mean prevalence
top_genes_by_mean <- function(bbs, n = 10) {
  p_mat <- get_p_matrix(bbs) # B x G
  gene_means <- colMeans(p_mat)
  top_genes <- names(sort(gene_means, decreasing = TRUE))[1:min(n, length(gene_means))]
  tibble(gene = top_genes, mean_prev = gene_means[top_genes])
}

# Detection curve for many genes at once
# Returns a tidy data.frame with columns: gene, m, detect_mean, detect_lo, detect_hi
compute_detection_curves_multi <- function(bbs, genes = NULL, top_n = 10, M = 1500, m_values = NULL) {
  p_mat <- get_p_matrix(bbs) # B x G
  B <- nrow(p_mat)
  gene_names <- colnames(p_mat)
  if (is.null(genes)) {
    # choose top_n
    genes <- top_genes_by_mean(bbs, n = top_n)$gene
  } else {
    # validate genes exist
    missing <- setdiff(genes, gene_names)
    if (length(missing) > 0) stop("The following genes are not present in bbs: ", paste(missing, collapse = ", "))
  }
  if (is.null(m_values)) m_values <- seq_len(M)
  
  # For each gene, compute detection probabilities across B draws for all m,
  # then summarise across draws to mean + 95% CI.
  res_list <- map(genes, function(g) {
    p_vec <- p_mat[, g]  # length B
    # matrix of detection probs: B x length(m_values)
    # detect_prob[b, j] = 1 - (1 - p_vec[b])^m_j
    mvals <- m_values
    # Vectorized: compute (1 - p)^m for each p and each m: use outer
    one_minus_p <- 1 - p_vec
    # outer(one_minus_p, mvals, `^`) gives B x length(mvals) of (1-p)^m
    det_mat <- 1 - outer(one_minus_p, mvals, `^`)  # B x J
    # summarise by column
    detect_mean <- colMeans(det_mat)
    detect_lo   <- apply(det_mat, 2, quantile, 0.025)
    detect_hi   <- apply(det_mat, 2, quantile, 0.975)
    tibble(gene = g,
           m = mvals,
           detect_mean = detect_mean,
           detect_lo = detect_lo,
           detect_hi = detect_hi)
  })
  bind_rows(res_list)
}

# Plotting helper: single panel (coloured lines) or faceted
plot_detection_curves <- function(curves_df, genes_to_plot = NULL, single_panel = TRUE) {
  df <- curves_df
  if (!is.null(genes_to_plot)) df <- df |> filter(gene %in% genes_to_plot)
  # choose colours
  ng <- length(unique(df$gene))
  pal <- scales::hue_pal()(ng)
  names(pal) <- unique(df$gene)
  
  if (single_panel) {
    ggplot(df, aes(x = m, y = detect_mean, colour = gene, fill = gene)) +
      geom_line(size = 0.9) +
      geom_ribbon(aes(ymin = detect_lo, ymax = detect_hi), alpha = 0.18, colour = NA) +
      scale_colour_manual(values = pal) +
      scale_fill_manual(values = pal) +
      scale_x_continuous(limits = c(1, max(df$m)), breaks = c(1,10,50,100,250,500,1000,1500)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Sample size (m)", y = "Detection probability", colour = "Gene", fill = "Gene",
           title = paste0("Detection probability curves (top ", ng, " genes)")) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "right")
  } else {
    ggplot(df, aes(x = m, y = detect_mean)) +
      geom_line() +
      geom_ribbon(aes(ymin = detect_lo, ymax = detect_hi), alpha = 0.2) +
      facet_wrap(~gene, ncol = 2, scales = "free_y") +
      scale_x_continuous(limits = c(1, max(df$m))) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Sample size (m)", y = "Detection probability") +
      theme_minimal(base_size = 12)
  }
}

# Plot
# assume 'bbs' is available (list length B with p_g for all genes)
# choose top 10 genes
top10 <- top_genes_by_mean(bbs, n = 10)$gene
top189 <- top_genes_by_mean(bbs, n = 189)$gene

# compute detection curves up to m = 1500
det_curves_top10 <- compute_detection_curves_multi(bbs, genes = top10, M = 500)
det_curves_top189 <- compute_detection_curves_multi(bbs, genes = top189, M = 1500)

# plot in single panel (coloured lines)
det_plot_top10 <- plot_detection_curves(det_curves_top10, single_panel = TRUE)
print(det_plot_top10)
# optionally save
ggsave("rarefaction/ecoli_bsi_detection_curve_ARG_top10.png", det_plot_top10, width = 10, height = 6, dpi = 300)


# plot all genes without a legend
df <- det_curves_top189
ng <- length(unique(df$gene))
#pal <- scales::hue_pal()(ng) # rainbow colour
pal <- scales::seq_gradient_pal(high = "#d9f0a3", low = "#00441b",space = "Lab")(seq(0, 1, length.out = ng))
#pal <- viridisLite::viridis(ng, option = "D", end = 0.9)
#pal <- viridisLite::viridis(ng, option = "E") dark blue to yellow

names(pal) <- unique(df$gene)
all_genes_det_prob_plot <- ggplot(df, aes(x = m, y = detect_mean, colour = gene, fill = gene)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = detect_lo, ymax = detect_hi), alpha = 0.02, colour = NA) +
  scale_colour_manual(values = pal ) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(limits = c(1, max(df$m))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_vline(xintercept = 1471, linetype = "dashed", colour = "black" ) +
  labs(x = "Sample size (m)", y = "Detection probability", colour = "Gene", fill = "Gene",
       title = paste0("Detection probability curves (top ", ng, " genes)")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
all_genes_det_prob_plot
# optionally save
ggsave("rarefaction/ecoli_bsi_detection_curve_ARG_top189.png", all_genes_det_prob_plot, width = 9, height = 6, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Per-isolate r distribution posterior mean ####
# combine r_tab across draws by summing probabilities for each r weighted equally
r_values <- unique(unlist(lapply(ecoli_bsi_bbs, function(x) x$r_tab$r)))
r_matrix <- matrix(0, nrow = B, ncol = length(r_values))
colnames(r_matrix) <- as.character(r_values)
for (b in seq_len(B)) {
  rt <- ecoli_bsi_bbs[[b]]$r_tab
  r_matrix[b, as.character(rt$r)] <- rt$prob
}
ecoli_bsi_r_summary <- tibble(
  r = as.integer(colnames(r_matrix)),
  mean_prob = colMeans(r_matrix),
  lo = apply(r_matrix, 2, quantile, 0.025),
  hi = apply(r_matrix, 2, quantile, 0.975)) |> 
  arrange(r) |>
  mutate(Genus = "Escherichia")
#print(ecoli_bsi_r_summary, n=25)

# plot distributon of number of genes per isolate
args_per_isolate_plot <- ggplot(data = ecoli_bsi_r_summary) +
  geom_line(aes(x = r, y = mean_prob), colour = "seagreen3") +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi), alpha = 0.3, fill = "seagreen3") +
  labs(x = "number of AMR genes per isolate",
       y = "probability",
       title = "Posterior distribution of number of AMR genes per isolate") +
  theme_minimal()
args_per_isolate_plot
# save
ggsave("rarefaction/ecoli_bsi_args_per_isolate_plot.png", args_per_isolate_plot, width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Co-occurrence: compute posterior means for pairwise co-occurrence probabilities ####
#average p_pair over draws
pair_sum <- Reduce("+", lapply(ecoli_bsi_bbs, function(x) x$p_pair)) / B
ecoli_bsi_pair_sum_mat <- pair_sum
colnames(ecoli_bsi_pair_sum_mat) <- rownames(ecoli_bsi_pair_sum_mat) <- gene_cols
ecoli_bsi_pair_sum_mat


genes <- gene_cols
B <- length(ecoli_bsi_bbs)

# compute posterior O/E draws
oe_draws <- map(ecoli_bsi_bbs, function(draw) {
  p_g <- draw$p_g
  p_pair <- draw$p_pair
  expand.grid(gene1 = genes, gene2 = genes) |>
    mutate(pA = p_g[gene1],
           pB = p_g[gene2],
           pAB = p_pair[cbind(gene1, gene2)],
           oe = pAB / (pA * pB))
  })

oe_draws_df <- bind_rows(oe_draws, .id = "draw") |>
  mutate(draw = as.integer(draw))

# summarise posterior OEs and CIs
oe_summary <- oe_draws_df |>
  group_by(gene1, gene2) |>
  summarise(oe_mean = mean(oe, na.rm = TRUE),
            oe_lo   = quantile(oe, 0.025, na.rm = TRUE),
            oe_hi   = quantile(oe, 0.975, na.rm = TRUE),
            .groups = "drop") |>
  mutate(log2_oe = log2(oe_mean),
         log2_lo = log2(oe_lo),
         log2_hi = log2(oe_hi) )

# mask where 95% CIs cross 0
oe_summary <- oe_summary |>
  mutate(significant = !(log2_lo <= 0 & log2_hi >= 0),
         log2_oe_masked = ifelse(significant, log2_oe, NA_real_))

# cluster genes for ordering
mat_for_clust <- oe_summary |>
  filter(gene1 == gene2) |>
  arrange(gene1)

gene_order <- mat_for_clust$gene1

oe_plot_df <- oe_summary |>
  mutate(gene1 = factor(gene1, levels = gene_order),
    gene2 = factor(gene2, levels = gene_order))
# remove same gene-gene
oe_plot_df <- oe_plot_df |>
  filter(gene1 != gene2)

# replace -inf with large negative
oe_plot_df <- oe_plot_df |>
   mutate(log2_oe_masked = case_when(log2_oe_masked == -Inf ~ -10,
                                     TRUE ~ log2_oe_masked))


arg_coocurr_heatmap <- ggplot(oe_plot_df, aes(x = gene2, y = gene1, fill = log2_oe_masked)) +
  geom_tile(color = "grey90") +
  scale_fill_viridis(  name = "log2(O / E)",   na.value = "white" , direction = -1) +
  labs( x = NULL, y = NULL,
    title = "Normalized ARG co-occurrence (Observed / Expected)",
    subtitle = "log2 scale; white cells = 95% CI includes independence; yellow = no co-occurrence"
    ) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
arg_coocurr_heatmap
 #save
ggsave("rarefaction/ecoli_bsi_arg_cooccurrence_heatmap.png", arg_coocurr_heatmap, width = 13, height = 12, dpi = 500)


View(oe_plot_df)
# interpretation:
#log2(O/E) = 0 -> no association
#log2(O/E) > 0 -> enriched co-occurence
#log2(O/E) < 0 -> avoidance
# white tile = insufficient evidence (95% CIs cross 0)
#log2(O/E) = 0 -> no association

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Klebsiella ARGs ####
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

# * * Prevalence summary ####
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

# * * Combined Prevalence Histogram with E. coli ####
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

# * * Cumulative population mass curves with increasing frequency ####
# * * * 1. get posterior prevalence matrix (rows = draws, cols = genes/features)
# assumes kleb_bsi_arg_bbs is a list of length B, each element has $p_g (named vector)
p_mat <- get_p_matrix(kleb_bsi_bbs)   # B x G
B <- nrow(p_mat)
G <- ncol(p_mat)

# * * 2. define frequency grid
# avoid 0 because log scale; include a tiny positive near-zero if you like
f_grid <- c(0, seq(0.00001, 0.0001, by = 0.000001), seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
f_grid <- unique(sort(f_grid))

# compute total mass per draw once 
total_mass_vec <- rowSums(p_mat)   # length B
# guard: if any total_mass == 0, set to NA to avoid division by zero later
total_mass_vec[total_mass_vec == 0] <- NA_real_

# corrected: compute normalised mass proportion for each f ----
# returns matrix: rows = draws (B), cols = length(f_grid)
prop_mat <- sapply(f_grid, function(f) {
  mask <- p_mat >= f                 # B x G logical
  masked_mass <- rowSums(p_mat * mask) # B-length numeric
  prop <- masked_mass / total_mass_vec
  prop
}) # result is B x length(f_grid) matrix (as columns)
#View(prop_mat)

# summarise across draws for each f ----
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
# (these are unchanged; they depend on f not on normalization)
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

# check min and max mass
signif(mass_df$mean_mass[1], 4)  # should be ≈ 1
signif(mass_df$mean_mass[nrow(mass_df)], 4)  # should be ≈ 0

# plot (mean + 95% CI) 
cumulative_mass_plot <- ggplot(kleb_bsi_mass_df, aes(x = f)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, fill = "darkorange", colour = NA) +
  geom_line(aes(y = mean_mass), size = 1.05, colour = "darkorange") +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(min(mass_df$f), 1)) +
  labs(
    x = "Feature frequency f (log scale)",
    y = "Proportion of population with features of frequency ≥ f",
    title = "Cumulative mass curve: proportion of population at frequency ≥ f",
    subtitle = "Mean (line) and 95% posterior interval (ribbon)"
  ) +
  theme_minimal(base_size = 13)

print(cumulative_mass_plot)
#save
ggsave("rarefaction/kleb_bsi_arg_cumulative_mass.png", cumulative_mass_plot, width =16, height = 8, dpi = 300)

# * * Combined E. coli and Klebsiella cumulative mass curve (normalised to 1) ####
# load data if required
#kleb_bsi_mass_df <- read.csv("rarefaction/kleb_bsi_arg_bayesboot_mass_df.csv")
#ecoli_bsi_mass_df <- read.csv("rarefaction/ecoli_bsi_arg_bayesboot_mass_df.csv")

# rbind
combined_bsi_arg_mass_df <- rbind(ecoli_bsi_mass_df, kleb_bsi_mass_df)
# save
write.csv(combined_bsi_arg_mass_df, "rarefaction/combined_bsi_arg_bayesboot_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~#
# read in saved df ####
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

# * * Cumulative sam_cov arg fraction with sample size plot ####
bootbayes_sample_coverage_plot <- ggplot(kleb_bsi_mass_df, aes(x = mean_mass, colour = "darkorange")) +
  geom_line(aes(y = min_sample_90), colour = "darkorange") +
  geom_line(aes(y = min_sample_95), colour = "darkorange3") +
  geom_line(aes(y = min_sample_99), colour = "darkorange4") +
  geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "darkorange1", colour = NA) +
  geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "darkorange3", colour = NA) +
  geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "darkorange4", colour = NA) +
  # add line for goods coverage estimator from 
  # add line for iNEXT and CIs coverage estimate
  # add line and error bars for preseqR coverage estimate and 
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(0,5000)) +
  labs(x = "Sample coverage",
       y = "Sample size" #,
       #title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bootbayes_sample_coverage_plot
ggsave("rarefaction/kleb_bsi_arg_bootbayes_sample_coverage_plot.png", plot = bootbayes_sample_coverage_plot, width = 6, height = 4, units = "in", dpi = 300)

# * * Combined E. coli and Klebsiella sample coverage plot
# * * Cumulative samarg fraction with sample size plot ####
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




# * * Detection curves per gene ####
top10 <- top_genes_by_mean(kleb_bsi_bbs, n = 10)$gene
top189 <- top_genes_by_mean(kleb_bsi_bbs, n = 189)$gene

# compute detection curves up to m = 1500
det_curves_top10 <- compute_detection_curves_multi(kleb_bsi_bbs, genes = top10, M = 500)
det_curves_top189 <- compute_detection_curves_multi(kleb_bsi_bbs, genes = top189, M = 1500)

# plot in single panel (coloured lines)
det_plot_top10 <- plot_detection_curves(det_curves_top10, single_panel = TRUE)
print(det_plot_top10)
# optionally save
ggsave("rarefaction/kleb_bsi_detection_curve_ARG_top10.png", det_plot_top10, width = 10, height = 6, dpi = 300)


# plot all genes without a legend
df <- det_curves_top189
ng <- length(unique(df$gene))
#pal <- scales::hue_pal()(ng) # rainbow colour
pal <- scales::seq_gradient_pal(high = "lightyellow", low = "darkorange3",space = "Lab")(seq(0, 1, length.out = ng))
#pal <- viridisLite::viridis(ng, option = "D", end = 0.9)
#pal <- viridisLite::viridis(ng, option = "E") dark blue to yellow

names(pal) <- unique(df$gene)
all_genes_det_prob_plot <- ggplot(df, aes(x = m, y = detect_mean, colour = gene, fill = gene)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = detect_lo, ymax = detect_hi), alpha = 0.02, colour = NA) +
  scale_colour_manual(values = pal ) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(limits = c(1, max(df$m))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_vline(xintercept = 468, linetype = "dashed", colour = "black" ) +
  labs(x = "Sample size (m)", y = "Detection probability", colour = "Gene", fill = "Gene",
       title = paste0("Detection probability curves (top ", ng, " genes)")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
all_genes_det_prob_plot
# optionally save
ggsave("rarefaction/kleb_bsi_detection_curve_ARG_top189.png", all_genes_det_prob_plot, width = 9, height = 6, dpi = 300)

# * * Per Isolate number of genes ####
r_values <- unique(unlist(lapply(kleb_bsi_bbs, function(x) x$r_tab$r)))
r_matrix <- matrix(0, nrow = B, ncol = length(r_values))
colnames(r_matrix) <- as.character(r_values)
for (b in seq_len(B)) {
  rt <- kleb_bsi_bbs[[b]]$r_tab
  r_matrix[b, as.character(rt$r)] <- rt$prob
}
kleb_bsi_r_summary <- tibble(
  r = as.integer(colnames(r_matrix)),
  mean_prob = colMeans(r_matrix),
  lo = apply(r_matrix, 2, quantile, 0.025),
  hi = apply(r_matrix, 2, quantile, 0.975) )|> 
  arrange(r) |>
  mutate(Genus = "Klebsiella")
#print(kleb_bsi_r_summary, n=25)

# plot distributon of number of genes per isolate
args_per_isolate_plot <- ggplot(data = kleb_bsi_r_summary) +
  geom_line(aes(x = r, y = mean_prob), colour = "darkorange") +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi), alpha = 0.3, fill = "darkorange") +
  labs(x = "number of AMR genes per isolate",
       y = "probability",
       title = "Posterior distribution of number of AMR genes per isolate") +
  theme_minimal()
args_per_isolate_plot
# save
ggsave("rarefaction/kleb_bsi_args_per_isolate_plot.png", args_per_isolate_plot, width = 6, height = 4, dpi = 300)

# * * * Combined E. coli and Klebsiella posterior distirbution of genes per isolate ####
combined_bsi_r_summary <- rbind(ecoli_bsi_r_summary, kleb_bsi_r_summary)

# plot distributon of number of genes per isolate
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
combined_args_per_isolate_plot <- ggplot(data = combined_bsi_r_summary, aes(colour = Genus, fill = Genus)) +
  geom_line(aes(x = r, y = mean_prob)) +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi), alpha = 0.3, colour = NA) +
  scale_x_continuous(limits = c(0,24)) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = "number of AMR genes per isolate",
       y = "probability",
       title = "Posterior distribution of number of AMR genes per isolate") +
  theme_minimal()
combined_args_per_isolate_plot
# save
ggsave("rarefaction/combined_bsi_args_per_isolate_plot.png", combined_args_per_isolate_plot, width = 8, height = 4, dpi = 300)


# * * Heatmap of gene co-occurrences/ avoidances ####
#average p_pair over draws
pair_sum <- Reduce("+", lapply(kleb_bsi_bbs, function(x) x$p_pair)) / B
kleb_bsi_pair_sum_mat <- pair_sum
colnames(kleb_bsi_pair_sum_mat) <- rownames(kleb_bsi_pair_sum_mat) <- gene_cols
kleb_bsi_pair_sum_mat

genes <- gene_cols
B <- length(kleb_bsi_bbs)

# compute posterior O/E draws
oe_draws <- map(kleb_bsi_bbs, function(draw) {
  p_g <- draw$p_g
  p_pair <- draw$p_pair
  expand.grid(gene1 = genes, gene2 = genes) |>
    mutate(pA = p_g[gene1],
           pB = p_g[gene2],
           pAB = p_pair[cbind(gene1, gene2)],
           oe = pAB / (pA * pB))
})

oe_draws_df <- bind_rows(oe_draws, .id = "draw") |>
  mutate(draw = as.integer(draw))

# summarise posterior OEs and CIs
oe_summary <- oe_draws_df |>
  group_by(gene1, gene2) |>
  summarise(oe_mean = mean(oe, na.rm = TRUE),
            oe_lo   = quantile(oe, 0.025, na.rm = TRUE),
            oe_hi   = quantile(oe, 0.975, na.rm = TRUE),
            .groups = "drop") |>
  mutate(log2_oe = log2(oe_mean),
         log2_lo = log2(oe_lo),
         log2_hi = log2(oe_hi) )

# mask where 95% CIs cross 0
oe_summary <- oe_summary |>
  mutate(significant = !(log2_lo <= 0 & log2_hi >= 0),
         log2_oe_masked = ifelse(significant, log2_oe, NA_real_))

# cluster genes for ordering
mat_for_clust <- oe_summary |>
  filter(gene1 == gene2) |>
  arrange(gene1)

gene_order <- mat_for_clust$gene1

oe_plot_df <- oe_summary |>
  mutate(gene1 = factor(gene1, levels = gene_order),
         gene2 = factor(gene2, levels = gene_order))
# remove same gene-gene
oe_plot_df <- oe_plot_df |>
  filter(gene1 != gene2)

# replace -inf with large negative
oe_plot_df <- oe_plot_df |>
  mutate(log2_oe_masked = case_when(log2_oe_masked == -Inf ~ -10,
                                    TRUE ~ log2_oe_masked))


arg_coocurr_heatmap <- ggplot(oe_plot_df, aes(x = gene2, y = gene1, fill = log2_oe_masked)) +
  geom_tile(color = "grey90") +
  scale_fill_viridis(  name = "log2(O / E)",   na.value = "white" , direction = -1) +
  labs( x = NULL, y = NULL,
        title = "Normalized ARG co-occurrence (Observed / Expected)",
        subtitle = "log2 scale; white cells = 95% CI includes independence; yellow = no co-occurrence"
  ) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
arg_coocurr_heatmap
#save
ggsave("rarefaction/kleb_bsi_arg_cooccurrence_heatmap.png", arg_coocurr_heatmap, width = 13, height = 12, dpi = 500)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Summary table for E. coli and Klebsiella ARGs ####

#library(dplyr)
#library(tidyr)
#library(purrr)
#library(stringr)

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
View(combined_bsi_arg_sample_coverage_summary_table)

write.csv(combined_bsi_arg_sample_coverage_summary_table, "rarefaction/combined_bsi_arg_sample_coverage_summary_table.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * E. coli plasmids ####
# prepare E. coli plasmids df
# master sample list (one row per sample)
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

# * * Run isolate-level bayesian bootstrapping for E.coli pling plasmids #### 
ecoli_bsi_pling_bbs <- isolate_bayesboot(iso_df, gene_cols, B = B, alpha_isolate = 1)

# save
saveRDS(ecoli_bsi_pling_bbs, "rarefaction/ecoli_bsi_PLING_bayesboot.rds")
#~~~~~~~~~~~~#
# read in saved data
#ecoli_bsi_pling_bbs <- read_rds("rarefaction/ecoli_bsi_PLING_bayesboot.rds")
#~~~~~~~~~~~~#

# * * Prevalence summary ####
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

# * * Cumulative population mass curves with increasing frequency ####
# * * * 1. get posterior prevalence matrix (rows = draws, cols = genes/features)
# assumes ecoli_bsi_pling_bbs is a list of length B, each element has $p_g (named vector)
get_p_matrix <- function(bbs) {
  p_list <- lapply(bbs, function(x) x$p_g)
  p_mat <- do.call(rbind, p_list)  # B x G
  colnames(p_mat) <- names(p_list[[1]])
  return(as.matrix(p_mat))
}

p_mat <- get_p_matrix(ecoli_bsi_pling_bbs)   # B x G
B <- nrow(p_mat)
G <- ncol(p_mat)
min(p_mat)

# * * 2. define frequency grid
# avoid 0 because log scale; include a tiny positive near-zero if you like
f_grid <- c(0, seq(0.00001, 0.0001, by = 0.000001), seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
f_grid <- unique(sort(f_grid))

# compute total mass per draw once 
total_mass_vec <- rowSums(p_mat)   # length B
# guard: if any total_mass == 0, set to NA to avoid division by zero later
total_mass_vec[total_mass_vec == 0] <- NA_real_

# corrected: compute normalised mass proportion for each f ----
# returns matrix: rows = draws (B), cols = length(f_grid)
prop_mat <- sapply(f_grid, function(f) {
  mask <- p_mat >= f                 # B x G logical
  masked_mass <- rowSums(p_mat * mask) # B-length numeric
  prop <- masked_mass / total_mass_vec
  prop
}) # result is B x length(f_grid) matrix (as columns)
#View(prop_mat)

# summarise across draws for each f ----
mean_mass <- colMeans(prop_mat, na.rm = TRUE)
q2.5 <- apply(prop_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
q97.5 <- apply(prop_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

ecoli_bsi_mass_df <- tibble(
  f = f_grid,
  mean_mass = mean_mass,
  q2.5 = q2.5,
  q97.5 = q97.5
) |> filter(!is.na(mean_mass))

# ---- optional: min sample sizes using f as feature frequency ----
# (these are unchanged; they depend on f not on normalization)
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


# check min and max mass
signif(mass_df$mean_mass[1], 4)  # should be ≈ 1
signif(mass_df$mean_mass[nrow(mass_df)], 4)  # should be ≈ 0

# plot (mean + 95% CI) 
cumulative_mass_plot <- ggplot(ecoli_bsi_mass_df, aes(x = f)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, fill = "seagreen3", colour = NA) +
  geom_line(aes(y = mean_mass), size = 1.05, colour = "seagreen3") +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
                ) +
  coord_cartesian(xlim = c(min(mass_df$f), 1)) +
  labs(
    x = "Feature frequency f (log scale)",
    y = "Proportion of population with features of frequency ≥ f",
    title = "Cumulative mass curve: proportion of population at frequency ≥ f",
    subtitle = "Mean (line) and 95% posterior interval (ribbon)"
  ) +
  theme_minimal(base_size = 13)
cumulative_mass_plot
#save
ggsave("rarefaction/ecoli_bsi_pling_cumulative_mass.png", cumulative_mass_plot, width =16, height = 8, dpi = 300)


# * * Cumulative sampling fraction with sample size plot ####
bootbayes_sample_coverage_plot <- ggplot(ecoli_bsi_mass_df, aes(x = mean_mass, colour = "seagreen3", fill = "seagreen3")) +
  geom_line(aes(y = min_sample_90), colour = "seagreen1") +
  geom_line(aes(y = min_sample_95), colour = "seagreen3") +
  geom_line(aes(y = min_sample_99), colour = "seagreen4") +
  geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "seagreen1", colour = NA) +
  geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "seagreen3", colour = NA) +
  geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "seagreen4", colour = NA) +
  # add line for goods coverage estimator from 
  # add line for iNEXT and CIs coverage estimate
  # add line and error bars for preseqR coverage estimate and 
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(0,10000)) +
  labs(x = "Sample coverage",
       y = "Minimum sample size",
       title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bootbayes_sample_coverage_plot
ggsave("rarefaction/ecoli_bsi_PLING_bootbayes_sample_coverage_plot.png", plot = bootbayes_sample_coverage_plot, width = 8, height = 6, units = "in", dpi = 300)




# * * Detection curves per gene ####
#View(ecoli_bsi_pling_prevalence_summary)
length(unique(ecoli_bsi_pling_prevalence_summary$gene)) #825
top10 <- top_genes_by_mean(ecoli_bsi_pling_bbs, n = 10)$gene
top825 <- top_genes_by_mean(ecoli_bsi_pling_bbs, n = 825)$gene

# compute detection curves up to m = 1500
det_curves_top10 <- compute_detection_curves_multi(ecoli_bsi_pling_bbs, genes = top10, M = 500)
det_curves_top825 <- compute_detection_curves_multi(ecoli_bsi_pling_bbs, genes = top825, M = 1500)

# plot in single panel (coloured lines)
det_plot_top10 <- plot_detection_curves(det_curves_top10, single_panel = TRUE)
print(det_plot_top10)
# optionally save
ggsave("rarefaction/ecoli_bsi_pling_detection_curve_top10.png", det_plot_top10, width = 10, height = 6, dpi = 300)


# plot all genes without a legend
df <- det_curves_top825
ng <- length(unique(df$gene))
pal <- scales::seq_gradient_pal(high = "#d9f0a3", low = "#00441b",space = "Lab")(seq(0, 1, length.out = ng))

names(pal) <- unique(df$gene)
all_genes_det_prob_plot <- ggplot(df, aes(x = m, y = detect_mean, colour = gene, fill = gene)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = detect_lo, ymax = detect_hi), alpha = 0.02, colour = NA) +
  scale_colour_manual(values = pal ) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(limits = c(1, max(df$m))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_vline(xintercept = 1471, linetype = "dashed", colour = "black" ) +
  labs(x = "Sample size (m)", y = "Detection probability", colour = "Gene", fill = "Gene",
       title = paste0("Detection probability curves (top ", ng, " genes)")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
all_genes_det_prob_plot
# optionally save
ggsave("rarefaction/ecoli_bsi_pling_detection_curve_top825.png", all_genes_det_prob_plot, width = 9, height = 6, dpi = 300)



# * * Per Isolate number of genes ####
r_values <- unique(unlist(lapply(ecoli_bsi_pling_bbs, function(x) x$r_tab$r)))
r_matrix <- matrix(0, nrow = B, ncol = length(r_values))
colnames(r_matrix) <- as.character(r_values)
for (b in seq_len(B)) {
  rt <- ecoli_bsi_pling_bbs[[b]]$r_tab
  r_matrix[b, as.character(rt$r)] <- rt$prob
}
ecoli_bsi_pling_r_summary <- tibble(
  r = as.integer(colnames(r_matrix)),
  mean_prob = colMeans(r_matrix),
  lo = apply(r_matrix, 2, quantile, 0.025),
  hi = apply(r_matrix, 2, quantile, 0.975) )|> 
  arrange(r) |>
  mutate(Genus = "Escherichia")
#print(ecoli_bsi_pling_r_summary, n=25)

# plot distributon of number of genes per isolate
args_per_isolate_plot <- ggplot(data = ecoli_bsi_pling_r_summary) +
  geom_line(aes(x = r, y = mean_prob), colour = "seagreen3") +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi), alpha = 0.3, fill = "seagreen3") +
  scale_x_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) +
  labs(x = "number of plasmids per isolate",
       y = "probability",
       title = "Posterior distribution of number of AMR genes per isolate") +
  theme_minimal()
args_per_isolate_plot
# save
ggsave("rarefaction/ecoli_bsi_pling_per_isolate_plot.png", args_per_isolate_plot, width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Heatmap of gene co-occurrences/ avoidances ####
# just with top plasmids




#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Heatmap of gene co-occurrences/ avoidances ####
install.packages("tdigest")
library(tdigest)

# computatioanlly efficient low memory use version:
# genes to include (choose top_n or provide your list)
genes_all <- gene_cols               # all genes if small
top_n <-  length(genes_all)                  # choose based on memory
# compute top_n by mean prevalence (uses get_p_matrix or your bbs)
top_genes <- top_genes_by_mean(ecoli_bsi_pling_bbs, n = top_n)$gene

genes_to_use <- top_genes  # set earlier (e.g. top 30 or 50)
G <- length(genes_to_use)
pairs_idx <- which(upper.tri(matrix(1, nrow=G, ncol=G)), arr.ind = TRUE)

welford_list <- list()
td_list <- list()

for (b in seq_len(B)) {
  draw <- ecoli_bsi_pling_bbs[[b]]
  p_g_sub <- draw$p_g[genes_to_use]
  p_pair_sub <- draw$p_pair[genes_to_use, genes_to_use, drop = FALSE]
  
  denom <- outer(p_g_sub, p_g_sub, "*") + 1e-12
  log2_oe <- log2(p_pair_sub / denom)
  
  for (r in seq_len(nrow(pairs_idx))) {
    i <- pairs_idx[r,1]; j <- pairs_idx[r,2]
    v <- log2_oe[i, j]
    if (!is.finite(v)) next
    
    key <- paste(genes_to_use[i], genes_to_use[j], sep = ":::")
    # initialize on first observation
    if (is.null(welford_list[[key]])) welford_list[[key]] <- welford_init()
    if (is.null(td_list[[key]])) td_list[[key]] <- tdigest::tdigest(vec = numeric(0), compression = 100)
    
    # update accumulators
    welford_list[[key]] <- welford_update(welford_list[[key]], v)
    tdigest::td_add(td_list[[key]], v, count = 1)
  }
  
  if (b %% 500 == 0) message("Processed draw ", b, " / ", B)
}

# finalize: build result tibble
res <- purrr::map_dfr(names(welford_list), function(k) {
  wf <- welford_finalize(welford_list[[k]])
  tdq <- quantile(td_list[[k]], c(0.025, 0.975))
  parts <- strsplit(k, ":::")[[1]]
  tibble(gene1 = parts[1], gene2 = parts[2],
         log2_mean = wf$mean, log2_sd = wf$sd,
         q2.5 = tdq[1], q97.5 = tdq[2])
})
#View(res)

res_summary <- res |>
  mutate(oe_mean = 2^ log2_mean ,
         oe_lo = 2^q2.5,
         oe_hi = 2^ q97.5) |>
  rename(log2_oe = log2_mean,
         log2_lo = q2.5, 
         log2_hi = q97.5) |>
  #masking
  mutate(significant = !(log2_lo <= 0 & log2_hi >= 0),
         log2_oe_masked = ifelse(significant, log2_oe, NA_real_))

View(res_summary)
# cluster genes for ordering
mat_for_clust <- res_summary |>
  arrange(gene1)

gene_order1 <- unique(mat_for_clust$gene1)
gene_order2 <- unique(mat_for_clust$gene2)
gene_order <- unique(c(gene_order1, gene_order2))

oe_plot_df <- res_summary |>
  mutate(gene1 = factor(gene1, levels = gene_order),
        gene2 = factor(gene2, levels = gene_order)) |>
  filter(gene1 != gene2) |>
  # replace -inf with large negative
  mutate(log2_oe_masked = case_when(log2_oe_masked == -Inf ~ -10,
                                    TRUE ~ log2_oe_masked))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# filter top hits
ecoli_bsi_pling_associations <- oe_plot_df |>
  filter(log2_oe_masked >10)
View(ecoli_bsi_pling_associations)
# save
write.csv(ecoli_bsi_pling_associations, "rarefaction/ecoli_bsi_pling_associations.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# fix table
# Create a lookup table from your summarised results
lookup <- res_summary |>
  select(gene1, gene2, log2_oe, log2_lo, log2_hi, log2_oe_masked, oe_mean)

# Build full ordered grid (exclude self-pairs if desired)
full_grid <- tidyr::expand_grid(
  gene1 = gene_order,
  gene2 = gene_order
) |>
  filter(gene1 != gene2)   # keep only off-diagonal ordered pairs

# Left-join direct matches, and also left-join the swapped match (B,A) as fallback
full_aug <- full_grid |>
  left_join(lookup, by = c("gene1", "gene2")) |>
  left_join(
    lookup |>
      rename(gene1_swap = gene1,
             gene2_swap = gene2,
             log2_oe_swap = log2_oe,
             log2_lo_swap = log2_lo,
             log2_hi_swap = log2_hi,
             log2_oe_masked_swap = log2_oe_masked,
             oe_mean_swap = oe_mean),
    by = c("gene1" = "gene2_swap", "gene2" = "gene1_swap")
  )

# Now resolve values:
# prefer the direct value; if missing, fall back to swapped; if both missing, set masked = -10
oe_plot_df_amended <- full_aug |>
  mutate(
    # numeric fields: prefer direct then swapped, keep NA if both missing
    log2_oe_final      = coalesce(log2_oe, log2_oe_swap),
    log2_lo_final      = coalesce(log2_lo, log2_lo_swap),
    log2_hi_final      = coalesce(log2_hi, log2_hi_swap),
    oe_mean_final      = coalesce(oe_mean, oe_mean_swap),
    
    # masked field: if present use it, else use swapped masked; if still missing -> -10
    log2_oe_masked_final = case_when(
      !is.na(log2_oe_masked) ~ log2_oe_masked,
      !is.na(log2_oe_masked_swap) ~ log2_oe_masked_swap,
      TRUE ~ -10
    )
  ) |>
  # bring names back to what you used previously
  transmute(
    gene1 = factor(gene1, levels = gene_order),
    gene2 = factor(gene2, levels = gene_order),
    log2_oe = log2_oe_final,
    log2_lo = log2_lo_final,
    log2_hi = log2_hi_final,
    oe_mean = oe_mean_final,
    # the user-requested masked/log2 column with -10 for missing combinations
    log2_oe_masked = log2_oe_masked_final
  ) |>
  # fix masking
  mutate(significant = !(log2_lo <= 0 & log2_hi >= 0),
         log2_oe_masked = case_when(significant == TRUE ~ log2_oe, 
                                    significant == FALSE ~ NA_real_,
                                    is.na(significant) ~ -10)) 



# Quick checks
nrow(oe_plot_df_amended)
#View(oe_plot_df_amended)
length(unique(c(as.character(oe_plot_df_amended$gene1), as.character(oe_plot_df_amended$gene2))))

# cluster similar association patterns for neat plotting:
library(stats)
# 1. make wide matrix: rows = gene, cols = gene, values = log2_oe (use masked or mean)
mat_wide <- oe_plot_df_amended |>
  select(gene1, gene2, log2_oe_masked) |>
  pivot_wider(names_from = gene2, values_from = log2_oe_masked) |>
  column_to_rownames("gene1") |>
  as.matrix()

# replace sentinel values first
mat_clean <- mat_wide
mat_clean[mat_clean == -10] <- NA

# remove genes with too few observations (optional but recommended)
min_non_na <- 5
keep <- rowSums(!is.na(mat_clean)) >= min_non_na
mat_clean <- mat_clean[keep, keep]

# correlation-based distance
cor_mat <- cor(t(mat_clean), use = "pairwise.complete.obs")

# remove any remaining NA correlations
cor_mat[is.na(cor_mat)] <- 0   # neutral similarity

dist_mat <- as.dist(1 - cor_mat)

# USE average linkage (not ward)
hc <- hclust(dist_mat, method = "average")
gene_order_clustered <- hc$labels[hc$order]
#View(oe_plot_df_amended)
# all genes you want to plot (from your oe_plot_df_amended)
all_genes <- sort(unique(c(as.character(oe_plot_df_amended$gene1), as.character(oe_plot_df_amended$gene2))))
# make full order by appending genes not included in the clustered list
missing_genes <- setdiff(all_genes, gene_order_clustered)
gene_order_full <- c(gene_order_clustered, sort(missing_genes))

# check
length(gene_order_full)  # should equal 
length(all_genes)

# apply to your plotting df
oe_plot_df_clustered_full <- oe_plot_df_amended |>
  mutate(
    gene1 = factor(as.character(gene1), levels = gene_order_full),
    gene2 = factor(as.character(gene2), levels = gene_order_full)
  )
View(oe_plot_df_clustered_full)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# filter top hits
ecoli_bsi_pling_associations <- oe_plot_df_clustered_full |>
  filter(log2_oe_masked >10 & !is.na(gene1) & !is.na(gene2)) |>
  filter(gene1 != gene2) |>
  mutate(gene_low  = pmin(as.character(gene1), as.character(gene2)),
         gene_high = pmax(as.character(gene1), as.character(gene2))) |>
  distinct(gene_low, gene_high, .keep_all = TRUE) |>
  select(-gene_low, -gene_high)
  
#View(ecoli_bsi_pling_associations)
nrow(ecoli_bsi_pling_associations) # 183
# save
write.csv(ecoli_bsi_pling_associations, "rarefaction/ecoli_bsi_pling_associations.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# plot
arg_coocurr_heatmap <- ggplot(oe_plot_df_clustered_full, aes(x = gene2, y = gene1, fill = log2_oe_masked)) +
  geom_tile(color = "grey90") +
  scale_fill_viridis(  name = "log2(O / E)",   na.value = "white" , direction = -1) +
  labs( x = NULL, y = NULL,
        title = "Normalized ARG co-occurrence (Observed / Expected)",
        subtitle = "log2 scale; white cells = 95% CI includes independence; yellow = no co-occurrence"
  ) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
arg_coocurr_heatmap
#save
ggsave("rarefaction/ecoli_bsi_pling_cooccurrence_heatmap.png", arg_coocurr_heatmap, width = 13, height = 12, dpi = 500)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Klebsiella plasmids ####
# prepare Klebsiella plasmids df
# master sample list (one row per sample)
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

# * * Run isolate-level bayesian bootstrapping for E.coli pling plasmids #### 
kleb_bsi_pling_bbs <- isolate_bayesboot(iso_df, gene_cols, B = B, alpha_isolate = 1)
# save
saveRDS(kleb_bsi_pling_bbs, "rarefaction/kleb_bsi_PLING_bayesboot.rds")
#~~~~~~~~~~~~~#
# read in saved data
kleb_bsi_pling_bbs <- read_rds("rarefaction/kleb_bsi_PLING_bayesboot.rds")
#~~~~~~~~~~~~~#

# * * Prevalence summary ####
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

# * * Combined Prevalence Histogram with E. coli ####
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



# * * Cumulative population mass curves with increasing frequency ####
# * * * 1. get posterior prevalence matrix (rows = draws, cols = genes/features)
# assumes kleb_bsi_pling_bbs is a list of length B, each element has $p_g (named vector)
p_mat <- get_p_matrix(kleb_bsi_pling_bbs)   # B x G
B <- nrow(p_mat)
G <- ncol(p_mat)

# * * 2. define frequency grid
# avoid 0 because log scale; include a tiny positive near-zero if you like
f_grid <- c(0, seq(0.00001, 0.0001, by = 0.000001), seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 1, by = 0.001))
f_grid <- unique(sort(f_grid))

# compute total mass per draw once 
total_mass_vec <- rowSums(p_mat)   # length B
# guard: if any total_mass == 0, set to NA to avoid division by zero later
total_mass_vec[total_mass_vec == 0] <- NA_real_

# corrected: compute normalised mass proportion for each f ----
# returns matrix: rows = draws (B), cols = length(f_grid)
prop_mat <- sapply(f_grid, function(f) {
  mask <- p_mat >= f                 # B x G logical
  masked_mass <- rowSums(p_mat * mask) # B-length numeric
  prop <- masked_mass / total_mass_vec
  prop
}) # result is B x length(f_grid) matrix (as columns)
#View(prop_mat)

# summarise across draws for each f ----
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
# (these are unchanged; they depend on f not on normalization)
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

# check min and max mass
signif(mass_df$mean_mass[1], 4)  # should be ≈ 1
signif(mass_df$mean_mass[nrow(mass_df)], 4)  # should be ≈ 0

# plot (mean + 95% CI) 
cumulative_mass_plot <- ggplot(kleb_bsi_mass_df, aes(x = f)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, fill = "darkorange", colour = NA) +
  geom_line(aes(y = mean_mass), size = 1.05, colour = "darkorange") +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(min(mass_df$f), 1)) +
  labs(
    x = "Feature frequency f (log scale)",
    y = "Proportion of population with features of frequency ≥ f",
    title = "Cumulative mass curve: proportion of population at frequency ≥ f",
    subtitle = "Mean (line) and 95% posterior interval (ribbon)"
  ) +
  theme_minimal(base_size = 13)

print(cumulative_mass_plot)
#save
ggsave("rarefaction/kleb_bsi_pling_cumulative_mass.png", cumulative_mass_plot, width =16, height = 8, dpi = 300)

# * * Comboined E. coli and Klebsiella cumulative mass curve (normalised to 1) ####
# load data if required
#kleb_bsi_mass_df <- read.csv("rarefaction/kleb_bsi_pling_bayesboot_mass_df.csv")
#ecoli_bsi_mass_df <- read.csv("rarefaction/ecoli_bsi_pling_bayesboot_mass_df.csv")

# rbind
combined_bsi_pling_mass_df <- rbind(ecoli_bsi_mass_df, kleb_bsi_mass_df)
# save
write.csv(combined_bsi_pling_mass_df, "rarefaction/combined_bsi_pling_bayesboot_mass_df.csv", row.names = FALSE)
#~~~~~~~~~~~~~~~~~~#
# * * read-in saved df ####
#combined_bsi_pling_mass_df <- read.csv("rarefaction/combined_bsi_pling_bayesboot_mass_df.csv")
#~~~~~~~~~~~~~~~~~~#


# plot (mean + 95% CI) 
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
combined_cumulative_mass_plot <- ggplot(combined_bsi_pling_mass_df, aes(x = f, colour = Genus, fill = Genus)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.25, colour = NA) +
  geom_line(aes(y = mean_mass), size = 1) +
  scale_x_log10(
    #breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(min(mass_df$f), 1)) +
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

# * * Cumulative sampling fraction with sample size plot ####
bootbayes_sample_coverage_plot <- ggplot(kleb_bsi_mass_df, aes(x = mean_mass, colour = "darkorange")) +
  geom_line(aes(y = min_sample_90), colour = "darkorange") +
  geom_line(aes(y = min_sample_95), colour = "darkorange3") +
  geom_line(aes(y = min_sample_99), colour = "darkorange4") +
  geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "darkorange1", colour = NA) +
  geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "darkorange3", colour = NA) +
  geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.2, fill = "darkorange4", colour = NA) +
  # add line for goods coverage estimator from 
  # add line for iNEXT and CIs coverage estimate
  # add line and error bars for preseqR coverage estimate and 
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(0,10000)) +
  labs(x = "Sample coverage",
       y = "Minimum sample size",
       title = "Minimum sample size to capture a certain proportion of the population",
  ) +
  theme_minimal()
bootbayes_sample_coverage_plot
ggsave("rarefaction/kleb_bsi_PLING_bootbayes_sample_coverage_plot.png", plot = bootbayes_sample_coverage_plot, width = 8, height = 6, units = "in", dpi = 300)

# * * Combined E. coli and Klebsiella sample coverage plot
# * * Cumulative sampling fraction with sample size plot ####
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
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



# * * Detection curves per gene ####
top10 <- top_genes_by_mean(kleb_bsi_pling_bbs, n = 10)$gene
n_pling_subcommunities <- (length(unique(kleb_bsi_amrfinder_metadata$community_subcommunity)) -1)
top417 <- top_genes_by_mean(kleb_bsi_pling_bbs, n = n_pling_subcommunities)$gene


# compute detection curves up to m = 1500
det_curves_top10 <- compute_detection_curves_multi(kleb_bsi_pling_bbs, genes = top10, M = 500)
det_curves_top417 <- compute_detection_curves_multi(kleb_bsi_pling_bbs, genes = top417, M = 1500)

# plot in single panel (coloured lines)
det_plot_top10 <- plot_detection_curves(det_curves_top10, single_panel = TRUE)
print(det_plot_top10)
# optionally save
ggsave("rarefaction/kleb_bsi_detection_curve_pling_top10.png", det_plot_top10, width = 10, height = 6, dpi = 300)


# plot all genes without a legend
df <- det_curves_top189
ng <- length(unique(df$gene))
pal <- scales::seq_gradient_pal(high = "lightyellow", low = "darkorange3",space = "Lab")(seq(0, 1, length.out = ng))

names(pal) <- unique(df$gene)
all_genes_det_prob_plot <- ggplot(df, aes(x = m, y = detect_mean, colour = gene, fill = gene)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = detect_lo, ymax = detect_hi), alpha = 0.02, colour = NA) +
  scale_colour_manual(values = pal ) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(limits = c(1, max(df$m))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  geom_vline(xintercept = 468, linetype = "dashed", colour = "black" ) +
  labs(x = "Sample size (m)", y = "Detection probability", colour = "Gene", fill = "Gene",
       title = paste0("Detection probability curves (top ", ng, " genes)")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
all_genes_det_prob_plot
# optionally save
ggsave("rarefaction/kleb_bsi_detection_curve_pling_top417.png", all_genes_det_prob_plot, width = 9, height = 6, dpi = 300)

# * * Per Isolate number of genes ####
r_values <- unique(unlist(lapply(kleb_bsi_pling_bbs, function(x) x$r_tab$r)))
r_matrix <- matrix(0, nrow = B, ncol = length(r_values))
colnames(r_matrix) <- as.character(r_values)
for (b in seq_len(B)) {
  rt <- kleb_bsi_pling_bbs[[b]]$r_tab
  r_matrix[b, as.character(rt$r)] <- rt$prob
}
kleb_bsi_pling_r_summary <- tibble(
  r = as.integer(colnames(r_matrix)),
  mean_prob = colMeans(r_matrix),
  lo = apply(r_matrix, 2, quantile, 0.025),
  hi = apply(r_matrix, 2, quantile, 0.975) )|> 
  arrange(r) |>
  mutate(Genus = "Klebsiella")
#print(kleb_bsi_pling_r_summary, n=25)

# plot distributon of number of genes per isolate
plings_per_isolate_plot <- ggplot(data = kleb_bsi_pling_r_summary) +
  geom_line(aes(x = r, y = mean_prob), colour = "darkorange") +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi), alpha = 0.3, fill = "darkorange") +
  scale_x_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) +
  labs(x = "number of AMR genes per isolate",
       y = "probability",
       title = "Posterior distribution of number of AMR genes per isolate") +
  theme_minimal()
plings_per_isolate_plot
# save
ggsave("rarefaction/kleb_bsi_plings_per_isolate_plot.png", plings_per_isolate_plot, width = 6, height = 4, dpi = 300)

# * * * Combined E. coli and Klebsiella posterior distirbution of genes per isolate ####
combined_bsi_pling_r_summary <- rbind(ecoli_bsi_pling_r_summary, kleb_bsi_pling_r_summary)

# plot distributon of number of genes per isolate
pal <- c("Escherichia"   = "seagreen3", "Klebsiella"   = "darkorange")
combined_plings_per_isolate_plot <- ggplot(data = combined_bsi_pling_r_summary, aes(colour = Genus, fill = Genus)) +
  geom_line(aes(x = r, y = mean_prob)) +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi), alpha = 0.3, colour = NA) +
  scale_x_continuous(limits = c(0,24)) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = "number of AMR genes per isolate",
       y = "probability",
       title = "Posterior distribution of number of AMR genes per isolate") +
  theme_minimal()
combined_plings_per_isolate_plot
# save
ggsave("rarefaction/combined_bsi_plings_per_isolate_plot.png", combined_plings_per_isolate_plot, width = 8, height = 4, dpi = 300)


# * * Heatmap of gene co-occurrences/ avoidances ####
#average p_pair over draws
pair_sum <- Reduce("+", lapply(kleb_bsi_pling_bbs, function(x) x$p_pair)) / B
kleb_bsi_pling_pair_sum_mat <- pair_sum
colnames(kleb_bsi_pling_pair_sum_mat) <- rownames(kleb_bsi_pling_pair_sum_mat) <- gene_cols
#kleb_bsi_pling_pair_sum_mat

genes <- gene_cols
B <- length(kleb_bsi_pling_bbs)

# compute posterior O/E draws
oe_draws <- map(kleb_bsi_pling_bbs, function(draw) {
  p_g <- draw$p_g
  p_pair <- draw$p_pair
  expand.grid(gene1 = genes, gene2 = genes) |>
    mutate(pA = p_g[gene1],
           pB = p_g[gene2],
           pAB = p_pair[cbind(gene1, gene2)],
           oe = pAB / (pA * pB))
})

oe_draws_df <- bind_rows(oe_draws, .id = "draw") |>
  mutate(draw = as.integer(draw))

# summarise posterior OEs and CIs
oe_summary <- oe_draws_df |>
  group_by(gene1, gene2) |>
  summarise(oe_mean = mean(oe, na.rm = TRUE),
            oe_lo   = quantile(oe, 0.025, na.rm = TRUE),
            oe_hi   = quantile(oe, 0.975, na.rm = TRUE),
            .groups = "drop") |>
  mutate(log2_oe = log2(oe_mean),
         log2_lo = log2(oe_lo),
         log2_hi = log2(oe_hi) )

# mask where 95% CIs cross 0
oe_summary <- oe_summary |>
  mutate(significant = !(log2_lo <= 0 & log2_hi >= 0),
         log2_oe_masked = ifelse(significant, log2_oe, NA_real_))

# cluster genes for ordering
mat_for_clust <- oe_summary |>
  filter(gene1 == gene2) |>
  arrange(gene1)

gene_order <- mat_for_clust$gene1

oe_plot_df <- oe_summary |>
  mutate(gene1 = factor(gene1, levels = gene_order),
         gene2 = factor(gene2, levels = gene_order))
# remove same gene-gene
oe_plot_df <- oe_plot_df |>
  filter(gene1 != gene2)

# replace -inf with negative
oe_plot_df <- oe_plot_df |>
  mutate(log2_oe_masked = case_when(log2_oe_masked == -Inf ~ -10,
                                    TRUE ~ log2_oe_masked))

# save dataset of top association pairs
nrow(oe_plot_df) # 173 472
kleb_bsi_pling_associations <- oe_plot_df |>
  filter(log2_oe_masked > 10) |>
  mutate(gene_low  = pmin(as.character(gene1), as.character(gene2)),
         gene_high = pmax(as.character(gene1), as.character(gene2))) |>
  distinct(gene_low, gene_high, .keep_all = TRUE) |>
  select(-gene_low, -gene_high)
nrow(kleb_bsi_pling_associations) # 116

# save
write.csv(kleb_bsi_pling_associations, "rarefaction/kleb_bsi_pling_associations.csv", row.names = FALSE)


# plot
pling_coocurr_heatmap <- ggplot(oe_plot_df, aes(x = gene2, y = gene1, fill = log2_oe_masked)) +
  geom_tile(color = "grey90") +
  scale_fill_viridis(  name = "log2(O / E)",   na.value = "white" , direction = -1) +
  labs( x = NULL, y = NULL,
        title = "Normalized pling co-occurrence (Observed / Expected)",
        subtitle = "log2 scale; white cells = 95% CI includes independence; yellow = no co-occurrence"
  ) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
pling_coocurr_heatmap
#save
ggsave("rarefaction/kleb_bsi_pling_cooccurrence_heatmap.png", pling_coocurr_heatmap, width = 13, height = 12, dpi = 500)


View(oe_plot_df)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Summary table for E. coli and Klebsiella plings ####

#library(dplyr)
#library(tidyr)
#library(purrr)
#library(stringr)

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
View(combined_bsi_pling_sample_coverage_summary_table)

write.csv(combined_bsi_pling_sample_coverage_summary_table, "rarefaction/combined_bsi_pling_sample_coverage_summary_table.csv", row.names = FALSE)

length(unique(kleb_bsi_amrfinder_metadata$community_subcommunity)) # 417 as one will be NA
length(unique(ecoli_bsi_amrfinder_metadata$community_subcommunity)) # 825


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Combined Panle plot for E. coli and Kleb posterior estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(ggplot2)
library(dplyr)
library(patchwork)

genus_colours <- c(
  "Escherichia" = "seagreen3",
  "Klebsiella"   = "darkorange"
)

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
    scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) +
    labs(x = xlab, y = ylab) +
    big_theme
  
  if (!is.null(xlim_min)) {
    p <- p + coord_cartesian(xlim = c(xlim_min, 1))
  }
  p
}

make_cov_plot <- function(df, xlab = "Sample coverage", ylab = "Sample size", y_limit = c(0, 5000)) {
  ggplot(df, aes(x = mean_mass, colour = Genus, fill = Genus)) +
    geom_line(aes(y = min_sample_90), linewidth = 0.8) +
    geom_line(aes(y = min_sample_95), linewidth = 0.8) +
    geom_line(aes(y = min_sample_99), linewidth = 0.8) +
    geom_ribbon(aes(y = min_sample_90, xmin = q2.5, xmax = q97.5), alpha = 0.15, colour = NA) +
    geom_ribbon(aes(y = min_sample_95, xmin = q2.5, xmax = q97.5), alpha = 0.425, colour = NA) +
    geom_ribbon(aes(y = min_sample_99, xmin = q2.5, xmax = q97.5), alpha = 0.70, colour = NA) +
    geom_vline(xintercept = 0.80, linetype = "dashed") +
    scale_fill_manual(values = genus_colours) +
    scale_colour_manual(values = genus_colours) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = y_limit) +
    labs(x = xlab, y = ylab) +
    big_theme
}

# make individual plots:
p1 <- make_mass_plot(
  mass_df_mlst,
  xlab = "MLST frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(mass_df_mlst$f[mass_df_mlst$f > 0], na.rm = TRUE)
)
p1

p2 <- make_cov_plot(mass_df_mlst)
p2

p3 <- make_mass_plot(
  mass_df_fastbaps_L3,
  xlab = "FastBAPS cluster frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(mass_df_fastbaps_L3$f[mass_df_fastbaps_L3$f > 0], na.rm = TRUE)
)
p3

p4 <- make_cov_plot(mass_df_fastbaps_L3)
p4

p5 <- make_mass_plot(
  combined_bsi_arg_mass_df,
  xlab = "ARG frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(combined_bsi_arg_mass_df$f[combined_bsi_arg_mass_df$f > 0], na.rm = TRUE)
)
p5

p6 <- make_cov_plot(combined_bsi_arg_mass_df)
p6

p7 <- make_mass_plot(
  combined_bsi_pling_mass_df,
  xlab = "Plasmid subcommunity frequency (f)",
  ylab = "Proportion of population\nwith feature of frequency ≥ f",
  xlim_min = min(combined_bsi_pling_mass_df$f[combined_bsi_pling_mass_df$f > 0], na.rm = TRUE)
)
p7

p8 <- make_cov_plot(combined_bsi_pling_mass_df)
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
ggsave("rarefaction/combined_4x2_mass_and_coverage_panel.png",
  plot = combined_figure, width = 10, height = 15, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Summary plot for 80% sample coverage for various metrics  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# neaten this up to pull estimates from summary tables directly.
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

# ---- plot ----
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

# ---- plot pct----
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
  labs(x = NULL, y = "Percent (%) BSIs", fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ss_80_barplot
ggsave("sample_size_80_cov_barplot_pct.png", ss_80_barplot, width = 5, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# repeat for 95% coverage
# summary plot for 80% sample coverage for various metrics 
df <- tribble(
  ~metric, ~genus, ~mean, ~lower, ~upper,
  "MLST", "Escherichia", 3364 , 3086, 3651,
  "MLST", "Klebsiella", 2302, 2138, 2721,
  "fastBAPS cluster", "Escherichia", 1760, 1575, 1995, # L3 cluster
  "fastBAPS cluster", "Klebsiella", 307, 228, 414, # L3 cluster
  "AMR gene", "Escherichia", 167, 143, 192,
  "AMR gene", "Klebsiella", 373, 305, 433,
  "Plasmid", "Escherichia", 4831, 4341 , 5349,
  "Plasmid", "Klebsiella", 2495, 2139, 2995
)

# order metrics by overall size (largest first)
metric_order <- df |>
  group_by(metric) |>
  summarise(m = max(mean)) |>
  arrange(desc(m)) |>
  pull(metric)

df$metric <- factor(df$metric, levels = metric_order)

# ---- plot ----
ss_95_barplot <- ggplot(df, aes(x = metric, y = mean, fill = genus)) +
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

ss_95_barplot
ggsave("sample_size_95_cov_barplot.png", ss_95_barplot, width = 8, height = 4, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# HIERARCHICAL BAYESIAN DIRICHLET MODEL - REGIONAL BAYESBOOT ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
install.packages("loo")

library(cmdstanr)
library(loo)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(gtools)
library(posterior)


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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Stan code: shared tau
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Stan code: region-specific taus (hierarchical over log tau_r)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit wrappers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
# * * Helpers to extract posterior draws ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# helper function to extract posterior frequencies of features (MLSTs)
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

# extract the estimates opsterior frequencies of MLSTs (features)
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Call model for E. coli ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
# this is not really meaningful as only 10 regions. Better to do k-fold cross validation (see code below)
# this is the reason fo pareto-k-diagnostics being too high
#loo_shared <- fit_u$fit$loo()
#loo_region <- fit_r$fit$loo()
#cmp <- compare_models_loo(fit_u, fit_r)
#print(cmp$comparison)

# run k-fold cross-validation
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
 

# run: 
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Compute mass curves ####
library(dplyr)
library(tibble)
library(gtools)
library(posterior)

compute_mass_curve_freq <- function(draws,
                                    f_grid,
                                    include_cols = NULL) {
  stopifnot(is.matrix(draws) || is.data.frame(draws))
  draws <- as.matrix(draws)
  
  if (is.null(include_cols)) {
    include_cols <- colnames(draws)
  }
  include_cols <- intersect(include_cols, colnames(draws))
  if (length(include_cols) == 0) stop("No matching columns found in draws.")
  
  d <- draws[, include_cols, drop = FALSE]
  
  masses <- sapply(f_grid, function(f) {
    rowSums(d * (d >= f))
  })
  
  tibble::tibble(
    f = f_grid,
    median = apply(masses, 2, stats::median),
    q2.5 = apply(masses, 2, stats::quantile, probs = 0.025, names = FALSE),
    q97.5 = apply(masses, 2, stats::quantile, probs = 0.975, names = FALSE),
    mean = apply(masses, 2, mean),
    sd = apply(masses, 2, stats::sd)) |>
    dplyr::mutate(
      min_sample_90 = ifelse(f > 0, log(1 - 0.90) / log(1 - f), NA_real_),
      min_sample_95 = ifelse(f > 0, log(1 - 0.95) / log(1 - f), NA_real_),
      min_sample_99 = ifelse(f > 0, log(1 - 0.99) / log(1 - f), NA_real_)
    )
}

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
  
  pi_draws <- extract_pi_draws(draws_df, feature_names)
  
  B <- nrow(pi_draws)
  draw_idx <- sample.int(B, size = min(n_post_draws, B), replace = FALSE)
  
  has_shared_tau <- "tau" %in% names(draws_df)
  has_region_tau <- any(grepl("^tau_r\\[", names(draws_df)))
  
  # Extract tau draws if present
  tau_shared <- if (has_shared_tau) draws_df$tau else NULL
  tau_r_mat <- NULL
  if (has_region_tau) {
    tau_r_cols <- grep("^tau_r\\[", names(draws_df), value = TRUE)
    tau_r_mat <- as.matrix(draws_df[, tau_r_cols, drop = FALSE])
    colnames(tau_r_mat) <- prep$regions
  }
  
  # Global curve from pi
  global_draws <- pi_draws[draw_idx, , drop = FALSE]
  global_curve <- compute_mass_curve_freq(global_draws, f_grid)
  
  region_curves <- vector("list", length(prep$regions))
  names(region_curves) <- prep$regions
  
  for (r in seq_along(prep$regions)) {
    region_name <- prep$regions[r]
    y_r <- as.numeric(prep$y[r, ])
    
    theta_draws <- matrix(NA_real_, nrow = length(draw_idx), ncol = length(feature_names))
    colnames(theta_draws) <- feature_names
    
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
      theta_draws[i, ] <- as.numeric(gtools::rdirichlet(1, alpha_post))
    }
    
    region_curves[[r]] <- compute_mass_curve_freq(theta_draws, f_grid) |>
      mutate(region = region_name, .before = 1)
  }
  
  list(
    global_curve = global_curve,
    region_curves = bind_rows(region_curves)
  )
}

f_grid <- c(
  0,
  seq(0.0001, 0.001, by = 0.00001),
  seq(0.0011, 0.01, by = 0.0001),
  seq(0.011, 1, by = 0.001)
)

curves <- compute_region_global_mass_curves(
  fit_obj = fit_r,   # or fit_r
  f_grid = f_grid,
  n_post_draws = 400
)

curves$global_curve
curves$region_curves

# plot global mass curve
plot <- ggplot(data = curves$global, aes(x = f, y = median)) +
  geom_line() +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.7, colour = NULL) +
  scale_x_log10() +
  theme_minimal()

plot


# plot regional mass curve
plot <- ggplot(data = curves$region_curves, aes(x = f, y = median, colour = region, fill = region)) +
  geom_line() +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.7, colour = NULL) +
  scale_x_log10() +
  theme_minimal()

plot

# plot global sample size curve


# plot regional sample size curve


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * Mass curves with global and regional masses  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
extract_pi_draws <- function(draws_df, feature_names) {
  cols <- paste0("pi[", seq_along(feature_names), "]")
  if (!all(cols %in% names(draws_df))) {
    stop("Missing pi columns in draws.")
  }
  pi <- as.matrix(draws_df[, cols, drop = FALSE])
  colnames(pi) <- feature_names
  pi
}


fit_obj <- fit_r
f_grid <- f_grid

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

# call function
f_grid <- c(seq(0.0001, 0.001, by = 0.00001),
            seq(0.0011, 0.01, by = 0.0001),
            seq(0.011, 1, by = 0.001))

curve_df <- compute_region_global_mass_curves(
  fit_obj = fit_r,     # or fit_r
  f_grid = f_grid,
  n_post_draws = 1000)
#View(curve_df)
colnames(curve_df)
 
unique(curve_df$region)
# save
write.csv(curve_df, "ecoli_bsi_mlst_mass_curves.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  * * colour-blind friendly plot ####
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


make_cov_plot_region <- function(df, xlab = "Sample coverage", ylab = "Sample size", y_limit = c(0, 5000), pal = region_pal) {
  ggplot(df |> dplyr::filter(f > 0), aes(x = global_median, colour = region, fill = region)) +
    #geom_line(aes(y = min_sample_90), linewidth = 0.8) +
    geom_line(aes(y = min_sample_95), linewidth = 0.8) +
    #geom_line(aes(y = min_sample_99), linewidth = 0.8) +
    #geom_ribbon(aes(y = min_sample_90, xmin = global_q2.5, xmax = global_q97.5), alpha = 0.15, colour = NA) +
    geom_ribbon(aes(y = min_sample_95, xmin = global_q2.5, xmax = global_q97.5), alpha = 0.2, colour = NA) +
    #geom_ribbon(aes(y = min_sample_99, xmin = global_q2.5, xmax = global_q97.5), alpha = 0.70, colour = NA) +
    geom_vline(xintercept = 0.80, linetype = "dashed") +
    scale_fill_manual(values = region_pal) +
    scale_colour_manual(values = region_pal) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = y_limit) +
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * summary tables for Ecoli BSI MLST by region ####
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Call HIERARCHICAL BAYESBOOT model for Klebsiella ####
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
# * * * Mass curves with global and regional masses  ####
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
write.csv(curve_df, "kleb_bsi_mlst_mass_curves.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# colourblind friendly plot:

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
  xlim_min = min(mass_df_mlst$f[mass_df_mlst$f > 0], na.rm = TRUE), 
  pal = region_pal
)
p1
ggsave("kleb_bsi_mlst_regional_cumulative_mass_curve_regionl_tau_bayesboot.png", p1, units = "in", width = 6, height = 4, dpi = 300)


p2 <- make_cov_plot_region(curve_df, pal = region_pal)
p2
# save
ggsave("kleb_bsi_mlst_regional_samplecoverage_vs_ss_regionl_tau_bayesboot.png", p2, units = "in", width = 6, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * summary tables for kleb BSI MLST by region ####
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * summary tables ####
# thresholds to evaluate
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * HIERARCHICAL BAYESBOOT FOR ECOLI FASTBAPS ####
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
# overall delta 54 +/- 25 so still some evidence of overall improvement of using regional tau model for ecolis


# save
saveRDS(kfold_compare, "ecoli_bsi_fastbaps_L3_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#kfold_compare <- readRDS("ecoli_bsi_fastbaps_L3_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Mass curves with global and regional masses  ####
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
write.csv(curve_df, "ecoli_bsi_fastbaps_L3_mass_curves.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# colourblind friendly plot:

# colourBrewerSet2 + Dark2 themed
region_pal <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
)
region_pal <- rev(region_pal)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * summary tables for ecoli BSI fastbaps_L3 by region ####
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * summary tables ####
# thresholds to evaluate
#thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

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
#View(ecoli_bsi_fastbaps_L3_regional_bhm_summary_cells_only)

write.csv(ecoli_bsi_fastbaps_L3_regional_bhm_summary_cells_only, "rarefaction/ecoli_bsi_fastbaps_L3_regional_bhm_summary_cells_only.csv", row.names = FALSE)
write.csv(ecoli_bsi_fastbaps_L3_regional_bhm_summary, "rarefaction/ecoli_bsi_fastbaps_L3_regional_bhm_summary.csv", row.names = FALSE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * HIERARCHICAL BAYESBOOT FOR KLEB FASTBAPS ####
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
saveRDS(kfold_compare, "kleb_bsi_fastbaps_L3_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
# read back in
#kfold_compare <- readRDS("kleb_bsi_fastbaps_L3_kfold_compare.rds")
#~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * Mass curves with global and regional masses  ####
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# colourblind friendly plot:

# colourBrewerSet2 + Dark2 themed
region_pal <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
)
region_pal <- rev(region_pal)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * summary tables for kleb BSI fastbaps_L3 by region ####
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * * summary tables ####
# thresholds to evaluate
#thresh <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * HIERARCHICAL BAYESBOOT FOR ECOLI ARG - BETA_BINOMIAL ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# write feature-level bayesboot functions, using Beta-Binomial hierarchical distributions
# prepare data
library(cmdstanr)
library(loo)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(gtools)
library(posterior)

# load and prepare data
amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
# fix regions
amrfinder_metadata_updated <- amrfinder_metadata_updated |>
  mutate(region = case_when(
    region == "North East & Yorkshire A" ~ "North East A",
    region == "North East & Yorkshire B" ~ "North East B",
    sample == "301476L" ~ "South West", # unknown in amrfinder metadata but annotated as north bristol hospital in samples metadata table, based on NB isolateID
    sample == "111369" ~ "South West", # based on isolateid starting with NB
    TRUE ~ region)) 
# resave
#write.csv(amrfinder_metadata_updated, "amrfinder_metadata_with_NAs_updated.csv", row.names = FALSE)


# filter only ecoli BSI from amrfinder df
ecoli_bsi_amrfinder_metadata <- amrfinder_metadata_updated |>
  filter(grepl("Escherichia", kraken2_species) & (sampletype == "blood" | sampletype == "unknown" | is.na(sampletype))) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
#length(unique(ecoli_bsi_amrfinder_metadata$sample)) #1471 non-duplicates
# prepare E.coli and Klebsiella data into format where 1 isolate per row, and genes are columns
ecoli_bsi_arg_presence <- ecoli_bsi_amrfinder_metadata |>
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


# filter only Kleb BSI from amrfinder df
kleb_bsi_amrfinder_metadata <- amrfinder_metadata_updated |>
  filter(grepl("Klebsiella", kraken2_species) & (sampletype == "blood" | sampletype == "unknown" | is.na(sampletype))) |>
  filter(duplicate_assembly_qc_pass != "duplicate")
#length(unique(kleb_bsi_amrfinder_metadata$sample)) #468 non-duplicates
# prepare Klebsiella data into format where 1 isolate per row, and genes are columns
kleb_bsi_arg_presence <- kleb_bsi_amrfinder_metadata |>
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


gene_cols <- setdiff(names(ecoli_bsi_arg_presence), c("sample", "region"))
prep_arg <- prep_subiso_data(
  df = ecoli_bsi_arg_presence,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols
)

fit_u_arg <- fit_shared_tau_subiso(prep_arg, iter_warmup = 1000, iter_sampling = 1000, chains = 4)
fit_r_arg <- fit_region_tau_subiso(prep_arg, iter_warmup = 1000, iter_sampling = 1000, chains = 4)

# save
saveRDS(prep_arg, "prep_arg_subiso.rds")
saveRDS(fit_u_arg, "fit_u_arg_subiso.rds")
saveRDS(fit_r_arg, "fit_r_arg_subiso.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * K-fold cross-validation ####
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

kfold_compare_subiso$summary
kfold_compare_subiso$per_region
kfold_compare_subiso$fold_id
# check model diagnostics as a few unterminated chains!


# check which is better

# cal cumulative mass grid
f_grid <- c(  seq(0.0001, 0.001, by = 0.00001),   seq(0.0011, 0.01, by = 0.0001),   seq(0.011, 1, by = 0.001))

curve_arg <- compute_subiso_mass_curves(
  fit_obj = fit_r_arg,
  f_grid = f_grid,
  n_post_draws = 1000)

#View(curve_arg)
colnames(curve_arg)
table(curve_arg$region)



write.csv(curve_arg, "ecoli_arg_subiso_mass_curves.csv", row.names = FALSE)

# plot mass curve
# colourBrewerSet2 + Dark2 themed
region_pal <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02")
region_pal <- rev(region_pal)

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


# save summary table


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * HIERARCHICAL BAYESBOOT FOR KLEB ARG ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gene_cols <- setdiff(names(kleb_bsi_arg_presence), c("sample", "region"))
prep_arg <- prep_subiso_data(
  df = kleb_bsi_arg_presence,
  region_col = "region",
  isolate_id_col = "sample",
  feature_cols = gene_cols
)

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

kfold_compare_subiso$summary
kfold_compare_subiso$per_region
kfold_compare_subiso$fold_id
# check model diagnostics as a few unterminated chains!

# check which is better

# cal cumulative mass grid
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

# plot mass curve
# colourBrewerSet2 + Dark2 themed
region_pal <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02")
region_pal <- rev(region_pal)

# make individual plots:
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * HIERARCHICAL BAYESBOOT FOR ECOLI PLASMIDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare E. coli plasmids df
# master sample list (one row per sample)
all_samples <- ecoli_bsi_samples_metadata |> distinct(sequencing_id)   

# detected pairs from AMRFinder (1 if detected)
detected <- ecoli_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)

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

kfold_compare_subiso$summary
kfold_compare_subiso$per_region
kfold_compare_subiso$fold_id
# check model diagnostics as a few unterminated chains!

# check which is better

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

# plot mass curve
# colourBrewerSet2 + Dark2 themed
region_pal <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02")
region_pal <- rev(region_pal)

# make individual plots:
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * HIERARCHICAL BAYESBOOT FOR KLEB PLASMIDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# prepare Kleb plasmids df
# master sample list (one row per sample)
all_samples <- kleb_bsi_samples_metadata |> distinct(sequencing_id)   

# detected pairs from AMRFinder (1 if detected)
detected <- kleb_bsi_amrfinder_metadata |>
  filter(!is.na(community_subcommunity)) |>
  distinct(sample, community_subcommunity) |>
  mutate(presence = 1L)

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

kfold_compare_subiso$summary
kfold_compare_subiso$per_region
kfold_compare_subiso$fold_id
# check model diagnostics as a few unterminated chains!

# check which is better

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

# plot mass curve
# colourBrewerSet2 + Dark2 themed
region_pal <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02")
region_pal <- rev(region_pal)

# make individual plots:
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * * PANEL PLOT OF ALL BY-REGION HIERARCHICAL MODELS ####
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

# plot functions defined above
# mass plots for ecoli and kleb

# colourBrewerSet2 + Dark2 themed
region_pal <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
)
region_pal <- rev(region_pal)

# make individual plots:
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# sample size plots for 95% confidence level for ecoli and kleb
# make individual plots:
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















#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PHYLOSAMP ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PHYLOSAMP SAMPLE SIZE ESTIMATION ####
#install.packages("phylosamp")
library(phylosamp)
#?phylosamp

# work through Vignetes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * V1 Estimating bias in observed variant prevalence ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# - not that relevant?
# higher coefficient of detection not really relevant to BSIs?
# estimate ratio of VOC vs rest of pop
coeff_of_detection <- vartrack_cod_ratio(phi_v1=0.975, phi_v2=0.95, gamma_v1=0.8, gamma_v2=0.6)
print(coeff_of_detection) # 1.356 # as higher test sensitivity (phi) and higher sequencing success, gamma, for V1

# test how observed proportion expected to change
v1_obs_freq <- varfreq_obs_freq(p_v1=0.2, c_ratio=coeff_of_detection)
v1_obs_freq # 0.25484

# so easier to detect prevalence, as infalted apparent prevalence, but harder to estimate prevalence with precision

# estimate multiplicative bias in observed variant prevalence 
bias <- varfreq_expected_mbias(p_v1=0.2, c_ratio=coeff_of_detection)
bias # 1.2745 = observed prev / actual prev

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * V2 Presence/absence or prevalence detection (Variant monitoring) in cross-sectional sample ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# required knowledge of the coefficient of detection of VOC vs rest of population
# ratio of any variable not provided assumed to be 1

# detect presence vs absence
ss_detect <- vartrack_samplesize_detect(p_v1=0.02, # minimum variant prevalence to detect (e.g.: 2%)
                                        prob=0.95,  # desired probability of detection
                                        omega=0.82, # sequencing success rate # parameterised from the NEKSUS study
                                        c_ratio=1,  # coefficient of detection ratio 
                                        sampling_freq="xsect")
ss_detect #180
# this is a simple power calculation
# apply bayesian bootstrapped posterior distribution of MLST/ lineage frequency to obtain what sample coverage would be achieved by min 2% frequency threshold

# detect prevalence with a certain precision
ss_prev <- vartrack_samplesize_prev(p_v1=0.02, # minimum variant prevalence to track (e.g.: 2%)
                                    prob=0.95,  # desired probability of detection
                                    precision = 0.25,
                                    omega=0.82, # sequencing success rate # parameterised from the NEKSUS study
                                    c_ratio=1,  # coefficient of detection ratio 
                                    sampling_freq="xsect")
ss_prev # 3672.8 <- much higher!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * V3 (post-hoc) Confidence of detection given fixed sample in cross-sectional sample ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

prob_detect <- vartrack_prob_detect(p_v1=0.005, # min prevalence to detect
                                    n= 1500, # actual number of samples selected for sequencing
                                    omega=0.82,  #sequencing success
                                    c_ratio=1, 
                                    sampling_freq="xsect")
prob_detect # 0.997 prob of detecting something at 0.5% prevalence

# precision of estimation
prob_prev <- vartrack_prob_prev(p_v1=0.1, 
                                n=468, 
                                omega=0.82, 
                                precision=0.25,
                                c_ratio=1, 
                                sampling_freq="xsect")

prob_prev # 0.897 so 90% confidence that any estimate of prevalence of a variant of at lest 10% prevalence is within 25% of true value, with this sample for Klebsiella
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * V4 Sample Sizes from periodic sampling ####
# detecting presence/absence and prevalence estimation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# calculate sample size needed to detect a particular variant in the population. Provide either:
# the number of days after introduction within which it should be detected
# OR  by the time a variant reaches a specific frequency
preiodic_ss_detect <- vartrack_samplesize_detect(prob=0.95, # desired probability of detection
                           p_v1=0.01, # desired prevalence to detect
                           p0_v1=1/10000,  # initial variant prevalence
                           r_v1=0.01, # estimated logistic growth rate per day
                           omega=0.82, # sequencing success
                           c_ratio=1, 
                           sampling_freq="cont")

preiodic_ss_detect # 34 samples per day to ensure detection of a variant by the time it reaches 1%


# instead, to detect a variant within the first month of its introduction into the population:
periodic_ss_detect_30d <- vartrack_samplesize_detect(prob=0.95,
                                                     t = 30, 
                                                     p0_v1=3/10000, 
                                                     r_v1=0.01, 
                                                     omega=0.82,
                                                     c_ratio=1, 
                                                     sampling_freq="cont")

periodic_ss_detect_30d # 64 samples per day, assuming 82% sequencing success,
#to detect a variant within 30 days of introduction, assuming it grows at 10% per day 
#and starts at 3 in 10,000 prevalence

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * V5 (Post-hoc) Probability of detection in periodic sampling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# probability of detecting a variant by the time it reaches 1% 
periodic_prob_detect <- vartrack_prob_detect(n=20, # perday samples
                                             p_v1=0.01, 
                                             omega=0.82, 
                                            p0_v1=3/10000, # initial prevalence
                                            r_v1=0.1, # exponential growth rate
                                             c_ratio=1, 
                                             sampling_freq="cont")
periodic_prob_detect # 0.825

# probability of detecting variant within first month of introductino into the population
periodic_prob_detect_30d <- vartrack_prob_detect(n=20,
                                                 t=30, 
                                                 omega=0.82, 
                                                 p0_v1=3/10000,
                                                 r_v1=0.1,
                                                 c_ratio=1, 
                                                 sampling_freq="cont")
periodic_prob_detect_30d # 61% probability of detecting a variant within 30d given 20 samples collected per day

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Sensitivity analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sample size to detect a certain prevalence after emergence with periodic sampling, with different growth rates and initial prevalences
# load packages
library(RColorBrewer)
library(scales)

#make a per year to per day conversion function for exponential growth rate
convert_to_daily_rate <- function(input_rate, time_in_days = 365) {
  rate_daily <- ((1 + input_rate)^ (1/time_in_days)) -1
  return(rate_daily)
}

# calculate reasonable exponential growth rates:
# from Gladstone et al. the following exponential growth rates (per year) on emergence of E. coli CC131 
# clade A in 2005 Ne = 3107, exp phase within 9 months (0.7 years)
# clade B in 1990s Ne = 650 # exp phase in 4 years
# clade C1 in 2003 Ne = 347 # exp phase in 3 years
# clade C2 in 2005 Ne = 131 per year # exp phase in 3 years

clade_a_daily_rate <- convert_to_daily_rate(3107) # 0.0222
clade_b_daily_rate <- convert_to_daily_rate(650)
clade_c1_daily_rate <- convert_to_daily_rate(131) # 0.01346
clade_c2_daily_rate <- convert_to_daily_rate(347)

# explore daily rates between 0.01 and 0.03
# input vectors
p_v1 <- c(seq(0.0001, 0.001, by = 0.00001), seq(0.0011, 0.01, by = 0.0001), seq(0.011, 0.2, by = 0.001))
prob <- 0.95
p0_v1 <- c(0.0005, 0.001, 0.002, 0.003) # initial prevalence
r_v1 <- c(0.01, 0.02, 0.03) # logistic growth rate per day
omega <- c(1)
c_ratio <- 1

# make df with all of these combos
df_grid <- tidyr::crossing(
  p_v1 = p_v1,
  prob = prob,
  p0_v1 = p0_v1,
  r_v1 = r_v1,
  omega = omega,
  c_ratio = c_ratio
) |>
  arrange(p0_v1, r_v1, omega, p_v1, prob, c_ratio)
View(df_grid)


# helper that calls the function and coerces to a single numeric
call_ss_single <- function(prob, p_v1_target, p0_v1, r_v1, omega, c_ratio, sampling_freq = "cont") {
  # call the underlying fn
  res <- vartrack_samplesize_detect(
    prob   = prob,
    p_v1   = p_v1_target,
    p0_v1  = p0_v1,
    r_v1   = r_v1,
    omega  = omega,
    c_ratio = c_ratio,
    sampling_freq = sampling_freq
  )
  
  # If result is a list, numeric vector, or something else:
  if (is.list(res) && length(res) == 1) {
    res_val <- res[[1]]
  } else {
    res_val <- res
  }
  
  # If it returned a vector (length > 1), try to pick the best matching element
  if (length(res_val) > 1) {
    # If result has names that correspond to prevalences, try to use them
    nms <- names(res_val)
    
    if (!is.null(nms)) {
      # try to find exact or closest name match to p_v1_target
      # numeric names sometimes appear as "0.001" etc.
      nms_num <- suppressWarnings(as.numeric(nms))
      if (!all(is.na(nms_num))) {
        # pick index of closest prevalence
        idx <- which.min(abs(nms_num - p_v1_target))
        return(as.numeric(res_val[idx]))
      }
    }
    
    # otherwise just return the first element (fallback)
    return(as.numeric(res_val[1]))
  }
  
  # otherwise length 1 -> coerce to numeric and return
  return(as.numeric(res_val))
}


# Now build df_grid (as before) and call the helper rowwise with pmap_dbl
df_grid <- tidyr::crossing(
  p_v1 = p_v1,
  p0_v1 = p0_v1,
  r_v1 = r_v1,
  omega = omega
) |>
  arrange(p0_v1, r_v1, omega, p_v1) |>
  filter(p_v1 > p0_v1) # filter out cases where target prevalence is less than initial prevalence

df_grid <- df_grid |>
  mutate(
    sample_size_daily = pmap_dbl(
      list(p_v1, p0_v1, r_v1, omega),
      function(p_v1_val, p0_v1_val, r_v1_val, omega_val) {
        call_ss_single(
          prob = prob,
          p_v1_target = p_v1_val,
          p0_v1 = p0_v1_val,
          r_v1 = r_v1_val,
          omega = omega_val,
          c_ratio = c_ratio,
          sampling_freq = "cont"
        )
      }
    ),
    sample_size_28d   = sample_size_daily * 28,
    sample_size_30.4d = sample_size_daily * 30.4
  )

# filter
plot_df <- df_grid |> filter(omega == 1 ) 

# prepare aesthetics:
# colourblind-friendly palette for r_v1 (3 values)
base_cols <- brewer.pal(n = length(unique(r_v1)), name = "Dark2") # colourblind friendly
names(base_cols) <- as.character(sort(unique(r_v1)))

alpha_map <- function(p0vals){
  # map numeric p0 to alpha between 1 (largest) and 0.4 (smallest)
  res <- rescale(p0vals, to = c(0.2, 1), from = range(plot_df$p0_v1))
  return(res)
}
plot_df <- plot_df |> mutate(alpha_level = alpha_map(p0_v1))
table(plot_df$p0_v1, plot_df$alpha_level)

# create grouping label so each combination draws one line
plot_df <- plot_df |> mutate(combo = interaction(r_v1, p0_v1, sep = " | p0="))
#View(plot_df)

# draw the plot
g <- ggplot(plot_df, aes(x = p_v1, y = sample_size_28d,
                         group = combo,
                         colour = factor(r_v1),
                         alpha = alpha_level)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(limits = c(0, 10000)) +
  scale_x_continuous(trans = "log10") + # variant prevalences often plotted on log scale
  scale_colour_manual(
    name = "r_v1 (growth rate per day)",
    values = base_cols
  ) +
  scale_alpha_continuous(
    name = "initial prevalence (p0_v1)",
    range = c(0.4, 1),
    breaks = alpha_map(sort(unique(p0_v1), decreasing = FALSE)),
    labels = paste0("p0 = ", format(sort(unique(p0_v1)), scientific = TRUE))
  ) +
  guides(
    colour = guide_legend(order = 1),
    alpha = guide_legend(order = 2, override.aes = list(size = 1.5))
  ) +
  labs(
    title = "Required samples (28 days) vs Variant prevalence (p_v1)",
    subtitle = "Only omega == 1; lines show each (r_v1, p0_v1) combination",
    x = "Variant prevalence to detect (p_v1)",
    y = "Samples required over 28 days",
    caption = "Colour = growth rate (r_v1). Lighter lines = smaller initial p0_v1"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

print(g)

one_percent <- plot_df |>
  filter(p_v1 == 0.01) |>
  pivot_wider(names_from = p0_v1, values_from = sample_size_28d, id_cols = r_v1)
View(one_percent)

# 4-weekly samples required to detect a variant just above 0.1% frequency (0.001), with low growth rate 
max(one_percent[one_percent$r_v1 == 0.01 & one_percent$p0_v1 == 0.0005,]$sample_size_28d) # 87
max(one_percent[one_percent$r_v1 == 0.01 & one_percent$p0_v1 == 0.001,]$sample_size_28d) # 91
max(one_percent[one_percent$r_v1 == 0.01 & one_percent$p0_v1 == 0.002,]$sample_size_28d) # 104
max(one_percent[one_percent$r_v1 == 0.01 & one_percent$p0_v1 == 0.003,]$sample_size_28d) # 118

# 4-weekly samples required to detect a variant just above 0.1% frequency, 
max(one_percent[one_percent$r_v1 == 0.02 & one_percent$p0_v1 == 0.0005,]$sample_size_28d) # 173
max(one_percent[one_percent$r_v1 == 0.02 & one_percent$p0_v1 == 0.001,]$sample_size_28d) # 184
max(one_percent[one_percent$r_v1 == 0.02 & one_percent$p0_v1 == 0.002,]$sample_size_28d) # 208
max(one_percent[one_percent$r_v1 == 0.02 & one_percent$p0_v1 == 0.003,]$sample_size_28d) # 235

# 4-weekly samples required to detect a variant just above 0.1% frequency, 
max(one_percent[one_percent$r_v1 == 0.03 & one_percent$p0_v1 == 0.0005,]$sample_size_28d) # 257
max(one_percent[one_percent$r_v1 == 0.03 & one_percent$p0_v1 == 0.001,]$sample_size_28d) # 269
max(one_percent[one_percent$r_v1 == 0.03 & one_percent$p0_v1 == 0.002,]$sample_size_28d) # 312
max(one_percent[one_percent$r_v1 == 0.03 & one_percent$p0_v1 == 0.003,]$sample_size_28d) # 345
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make summary table of sample sizes for varios frequencies fo detection


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Now varying by time to detection ####
# instead, to detect a variant within the first month of its introduction into the population:
periodic_ss_detect_30d <- vartrack_samplesize_detect(prob=0.95,
                                                     t = 30, 
                                                     p0_v1=3/10000, 
                                                     r_v1=0.01, 
                                                     omega=0.82,
                                                     c_ratio=1, 
                                                     sampling_freq="cont")

periodic_ss_detect_30d # 64 samples per day, assuming 82% sequencing success,
#to detect a variant within 30 days of introduction, assuming it grows at 10% per day 

# input vectors
t <- c(30, 61, 91, 122, 152, 183, 213, 244, 274, 305, 335, 365)
prob <- 0.95
p0_v1 <- c(0.0005, 0.001, 0.002, 0.003) # initial prevalence
r_v1 <- c(0.01, 0.02, 0.03) # logistic growth rate per day
omega <- c(1)
c_ratio <- 1

# make df with all of these combos
df_grid <- tidyr::crossing(
  t = t,
  prob = prob,
  p0_v1 = p0_v1,
  r_v1 = r_v1,
  omega = omega,
  c_ratio = c_ratio
) |>
  arrange(p0_v1, r_v1, omega, t, prob, c_ratio)
View(df_grid)

#~~~~~~~~~#
# helper that calls the function and coerces to a single numeric
call_ss_single_time <- function(prob, t_target, p0_v1, r_v1, omega, c_ratio, sampling_freq = "cont") {
  # call the underlying fn
  res <- vartrack_samplesize_detect(
    prob   = prob,
    t   = t_target,
    p0_v1  = p0_v1,
    r_v1   = r_v1,
    omega  = omega,
    c_ratio = c_ratio,
    sampling_freq = sampling_freq
  )
  
  # If result is a list, numeric vector, or something else:
  if (is.list(res) && length(res) == 1) {
    res_val <- res[[1]]
  } else {
    res_val <- res
  }
  
  # If it returned a vector (length > 1), try to pick the best matching element
  if (length(res_val) > 1) {
    # If result has names that correspond to prevalences, try to use them
    nms <- names(res_val)
    
    if (!is.null(nms)) {
      # try to find exact or closest name match to p_v1_target
      # numeric names sometimes appear as "0.001" etc.
      nms_num <- suppressWarnings(as.numeric(nms))
      if (!all(is.na(nms_num))) {
        # pick index of closest prevalence
        idx <- which.min(abs(nms_num - p_v1_target))
        return(as.numeric(res_val[idx]))
      }
    }
    
    # otherwise just return the first element (fallback)
    return(as.numeric(res_val[1]))
  }
  
  # otherwise length 1 -> coerce to numeric and return
  return(as.numeric(res_val))
}


# Now build df_grid (as before) and call the helper rowwise with pmap_dbl
df_grid <- tidyr::crossing(
  t = t,
  p0_v1 = p0_v1,
  r_v1 = r_v1,
  omega = omega
) |>
  arrange(p0_v1, r_v1, omega, t) 

df_grid <- df_grid |>
  mutate(
    sample_size_daily = pmap_dbl(
      list(t, p0_v1, r_v1, omega),
      function(t_val, p0_v1_val, r_v1_val, omega_val) {
        call_ss_single_time(
          prob = prob,
          t = t_val,
          p0_v1 = p0_v1_val,
          r_v1 = r_v1_val,
          omega = omega_val,
          c_ratio = c_ratio,
          sampling_freq = "cont"
        )
      }
    ),
    sample_size_28d   = sample_size_daily * 28,
    sample_size_30.4d = sample_size_daily * 30.4
  )

# filter
plot_df_time <- df_grid |> filter(omega == 1 ) 

# prepare aesthetics:
# colourblind-friendly palette for r_v1 (3 values)
base_cols <- brewer.pal(n = length(unique(r_v1)), name = "Dark2") # colourblind friendly
names(base_cols) <- as.character(sort(unique(r_v1)))

plot_df_time <- plot_df_time |> mutate(alpha_level = alpha_map(p0_v1))
table(plot_df_time$p0_v1, plot_df_time$alpha_level)

# create grouping label so each combination draws one line
plot_df_time <- plot_df_time |> 
  mutate(combo = interaction(r_v1, p0_v1, sep = " | p0="),
         time_months = round(t/30.4))
View(plot_df_time)

# draw the plot
g <- ggplot(plot_df_time, aes(x = time_months, y = sample_size_28d,
                         group = combo,
                         colour = factor(r_v1),
                         alpha = alpha_level)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(limits = c(0, 5000)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12)) + # variant prevalences often plotted on log scale
  scale_colour_manual(
    name = "r_v1 (growth rate per day)",
    values = base_cols
  ) +
  scale_alpha_continuous(
    name = "initial prevalence (p0_v1)",
    range = c(0.4, 1),
    breaks = alpha_map(sort(unique(p0_v1), decreasing = FALSE)),
    labels = paste0("p0 = ", format(sort(unique(p0_v1)), scientific = TRUE))
  ) +
  guides(
    colour = guide_legend(order = 1),
    alpha = guide_legend(order = 2, override.aes = list(size = 1.5))
  ) +
  labs(
    title = "Required samples (28 days) vs Variant prevalence (p_v1)",
    subtitle = "Only omega == 1; lines show each (r_v1, p0_v1) combination",
    x = "Variant prevalence to detect (p_v1)",
    y = "Samples required over 28 days",
    caption = "Colour = growth rate (r_v1). Lighter lines = smaller initial p0_v1"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

print(g)

six_months <- plot_df_time |>
  filter(time_months == 6) |>
  pivot_wider(names_from = p0_v1, values_from = sample_size_28d, id_cols = r_v1)
View(six_months)

# generate estimate for within 6-month detection
# 4-weekly samples required to detect a variant just above 0.1% frequency (0.001), with low growth rate 
max(one_percent[one_percent$r_v1 == 0.01 & one_percent$p0_v1 == 0.0005,]$sample_size_28d) # 87
max(one_percent[one_percent$r_v1 == 0.01 & one_percent$p0_v1 == 0.001,]$sample_size_28d) # 91
max(one_percent[one_percent$r_v1 == 0.01 & one_percent$p0_v1 == 0.002,]$sample_size_28d) # 104
max(one_percent[one_percent$r_v1 == 0.01 & one_percent$p0_v1 == 0.003,]$sample_size_28d) # 118

# 4-weekly samples required to detect a variant just above 0.1% frequency, 
max(one_percent[one_percent$r_v1 == 0.02 & one_percent$p0_v1 == 0.0005,]$sample_size_28d) # 173
max(one_percent[one_percent$r_v1 == 0.02 & one_percent$p0_v1 == 0.001,]$sample_size_28d) # 184
max(one_percent[one_percent$r_v1 == 0.02 & one_percent$p0_v1 == 0.002,]$sample_size_28d) # 208
max(one_percent[one_percent$r_v1 == 0.02 & one_percent$p0_v1 == 0.003,]$sample_size_28d) # 235

# 4-weekly samples required to detect a variant just above 0.1% frequency, 
max(one_percent[one_percent$r_v1 == 0.03 & one_percent$p0_v1 == 0.0005,]$sample_size_28d) # 257
max(one_percent[one_percent$r_v1 == 0.03 & one_percent$p0_v1 == 0.001,]$sample_size_28d) # 269
max(one_percent[one_percent$r_v1 == 0.03 & one_percent$p0_v1 == 0.002,]$sample_size_28d) # 312
max(one_percent[one_percent$r_v1 == 0.03 & one_percent$p0_v1 == 0.003,]$sample_size_28d) # 345

