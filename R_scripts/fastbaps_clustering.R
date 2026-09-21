#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# fastbaps Clustering ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Install and load other packages
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
# Install BiocManager (if not already installed)
install.packages("BiocManager")
# Install ggtree via Bioconductor
BiocManager::install("ggtree")
library(ggtree)
library(ape)
library(phytools)


# install and load fastbaps
install.packages("devtools")
devtools::install_github("gtonkinhill/fastbaps")
library(fastbaps)

# set workign directory
setwd("~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Make input files for Panaroo multiple sequence alignment ####

# load data if not already loaded
amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")

# all ecolis
all_ecoli_deduplicated <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(grepl("Escherichia", kraken2_species)) |>
  group_by(sample, run) |>
  slice_head()  |>
  ungroup() |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  mutate(gff3_filepath = paste0("/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/",run, "/bakta/autocycler_medaka_full/", sample_corrected ,"_bakta/", sample_corrected, ".gff3")) |>
  select(sample, sample_corrected, run, gff3_filepath, hospital, sampletype_cat, kraken2_species, species.x)

#View(all_ecoli_deduplicated)
nrow(all_ecoli_deduplicated) # 1559
length(unique(all_ecoli_deduplicated$sample))

# save txt of filepaths
ecoli_gff3_filepaths <- all_ecoli_deduplicated |>  pull(gff3_filepath) 
# save as txt file
writeLines(ecoli_gff3_filepaths, "ecoli_gff3_filepaths.txt")

#~~~~~~~~~~~~~~~~~~#
# all klebsiellas
all_kleb_deduplicated <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(grepl("Klebsiella", kraken2_species)) |>
  group_by(sample, run) |>
  slice_head()  |>
  ungroup() |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  mutate(gff3_filepath = paste0("/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/",run, "/bakta/autocycler_medaka_full/", sample_corrected ,"_bakta/", sample_corrected, ".gff3")) |>
  select(sample, sample_corrected, run, gff3_filepath, hospital, sampletype_cat, kraken2_species, species.x)

#View(all_kleb_deduplicated)
nrow(all_kleb_deduplicated) # 518
length(unique(all_kleb_deduplicated$sample)) # 518

# save txt of filepaths
kleb_gff3_filepaths <- all_kleb_deduplicated |>   pull(gff3_filepath) 
#View(kleb_gff3_filepaths)
# save as txt file
writeLines(kleb_gff3_filepaths, "kleb_gff3_filepaths.txt")

#~~~~~~~~~~~~~~~~#
# Kleb pneumo species complex
#amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
colnames(amrfinder_metadata_updated)
kpsc_deduplicated <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(typing_scheme == "klebsiella") |>
  group_by(sample, run) |>
  slice_head()  |>
  ungroup() |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  mutate(gff3_filepath = paste0("/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/",run, "/bakta/autocycler_medaka_full/", sample_corrected ,"_bakta/", sample_corrected, ".gff3")) |>
  select(sample, sample_corrected, run, gff3_filepath, hospital, sampletype_cat, kraken2_species, species.x)

View(kpsc_deduplicated)
nrow(kpsc_deduplicated) # 406
length(unique(kpsc_deduplicated$sample)) # 406
table(kpsc_deduplicated$sampletype_cat, useNA = "ifany") # 353 blood, 39 screening
table(kpsc_deduplicated$hospital,  kpsc_deduplicated$sampletype_cat, useNA = "ifany")
table(kpsc_deduplicated$hospital,  kpsc_deduplicated$run, useNA = "ifany")

# save txt of filepaths
kpsc_gff3_filepaths <- kpsc_deduplicated |>
  pull(gff3_filepath) 
kpsc_gff3_filepaths

# save as txt file
writeLines(kpsc_gff3_filepaths, "kpsc_gff3_filepaths.txt")

#~~~~~~~~~~~~~~~~#
# Kleb oxy species complex
#amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
kosc_deduplicated <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(typing_scheme == "koxytoca") |>
  group_by(sample, run) |>
  slice_head()  |>
  ungroup() |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  mutate(gff3_filepath = paste0("/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/",run, "/bakta/autocycler_medaka_full/", sample_corrected ,"_bakta/", sample_corrected, ".gff3")) |>
  select(sample, sample_corrected, run, gff3_filepath, hospital, sampletype_cat, kraken2_species, species.x)

#View(kosc_deduplicated)
nrow(kosc_deduplicated) # 93
length(unique(kosc_deduplicated$sample)) # 93
table(kosc_deduplicated$sampletype_cat, useNA = "ifany") # 81 blood
table(kosc_deduplicated$hospital,  kosc_deduplicated$sampletype_cat, useNA = "ifany")
table(kosc_deduplicated$hospital,  kosc_deduplicated$run, useNA = "ifany")

# save txt of filepaths
kosc_gff3_filepaths <- kosc_deduplicated |>
  pull(gff3_filepath) 
kosc_gff3_filepaths

# save as txt file
writeLines(kosc_gff3_filepaths, "kosc_gff3_filepaths.txt")

#~~~~~~~~~~~~~~~~#
# Kleb aerogenes complex
#amrfinder_metadata_updated <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
kaerogenes_deduplicated <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(typing_scheme == "kaerogenes") |>
  group_by(sample, run) |>
  slice_head()  |>
  ungroup() |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  mutate(gff3_filepath = paste0("/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/",run, "/bakta/autocycler_medaka_full/", sample_corrected ,"_bakta/", sample_corrected, ".gff3")) |>
  select(sample, sample_corrected, run, gff3_filepath, hospital, sampletype_cat, kraken2_species, species.x)

nrow(kaerogenes_deduplicated) # 21
length(unique(kaerogenes_deduplicated$sample)) # 21
table(kaerogenes_deduplicated$sampletype_cat, useNA = "ifany") # 18 blood
table(kaerogenes_deduplicated$hospital,  kaerogenes_deduplicated$sampletype_cat, useNA = "ifany")
table(kaerogenes_deduplicated$hospital,  kaerogenes_deduplicated$run, useNA = "ifany")

# save txt of filepaths
kaerogenes_gff3_filepaths <- kaerogenes_deduplicated |>
  pull(gff3_filepath) 
kaerogenes_gff3_filepaths

# save as txt file
writeLines(kaerogenes_gff3_filepaths, "kaerogenes_gff3_filepaths.txt")

#~~~~~~~~~~~~~~~~#
# just kleb pneumo (no species complex)
# Kleb pneumo species 
kp_deduplicated <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(typing_scheme == "klebsiella") |>
  filter(kraken2_species == "Klebsiella pneumoniae") |>
  group_by(sample, run) |>
  slice_head()  |>
  ungroup() |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  mutate(gff3_filepath = paste0("/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/",run, "/bakta/autocycler_medaka_full/", sample_corrected ,"_bakta/", sample_corrected, ".gff3")) |>
  select(sample, sample_corrected, run, gff3_filepath, hospital, sampletype_cat, kraken2_species, species.x)

View(kp_deduplicated)
nrow(kp_deduplicated) # 352
length(unique(kp_deduplicated$sample)) # 352
table(kp_deduplicated$sampletype_cat, useNA = "ifany") # 3502 blood, 36 screening
#table(kp_deduplicated$hospital,  kp_deduplicated$sampletype_cat, useNA = "ifany")
#table(kp_deduplicated$hospital,  kp_deduplicated$run, useNA = "ifany")

# save txt of filepaths
kp_gff3_filepaths <- kp_deduplicated |>
  pull(gff3_filepath) 
kp_gff3_filepaths

# save as txt file
writeLines(kp_gff3_filepaths, "kp_gff3_filepaths.txt")

#~~~~~~~~~~~~~~~~~~~#
# filter kp deduplicated to include only those in largest pastbaps cluster

#kp_multi_3_fa <- read.csv("kp_multi_3_fa.csv")
#view(kp_multi_3_fa)
kp_only_main_partition <- kp_multi_3_fa |>
  filter(`Level.3` == 29) |>
  pull(Isolates)

kp_deduplicated_cluster_29 <- kp_deduplicated |>
  filter(sample_corrected %in% kp_only_main_partition)
nrow(kp_deduplicated_cluster_29) # 259
# save txt of filepaths
kp_cluster29_gff3_filepaths <- kp_deduplicated_cluster_29 |>
  pull(gff3_filepath) 
kp_cluster29_gff3_filepaths

# save as txt file
writeLines(kp_cluster29_gff3_filepaths, "kp_cluster29_gff3_filepaths.txt")

if FALSE {
# check these are the ones in the large fastbaps cluster
#read in
all_kleb_partition <- read.csv("kleb_by_sc_fastbaps_clusters.csv")
View(all_kleb_partition)
table(all_kleb_partition$Level.3, all_kleb_partition$kraken2_species)

all_kleb_partition_main_cluster_isolates <- all_kleb_partition |>
  filter(Level.3 == "kpsc_14") |>
  pull(sample_corrected)

length(all_kleb_partition_main_cluster_isolates) # 269 with whole KpSC alignment
length(kp_only_main_partition) # 259 with kp only alignment
setdiff(all_kleb_partition_main_cluster_isolates, kp_only_main_partition) # 11 in pan-Kpsc cluster but not in others.
setdiff(kp_only_main_partition, all_kleb_partition_main_cluster_isolates ) # 1 in kp main cluster but not in other one
# try again with core genome alignment of just main Kp cluster, then if still no good, modify parameters of core genome calc. (use lower than 95% threshold to define core genome).
}

#~~~~~~~~~~~~~~~~#
# all citrobacter
all_citrobacter_deduplicated <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(grepl("Citrobacter", kraken2_species)) |>
  group_by(sample, run) |>
  slice_head()  |>
  ungroup() |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  mutate(gff3_filepath = paste0("/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/",run, "/bakta/autocycler_medaka_full/", sample_corrected ,"_bakta/", sample_corrected, ".gff3")) |>
  select(sample, sample_corrected, run, gff3_filepath, hospital, sampletype_cat, kraken2_species, species.x)

#View(all_citrobacter_deduplicated)
nrow(all_citrobacter_deduplicated) # 20 deduplicated (26 otherwise)
length(unique(all_citrobacter_deduplicated$sample)) # 518

# save txt of filepaths
citrobacter_gff3_filepaths <- all_citrobacter_deduplicated |>   pull(gff3_filepath) 
#View(citrobacter_gff3_filepaths)
# save as txt file
writeLines(citrobacter_gff3_filepaths, "citrobacter_gff3_filepaths.txt")

#~~~~~~~~~~~~~~~~~~~#
#all enterobacters
all_enterobacter_deduplicated <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(grepl("Enterobacter", kraken2_species)) |>
  group_by(sample, run) |>
  slice_head()  |>
  ungroup() |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  mutate(gff3_filepath = paste0("/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/",run, "/bakta/autocycler_medaka_full/", sample_corrected ,"_bakta/", sample_corrected, ".gff3")) |>
  select(sample, sample_corrected, run, gff3_filepath, hospital, sampletype_cat, kraken2_species, species.x)

#View(all_enterobacter_deduplicated)
nrow(all_enterobacter_deduplicated) # 21 deduplicated (24 without deduplication)
length(unique(all_enterobacter_deduplicated$sample)) # 21

# save txt of filepaths
enterobacter_gff3_filepaths <- all_enterobacter_deduplicated |>   pull(gff3_filepath) 
#View(enterobacter_gff3_filepaths)
# save as txt file
writeLines(enterobacter_gff3_filepaths, "enterobacter_gff3_filepaths.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# E. coli Overall ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# try with just snps as the pile is too large otherwise
#fasta.file.name <- "/mnt/nfs-gram-negative-study/assembly/assembly_main_pipeline_v2/pangenomes/all_ecoli_gff/panaroo_results2/core_gene_alignment.aln"
fasta.file.name <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/all_ecoli_core_snps.fasta"
sparse.data <- import_fasta_sparse_nt(fasta.file.name)
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric") #> [1] "Optimised hyperparameter: 0.011"
ecoli.baps.hc <- fast_baps(sparse.data)

ecoli_multi_3 <- multi_res_baps(sparse.data, levels = 3)
table(ecoli_multi_3$`Level 1`)
table(ecoli_multi_3$`Level 2`)
table(ecoli_multi_3$`Level 3`)

# merge
all_ecoli_partition <- all_ecoli_deduplicated |>
  left_join(ecoli_multi_3, by = c("sample_corrected" = "Isolates"))
nrow(all_ecoli_partition)
View(all_ecoli_partition)
sum(is.na(all_ecoli_partition$`Level 3`))

# save
write.csv(all_ecoli_partition, "ecoli_fastbaps_clusters.csv", row.names = FALSE)
#read in
#all_ecoli_partition <- read.csv("ecoli_fastbaps_clusters.csv")
#View(all_ecoli_partition)
#~~~~~~~~~~~~~~~#
# look at stability of baps clusters
boot.result <- boot_fast_baps(sparse.data)
dendro <- as.dendrogram(fast_baps(sparse.data))
#> [1] "Calculating initial clustering..."
#> [1] "Calculating initial dk values..."
#> [1] "Clustering using hierarchical Bayesian clustering..."
gplots::heatmap.2(boot.result, dendro, dendro, tracecol = NA)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Klebsiella Overall ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#fasta.file.name <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/all_kleb_core_gene_alignment.aln"
fasta.file.name <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/all_kleb_core_snps.fasta"
sparse.data <- import_fasta_sparse_nt(fasta.file.name)
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric") # 0.011
kleb_multi_3 <- multi_res_baps(sparse.data, levels = 3)
table(kleb_multi_3$`Level 1`)
table(kleb_multi_3$`Level 2`)
table(kleb_multi_3$`Level 3`)


# merge
all_kleb_partition <- all_kleb_deduplicated |>
  left_join(kleb_multi_3, by = c("sample_corrected" = "Isolates"))
nrow(all_kleb_partition)
#View(all_kleb_partition)
sum(is.na(all_kleb_partition$`Level 3`))
sum(is.na(all_kleb_partition$hospital))

# save
write.csv(all_kleb_partition, "kleb_fastbaps_clusters.csv", row.names = FALSE)
#read in
#all_kleb_partition <- read.csv("kleb_fastbaps_clusters.csv")
#View(all_kleb_partition)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Separately for each Kleb species complex
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Klebsiella pneumoniae ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fasta.file.name.fa <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/kpsc_core_snps.fasta"
sparse.data.fa <- import_fasta_sparse_nt(fasta.file.name.fa)
sparse.data.fa <- optimise_prior(sparse.data.fa, type = "optimise.symmetric") #> [1] "Optimised hyperparameter: 0.025"
kpsc_multi_3_fa <- multi_res_baps(sparse.data.fa, levels = 3)
# try clustering to more levels, as 3 still groups 244 Kpsc together
kpsc_multi_6_fa <- multi_res_baps(sparse.data.fa, levels = 6)
kpsc_multi_20_fa <- multi_res_baps(sparse.data.fa, levels = 20)

nrow(kpsc_multi_3_fa) #406
max(kpsc_multi_3_fa$`Level 1`) # 3
max(kpsc_multi_3_fa$`Level 2`) # 7
max(kpsc_multi_3_fa$`Level 3`) # 15
table(kpsc_multi_3_fa$`Level 3`) # 15


View(kpsc_multi_6_fa) # 3
max(kpsc_multi_6_fa$`Level 1`) # 3
max(kpsc_multi_6_fa$`Level 2`) # 7
max(kpsc_multi_6_fa$`Level 3`) # 15
max(kpsc_multi_6_fa$`Level 4`) # 26
max(kpsc_multi_6_fa$`Level 5`) # 36
max(kpsc_multi_6_fa$`Level 6`) # 45
table(kpsc_multi_6_fa$`Level 6`)
table(kpsc_multi_6_fa$`Level 5`)
table(kpsc_multi_6_fa$`Level 4`)
table(kpsc_multi_6_fa$`Level 3`)
# the 269 in cluster 14 (L3), cluster 23 (L4), cluster 27 (L5), or cluster 29 (L6) are all still together 

max(kpsc_multi_20_fa$`Level 1`) # 3
max(kpsc_multi_20_fa$`Level 2`) # 7
max(kpsc_multi_20_fa$`Level 3`) # 15
max(kpsc_multi_20_fa$`Level 4`) # 26
max(kpsc_multi_20_fa$`Level 5`) # 36
max(kpsc_multi_20_fa$`Level 6`) # 43
max(kpsc_multi_20_fa$`Level 7`) # 45
max(kpsc_multi_20_fa$`Level 8`) # 45
max(kpsc_multi_20_fa$`Level 9`) # 45
max(kpsc_multi_20_fa$`Level 10`) # 36
max(kpsc_multi_20_fa$`Level 11`) # 36
max(kpsc_multi_20_fa$`Level 12`) # 45
max(kpsc_multi_20_fa$`Level 13`) # 45
max(kpsc_multi_20_fa$`Level 14`) # 45
max(kpsc_multi_20_fa$`Level 15`) # 45
max(kpsc_multi_20_fa$`Level 16`) # 45
max(kpsc_multi_20_fa$`Level 17`) # 45
max(kpsc_multi_20_fa$`Level 18`) # 45
max(kpsc_multi_20_fa$`Level 19`) # 45
max(kpsc_multi_20_fa$`Level 20`) # 45
table(kpsc_multi_20_fa$`Level 20`) # 45 is may cluster number - 269 isolates still all together in main Kleb cluster



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Klebsiella oxytoca ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fasta.file.name.fa <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/kosc_core_snps.fasta"
sparse.data.fa <- import_fasta_sparse_nt(fasta.file.name.fa)
sparse.data.fa <- optimise_prior(sparse.data.fa, type = "optimise.symmetric") #> [1] "Optimised hyperparameter: 0.008"
kosc_multi_3_fa <- multi_res_baps(sparse.data.fa, levels = 3)
nrow(kosc_multi_3_fa) # 93
max(kosc_multi_3_fa$`Level 1`) # 5
max(kosc_multi_3_fa$`Level 2`) # 8
max(kosc_multi_3_fa$`Level 3`) # 15

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Klebsiella aerogenes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fasta.file.name.fa <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/kaerogenes_core_snps.fasta"
sparse.data.fa <- import_fasta_sparse_nt(fasta.file.name.fa)
sparse.data.fa <- optimise_prior(sparse.data.fa, type = "optimise.symmetric") #> [1] "Optimised hyperparameter: 0.062"
kaerogenes_multi_3_fa <- multi_res_baps(sparse.data.fa, levels = 3)
nrow(kaerogenes_multi_3_fa) #21
max(kaerogenes_multi_3_fa$`Level 1`) # 2
max(kaerogenes_multi_3_fa$`Level 2`) # 2
max(kaerogenes_multi_3_fa$`Level 3`) # 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# only Klebsiella pneumoniae senso stricto ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fasta.file.name.fa <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/kp_core_snps.fasta"
sparse.data.fa <- import_fasta_sparse_nt(fasta.file.name.fa)
sparse.data.fa <- optimise_prior(sparse.data.fa, type = "optimise.symmetric") #> [1] "Optimised hyperparameter: 0.008"
kp_multi_3_fa <- multi_res_baps(sparse.data.fa, levels = 3)
nrow(kp_multi_3_fa) # 352
max(kp_multi_3_fa$`Level 1`) # 4
max(kp_multi_3_fa$`Level 2`) # 12
max(kp_multi_3_fa$`Level 3`) # 29
table(kp_multi_3_fa$`Level 1`) # 4
table(kp_multi_3_fa$`Level 2`) # 12
table(kp_multi_3_fa$`Level 3`) # 29


kp_multi_30_fa <- multi_res_baps(sparse.data.fa, levels = 30)
max(kp_multi_3_fa$`Level 2`) # 12
max(kp_multi_3_fa$`Level 3`) # 29
max(kp_multi_30_fa$`Level 30`) # 41
max(kp_multi_30_fa$`Level 28`) # 41
max(kp_multi_30_fa$`Level 29`) # 41

table(kp_multi_3_fa$`Level 1`) # 
table(kp_multi_3_fa$`Level 2`) # 
table(kp_multi_30_fa$`Level 30`) # still 259 in largest group

# save 
write.csv(kp_multi_3_fa, "kp_multi_3_fa.csv", row.names = FALSE)
# read in
#kp_multi_3_fa <- read.csv("kp_multi_3_fa.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# only Klebsiella pneumoniae  cluster 29 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fasta.file.name.fa <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/kp_cluster29_core_snps.fasta"
sparse.data.fa <- import_fasta_sparse_nt(fasta.file.name.fa)
sparse.data.fa <- optimise_prior(sparse.data.fa, type = "optimise.symmetric") #> [1] "Optimised hyperparameter: 0.051"
kp_cl29_multi_3_fa <- multi_res_baps(sparse.data.fa, levels = 3)

nrow(kp_cl29_multi_3_fa) # 
max(kp_cl29_multi_3_fa$`Level 1`) # 
max(kp_cl29_multi_3_fa$`Level 2`) # 
max(kp_cl29_multi_3_fa$`Level 3`) # 
table(kp_cl29_multi_3_fa$`Level 1`) # 
table(kp_cl29_multi_3_fa$`Level 2`) # 
table(kp_cl29_multi_3_fa$`Level 3`) # 


kp_cl29_multi_30_fa <- multi_res_baps(sparse.data.fa, levels = 30)
max(kp_cl29_multi_30_fa$`Level 3`) # 
max(kp_cl29_multi_30_fa$`Level 30`) # 
max(kp_cl29_multi_30_fa$`Level 28`) # 
max(kp_cl29_multi_30_fa$`Level 29`) # 

table(kp_cl29_multi_30_fa$`Level 1`) # 
table(kp_cl29_multi_30_fa$`Level 2`) # 
table(kp_cl29_multi_30_fa$`Level 30`) # 
# all single cluster??? why is this? 
# save 
write.csv(kp_cl29_multi_3_fa, "kp_cl29_multi_3_fa.csv", row.names = FALSE)
# read in
#kp_cl29_multi_3_fa <- read.csv("kp_cl29_multi_3_fa.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# only Klebsiella pneumoniae  cluster 29 - at 90% threshold for core genome ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fasta.file.name.fa <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/kp_cl29_90th_core_snps.fasta"
sparse.data.fa <- import_fasta_sparse_nt(fasta.file.name.fa)
sparse.data.fa <- optimise_prior(sparse.data.fa, type = "optimise.baps") #> [1] "Optimised hyperparameter: 0.051" # 0.151 with optimised fastbaps
kp_cl29_90th_multi_3_fa <- multi_res_baps(sparse.data.fa, levels = 3)

nrow(kp_cl29_90th_multi_3_fa) # 
max(kp_cl29_90th_multi_3_fa$`Level 1`) # 
max(kp_cl29_90th_multi_3_fa$`Level 2`) # 
max(kp_cl29_90th_multi_3_fa$`Level 3`) # 
table(kp_cl29_90th_multi_3_fa$`Level 1`) # 
table(kp_cl29_90th_multi_3_fa$`Level 2`) # 
table(kp_cl29_90th_multi_3_fa$`Level 3`) # 


kp_cl29_multi_30_fa <- multi_res_baps(sparse.data.fa, levels = 30)
max(kp_cl29_multi_30_fa$`Level 3`) # 
max(kp_cl29_multi_30_fa$`Level 30`) # 
max(kp_cl29_multi_30_fa$`Level 28`) # 
max(kp_cl29_multi_30_fa$`Level 29`) # 

table(kp_cl29_multi_30_fa$`Level 1`) # 
table(kp_cl29_multi_30_fa$`Level 2`) # 
table(kp_cl29_multi_30_fa$`Level 30`) # 
# all single cluster??? why is this? 
# save 
write.csv(kp_cl29_multi_3_fa, "kp_cl29_multi_3_fa.csv", row.names = FALSE)
# read in
#kp_cl29_multi_3_fa <- read.csv("kp_cl29_multi_3_fa.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  combine all 3 tables 
# add a species complex column
kpsc_multi_3_fa <- kpsc_multi_3_fa |>
  mutate(species_complex = "kpsc",
         `Level 1` = paste0(species_complex, "_", `Level 1`),
         `Level 2` = paste0(species_complex, "_", `Level 2`),
         `Level 3` = paste0(species_complex, "_", `Level 3`)
                            )
kosc_multi_3_fa <- kosc_multi_3_fa |>
  mutate(species_complex = "kosc",
         `Level 1` = paste0(species_complex, "_", `Level 1`),
         `Level 2` = paste0(species_complex, "_", `Level 2`),
         `Level 3` = paste0(species_complex, "_", `Level 3`)
  )
kaerogenes_multi_3_fa <- kaerogenes_multi_3_fa |>
  mutate(species_complex = "kaerogenes",
         `Level 1` = paste0(species_complex, "_", `Level 1`),
         `Level 2` = paste0(species_complex, "_", `Level 2`),
         `Level 3` = paste0(species_complex, "_", `Level 3`)
  )

all_kleb_fastbaps_3l <- rbind(kpsc_multi_3_fa, kosc_multi_3_fa, kaerogenes_multi_3_fa)
nrow(all_kleb_fastbaps_3l) # 520
table(all_kleb_fastbaps_3l$`Level 1`, all_kleb_fastbaps_3l$`Level 2`)
length(unique(all_kleb_fastbaps_3l$`Level 1`)) # 10
length(unique(all_kleb_fastbaps_3l$`Level 2`)) # 17
length(unique(all_kleb_fastbaps_3l$`Level 3`)) #32
# merge with all_kleb deduplicated table


# merge
all_kleb_partition <- all_kleb_deduplicated |>
  left_join(all_kleb_fastbaps_3l, by = c("sample_corrected" = "Isolates"))
nrow(all_kleb_partition) # 518
#View(all_kleb_partition)
sum(is.na(all_kleb_partition$`Level 3`))
sum(is.na(all_kleb_partition$hospital))

# save
write.csv(all_kleb_partition, "kleb_by_sc_fastbaps_clusters.csv", row.names = FALSE)
#read in
#all_kleb_partition <- read.csv("kleb_by_sc_fastbaps_clusters.csv")
#View(all_kleb_partition)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Enterobacter Overall ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fasta.file.name <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/all_enterobacter_core_gene_alignment.aln"
sparse.data <- import_fasta_sparse_nt(fasta.file.name)
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric") #> [1] "Optimised hyperparameter: 0.016"
entero_multi_3 <- multi_res_baps(sparse.data, levels = 3)
table(entero_multi_3$`Level 1`)
table(entero_multi_3$`Level 2`)
table(entero_multi_3$`Level 3`)

# merge
all_entero_partition <- all_enterobacter_deduplicated |>
  left_join(entero_multi_3, by = c("sample_corrected" = "Isolates"))
nrow(all_entero_partition)
#View(all_entero_partition)
sum(is.na(all_entero_partition$`Level 3`))
sum(is.na(all_entero_partition$hospital))

# save
write.csv(all_entero_partition, "entero_fastbaps_clusters.csv", row.names = FALSE)
#read in
#all_entero_partition <- read.csv("entero_fastbaps_clusters.csv")
#View(all_entero_partition)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Citrobacter Overall ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fasta.file.name <- "~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/all_citrobacter_core_gene_alignment.aln"
sparse.data <- import_fasta_sparse_nt(fasta.file.name)
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric") 
#"Optimised hyperparameter: 0.021" # second least conservative
# 0.211 for "optimised.baps"

citro_multi_3 <- multi_res_baps(sparse.data, levels = 3)
table(citro_multi_3$`Level 1`)
table(citro_multi_3$`Level 2`)
table(citro_multi_3$`Level 3`)

# merge
all_citro_partition <- all_citrobacter_deduplicated |>
  left_join(citro_multi_3, by = c("sample_corrected" = "Isolates"))
nrow(all_citro_partition)
#View(all_citro_partition)
sum(is.na(all_citro_partition$`Level 3`))
sum(is.na(all_citro_partition$hospital))

# save
write.csv(all_citro_partition, "citro_fastbaps_clusters.csv", row.names = FALSE)
#read in
#all_citro_partition <- read.csv("citro_fastbaps_clusters.csv")
#View(all_citro_partition)


#citro.baps.hc <- fast_baps(sparse.data)
#citro.clusters <- best_baps_partition(sparse.data, as.phylo(citro.baps.hc))
#citro.clusters
#table(citro.clusters) # same for optimise baps, optimise symmetric, hc, and baps
#1  2  3  4  5 
#3  2 13  1  1 
  

#alternatively
#d <- snp_dist(sparse.data)
#d <- as.dist(d/max(d))
#h <- hclust(d, method="ward.D2")
#partition <- best_baps_partition(sparse.data, h)
#partition
#table(partition)
#citro_cluster_df <- data.frame(
#  isolate = names(citro.clusters),
#  cluster = citro.clusters)
#save
#write.csv(citro_cluster_df, "citro_fastbaps_clusters.csv", row.names = FALSE)
#read in
#citro_clusters_df <- read.csv("citro_fastbaps_clusters.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot results with pre-calculated tree
library(ggnewscale)        # for new_scale_fill()
library(RColorBrewer)

#tree <- read.tree('mashtrees/ecoli_mashtree.bootstrap.dnd')
#tree$edge.length<-sqrt(tree$edge.length) #this just makes the tree look better
#tree<-midpoint_root(tree) #so does this

tree.file.name <- "mashtrees/ecoli_mashtree.bootstrap.dnd"
mashtree <- phytools::read.newick(tree.file.name)

# add names to all_ecoli_partition to match tree tip names
all_ecoli_partition <- all_ecoli_partition |>
  mutate(tip_names = paste0(run, "_", sample_corrected),
         tip_names = gsub("run0", "assembly_main_pipeline_no_db",  tip_names),
         tip_names = gsub("microbesng_failed", "microbesng_failed_assembly",  tip_names))

#check overalp between tip names and 
setdiff(mashtree$tip.label, all_ecoli_partition$tip_names) # 68 in mashtree tips, but not in all ecoli deduplicated
extra_tree_tips <- setdiff(mashtree$tip.label, all_ecoli_partition$tip_names) # 68 in mashtree tips, but not in all ecoli deduplicated
setdiff( all_ecoli_partition$tip_names, mashtree$tip.label) # 1 among the gff3s, but not in tree tips??
common_tips <- intersect(mashtree$tip.label, all_ecoli_partition$tip_names)
length(intersect(mashtree$tip.label, all_ecoli_partition$tip_names)) # 1558 - so almost all.

plot.df <- all_ecoli_partition |>
  select(tip_names, `Level 1`)|>
  rename(id = tip_names,
         fastbaps = `Level 1`) 

#View(plot.df)
plot.df$id
# filter plot tips to only include tips in both tree and clustering
plot.df <- plot.df |>
  filter(id %in% common_tips )
nrow(plot.df) # 1558

# filter tree tips to inly include 
mashtree <- drop.tip(mashtree, extra_tree_tips)
length(mashtree$tip.label) # 1558

# root tree, make neater
mashtree$edge.length<-sqrt(mashtree$edge.length) #this just makes the tree look better
mashtree<-midpoint_root(mashtree) #so does this

gg <- ggtree(mashtree)
print(gg)

f2 <- facet_plot(gg, panel = "fastbaps", data = plot.df, geom = geom_tile, aes(x = fastbaps), 
                 color = "blue")
f2


# create metadata panels data subsets
amrfinder_ecoli <- amrfinder_metadata_updated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(grepl("Escherichia", kraken2_species)) |>
  mutate(sample_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                      run != "run0" ~ gsub("AH", "AH-", sample)),
         sample_corrected = gsub("AB90", "AB90-A", sample_corrected)) |>
  
  mutate(tip_names = paste0(run, "_", sample_corrected),
         tip_names = gsub("run0", "assembly_main_pipeline_no_db",  tip_names),
         tip_names = gsub("microbesng_failed", "microbesng_failed_assembly",  tip_names)) |>
  select(c(tip_names, everything()))

data_full <- amrfinder_ecoli |>
  group_by(sample, run) |>
  slice_head() |>
  ungroup() |>
  group_by(typing_scheme, mlst_profile) |>
  mutate(mlst_count = n()) |>
  mutate(mlst_grouped = case_when(mlst_count == 1 ~ "single",
                                  mlst_count >=2 & mlst_count <=5 ~ "minor (2-5 occ.)",
                                  mlst_count >=6 & mlst_count <=20 ~ "minor (6-20 occ.)",
                                  # mlst_count >20 & mlst_count <=50 ~ "minor (21-50 occ.)",
                                  mlst_count >20 ~ paste0("ST", mlst),
                                  mlst_count == "" ~ NA,
                                  is.na(mlst_count) ~ NA)) |>
  ungroup() 

colnames(data_full)
table(data_full$mlst_grouped)
  
  
time_place_hcai <- amrfinder_ecoli |>
  group_by(sample, run) |>
  slice_head() |>
  ungroup() |>
  select(tip_names, sampletype_cat, region, healthcare_exposure) |>
  rename(sampletype = sampletype_cat) |>
  filter(tip_names %in% common_tips) |>
pivot_longer(cols = c(sampletype, region, healthcare_exposure) , names_to = "var")
#View(time_place_hcai)
time_place_hcai$var <- factor(time_place_hcai$var, levels = c("sampletype", "region", "healthcare_exposure"))
time_place_hcai$value <- as.factor(time_place_hcai$value)

setdiff(mashtree$tip.label, unique(time_place_hcai$tip_names))
setdiff(unique(time_place_hcai$tip_names), mashtree$tip.label)
length(unique(time_place_hcai$tip_names)) # 1558

# panel data for clusters
fastbaps_partitions <- all_ecoli_partition |>
  select(c(tip_names, `Level 1`, `Level 2`, `Level 3`)) |>
  pivot_longer(cols = c(`Level 1`, `Level 2`, `Level 3`), names_to = "var")
#View(fastbaps_partitions)
fastbaps_partitions$value <- as.factor(fastbaps_partitions$value)

# gene counts
gene_counts <- data.frame(table(amrfinder_ecoli$Element.symbol)) #|> 
# filter(Freq>0) #e.g. filter out if occur <5 times in population

carb <- filter(amrfinder_ecoli, grepl('CARBAPENEM',Subclass)) |> 
  distinct(Element.symbol) #|> filter(Element.symbol %in% gene_counts$Var1) #do not filter
carb_selected <- amrfinder_ecoli |>
  filter(Element.symbol %in% carb$Element.symbol) |>
  select(tip_names, Element.symbol, molecule_type)|>
  complete( tip_names = mashtree$tip.label, Element.symbol = carb$Element.symbol, fill = list(molecule_type= "absent")) |>
  # remove multi-copy entries:
  group_by(tip_names, Element.symbol, molecule_type) |>
  slice_head() |>
  ungroup()

# remove a single sample that is present in both chromosome and plasmid
#failed_assembly_qc_AF235  
both_ch_and_p <- carb_selected  |>
  group_by(tip_names, Element.symbol) |>
  summarise(count = n()) |>
  filter(count >1) |>
  pull(tip_names)
carb_selected <- carb_selected |>
  mutate(molecule_type = case_when(tip_names %in% both_ch_and_p & Element.symbol == "blaNDM-5" ~ "both",
                                   TRUE ~ molecule_type)) |>
  group_by(tip_names, Element.symbol) |>
  slice_head() |>
  ungroup() |>
  rename(present = molecule_type,
         gene = Element.symbol)
carb_selected$present <- factor(carb_selected$present, levels = c("absent", "chromosome", "plasmid", "both"))
#View(carb_selected)



# set colours for different panels:
# sample type 
colours_type <- c("blood" = "darkred","screening" = "goldenrod1", "unknown" = "grey")
#region
colours_region <- c(
  # South & London → warm reds/oranges
  "London"              = "#E41A1C",  # strong red (warm, central)
  "South East A"        = "#FB8072",  # light coral
  "South East B"        = "#FCAE91",  # salmon
  "South East C"        = "#FEE0D2",  # pale rose
  "South West"          = "#FD8D3C",  # orange
  # Midlands → neutral yellow–green transition
  "Midlands"            = "#B3DE69",  # warm yellow-green
  # East → intermediate green (transitional)
  "East"                = "#66C2A5",  # turquoise-green
  # North regions → cool blues
  "North West"          = "#377EB8",  # deep blue
  "North East & Yorkshire A" = "#807DBA",  # purple-blue
  "North East & Yorkshire B" = "#B3B3E6"  # light lavender
  # Unknown → neutral
  
)

colours_hcai <- c(
  "community" = "#FFF3B0",
  "quasi-community" = "#FFD37F",
  "quasi-nosocomial" = "#FCA17D",
  "nosocomial" = "#C75197"
)

colours_combined <- c(colours_type, colours_region)
colours_combined <- c(colours_combined, c("-" = "white" ))
colours_combined <- c(colours_combined, colours_hcai)
# reorder:

# Define category order explicitly
sampletype_vars <- names(colours_type)
region_vars <- names(colours_region)
healthcare_vars <- names(colours_hcai)

# Combine in the desired order, plus the "-" category
ordered_levels <- c(sampletype_vars, "-", region_vars, healthcare_vars)

# Apply ordering to the factor
time_place_hcai$value <- factor(time_place_hcai$value, levels = ordered_levels)

# carbapenemase genes
colours <- c("absent" = "#ffffff","plasmid" = "#0d0c3b", "chromosome" = "grey", both = "#377EB8")


#set order
library(colorspace)
mlst_order <- c("ST131", "ST73", "ST69", "ST95", "ST12", "ST127", "ST1193", 
                "ST10", "ST38", "ST404", "ST141", "ST648", 
                "minor (6-20 occ.)", "minor (2-5 occ.)", "single")
data_full$mlst_grouped <- factor(data_full$mlst_grouped, levels = mlst_order)

# set colours:
# Paul's tol bright:
mlst_colours <- c(
  "#EE6677", "#228833", "#4477AA", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB", 
  "#2288AA", "#AA4499", "#44AA99", "#999933", "#882255",  
  "#DDCC77", "#88CCEE", lighten("goldenrod4", amount = 0.7)
)  
names(mlst_colours) <- mlst_order

# fastbaps levels
# get unique fastbaps levels (as characters)
fb_levels <- sort(unique(as.character(fastbaps_partitions$value)))
n_fb <- length(fb_levels)

# create expanded colour palette
base_set3 <- brewer.pal(min(12, 12), "Set3")     # base palette (12)
fastbaps_colours <- colorRampPalette(base_set3)(n_fb)
# name colours so mapping is explicit
names(fastbaps_colours) <- fb_levels


fastbaps_colours <- RColorBrewer::brewer.pal(165, name = "Set3")

#~~~~~~~~~~~~~~~#
# plot tree
tree_plot<-ggtree(mashtree)
p <- tree_plot %<+% data_full + 
  geom_tippoint(aes(color=mlst_grouped), size =3) +  
  scale_color_manual(values = mlst_colours) + 
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_treescale(x = 0, y = length(mashtree$tip.label)*0.95, offset = length(mashtree$tip.label)/100,
                 width = 0.1, color = "black") + theme_tree()

p
#~~~~~~~~~~~~~~#
p2 <- p + 
  #  panel for sampletype
  geom_facet( panel = "REGION, TYPE & HCAI", data = time_place_hcai ,geom = geom_tile,
              mapping = aes(x = var, fill = value),width = 0.8) +
  scale_fill_manual(name = "Sample Type, Region & \n Healthcare exposure", values = colours_combined, 
                    guide = guide_legend(title.position = "top")) +
  ggnewscale::new_scale_fill() + # reset
  geom_facet( panel = "fastBAPS partition", data = fastbaps_partitions ,geom = geom_tile,
              mapping = aes(x = var, fill = value),width = 0.8) +
  scale_fill_manual(name = "FastBAPS partitions", values = fastbaps_colours, 
                    guide = guide_legend(title.position = "top", ncol = 10)) +
  new_scale_fill() + # reset
  # carb res genes
  geom_facet(panel = "CARBAPENEMASE GENES", data = carb_selected, geom = geom_tile, 
            aes(x = gene, fill = present),  width = 0.8) + 
  scale_fill_manual(name = "Gene Location", values = colours, guide = guide_legend(title.position = "top")) +
  # theme 
  theme(axis.text.x = element_text(angle=90),strip.text = element_text(size = 8),strip.background = element_rect(fill = "white", color = "black")) +
  scale_x_discrete(expand = c(0, 0)) +   # make tree take up more of lateral space
  coord_cartesian(clip = "off") +
  labs(color="MLST")

p2
p3 <-facet_widths(p2,widths = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) 
ggsave(file ="all_ecoli_mashtree_carbapenemases_plot_fastbaps_partitions.png", plot =p2, width =17, height = 11, dpi = 300)

p3


