# How much is enough for genomic surveillance of Escherichia and Klebsiella BSIs in England?

This repository contains code to estiamte the same size of bloostream infection (BSI)-associated E. coli and Klebsiella isolates required for genomic sequencing to achieve a specified sampling coverage using Bayesian bootstrapping combined with power calculation. These analyses use data from the NEKSUS (National E. coli and Klebsiella bloodstream infection and CPE UK Surveillance) Study. Data for this study are underembargo until April 2027. The results of this analysis are presented as a poster at ESCMID 2026. 

# Extended methods
## Isolate collection
Ten NHS microbiology labs from nine English NHS Trusts (groups of hospitals under the same administration) representing the largest in terms of number of emergency admissions across all seven NHS England regions were recruited to the NEKSUS Consortium. Consecutive, unselected BSI and CPE-positive rectal screening isolates were collected between October 2023 and March 2024 as part of routine clinical practice. Isolates were stored in brain–heart infusion broth with 10% glycerol at −70 °C and then grown on blood agar for 24 h at 37 °C, following which a colony sweep of the pure bacterial culture was suspended in 1 ml phosphate buffer saline, pelleted and cold-packed (for GENEWIZ sequencing) or suspended in 1 ml DNA/RNA Shield (Zymo Research, USA) and shipped at ambient temperature following MicrobesNG strain submission procedures (for MicrobesNG sequencing). Bacteria were subcultured for a further 24 h at 37 °C where there was insufficient growth after 24 h.

## DNA extraction and sequencing
DNA extraction, library preparation and sequencing were conducted at GENEWIZ Germany GmbH (Leipzig, Germany) and MicrobesNG (Birmingham, United Kingdom). For isolates processed by GENEWIZ, DNA was extracted using the MagMAX Microbiome Ultra Nucleic Acid Isolation Kit with bead plate (Life Technologies, Carlsbad, CA, USA). Genomic DNA was quantified using the Qubit 4.0 Fluorometer and qualified using the Agilent 5600 Fragment Analyzer. For isolates processed by MicrobesNG, cells were lysed using TE buffer containing lysozyme (MPBio, USA), metapolyzyme (Sigma-Aldrich, USA) and RNase A (ITW Reagents, Spain), followed by treatment with proteinase K (VWR Chemicals, Ohio, USA) and SDS (Sigma-Aldrich, Missouri, USA). Genomic DNA was purified using SPRI beads, resuspended in EB buffer (10 mM Tris-HCl, pH 8.0), and quantified with the Quant-iT dsDNA HS (ThermoFisher Scientific) assay in an Eppendorf AF2200 plate reader (Eppendorf UK Ltd, United Kingdom) and diluted as appropriate. 

For both sequencing providers, DNA libraries were prepared for sequencing using Rapid Barcoding Kit 96 V14 (Oxford Nanopore Technologies, Oxford, UK) according to the manufacturer’s recommendations. Briefly, sequencing libraries were generated using a transposase, which simultaneously cleaves template molecules and attaches barcoded tags to the cleaved ends. The barcoded samples were then pooled (96-plexed) before solid-phase reversible immobilisation-cleaning and addition of Rapid Adapters to the tagged ends. The library pools were loaded onto ONT PromethION flow cells (R10 [M Version]) – one 96-plex pool per flow cell – and sequenced on a PromethION P2 Solo for 72 h according to the manufacturer’s instructions (GENEWIZ), or were loaded onto a FLO-PRO114M (R.10.4.1) flow cell.

## Bioinformatic analysis
Computational analysis was performed on a virtual machine in the Oracle Cloud Infrastructure. POD5 files were basecalled and demultiplexed using Dorado [27] v5.0.0 (super-high accuracy 5mCG, 5hmCG and 6mA methylation aware simplex DNA model). All bioinformatic tools were run using default settings unless otherwise specified. Raw-read quality was evaluated with SeqKit [28] v2.9.0. Long-reads were subsampled to generate 4 independent sets of 60× reads using the built-in subsampling and genome size estimation scripts from Autocycler [20, 21] v0.2.1. Genome assembly was done using Autocycler v0.2.1 using 20 input assemblies (4x each of Canu [31] v2.2, Flye [30], Raven [32] v1.8.3, Miniasm [33] v0.3 and Plassembler). A single 60× subsampling depth was used based on findings of previous benchmarking work [1, 15, 34] that the assembly performance of Flye is comparable in the range 60–100× read depth. The lower end of this range, 60×, was selected to minimise compute use. Assemblies were polished using Medaka v2.0.1 using un-subsampled long-reads ( `--bacteria` flag). 

QC

##AMR Annotation
bakta -> panaroo -> snp-sites -> fastbaps
AMRFInderPLUS
MLST/ kleborate
MOB-suite-> PLING


## Bayesian Bootstrapping

## Statistical analyses and visualisation
