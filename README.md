# How much is enough for genomic surveillance of *Escherichia* and *Klebsiella* BSIs in England?

This repository contains code to estiamte the same size of bloostream infection (BSI)-associated *E. coli* and *Klebsiella* isolates required for genomic sequencing to achieve a specified sampling coverage using Bayesian bootstrapping combined with power calculation. These analyses use data from the NEKSUS (National *E. coli* and *Klebsiella* bloodstream infection and CPE UK Surveillance) Study. Data for this study are underembargo until April 2027. The results of this analysis are presented as a poster at ESCMID 2026. 

# Extended methods
## Isolate collection
Ten NHS microbiology labs from nine English NHS Trusts (groups of hospitals under the same administration) representing the largest in terms of number of emergency admissions across all seven NHS England regions were recruited to the NEKSUS Consortium. Consecutive, unselected *E. coli* and *Klebsiella* BSI and CPE-positive rectal screening isolates were collected between October 2023 and March 2024 as part of routine clinical practice. Isolates of the same bacterial species from the same patient were deduplicated using a standard 14-day rolling window. Isolates were stored in brain–heart infusion broth with 10% glycerol at −70 °C and then grown on blood agar for 24 h at 37 °C, following which a colony sweep of the pure bacterial culture was suspended in 1 ml phosphate buffer saline, pelleted and cold-packed (for GENEWIZ sequencing) or suspended in 1 ml DNA/RNA Shield (Zymo Research, USA) and shipped at ambient temperature following MicrobesNG strain submission procedures (for MicrobesNG sequencing). Bacteria were subcultured for a further 24 h at 37 °C where there was insufficient growth after 24 h.

## DNA extraction and sequencing
DNA extraction, library preparation and sequencing were conducted at GENEWIZ Germany GmbH (Leipzig, Germany) and MicrobesNG (Birmingham, United Kingdom). For isolates processed by GENEWIZ, DNA was extracted using the MagMAX Microbiome Ultra Nucleic Acid Isolation Kit with bead plate (Life Technologies, Carlsbad, CA, USA). Genomic DNA was quantified using the Qubit 4.0 Fluorometer and qualified using the Agilent 5600 Fragment Analyzer. For isolates processed by MicrobesNG, cells were lysed using TE buffer containing lysozyme (MPBio, USA), metapolyzyme (Sigma-Aldrich, USA) and RNase A (ITW Reagents, Spain), followed by treatment with proteinase K (VWR Chemicals, Ohio, USA) and SDS (Sigma-Aldrich, Missouri, USA). Genomic DNA was purified using SPRI beads, resuspended in EB buffer (10 mM Tris-HCl, pH 8.0), and quantified with the Quant-iT dsDNA HS (ThermoFisher Scientific) assay in an Eppendorf AF2200 plate reader (Eppendorf UK Ltd, United Kingdom) and diluted as appropriate. 

For both sequencing providers, DNA libraries were prepared for sequencing using Rapid Barcoding Kit 96 V14 (Oxford Nanopore Technologies, Oxford, UK) according to the manufacturer’s recommendations. Briefly, sequencing libraries were generated using a transposase, which simultaneously cleaves template molecules and attaches barcoded tags to the cleaved ends. The barcoded samples were then pooled (96-plexed) before solid-phase reversible immobilisation-cleaning and addition of Rapid Adapters to the tagged ends. The library pools were loaded onto ONT PromethION flow cells (R10 [M Version]) – one 96-plex pool per flow cell – and sequenced on a PromethION P2 Solo for 72 h according to the manufacturer’s instructions (GENEWIZ), or were loaded onto a FLO-PRO114M (R.10.4.1) flow cell and sequenced according to the manufacturer’s instructions (MicrobesNG).

## Bioinformatic analysis
Computational analysis was performed on a virtual machine in the Oracle Cloud Infrastructure. POD5 files were basecalled and demultiplexed using Dorado<sup>1</sup> v5.0.0 (super-high accuracy 5mCG, 5hmCG and 6mA methylation aware simplex DNA model).  long-read-only assembly and annotation was implemented in a Nextflow<sup>2</sup> v24.04.3.5916 pipeline<sup>3</sup>, that has previously been shown to produce Enterobacterales assemblies of superior completeness and comparable accuracy to hybrid assembly<sup>4</sup>. All bioinformatic tools were run using default settings unless otherwise specified. Raw-read quality was evaluated with SeqKit<sup>5</sup> v2.9.0. Long-reads were subsampled to generate 4 independent sets of 60× reads using the built-in subsampling and genome size estimation scripts from Autocycler<sup>6,7</sup> v0.2.1. Genome assembly was done using Autocycler<sup>6,7</sup> v0.2.1 using 20 input assemblies (4x each of Canu<sup>8</sup> v2.2, Flye<sup>9</sup>, Raven<sup>10</sup> v1.8.3, Miniasm<sup>11</sup> v0.3 and Plassembler<sup>12</sup>). A single 60× subsampling depth was used based on findings of previous benchmarking work<sup>13-15</sup> that the assembly performance of Flye is comparable in the range 60–100× read depth. The lower end of this range, 60×, was selected to minimise compute use. Assemblies were polished using Medaka<sup>16</sup> v2.0.1 (`--bacteria` flag) using un-subsampled long-reads<sup>4</sup>. Assembly quality was assessed for completeness (>99%; longest contig at least 4Mb) and contamination (<5%) using CheckM2<sup>17</sup> and Autocycler summary outputs.


## Annotation
### Bacterial lineages
Polished assemblies passing quality control were grouped using two methods: 1) multi-locus sequence type (MLST) and 2) fastBAPS clusters. MLSTs, and clonal complex/clonal group where appropriate, were annotated using Kleborate<sup>18</sup> v3.1.2. fastBAPS<sup>19</sup> v1.0.8 clusters were generated by annotating assemblies with Bakta<sup>20</sup> v1.10.4, then generating core genome alignments using Panaroo<sup>21</sup> v1.6.0 (95% presence threshold used to define 'core' genes) and extracting variable sites from the core genome alignment using snp-sites<sup>22</sup> v2.5.1. Separate core genome alignments were generated for each species complex (*E. coli*, *Klebsiella penumoniae*, *K. oxyotca*, *K. aerogenes*). Symmetric optimisation of hyperparameters was used, and the third level of clustering (level 3) was selected for each species complex, as this gave the most appropriate lineage resolution by partitionings large MLSTs (e.g. *E. coli* ST131) while minimising the number of singletons.

### Plasmids
Plasmids were reconstructed and annotated using MOB-suite<sup>23</sup> v3.1.9 (`mob_recon` and `mob_typer`). Plasmids were clustered into communities and subcommunities using pling<sup>24</sup> v2.0.0 using a default containment threshold of 0.5 for comunities and a default DCJ-Indel distance of 4 for subcommunities.

### Antimicrobial resistance (AMR) genes
AMR genes were annotated using AMRFinderPlus<sup>25</sup> v4.2.27, using the ```Echerichia``` or ```Klebsiella``` ```--species``` flag based on the species assigned by Kraken2<sup>26</sup> v2.1.3.

## Bayesian Bootstrapping
Bayesian bootstrapping<sup>27</sup> was used to estimate the posterior frequency distributions of genetic features (MLSTs, fastBAPS clusters, plasmid subcommunities, and AMR genes). with uncertainty and allowing for novel (unseen) features. Briefly, 

A power calculaiton (Wohl et al., 2023) was applied to posterior Bayesian estiamtes of frequency distributions, using equation ...m to estimate the sample size required to detect isolates with a feature at least as frequent as f, which would represent a certain sample coverage at a certain certainty level, p.

Regional estimates were obtained using hierarhical bayesian modelling, where each region varied from the global parameters by a region-specific shrinkage rameter, tau. Region specific shrinkage was selected as this outperformed universal shrinkage for all genomic units using k-fold cross-validation[]. Dirichlet distributions with multinomial priors were used for MLSTs and fastBAPS clusters, and Binomial distributions with Beta distribution priors were used for sub-isolate-level features such as plasmids and AMR genes.

## Statistical analyses and visualisation
Statistical analyses and visulaisations were conducted in R[] v4.5.2.

# References
1. Oxford Nanopore Technologies. Dorado v0.9 2024 [Available from: https://github.com/nanoporetech/dorado?tab=readme-ov-file#alignment].
2. Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nature Biotechnology. 2017;35(4):316-9.10.1038/nbt.3820
3. Nagy D. Autocycler ONT long-read assembly pipeline. 2025.
4. Nagy D, Pennetta V, Rodger G, Hopkins K, Jones CR, Consortium TN, et al. Nanopore long-read-only genome assembly of clinical Enterobacterales isolates is complete and accurate. Microbial Genomics. 2026;12(2).https://doi.org/10.1099/mgen.0.001631
5. Shen W, Le S, Li Y, Hu F. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE. 2016;11(10):e0163962.10.1371/journal.pone.0163962
6. Wick RR, Howden BP, Stinear TP. Autocycler: long-read consensus assembly for bacterial genomes. Bioinformatics. 2025;41(9).10.1093/bioinformatics/btaf474
7. Wick RR. Autocycler. 2025.
8. Koren S, Walenz BP, Berlin K, Miller JR, Bergman NH, Phillippy AM. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome Res. 2017;27(5):722-36.10.1101/gr.215087.116
9. Kolmogorov M, Yuan J, Lin Y, Pevzner P. Assembly of Long Error-Prone Reads Using Repeat Graphs. Nature Biotechnology. 2019.doi:10.1038/s41587-019-0072-8
10. Vaser R, Šikić M. Time- and memory-efficient genome assembly with Raven. Nature Computational Science. 2021;1(5):332-6.10.1038/s43588-021-00073-4
11. Li H. Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences. Bioinformatics. 2016;32(14):2103-10.10.1093/bioinformatics/btw152
12. Bouras G, Sheppard AE, Mallawaarachchi V, Vreugde S. Plassembler: an automated bacterial plasmid assembly tool. Bioinformatics. 2023;39(7).10.1093/bioinformatics/btad409
13. Sanderson ND, Kapel N, Rodger G, Webster H, Lipworth S, Street TL, et al. Comparison of R9.4.1/Kit10 and R10/Kit12 Oxford Nanopore flowcells and chemistries in bacterial genome reconstruction. Microbial Genomics. 2023;9(1).https://doi.org/10.1099/mgen.0.000910
14.	Sanderson ND, Hopkins KMV, Colpus M, Parker M, Lipworth S, Crook D, et al. Evaluation of the accuracy of bacterial genome reconstruction with Oxford Nanopore R10.4.1 long-read-only sequencing. Microb Genom. 2024;10(5).10.1099/mgen.0.001246
15. Wick RR. Oxford Nanopore accuracy vs depth (Ryan Wick's bioinformatics blog) [Internet]2021. [cited 2025]. Available from: https://rrwick.github.io/2021/08/10/accuracy-vs-depth.html.
16. Oxford Nanopore Technologies Ltd. Medaka. 2024.
18. Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nature Methods. 2023;20(8):1203-12.10.1038/s41592-023-01940-w
19. Lam MMC, Wick RR, Watts SC, Cerdeira LT, Wyres KL, Holt KE. A genomic surveillance framework and genotyping tool for Klebsiella pneumoniae and its related species complex. Nature Communications. 2021;12(1):4188.10.1038/s41467-021-24448-3
20. Tonkin-Hill G, Lees JA, Bentley SD, Frost SDW, Corander J. Fast hierarchical Bayesian analysis of population structure. Nucleic Acids Res. 2019;47(11):5539-49.10.1093/nar/gkz361
21. Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics. 2021;7(11).https://doi.org/10.1099/mgen.0.000685
22. Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, Lees JA, et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol. 2020;21(1):180.10.1186/s13059-020-02090-4
23. Page AJ, Taylor B, Delaney AJ, Soares J, Seemann T, Keane JA, et al. SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments. Microbial Genomics. 2016;2(4).https://doi.org/10.1099/mgen.0.000056
24. Robertson J, Nash JHE. MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. Microbial Genomics. 2018;4(8).https://doi.org/10.1099/mgen.0.000206
25. Frolova D, Lima L, Roberts LW, Bohnenkamper L, Wittler R, Stoye J, et al. Applying rearrangement distances to enable plasmid epidemiology with pling. Microb Genom. 2024;10(10).10.1099/mgen.0.001300
26. Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728.10.1038/s41598-021-91456-0
27. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biology. 2019;20(1):257.10.1186/s13059-019-1891-0
28. Rubin DB. The Bayesian Bootstrap. The Annals of Statistics. 1981;9(1):130-4



