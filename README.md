**Methylation Array Analysis Pipeline**
An R-based workflow for analyzing Illumina 450K methylation array data, including preprocessing, quality control, differential methylation analysis, and visualization.

**The full Analysis report with codes can be found here:**

The https://parulkuls26.github.io/Methylation-array-analysis/methylation_array_analysis-2024.html

**Overview**
This repository provides tools and scripts for comprehensive analysis of DNA methylation data from Illumina HumanMethylation450 BeadChip arrays. The pipeline covers:

1.Data import and preprocessing
2.Quality control and normalization
3.Differential methylation analysis at probe and region levels
4.Visualization and annotation of results

**Requirements**
RStudio
R Version
R >= 4.0.0 recommended

**Dependencies**
Bioconductor packages:
limma
minfi
IlluminaHumanMethylation450kanno.ilmn12.hg19
IlluminaHumanMethylation450kmanifest
missMethyl
minfiData
Gviz
DMRcate

CRAN packages:

knitr
RColorBrewer
stringr
matrixStats

Package Descriptions

Methylation-Specific Packages

minfi — Core package for reading, preprocessing, and analyzing 450K array data

IlluminaHumanMethylation450kanno.ilmn12.hg19 — Annotation package for 450K probes (hg19 genome)

IlluminaHumanMethylation450kmanifest — Provides the Illumina manifest as an R object containing all annotation information for each CpG probe on the 450K array, useful for determining the genomic context of differentially methylated probes

missMethyl — Gene ontology and KEGG pathway analysis designed specifically for methylation data

minfiData — Example datasets for testing and tutorials

DMRcate — Identification of differentially methylated regions (DMRs)

Visualization Packages

RColorBrewer — Color palettes for visualizations
Gviz — Visualization of genomic data and methylation tracks

Statistical Testing

limma — Linear models for testing differential methylation

Utility Packages

knitr — Report generation and reproducible research
matrixStats — Fast row/column operations on matrices used throughout the workflow
stringr — String manipulation functions used in the workflow

Project Structure

R/ — R scripts and functions
data/ — Raw and processed data (IDAT files, sample sheets)
results/ — Analysis outputs (CSV files, DMPs, DMRs)
figures/ — Generated plots and visualizations
reports/ — Knitted HTML/PDF reports
README.md — Project documentation

Analysis Workflow
1. Data Import

Sample Sheet: Read sample metadata containing sample names, chip IDs, positions, and experimental groups
IDAT Files: Load raw red and green channel intensity data into an RGChannelSet object

2. Quality Control

Detection P-values: Calculate detection p-values to assess signal quality at each probe
Mean Detection P-values: Examine per-sample quality to identify failed samples
QC Report: Generate comprehensive PDF report with density plots, control probe performance, and sex prediction

3. Sample and Probe Filtering

Remove Failed Samples: Exclude samples with mean detection p-value > 0.05
Remove Failed Probes: Filter probes with detection p-value > 0.01 in any sample
Remove Sex Chromosome Probes: Exclude chrX and chrY probes when analyzing mixed-sex samples
Remove SNP-Affected Probes: Filter probes where SNPs may affect CpG measurement
Remove Cross-Reactive Probes: Exclude probes that map to multiple genomic locations (Chen et al., 2013)

Filtering Summary:

Initial (after normalization) — 484,535 probes
Failed detection p-value — 977 probes removed → 484,535 remaining
Sex chromosomes — 11,608 probes removed → 472,927 remaining
SNP-affected probes — ~17,000 probes removed → 455,924 remaining
Cross-reactive probes — 26,527 probes removed → 429,397 remaining

4. Preprocessing and Normalization

Quantile Normalization: Apply stratified quantile normalization using preprocessQuantile()
MDS Plots: Visualize sample relationships and identify batch effects or outliers
M-values and Beta Values: Extract M-values for statistical analysis and beta values for visualization

5. Data Exploration

MDS Plots: Examine principal components to identify sources of variation
Density Plots: Compare distributions of beta and M-values across samples
Higher Dimensions: Explore PC3, PC4 to identify additional variation sources

6. Differential Methylation Analysis (Probe-Level)

Design Matrix: Create model including cell type (factor of interest) and individual (blocking factor)
Contrast Matrix: Define pairwise comparisons between groups
Linear Model: Fit limma model to M-values
Empirical Bayes: Apply eBayes moderation for improved statistical power
Results Extraction: Get differentially methylated probes (DMPs) with annotation

DMP Results Include:

logFC: Log2 fold change in M-values
P.Value: Raw p-value
adj.P.Val: FDR-adjusted p-value
Annotation: Gene names, chromosome, CpG island relation

7. Differential Methylation Analysis (Region-Level)

CpG Annotation: Annotate probes with genomic coordinates and statistics using DMRcate
DMR Identification: Identify differentially methylated regions using kernel smoothing
Visualization: Plot DMRs with gene tracks, CpG islands, and sample methylation values

DMR Results Table Columns:

width — Size of the DMR in base pairs
no.cpgs — Number of CpG probes within the DMR
min_smoothed_fdr — Minimum smoothed FDR across CpGs in the region
Stouffer — Stouffer's combined p-value from all CpGs in the region
HMFDR — Harmonic mean FDR
Fisher — Fisher's combined p-value
maxdiff — Maximum beta value difference between groups
meandiff — Mean beta value difference across all CpGs in the DMR
overlapping.genes — Genes that overlap with the DMR

Interpreting meandiff:

Positive values: Hypermethylation in the first group compared to the second
Negative values: Hypomethylation in the first group compared to the second

Output Files

DMPs.csv — Differentially methylated probes with annotation
DMR_results.csv — Differentially methylated regions summary
qcReport.pdf — Quality control report

Key Findings
Differentially Methylated Regions (DMRs)
The analysis identified several significant DMRs between naive and rTreg cell types. Top DMRs include:

TIGIT — 5 CpGs, 597 bp, maxdiff = 54.3% methylation difference
PHF1, CUTA — 14 CpGs, 1903 bp, maxdiff = 54.7% methylation difference (shown in DMR plot on Chromosome 6)
TNFRSF8 — 3 CpGs, 136 bp, maxdiff = 37.4% methylation difference
ZC3H12D — 10 CpGs, 738 bp, maxdiff = 48.8% methylation difference
SLC9A1 — 3 CpGs, 259 bp, maxdiff = -45.3% (hypomethylated)
LAMA3 — 7 CpGs, 402 bp, maxdiff = -38.7% (hypomethylated)
ANKRD11 — 7 CpGs, 513 bp, maxdiff = 32.4% methylation difference
DUSP6 — 11 CpGs, 1405 bp, maxdiff = -29.2% (hypomethylated)
NTRK1, SH2D2A — 6 CpGs, 1859 bp, maxdiff = 38.7% methylation difference
LTBR, SCNN1A — 10 CpGs, 1544 bp, maxdiff = -34.2% (hypomethylated)

Key Observations

All DMRs showed extremely high statistical significance (p-values ranging from 10⁻¹⁰⁵ to 10⁻²⁰⁰)
Both hypermethylation (positive meandiff) and hypomethylation (negative meandiff) patterns were observed
The PHF1/CUTA region on Chromosome 6 showed clear separation between rTreg samples (higher methylation) and other cell types
Inter-individual variation was the largest source of variation but was effectively controlled by including individual as a blocking factor in the statistical model


