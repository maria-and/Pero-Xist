# Pero-Xist
This repository contains analysis scripts and processed data associated with the manuscript:  “Xist expression in male Peromyscus leucopus is associated with restricted chromatin repression and incomplete X-to-autosome dosage compensation.”

The project investigates Xist RNA  in male P. leucopus using bulk RNA-seq, CUT&Tag, RNA FISH, RT-qPCR, and single-nucleus Paired-Tag multiomic profiling. Analyses focus on X-linked gene regulation, H3K27me3 enrichment, X-to-autosome dosage compensation, and Xist-associated gene silencing across sexes and species.

Repository Contents
X:A analyses
CUT&Tag H3K27me3 enrichment analyses
Single-nucleus Paired-Tag RNA/chromatin analyses
Xist-positive versus Xist-negative cell comparisons
Candidate silenced gene identification
Figure generation and statistical analyses
Metadata and processed summary tables

Data Availability

Raw sequencing data are available through NCBI:

Paired-Tag: PRJNA1466601
Bulk CUT&Tag: PRJNA1466562

Bulk RNA-seq datasets used in this study are described in the manuscript Methods section.

Software

Analyses were performed primarily in R using packages including:

Seurat
dplyr
ggplot2
patchwork
scales

Additional processing was performed using:

CLC Genomics Workbench
CellRanger / CellRanger-ATAC
Bowtie2
MACS2
deepTools
Purpose

This repository is intended to provide reproducible code and analysis workflows supporting the findings of the manuscript, including:

persistent Xist expression in male P. leucopus
localized H3K27me3-associated repression
selective silencing of X-linked genes
incomplete X-to-autosome dosage compensation in P. leucopus
