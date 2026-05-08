# Essential Tremor GWAS-SMR Analysis Pipeline

## Overview
Integration of GWAS summary statistics with brain eQTL data 
using SMR to identify genes associated with Essential Tremor (ET).

## Data Sources
- **GWAS:** DECODE (Skuladottir et al. 2024) — 16,480 cases, 1,936,173 controls (GRCh38)
- **eQTL:** BrainMeta cis-eQTL (GRCh37)
- **LD Reference:** 1000 Genomes EUR (GRCh37, N=503)

## Pipeline Steps
1. GWAS QC (allele filtering, multi-allelic SNP removal)
2. Liftover hg38 → hg19 (CrossMap)
3. SMR analysis — all 22 chromosomes (for loop)
4. Multiple testing correction (Benjamini-Hochberg FDR)

## Key Findings
Significant genes (FDR < 0.05):
- SNX17, NRBP1, KRTCAP3, IFT172, AC107871.2, LINC00323

## Tools Used
- SMR v1.3.2
- PLINK v1.9
- CrossMap v0.7.3
- R (ggplot2, dplyr)
- TRUBA HPC cluster

## Repository Structure
scripts/   → Analysis scripts (R, bash)
results/   → SMR output files
data/      → Intermediate processed files (raw data not included)
