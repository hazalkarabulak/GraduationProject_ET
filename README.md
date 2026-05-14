# Essential Tremor — GWAS-SMR-TWAS Analysis Pipeline

## Overview
This repository contains the full bioinformatics pipeline for identifying
genes associated with Essential Tremor (ET) by integrating GWAS summary
statistics with brain eQTL, mQTL, and transcriptomic data.

## Data Sources
| Data | Source | N |
|---|---|---|
| GWAS | DECODE (Skuladottir et al. 2024) | 16,480 cases / 1,936,173 controls |
| eQTL | BrainMeta cis-eQTL | N=2,865 |
| mQTL | Brain-mMeta | N~1,160 |
| TWAS weights | GTEx v8 Brain Cortex | 1,049 genes |
| LD Reference | 1000 Genomes EUR | N=503 |

## Analysis Pipeline

### 1. GWAS Quality Control
- Raw SNPs: 10,301,333
- After QC: 5,661,054 SNPs
- Filters: allele frequency, multi-allelic SNP removal, strand ambiguity
- Liftover: hg38 → hg19 (CrossMap v0.7.3)

### 2. eQTL-SMR (Summary-based Mendelian Randomization)
- Tool: SMR v1.4.0
- Reference: 1000 Genomes EUR LD panel
- Input: 13,497 gene-SNP pairs (BrainMeta)
- Nominal significant (p < 0.05): 967 genes
- After FDR < 0.05 correction: 6 genes
- After HEIDI test (p > 0.05, removes pleiotropy): **5 genes**
  - KRTCAP3 (Chr2p23.3) — positive effect (b=+0.324)
  - IFT172 (Chr2p23.3) — positive effect (b=+0.265)
  - NRBP1 (Chr2p23.3) — negative effect (b=-0.424)
  - LINC00323 (Chr21) — negative effect (b=-0.172)
  - AC107871.2 (Chr15) — negative effect (b=-0.056)
- SNX17 excluded: HEIDI p=0.005 (pleiotropy detected)

### 3. mQTL-SMR
- Tool: SMR v1.4.0
- Reference: Brain-mMeta methylation QTL
- Significant CpG probes (FDR < 0.05): **5,353**
- Strongest signal: cg09231418 (Chr19, p_FDR=8.97e-142)

### 4. TWAS (Transcriptome-Wide Association Study)
- Tool: FUSION
- Weights: GTEx v8 Brain Cortex (1,049 genes)
- Significant genes (Bonferroni p < 4.77e-05): **858 genes**
- KRTCAP3 significant in both SMR and TWAS (TWAS.Z=-40.943)

### 5. Functional Enrichment
- Tool: FUMA GENE2FUNC
- Input: TWAS significant genes
- Analysis: Gene set enrichment, tissue expression, pathway analysis

## Key Findings
- **KRTCAP3** (Chr2p23.3): identified as the top candidate gene,
  significant in both eQTL-SMR and TWAS analyses
- Chr2p23.3 locus contains multiple co-localized genes (KRTCAP3, IFT172, NRBP1)
- 5,353 CpG probes show methylation-mediated effects on ET risk

## Tools & Software
| Tool | Version | Purpose |
|---|---|---|
| SMR | v1.4.0 | eQTL-SMR and mQTL-SMR |
| FUSION | latest | TWAS |
| CrossMap | v0.7.3 | Liftover hg38→hg19 |
| PLINK | v1.9 | LD calculations |
| R | v4.3.2 | Visualization |
| FUMA | web | Functional enrichment |

## Repository Structure
cd ~/EssentialTremor_Proje/
git add README.md
git commit -m "Update README: add full pipeline description, TWAS and mQTL results"
git push origin main
cat > ~/EssentialTremor_Proje/README.md << 'EOF'
# Essential Tremor — GWAS-SMR-TWAS Analysis Pipeline

## Overview
This repository contains the full bioinformatics pipeline for identifying
genes associated with Essential Tremor (ET) by integrating GWAS summary
statistics with brain eQTL, mQTL, and transcriptomic data.

## Data Sources
| Data | Source | N |
|---|---|---|
| GWAS | DECODE (Skuladottir et al. 2024) | 16,480 cases / 1,936,173 controls |
| eQTL | BrainMeta cis-eQTL | N=2,865 |
| mQTL | Brain-mMeta | N~1,160 |
| TWAS weights | GTEx v8 Brain Cortex | 1,049 genes |
| LD Reference | 1000 Genomes EUR | N=503 |

## Analysis Pipeline

### 1. GWAS Quality Control
- Raw SNPs: 10,301,333
- After QC: 5,661,054 SNPs
- Filters: allele frequency, multi-allelic SNP removal, strand ambiguity
- Liftover: hg38 to hg19 (CrossMap v0.7.3)

### 2. eQTL-SMR
- Tool: SMR v1.4.0
- Input: 13,497 gene-SNP pairs (BrainMeta)
- Nominal significant (p < 0.05): 967 genes
- After FDR < 0.05 and HEIDI > 0.05: 5 genes
  - KRTCAP3 (Chr2p23.3) b=+0.324
  - IFT172 (Chr2p23.3) b=+0.265
  - NRBP1 (Chr2p23.3) b=-0.424
  - LINC00323 (Chr21) b=-0.172
  - AC107871.2 (Chr15) b=-0.056
- SNX17 excluded: HEIDI p=0.005 (pleiotropy)

### 3. mQTL-SMR
- Tool: SMR v1.4.0
- Significant CpG probes (FDR < 0.05): 5,353
- Strongest: cg09231418 (Chr19, p_FDR=8.97e-142)

### 4. TWAS
- Tool: FUSION
- Weights: GTEx v8 Brain Cortex (1,049 genes)
- Significant (Bonferroni p < 4.77e-05): 858 genes
- KRTCAP3 significant in both SMR and TWAS (TWAS.Z=-40.943)

### 5. Functional Enrichment
- Tool: FUMA GENE2FUNC
- Input: TWAS significant genes

## Key Findings
- KRTCAP3 (Chr2p23.3): top candidate, significant in both eQTL-SMR and TWAS
- Chr2p23.3 locus contains co-localized genes: KRTCAP3, IFT172, NRBP1
- 5,353 CpG probes show methylation-mediated effects on ET risk

## Tools
| Tool | Version | Purpose |
|---|---|---|
| SMR | v1.4.0 | eQTL-SMR and mQTL-SMR |
| FUSION | latest | TWAS |
| CrossMap | v0.7.3 | Liftover |
| PLINK | v1.9 | LD calculations |
| R | v4.3.2 | Visualization |
| FUMA | web | Functional enrichment |

## Repository Structure
scripts/ - R visualization scripts
results/ - Key output files and figures
data/    - Intermediate processed files

## HPC
Analyses performed on TRUBA HPC cluster (ARF partition).
