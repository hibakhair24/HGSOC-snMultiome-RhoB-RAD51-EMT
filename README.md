HGSOC-snMultiome-RhoB-RAD51-EMT
Single-Nucleus Multi-Omic & Post-Transcriptional Modeling of the RhoB–RAD51–EMT Persister Axis in HGSOC

[![R Version](https://shields.io)](https://r-project.org)
[![Seurat v5](https://shields.io)](https://satijalab.org)
[![License: MIT](https://shields.io)](https://opensource.org)

 Overview
High-Grade Serous Ovarian Cancer (HGSOC) exhibits an ~80% relapse rate following platinum chemotherapy and PARP inhibitors (PARPi), driven by non-genetic cellular reprogramming into a Drug-Tolerant Persister (DTP) state. 

This repository contains a memory-optimized, disk-backed single-nucleus multi-omic (snRNA-seq + snATAC-seq) and post-transcriptional small non-coding RNA (sncRNA) target-inference pipeline investigating how chemotherapy pressure synchronizes RhoB-mediated nuclear RAD51 DNA double-strand break repair with Epithelial-to-Mesenchymal Transition (EMT) plasticity.

Biological Hypothesis & Analytical Architecture

[ Platinum / PARPi Genotoxic Stress ]
                  │
                  ▼
      ┌───────────────────────┐
      │   RhoB Upregulation   │  ◄── Derepression via miR-506 / miR-200c Loss
      └──────────┬─┬──────────┘
       ┌─────────┘ └─────────┐
       ▼ (Cytoplasm)         ▼ (Nucleus)
 Actin Remodeling      Scaffolding of RAD51
       │                     │
       ▼                     ▼
 EMT Activation        Hyper-Active Homologous
 (ZEB1 / Snail)        Recombination Repair (HRR)
       │                     │
       └──────────┬──────────┘
                  ▼
   [ Drug-Tolerant Persister Niche ]

Key Methodological Safeguards
1. Elimination of Enzymatic Dissociation Artifacts: Utilizes flash-frozen single-nucleus multi-omics (GSE293459) to prevent whole-cell digestion stress from artificially elevating baseline RAD51 and DNA damage transcripts.
2. Cell-Cycle Confounding Control: Stratifies all malignant epithelial nuclei (PAX8+ / WT1+) by cell-cycle phase (Seurat::CellCycleScoring) to decouple proliferation-driven RAD51 peaks from quiescent/slow-cycling persister states.
3. Epigenetic & Regulon Verification: Pairs transcript correlation with Signac promoter chromatin accessibility (RHOB, RAD51, ZEB1) and AUCell transcription factor target regulons to overcome post-translational mRNA disconnects.
4. Post-Transcriptional sncRNA Target Inference: Bridges poly(A)-capture limits by projecting validated multiMiR target networks (miR-506-3p, miR-200c-3p, miR-21-5p) onto single nuclei via AUCell.
5. 16 GB RAM Workstation Compatibility: Employs Seurat v5 + BPCells bit-packed binary disk streaming, keeping active memory overhead below 2 GB.



 Repository Structure

HGSOC-RhoB-RAD51-EMT/
├── data/                       Raw GSE293459 H5/fragment files & BPCells binary directories
├── scripts/
│   ├── 01_QC_Malignant.R        BPCells ingestion, QC, PAX8/WT1 malignant isolation & CellCycleScoring
│   ├── 02_Nuclear_Coexpression_and_Regulons.R  Phase-stratified Spearman correlation & EMT scoring
│   ├── 03_snATAC_Chromatin_Accessibility.R     Signac promoter accessibility & coverage tracks
│   └── 04_sncRNA_AUCell.R      multiMiR target retrieval & per-nucleus miRNA derepression scoring
├── results/                    Processed Seurat objects (.rds) & statistical summary tables (.csv)
├── figures/                     High-resolution publication figures (Fig1–Fig4, FigS1)
└── README.md

Data Sources
Single-Nucleus Multi-Omics (snRNA-seq + snATAC-seq): NCBI GEO Accession [GSE293459](https://nih.gov) (Treatment-Naive vs. Post-NACT HGSOC clinical tumors).
Paired Bulk miRNA/mRNA Reference:[TCGA-OV](https://cancer.gov) via Genomic Data Commons (GDC) and multiMiR (miRTarBase / TarBase validated interactions).

 Installation & Quick Start

 1. Hardware & Environment Requirements
OS: Linux, macOS, or Windows via WSL2 (Ubuntu)
RAM: 16 GB minimum (via `BPCells` disk-backed streaming)
Storage: ≥ 100 GB free NVMe SSD space
R Version: ≥ 4.4.0

2. Install Dependencies
install.packages(c("Seurat", "Signac", "tidyverse", "patchwork", "ggpubr", "remotes"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "AUCell", "multiMiR", "GSEABase", "GenomicRanges",
  "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38", "TFBSTools", "JASPAR2020"
))

remotes::install_github("bnprks/BPCells/r")
```

 3. Execute the Pipeline Sequentially
Place the downloaded GSE293459 count matrices, fragment files, and metadata CSV inside `data/`, then run:
```bash
Rscript scripts/01_QC_Malignant.R
Rscript scripts/02_Nuclear_Coexpression_and_Regulons.R
Rscript scripts/03_snATAC_Chromatin_Accessibility.R
Rscript scripts/04_sncRNA_AUCell.R
```

---

 Planned Output Figures
Fig1_Single_Nucleus_Landscape.pdf: UMAP projection of purified malignant HGSOC nuclei stratified by treatment status alongside RHOB, RAD51, and continuous EMT scores.
Fig2_Phase_Stratified_Correlation.pdf: Cell-cycle-stratified (`G1`, `S`, `G2M`) Spearman correlation plots comparing RHOB–RAD51 coupling in High-EMT persister nuclei (`Q4`) across Naive vs. Post-NACT cohorts.
Fig3_snATAC_Chromatin_Coverage.pdf: Signac chromatin accessibility tracks across RHOB, RAD51, and ZEB1 promoter regions (±2 kb TSS).
Fig4_sncRNA_Target_Derepression.pdf: Single-nucleus AUCell target enrichment distributions demonstrating miR-506-3p and miR-200c-3p loss of repression within post-treatment persister niches.

---

Citation & Contact
If you use this pipeline or adapt this analytical framework, please cite this repository and the underlying public data cohorts (`GSE293459`). For collaboration inquiries, please open an issue or contact the repository maintainer.
