# Group 7 – DNA Methylation Analysis Project

[![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Illumina 450K](https://img.shields.io/badge/Illumina_Infinium_450K-FF6600?logo=dna&logoColor=white)](https://emea.illumina.com/techniques/microarrays/methylation-arrays.html)
[![DNA Methylation](https://img.shields.io/badge/DNA_Methylation-800080?logo=helix&logoColor=white)](https://www.illumina.com/techniques/multiomics/epigenetics/dna-methylation-analysis.html)
[![minfi](https://img.shields.io/badge/minfi-B200ED?logo=bioconductor&logoColor=white)](https://bioconductor.org/packages/release/bioc/html/minfi.html)
[![preprocessQuantile](https://img.shields.io/badge/preprocessQuantile-4682B4?logo=bioconductor&logoColor=white)](https://bioconductor.org/packages/release/bioc/html/minfi.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Table of Contents

- [Project Overview](#project-overview)
- [Assigned Parameters](#assigned-parameters)
- [Tools and Technologies](#tools-and-technologies)
- [Repository Structure](#repository-structure)
- [Workflow Summary](#workflow-summary)
- [Key Results](#key-results)
- [How to Reproduce](#how-to-reproduce)
- [Resources and References](#resources-and-references)
- [License](#license)
- [Contact](#contact)

---

## Project Overview

This repository contains the DNA methylation analysis developed by **Group 7** for the *DNA/RNA Dynamics* course (Module 2), MSc in Bioinformatics – University of Bologna (Academic Year 2025–2026).

The project investigates genome-wide CpG methylation patterns using data from the Illumina HumanMethylation450K BeadChip, with the goal of identifying differentially methylated positions (DMPs) between a control group (CTRL) and a disease group (DIS). The analysis is performed entirely in R using Bioconductor packages and integrates raw-data import, quality control, normalisation, principal component analysis, statistical testing and biological interpretation.

Designed as both a scientific case study and an educational exercise, the repository provides a reproducible framework for exploring epigenetic variation and its potential role in disease mechanisms.

---

## Assigned Parameters

| Parameter                   | Value                |
| --------------------------- | -------------------- |
| Group ID                    | 7                    |
| Probe Address               | 54740449             |
| Detection p-value cut-off   | 0.05                 |
| Normalization method        | preprocessQuantile   |
| Differential analysis test  | t-test               |

> Address `54740449` corresponds to **AddressB of cg21802055** — a Type I probe measured in the Red channel.

---

## Tools and Technologies

- **Language**: R
- **Platform**: Illumina HumanMethylation450 BeadChip
- **Core packages**: `minfi`, `BiocManager`, `gplots`, `qqman`
- **Annotation**: `Illumina450Manifest_clean.RData`

---

## Repository Structure

### `/Input_data/` — Raw Data & Metadata

- `*.idat` (16 files): Red and green channel raw intensity files (2 per sample × 8 samples).
- `SampleSheet_Report_II.csv`: Sample metadata — IDs, experimental groups, Sentrix ID/position, and `.idat` file references.

> **Not tracked in git.** IDAT files and `.RData` objects are excluded via `.gitignore` (see [How to Reproduce](#how-to-reproduce)). Place them here locally before running the pipeline.

### `/scripts/` — Main Pipeline & Report

- `DRD_Group7_script.R`: Standalone R script with the full analysis code, runnable end to end.
- `DRD_Group7_report.docx` *(or `.html` / `.pdf`)*: The final written report with embedded figures, captions and references.

### `/outputs/` — Results & Visualizations

**QC & Intensity**
- `qc_plot.png`: Raw MSet quality-control plot (log median methylated vs unmethylated).
- `beta_m_values.png`: Distribution of Beta and M values (CTRL vs DIS).
- `Raw_normalised_beta.png`: Raw vs `preprocessQuantile`-normalised values (mean, SD, density, boxplot).
- `controlStripPlot.png`: Negative-control background assessment.

**PCA**
- `PCA_groups.png`: PCA coloured by experimental group (CTRL/DIS).
- `PCA_sex.png`: PCA coloured by sex.
- `scree_plot.png`: Variance explained per principal component.

**Clustering**
- `Complete_linkage_heatmap.png`: Heatmap, complete linkage.
- `Average_linkage_heatmap.png`: Heatmap, average linkage.
- `Single_linkage_heatmap.png`: Heatmap, single linkage (chaining effect).

**Statistics**
- `Histogram_pvalues.png`: Distribution of raw t-test p-values.
- `Boxplot_corrections.png`: Raw vs BH vs Bonferroni adjusted p-values.
- `manhattan_plot.png`: −log₁₀ p-values across genomic positions.
- `volcano_plot.png`: ΔBeta vs −log₁₀ p-value (effect size vs significance).

### `/diagram_workflow/` — Workflow Overview

- `workflow.png`: Diagram of the full pipeline, from raw IDAT import to final visualisations.

### `/report_pipeline/` — Report Guidelines

- `Report_pipeline_FINAL.pdf`: The professor's official instructions and required structure for the final report.

---

## Workflow Summary

```
Raw IDAT files
      │
      ▼
[1] Import  ──►  RGChannelSet (RGset)
      │
      ▼
[2] Extract Red / Green fluorescence
      │
      ▼
[3] Inspect assigned address 54740449  (cg21802055, Type I, Red)
      │
      ▼
[4] MSet.raw  ──►  QC plot
      │
      ▼
[5] Detection p-value filtering  (threshold 0.05)
      │
      ▼
[6] Beta & M values  ──►  group-wise distributions
      │
      ▼
[7] Normalisation: preprocessQuantile  ──►  raw vs normalised comparison
      │
      ▼
[8] PCA  (group / sex)
      │
      ▼
[9] Differential methylation: t-test  ──►  DMPs
      │
      ▼
[10] Multiple-testing correction (BH, Bonferroni)
      │
      ▼
[11] Manhattan + Volcano plots
      │
      ▼
[12] Heatmap of top 100 probes (hierarchical clustering)
```

---

## Key Results

- **Quality control:** The QC plot flagged 5 of 8 samples under the conservative 10.5 cutoff, but all were retained because their detection-p-value failure fractions stayed below 0.6%.
- **PCA:** CTRL and DIS samples do **not** separate globally, indicating that methylation differences are localised rather than genome-wide.
- **Differential methylation:** The t-test identified **16,107** nominally significant DMPs; after **Bonferroni** and **BH** correction, **no probe** remained significant at 0.05 — consistent with the limited power of n = 4 per group.
- **Candidate signals:** Volcano and Manhattan plots highlighted candidate DMPs with biologically meaningful effect sizes (|Δβ| > 0.1). The heatmap of the top 100 probes separated CTRL and DIS under hierarchical clustering — an illustrative result, since those probes were selected on the group contrast (selection bias).

> **Limitations.** The small cohort (n = 4 per group) and the unbalanced sex composition (CTRL: 2M/2F; DIS: 1M/3F) limit statistical power and confound the group contrast. Larger cohorts and models that account for covariates such as sex and batch would be needed to robustly identify and validate DMPs.

---

## How to Reproduce

1. Clone the repository and `cd` into it.
2. Place the raw data inside `Input_data/`:
   - the 16 `*.idat` files,
   - `SampleSheet_Report_II.csv`.
3. Put `Illumina450Manifest_clean.RData` alongside the script (or update the `load()` path in the script).
4. Open `scripts/DRD_Group7_script.R`, set the working directory at the top (`setwd(...)`) to the repo root, and run it.

Required packages (install once):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("minfi")
install.packages(c("gplots", "qqman"))
```

> IDAT files and `.RData` objects are intentionally excluded from version control via `.gitignore` (they are large binaries and re-downloadable). The repo tracks the code, report, figures and documentation only.

---

## Resources and References

- **BiocManager** — installs and manages Bioconductor packages. [Vignette](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
- **minfi** — core package for Illumina 450K/EPIC array analysis (preprocessing, QC, DMP analysis). Aryee MJ *et al.* (2014), *Bioinformatics* 30(10):1363–1369. DOI: [10.1093/bioinformatics/btu049](https://doi.org/10.1093/bioinformatics/btu049)
- **preprocessQuantile** — stratified quantile normalisation for Infinium arrays. Touleimat N, Tost J (2012), *Epigenomics* 4(3):325–341.
- **qqman** — Manhattan and Q-Q plots. Turner SD (2014), bioRxiv. DOI: [10.1101/005165](https://doi.org/10.1101/005165)
- **gplots** — plotting tools including `heatmap.2()`. [CRAN](https://cran.r-project.org/web/packages/gplots/index.html)
- **Beta vs M values** — Du P *et al.* (2010), *BMC Bioinformatics* 11:587.
- **450K platform** — Bibikova M *et al.* (2011), *Genomics* 98(4):288–295.
- **Illumina 450K Product Files** — [support.illumina.com](https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html)

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contact

For questions, feedback, or reproducibility concerns, reach out to the project members:

- Setayesh Khalili — setayesh.khalili@studenti.unibo.it
- Kiana Moradiani — kiana.moradiani@studenti.unibo.it
- Roksana Soheylian — name@studenti.unibo.it
- Sara Yari Mehmandoostsofla — name@studenti.unibo.it
- Sajjad Rezvani Khaledi — name@studenti.unibo.it
- Mozhdeh Asadimonfared — name@studenti.unibo.it
- Mohammad Reza Rezaei — name@studenti.unibo.it

---

> *This repository documents a reproducible methylation analysis workflow combining theoretical insights and practical bioinformatics skills.*
