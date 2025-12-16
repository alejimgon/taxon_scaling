# Taxon Scaling Analysis

Transcriptomic data normalization and differential expression analysis pipeline for microbial consortia.

## Overview

This repository contains the R script `taxon_scaling_analisys.r` for analyzing gene expression across multiple microbial taxa. The pipeline was used in:

**"Comprehensive analysis of the microbial consortium in the culture of flagellate Monocercomonoides exilis"**  
[Published in *Environmental Microbiome*](https://link.springer.com/article/10.1186/s40793-025-00758-7)

## Features

- **Normalization**: DESeq2-based count normalization
- **Differential Expression**: Multi-contrast analysis across time points (D2, D3, D5)
- **Visualization**: Automated volcano plot generation
- **Organisms analyzed**: 8 microbial taxa (Bacteroides spp., Citrobacter, Fusobacterium, Kerstersia, Monocercomonoides, Parabacteroides)

## Quick Start

### Prerequisites

- **R** ≥ 3.6
- **RStudio** (recommended)

### Installation

1. Clone or download this repository
2. Open the project in RStudio: `File → Open Project`
3. Install required packages:

```r
# Install BiocManager (if needed)
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "biobroom", "EnhancedVolcano"))

# Install CRAN packages
install.packages(c("tidyverse", "stringr", "readxl", "ggplot2", "gridExtra", "here", "plyr"))
```

### Running the Analysis

1. Ensure your data folder structure matches:
```
project_root/
├── featureCounts/
│   ├── samples.txt              # Sample metadata
│   ├── bacteroides_a_*.tsv
│   ├── bacteroides_b_*.tsv
│   ├── bacteroides_c_*.tsv
│   ├── citrobacter_*.tsv
│   ├── fusobacterium_*.tsv
│   ├── kerstersia_*.tsv
│   ├── monocercomonoides_*.tsv
│   └── parabacteroides_*.tsv
└── taxon_scaling_analisys.r
```

2. Open `taxon_scaling_analisys.r` in RStudio
3. Click **Source** (Cmd + Shift + S) or run: `source("taxon_scaling_analisys.r")`

## Input Data

### Required Files

- **featureCounts/samples.txt**: Tab-separated metadata with columns:
  - `Sample ID`: Unique identifier (e.g., RD2R1)
  - `Day`: Time point (D2, D3, D5)
  - Other metadata columns as needed

- **featureCounts/*.txt**: FeatureCounts output files with standardized naming patterns for each organism

### Sample Metadata Format

| Sample ID | Day | ... |
|-----------|-----|-----|
| RD2R1     | D2  | ... |
| RD2R2     | D2  | ... |
| RD3R1     | D3  | ... |

## Output

The script generates:

- **Normalized matrices**: Size factor-normalized count data
- **Differential expression tables**:
  - `comparisons_table.tsv`: All results with annotations
  - `comparisons_table_sorted.tsv`: Filtered DE genes (|log2FC| > 1.0, padj < 0.05)
- **Volcano plots** (PDF):
  - `volcano_plots_all.pdf`: All organisms, all contrasts
  - `volcano_plots_bacteria.pdf`: Bacterial species only
  - `volcano_plots_monocercomonoides.pdf`: Monocercomonoides species only

## Configuration

Edit these parameters in the script to customize the analysis:

```r
FC_THRESHOLD <- 1.0      # Log2 fold-change threshold
P_THRESHOLD <- 0.05      # Adjusted p-value threshold
```

## Methods

- **Normalization**: DESeq2 median-of-ratios method
- **Statistical testing**: Negative binomial generalized linear model (DESeq2)
- **P-value adjustment**: Benjamini-Hochberg FDR correction
- **Visualization**: ggplot2 with log2FC and -log10(p-value) coordinates

## References

Klingenberg H., Meinicke P. 2017
Christel S. et al. 2018

## Citation

If you use this script in your research, please cite:

```
Jiménez-González, A., Treitli, S.C., Peña-Diaz, P. et al. Comprehensive analysis of the microbial consortium in the culture of flagellate Monocercomonoides exilis. Environmental Microbiome 20, 97 (2025). https://doi.org/10.1186/s40793-025-00758-7
```

## Troubleshooting

**Error: "featureCounts/samples.txt not found"**
- Check that your data folder structure matches the expected layout
- Verify the working directory is set correctly

**Error: "there is no package called 'X'"**
- Install the missing package using `install.packages("X")` or `BiocManager::install("X")`

**Script runs slowly**
- DESeq2 normalization is computationally intensive; patience is needed for large datasets
- Consider running on a machine with ≥8GB RAM

## License

Please refer to the original publication for usage terms.

## Contact

For questions or issues, please open an issue on GitHub or contact the repository maintainer.