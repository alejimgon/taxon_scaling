
# Taxon Scaling Analysis

This repository contains the R script `taxon_scaling_analisys.r`, which was used for the normalization and differential expression analysis of transcriptomic data in the study:

[**"Comprehensive analysis of the microbial consortium in the culture of flagellate Monocercomonoides exilis"**](https://link.springer.com/article/10.1186/s40793-025-00758-7)

## Description

The script performs normalization of feature counts and differential expression analysis using DESeq2, generating volcano plots and summary tables for various microbial taxa. It is adapted from Klingenberg H., Meinicke P. (2017) and Christel S. et al. (2018).

## Usage

- Input: FeatureCounts tables and sample metadata.
- Output: Normalized count matrices, differential expression results, volcano plots.

## Requirements

- R (≥ 3.6)
- Packages: tidyverse, DESeq2, stringr, biobroom, readxl, ggplot2, EnhancedVolcano, gridExtra, grid

## Citation

If you use this script, please cite the original paper:

> [**"Comprehensive analysis of the microbial consortium in the culture of flagellate Monocercomonoides exilis"**](https://link.springer.com/article/10.1186/s40793-025-00758-7)
