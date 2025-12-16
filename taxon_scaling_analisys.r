# ============================================================================
# Taxon Scaling Analysis - Transcriptomic Data Normalization
# ============================================================================
# Purpose: Normalize transcriptomic count data using DESeq2 and identify 
#          differentially expressed genes across time points
# 
# Authors: Original analysis adapted from:
#   - Klingenberg H., Meinicke P. 2017
#   - Christel S. et al. 2018
# 
# Date created: 2022-09-7
# Last modified: 2025-12-16
# ============================================================================

library(tidyverse)
library(DESeq2)
library(stringr)
library(biobroom)
library(readxl)
library(ggplot2)
library(EnhancedVolcano)
library(gridExtra)
library(grid)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Define organisms to analyze
organisms <- list(
  BACA = "bacteroides_a",
  BACB = "bacteroides_b",
  BACC = "bacteroides_c",
  CITR = "citrobacter",
  FUS = "fusobacterium",
  KER = "kerstersia",
  MONO = "monocercomonoides",
  PARBA = "parabacteroides"
)

# Define contrasts (day comparisons)
contrasts_list <- list(
  c("condition", "D3", "D2", "D2_D3"),
  c("condition", "D5", "D3", "D3_D5"),
  c("condition", "D5", "D2", "D2_D5")
)

# DE gene thresholds
FC_THRESHOLD <- 1.0
P_THRESHOLD <- 0.05

# ============================================================================
# LOAD DATA
# ============================================================================

message("Loading sample metadata...")
setwd(here::here())

if (!file.exists("featureCounts/samples.txt")) {
  stop("Error: featureCounts/samples.txt not found")
}

samples_Mono_Meta <- read_tsv("featureCounts/samples.txt")
message(sprintf("Loaded %d samples", nrow(samples_Mono_Meta)))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Load and process feature counts for an organism
load_organism_counts <- function(pattern) {
  message(sprintf("  Loading files matching pattern: %s", pattern))
  
  input_files <- list.files("featureCounts/", 
                           pattern = pattern, 
                           full.names = TRUE)
  
  if (length(input_files) == 0) {
    warning(sprintf("No files found for pattern: %s", pattern))
    return(NULL)
  }
  
  input_files %>%
    map(read_tsv, comment = "#", col_types = cols()) %>%
    set_names(str_extract(., "RD\\dR\\d")) %>%
    plyr::join_all() -> featcounts
  
  return(featcounts)
}

#' Rename columns and reshape feature counts
process_organism_counts <- function(featcounts, org_name) {
  message(sprintf("  Processing %s...", org_name))
  
  colnames(featcounts) <- gsub(".+(RD\\dR\\d).+", "\\1", colnames(featcounts))
  colnames(featcounts) <- plyr::mapvalues(
    names(featcounts),
    from = samples_Mono_Meta$`Sample ID`,
    to = as.character(samples_Mono_Meta$Day)
  )
  
  featcounts %>%
    gather(sampleID, readcounts, -Geneid, -Chr, -Start, -End, -Strand, -Length) %>%
    mutate(day = gsub("(D\\d)R\\d", "\\1", sampleID)) -> processed
  
  return(processed)
}

#' Convert to numeric matrix for DESeq2
prepare_count_matrix <- function(featcounts_tidy) {
  count_matrix <- featcounts_tidy %>%
    select(Geneid, sampleID, readcounts) %>%
    spread(sampleID, readcounts)
  
  rownames(count_matrix) <- count_matrix$Geneid
  return(as.matrix(count_matrix[, -1]))
}

#' Normalize counts using DESeq2
DESeq2.norm.mat <- function(Xmat, cond, type) {
  Xmat.col <- ncol(Xmat)
  Xmat.row <- nrow(Xmat)
  colData <- data.frame(condition = cond, type = type)
  Xmat <- round(Xmat)
  storage.mode(Xmat) <- 'integer'
  
  dds <- DESeqDataSetFromMatrix(countData = Xmat, 
                                colData = colData, 
                                design = ~condition)
  colData(dds)$condition <- factor(colData(dds)$condition, 
                                   levels = unique(cond))
  dds <- DESeq(dds, quiet = TRUE)
  
  YMat <- Xmat / rep(dds@colData@listData$sizeFactor, each = Xmat.row)
  return(YMat)
}

#' Run DESeq2 differential expression analysis
DESeq2.result <- function(Xmat, cond, type, CONTRAST) {
  Xmat.col <- ncol(Xmat)
  Xmat.row <- nrow(Xmat)
  Xmat <- round(Xmat)
  storage.mode(Xmat) <- 'integer'
  colData <- data.frame(condition = cond, type = type)
  
  dds <- DESeqDataSetFromMatrix(countData = Xmat, 
                                colData = colData, 
                                design = ~condition)
  
  colData(dds)$condition <- factor(colData(dds)$condition, 
                                   levels = unique(cond))
  
  normFactors <- matrix(1, ncol = Xmat.col, nrow = Xmat.row)
  normalizationFactors(dds) <- normFactors
  dds <- DESeq(dds, quiet = FALSE)
  res <- results(dds, contrast = CONTRAST)
  
  return(res)
}

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

message("Loading and processing feature counts for all organisms...")

# Load all organisms
feature_counts_list <- list()
for (org_code in names(organisms)) {
  pattern <- organisms[[org_code]]
  message(sprintf("Processing %s (%s)", org_code, pattern))
  
  counts <- load_organism_counts(pattern)
  if (!is.null(counts)) {
    processed <- process_organism_counts(counts, org_code)
    feature_counts_list[[org_code]] <- processed
  }
}

# Prepare matrices
message("Preparing count matrices...")
count_matrices <- lapply(feature_counts_list, prepare_count_matrix)

# Setup conditions
cond.vec <- gsub("(D\\d)R\\d", "\\1", samples_Mono_Meta$Day)
type.vec <- rep("pe", length(cond.vec))

# Normalize all matrices
message("Normalizing matrices with DESeq2...")
normalized_matrices <- lapply(count_matrices, 
                              function(mat) DESeq2.norm.mat(mat, cond.vec, type.vec))

# Run DESeq2 for all contrasts
message("Running differential expression analysis...")
results_list <- list()

for (org_code in names(normalized_matrices)) {
  message(sprintf("  Analyzing %s...", org_code))
  
  for (i in seq_along(contrasts_list)) {
    contrast <- contrasts_list[[i]]
    condition_name <- contrast[4]
    
    res <- DESeq2.result(normalized_matrices[[org_code]], 
                        cond.vec, type.vec, 
                        CONTRAST = contrast[1:3]) %>%
      tidy() %>%
      mutate(term = condition_name, org = org_code)
    
    results_list[[paste0(org_code, "_", i)]] <- res
  }
}

# Combine results
message("Combining results...")
master <- bind_rows(results_list)

# Annotate differential expression
master.annot <- master %>%
  rename(estimate = "log2FoldChange") %>%
  mutate(
    diffexpressed = ifelse(
      log2FoldChange > FC_THRESHOLD & p.adjusted < P_THRESHOLD, "UP",
      ifelse(log2FoldChange < -FC_THRESHOLD & p.adjusted < P_THRESHOLD, "DOWN", "NO")
    )
  )

# ============================================================================
# OUTPUT RESULTS
# ============================================================================

message("Writing output tables...")
write.table(master.annot, 
           "featureCounts/comparisons_table.tsv",
           quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

master.annot.sort <- master.annot %>%
  filter(abs(log2FoldChange) > FC_THRESHOLD) %>%
  filter(p.adjusted < P_THRESHOLD) %>%
  arrange(desc(log2FoldChange))

write.table(master.annot.sort,
           "featureCounts/comparisons_table_sorted.tsv", 
           quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

message("Creating visualizations...")

# ============================================================================
# VISUALIZATION - VOLCANO PLOTS
# ============================================================================

# Color scheme for differential expression
volcano_colors <- c("UP" = "red2", "NO" = "grey30", "DOWN" = "royalblue")

#' Create a volcano plot for differential expression results
#'
#' @param results_table Data frame with log2FoldChange, p.adjusted, diffexpressed columns
#' @param title Plot title (e.g., "Organism D3 vs D2")
#' @param fc_threshold Fold change threshold for cutoff line
#' @param p_threshold P-value threshold for cutoff line
#'
#' @return ggplot object
create_volcano_plot <- function(results_table, title, 
                               fc_threshold = 1.0, 
                               p_threshold = 0.05) {
  
  plot <- ggplot(data = results_table, 
                 aes(x = log2FoldChange, y = -log10(p.adjusted), 
                     col = diffexpressed)) +
    labs(title = title, 
         x = "log2(Fold Change)", 
         y = "-log10(adjusted p-value)",
         color = "Expression") +
    geom_point(alpha = 0.7, size = 2.5) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "right") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               col = "red", linetype = 2, alpha = 0.7) +
    geom_hline(yintercept = -log10(p_threshold), 
               col = "red", linetype = 2, alpha = 0.7) +
    scale_color_manual(values = volcano_colors)
  
  return(plot)
}

#' Generate all volcano plots for an organism
#'
#' @param org_code Organism code (e.g., "BACA", "MONO")
#' @param org_name Full organism name for titles
#' @param master_data Master annotation table with all results
#'
#' @return List of 3 ggplot objects (D2_D3, D3_D5, D2_D5)
generate_organism_plots <- function(org_code, org_name, master_data) {
  
  contrasts <- list(
    list(name = "D2_D3", title = paste(org_name, "D3 vs D2")),
    list(name = "D3_D5", title = paste(org_name, "D5 vs D3")),
    list(name = "D2_D5", title = paste(org_name, "D5 vs D2"))
  )
  
  plots <- list()
  
  for (i in seq_along(contrasts)) {
    contrast_name <- contrasts[[i]]$name
    contrast_title <- contrasts[[i]]$title
    
    filtered_data <- master_data %>%
      filter(org == org_code) %>%
      filter(term == contrast_name)
    
    if (nrow(filtered_data) > 0) {
      plots[[contrast_name]] <- create_volcano_plot(
        results_table = filtered_data,
        title = contrast_title
      )
      message(sprintf("  Created plot: %s (%s)", org_code, contrast_name))
    } else {
      warning(sprintf("No data found for %s - %s", org_code, contrast_name))
    }
  }
  
  return(plots)
}

# Generate all plots programmatically
message("Generating volcano plots for all organisms...")

# Define organism display names
organism_names <- list(
  BACA = "Bacteroides sp. A",
  BACB = "Bacteroides sp. B",
  BACC = "Bacteroides sp. C",
  CITR = "Citrobacter sp.",
  FUS = "Fusobacterium sp.",
  KER = "Kerstersia sp.",
  MONO = "Monocercomonoides sp.",
  PARBA = "Parabacteroides sp."
)

# Generate all plots
all_plots <- list()

for (org_code in names(organism_names)) {
  message(sprintf("Processing %s...", org_code))
  
  org_plots <- generate_organism_plots(
    org_code = org_code,
    org_name = organism_names[[org_code]],
    master_data = master.annot
  )
  
  all_plots[[org_code]] <- org_plots
}

# ============================================================================
# SAVE VOLCANO PLOTS TO PDF
# ============================================================================

message("Saving volcano plots to PDF files...")

#' Save volcano plots to PDF
#'
#' @param plots List of ggplot objects
#' @param filename Output PDF filename
#' @param title_prefix Prefix for document title
#' @param width PDF width in inches
#' @param height PDF height in inches
save_volcano_pdf <- function(plots, filename, title_prefix = "", 
                            width = 18, height = 16) {
  
  if (length(plots) == 0) {
    warning(sprintf("No plots to save for %s", filename))
    return(invisible(NULL))
  }
  
  pdf(file = filename, width = width, height = height)
  
  if (!is.null(title_prefix)) {
    grid.newpage()
    grid.text(title_prefix, gp = gpar(fontsize = 16, fontface = "bold"))
  }
  
  do.call(grid.arrange, c(plots, ncol = 3))
  dev.off()
  
  message(sprintf("✓ Saved: %s", filename))
}

# Flatten plots for each category
all_volcano_plots <- unlist(all_plots, recursive = FALSE)
bacteria_plots <- unlist(all_plots[names(all_plots) != "MONO"], recursive = FALSE)
mono_plots <- all_plots$MONO

# Save complete volcano plot collection
save_volcano_pdf(
  plots = all_volcano_plots,
  filename = "volcano_plots_all.pdf",
  title_prefix = "Differential Expression Analysis - All Organisms"
)

# Save bacteria-only plots
save_volcano_pdf(
  plots = bacteria_plots,
  filename = "volcano_plots_bacteria.pdf",
  title_prefix = "Differential Expression Analysis - Bacterial Species"
)

# Save monocercomonoides-only plots
save_volcano_pdf(
  plots = mono_plots,
  filename = "volcano_plots_monocercomonoides.pdf",
  title_prefix = "Differential Expression Analysis - Monocercomonoides sp.",
  width = 10,
  height = 6
)

message("All volcano plots generated successfully!")