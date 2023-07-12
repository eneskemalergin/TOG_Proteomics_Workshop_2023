# Proteomics Data Analysis with protti 
#   an R package for proteomics data analysis
# TOG - Proteomics Workshop 2023
# Author: Enes K. Ergin
# Date: 2023-07-13

## Load the packages required

library(tidyverse) # Data manipulation
library(protti) # Main proteomics analysis package

## Prepare the Data for Analysis

# Metadata
metadata <- read_csv("data/meta_data.csv")
head(metadata)
# List unique values in the metadata
unique(metadata$SampleType)
unique(metadata$Disease)
unique(metadata$Stage)
# Find maximum number of Replicates
max(metadata$Replica)

# Peptide-Level Data (Removed to focus only on protein-level data)
# peptide_data <- read_csv("data/peptide_data.csv")

# Protein-Level Data
protein_data <- read_csv("data/protein_data.csv")
head(protein_data)

# Convert the peptide_data to long format
protein_data_long <- tidyr::pivot_longer(
    protein_data,
    cols = -c("PG.ProteinAccessions"),
    names_to = "ID",
    values_to = "Intensity"
)

# Merge the metadata and peptide_data_longer
protein_data_long <- dplyr::left_join(
    protein_data_long,
    metadata,
    by = "ID"
)
# Rename the columns
protein_data_long <- dplyr::rename(
    protein_data_long,
    "Protein" = "PG.ProteinAccessions",
    "Sample" = "ID"
)

# BONUS 1 - Synthetic Data Generation
# Creates 100 Proteins and their replicates with 5 replicates over 2 conditions
# data <- protti::create_synthetic_data(
#     n_proteins = 100,       # number of proteins 
#     frac_change = 0.05,
#     n_replicates = 5,
#     n_conditions = 2,
#     method = "effect_random",
#     additional_metadata = TRUE
# )
# # Preview the data
# head(data)


## Quality Checks and Data Filtering

# 1. Number of Identified Proteins per Samples
plot <- qc_ids(
    data = protein_data_long,
    sample = SampleName, 
    grouping = Protein,
    condition = SampleType,
    intensity = Intensity,
    plot = TRUE
)
# Visualize the plot
plot

# Sort the bars by the number of proteins
plot + ggplot2::coord_flip()

# Rotate the x-axis labels
plot + ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
)

# 2. Intensity Distribution and Median Intensity Plots
#   a. Intensity Distribution
qc_intensity_distribution(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity_log2 = Intensity,
    plot_style = "boxplot"
)
# Create a log2 transformed intensity column
protein_data_long$Intensity_log2 <- log2(protein_data_long$Intensity)
# Plot again with log2 transformed intensity
qc_intensity_distribution(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity_log2 = Intensity_log2,
    plot_style = "boxplot"
)

#   b. Median Intensity Plot
# DESCRIPTION: Plot the median intensity of each protein across samples as a lineplot
qc_median_intensities(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity = Intensity_log2
)

# 3. Co-efficients of Variation (CV)
# protti::qc_cvs() # Return DF or Plot (boxplot, violin, or density)

# Within Sample CVs
qc_cvs(
    data = protein_data_long,
    grouping = Protein,
    condition = SampleName,
    intensity = Intensity,
    plot = TRUE,
    plot_style="boxplot"
)

# Within SampleType CVs (Patient, Cell Line, PDX)
qc_cvs(
    data = protein_data_long,
    grouping = Protein,
    condition = SampleType,
    intensity = Intensity,
    plot = TRUE,
    plot_style="boxplot"
)

# CVs between Samples in Patient+BALL, Patient+ALL, Cell Line+BALL, Cell Line+ALL, PDX+BALL, PDX+ALL
# Combine CVs between SampleType+Disease
protein_data_long$SampleType_Disease <- paste(
    protein_data_long$SampleType,
    protein_data_long$Disease,
    sep = "_"
)

qc_cvs(
    data = protein_data_long,
    grouping = Protein,
    condition = SampleType_Disease,
    intensity = Intensity,
    plot = TRUE,
    plot_style="boxplot"
)

# 4. Data Completeness
qc_data_completeness(
    data = protein_data_long,
    sample = Sample,
    grouping = Protein,
    intensity = Intensity_log2,
    plot = TRUE
)

# 5. Sample Correlation
# NOTE: Requires pheatmap & seriation packages
qc_sample_correlation(
    data = protein_data_long,
    sample = Sample,
    grouping = Protein,
    intensity = Intensity_log2,
    condition = SampleType,
    interactive = FALSE, 
    method = "pearson",
)

# 6. Principal Component Analysis (PCA)

qc_pca(
  data = protein_data_long,
  sample = Sample,
  grouping = Protein,
  intensity = Intensity_log2,
  condition = SampleType,
  digestion = NULL,
  plot_style = "scree"
)

qc_pca(
  data = protein_data_long,
  sample = Sample,
  grouping = Protein,
  intensity = Intensity_log2,
  condition = SampleType,
  digestion = NULL,
  plot_style = "pca"
)

# 7. Ranked Intensity Distribution (Protein-Rank Plot)
qc_ranked_intensities(
  data = protein_data_long,
  sample = Sample,
  grouping = Protein,
  intensity_log2 = Intensity_log2,
  plot = TRUE,
  y_axis_transformation = "log2"
)

# 8. Removing Problematic Samples 
# Remove instances where SampleName == "NP21"
protein_data_long <- dplyr::filter(
    protein_data_long,
    SampleName != "NP21"
)

# 9. Removing Highly Missing Proteins
# Remove proteins with missing values in more than 75% of samples
# TODO: This is a custom solution, WIP


