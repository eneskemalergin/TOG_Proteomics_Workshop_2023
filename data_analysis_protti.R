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

# Co-efficients of Variation (CV)

