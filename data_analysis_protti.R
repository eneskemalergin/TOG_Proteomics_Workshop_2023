# Proteomics Data Analysis with protti
#   an R package for proteomics data analysis
# TOG - Proteomics Workshop 2023
# Author: Enes K. Ergin
# Date: 2023-07-13

## Load the packages required

# Install the packages if not installed
if (!requireNamespace("protti", quietly = TRUE))
    install.packages("protti")
if (!requireNamespace("tidyverse", quietly = TRUE))
    install.packages("tidyverse")
library(tidyverse) # Data manipulation
library(protti) # Main proteomics analysis package
# Load functions from utils.R
source("utils.R")

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

# Create a log2 transformed intensity column
protein_data_long$Intensity_log2 <- log2(protein_data_long$Intensity)

# BONUS 1 - Synthetic Data Generation
# Creates 100 Proteins and their replicates with 5 replicates over 2 conditions
data <- protti::create_synthetic_data(
    n_proteins = 100,       # number of proteins
    frac_change = 0.05,
    n_replicates = 5,
    n_conditions = 2,
    method = "effect_random",
    additional_metadata = TRUE
)
# Preview the data
head(data)


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

# Plot again with log2 transformed intensity
qc_intensity_distribution(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity_log2 = Intensity_log2,
    plot_style = "boxplot"
)

#   b. Median Intensity Plot
qc_median_intensities(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity = Intensity_log2
)

# 3. Co-efficients of Variation (CV)
# Within Sample CVs
qc_cvs(
    data = protein_data_long,
    grouping = Protein,
    condition = SampleName,
    intensity = Intensity,
    plot = TRUE,
    plot_style = "boxplot"
)

# Within SampleType CVs (Patient, Cell Line, PDX)
qc_cvs(
    data = protein_data_long,
    grouping = Protein,
    condition = SampleType,
    intensity = Intensity,
    plot = TRUE,
    plot_style = "boxplot"
)

# CVs between Samples in
# Patient+BALL, Patient+ALL,
# Cell Line+BALL, Cell Line+ALL
# PDX+BALL, PDX+ALL
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
    plot_style = "boxplot"
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
    method = "pearson"
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

# 10. Re-visiting Some QC Plots
qc_pca(
  data = protein_data_long,
  sample = Sample,
  grouping = Protein,
  intensity = Intensity_log2,
  condition = SampleType,
  digestion = NULL,
  plot_style = "pca"
)

## Data Normalization

# Median Normalization using protti
protein_data_long <- normalise(
    data = protein_data_long,
    sample = Sample,
    intensity = Intensity_log2,
    method = "median"
) # Adds normalized_"intensity" column

# Plot again with log2 transformed intensity
qc_intensity_distribution(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity_log2 = normalised_intensity_log2,
    plot_style = "boxplot"
)

b1 <- qc_median_intensities(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity = Intensity_log2
) + ggplot2::ggtitle("Before Normalization")

b2 <- qc_median_intensities(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity = normalised_intensity_log2
) + ggplot2::ggtitle("After Normalization")

# Plot the two plots side-by-side
# TODO: Check if the cowplot package is installed on WIN11
cowplot::plot_grid(
    b1, b2,
    ncol = 2
)

# Imputation of Missing Values

# Impute missing values
protein_data_long <- impute_with_downshifted_normal(
    data = protein_data_long,
    intensity_log2 = normalised_intensity_log2,
    prctl = 0.05,
    downshift_mag = 1.5,
    downshift_min = 0.1
)

p1 <- qc_intensity_distribution(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity_log2 = normalised_intensity_log2,
    plot_style = "violin"
) + ggplot2::theme(
    axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1
    )
) + ggplot2::ggtitle(
    "Non-Imputed Intensity Distribution"
)


p2 <- qc_intensity_distribution(
    data = protein_data_long,
    sample = SampleName,
    grouping = Protein,
    intensity_log2 = imputed_intensity_log2,
    plot_style = "violin"
) + ggplot2::theme(
    axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1
    )
) + ggplot2::ggtitle(
    "Imputed Intensity Distribution"
)

cowplot::plot_grid(
    p1, p2,
    ncol = 2,
    align = "v"
)

## Statistical Testing with Weighted Limma package

# Transform the long to wide format
# columns:Sample, rows: Protein, values: imputed_intensity_log2
protein_data_wide <- protein_data_long %>%
    dplyr::select(
        Sample,
        Protein,
        imputed_intensity_log2
    ) %>%
    tidyr::pivot_wider(
        names_from = Sample,
        values_from = imputed_intensity_log2
    ) %>%
    tibble::column_to_rownames(
        var = "Protein"
    )
protein_data_wide_non_imputed <- protein_data_long %>%
    dplyr::select(
        Sample,
        Protein,
        normalised_intensity_log2
    ) %>%
    tidyr::pivot_wider(
        names_from = Sample,
        values_from = normalised_intensity_log2
    ) %>%
    tibble::column_to_rownames(
        var = "Protein"
    )

# Find the indices of the missing values
na_index <- which(is.na(protein_data_wide_non_imputed))

# Initialize the weight matrix
weight_matrix <- matrix(
    data = 1,
    nrow = nrow(protein_data_wide),
    ncol = ncol(protein_data_wide)
)
# Weighting of Missing Values
na_weight <- 0.0001
# Replace the missing values with the na_weight
weight_matrix[na_index] <- na_weight

# Create the design matrix
cur_meta <- as.data.frame(metadata)
rownames(cur_meta) <- cur_meta$ID
cur_meta <- cur_meta[
    which(
        rownames(cur_meta) %in% colnames(protein_data_wide)
    ),
]
# Make sure the columns are in the same order
design_matrix <- model.matrix(
    ~ cur_meta[, "SampleType"]
)

# Variables
pval_thr <- 0.05
log2_fc_thr <- 1.0

# Fit a linear model with weights
fit <- limma::lmFit(
    protein_data_wide,
    design = design_matrix,
    weights = weight_matrix
)
# Run the model
fit_eb <- limma::eBayes(fit)
# Get the log2 fold changes
log2_fc <- fit_eb$coefficients[, 2]
# Get the average intensities
average <- fit_eb$Amean
# Get the pvalues
pvalues <- fit_eb$p.value[, 2]
# Get the adjusted pvalues
adj_pvalues <- p.adjust(
    pvalues,
    method = "fdr"
)

stat_res_df <- data.frame(
    Protein = rownames(protein_data_wide),
    log2_fc = log2_fc,
    average = average,
    pvalues = pvalues,
    adj_pvalues = adj_pvalues
)

# Initialize significance column with no significance
stat_res_df$significance <- "no significance"
# Find Up and Down regulated proteins
stat_res_df$significance[
    (stat_res_df$adj_pvalues < pval_thr) &
    (stat_res_df$log2_fc >= log2_fc_thr)
] <- "Up regulated"
stat_res_df$significance[
    (stat_res_df$adj_pvalues < pval_thr) &
    (stat_res_df$log2_fc <= -log2_fc_thr)
] <- "Down regulated"
# Make the significance column a factor column
stat_res_df$significance <- as.factor(stat_res_df$significance)

view(stat_res_df)

ggplot2::ggplot(
        stat_res_df,
        ggplot2::aes(
            x = log2_fc,
            y = -log10(adj_pvalues),
            color = significance,
            alpha = significance
        )
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(
        xintercept = log2_fc_thr,
        linetype = "dashed",
        color = "darkgrey"
    ) +
    ggplot2::geom_vline(
        xintercept = -log2_fc_thr,
        linetype = "dashed",
        color = "darkgrey"
    ) +
    ggplot2::geom_hline(
        yintercept = -log10(pval_thr),
        linetype = "dashed",
        color = "darkgrey"
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "Up regulated" = "#e63946",
            "Down regulated" = "#1d3557",
            "no significance" = "#b1a7a6"
        )
    ) +
    ggplot2::scale_alpha_manual(
        values = c(
            "Up regulated" = 1.0,
            "Down regulated" = 1.0,
            "no significance" = 0.2
        )
    ) +
    ggplot2::ggtitle("") +
    ggplot2::labs(
        x = "log2(Fold-Change)",
        y = "-log10(adjusted p-value)"
    )

# Plot MA
ggplot2::ggplot(
        stat_res_df,
        ggplot2::aes(
            x = average,
            y = log2_fc,
            color = significance,
            alpha = significance
        )
    ) +
    ggplot2::geom_point(
        size = 3.5,
    ) +
    ggplot2::geom_hline(
        yintercept = log2_fc_thr,
        linetype = "dashed",
        color = "darkgrey"
    ) +
    ggplot2::geom_hline(
        yintercept = -log2_fc_thr,
        linetype = "dashed",
        color = "darkgrey"
    ) +
    # Draw a vertical line at the average intensity of 8
    ggplot2::geom_vline(
        xintercept = 8,
        linetype = "dashed",
        color = "#201717"
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "Up regulated" = "#e63946",
            "Down regulated" = "#1d3557",
            "no significance" = "#b1a7a6"
        )
    ) +
    ggplot2::scale_alpha_manual(
        values = c(
            "Up regulated" = 1.0,
            "Down regulated" = 1.0,
            "no significance" = 0.2
        )
    ) +
    ggplot2::ggtitle("") +
    ggplot2::labs(
        x = "Average Intensity",
        y = "log2(Fold-Change)"
    )

# Enrichment Analysis with gProfiler2 package

# Get the list of proteins
# Background is the list of all proteins
background <- rownames(protein_data_wide)
# Initilize the query list
go_list <- list()
# Add the list of up regulated proteins
go_list[["up_regulated"]] <- rownames(
    stat_res_df[
        which(
            stat_res_df$significance == "Up regulated"
        ),
    ]
)
# Add the list of down regulated proteins
go_list[["down_regulated"]] <- rownames(
    stat_res_df[
        which(
            stat_res_df$significance == "Down regulated"
        ),
    ]
)

# Run the enrichment analysis with gProfiler2
gostres <- gprofiler2::gost(
    query = go_list[["up_regulated"]],
    organism = "hsapiens",  # Human
    user_threshold = 0.05,  # Enrichment p-value threshold
    custom_bg = background, # Custom Background
    domain_scope = "custom_annotated",
    sources = c(            # Source of annotations
        "GO:BP",
        "GO:MF",
        "GO:CC",
        "KEGG",
        "REAC"
        # "WP",
        # "HPA",
        # "CORUM",
        # "MIRNA",
        # "TF",
        # "HP"
    ),
    multi_query = FALSE,
    correction_method = "fdr",
)
# Calculate Gene Ratio
gostres$result$GeneRatio <- (
    gostres$result$intersection_size / gostres$result$term_size
)

gprofiler2::gostplot(
    gostres,
    capped = TRUE,
    interactive = TRUE,
)