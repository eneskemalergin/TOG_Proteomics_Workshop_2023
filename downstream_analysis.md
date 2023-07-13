Hands-on Downstream Analysis of Protein Level Data
================

## Notebook Setup

**Check and Install if Necessary**

``` r
if (!requireNamespace("protti", quietly = TRUE))
    install.packages("protti")
if (!requireNamespace("tidyverse", quietly = TRUE))
    install.packages("tidyverse")
if (!requireNamespace("pheatmap", quietly = TRUE))
    install.packages("pheatmap")
if (!requireNamespace("seriation", quietly = TRUE))
    install.packages("seriation")
if (!requireNamespace("gprofiler2", quietly = TRUE))
    install.packages("gprofiler2")
if (!requireNamespace("limma", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("limma")
}
```

**Load Packages Necessary**

``` r
library(limma) # for stat-test
library(protti) # for QC and process 
library(tidyverse) # data manipulation & visual
library(gprofiler2) # Enrichment of protein sets
source("utils.R")  # Custom functions
```
