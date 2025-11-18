---
title: "README"
author: "Zihao.Wang"
date: "2025-11-18"
output: html_document
---

## m6APrediction: Predict m6A RNA Modification Status

m6APrediction is an R package designed to predict the m6A (N6-methyladenosine) modification status of RNA sequences using a trained machine learning model. It supports both single sequence predictions and batch predictions for multiple sequences, leveraging features such as GC content, RNA type, RNA region, and DNA sequence motifs.

### Installation

To install the m6APrediction package from GitHub, use the devtools or remotes package:

```{r}
# Install using devtools
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("Ruza0501/m6APrediction")  

# Or install using remotes
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("Ruza0501/m6APrediction")  
```

### Usage Examples

1.  Predict m6A Status for Multiple Sequences (prediction_multiple) This function takes a data frame of features and returns predictions for all sequences.

```{r}

library(m6APrediction)

# Load the pre-trained model (included in the package)
ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))

# Load example input data (included in the package)
feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))

# Generate predictions with a custom threshold (e.g., 0.6)
predictions <- prediction_multiple(
  ml_fit = ml_fit,
  feature_df = feature_df,
  positive_threshold = 0.6
)

# View results
head(predictions[, c("DNA_5mer", "predicted_m6A_prob", "predicted_m6A_status")])
```

2.  Predict m6A Status for a Single Sequence (prediction_single) This function takes individual feature values for a single sequence and returns its prediction.

```{r}
library(m6APrediction)

# Load the pre-trained model
ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))

# Predict for a single sequence
result <- prediction_single(
  ml_fit = ml_fit,
  gc_content = 0.5,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 120,
  distance_to_junction = 50,
  evolutionary_conservation = 0.7,
  DNA_5mer = "ATGAT",
  positive_threshold = 0.4
)

# View results
print(result)
# Output:
# predicted_m6A_prob  predicted_m6A_status 
#           "0.62"               "Positive" 
```

### Model Performance: ROC and Precision-Recall Curves

The following curves demonstrate the strong performance of the m6A prediction model. The ROC curve illustrates the trade-off between sensitivity and specificity, while the Precision-Recall curve highlights precision across different recall levels. ![Figure 1: ROC curve of the m6A prediction model (AUC = 0.8854). A higher AUC indicates superior ability to distinguish between positive and negative samples.](%22C:/Users/12295/Documents/Practical7/m6APrediction/ROC_PRC.png%22)
