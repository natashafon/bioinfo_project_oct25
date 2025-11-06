# TASK 8: Multi-column operations per group

#' TASK 8 - DATA.TABLE VERSION
#'
#' @description
#' This function performs multi-column summary operations per group using `data.table`:
#' it merges counts and metadata by `sample_id`, computes mean, median, Q1, and Q3
#' per gene and condition, reshapes the result to wide format, and filters genes
#' where the treated mean is at least twice the control mean.
#'
#' @param counts_file Path to CSV with columns containing gene, sample_id, count
#' @param meta_file Path to CSV with columns containing sample_id, condition, batch, patient_id, and timepoint
#' @return list(statsum_dt8, wide_dt8, filtered_dt8)
#' @export

task8_dt <- function(counts_file = "bulk_counts_long.csv",
                     meta_file   = "sample_metadata.csv") {

# Load library
library(data.table)

# Read the CSV files
counts_dt8 <- fread(counts_file)
meta_dt8   <- fread(meta_file)

# Merge counts with metadata by sample_id
merged_dt8 <- merge(counts_dt8, meta_dt8, by = "sample_id", sort = FALSE)

# Compute per-condition robust summary stats for each gene (mean, median, Q1, Q3)
statsum_dt8 <- merged_dt8[, .(
mean_count   = mean(count, na.rm = TRUE),
median_count = median(count, na.rm = TRUE),
    Q1           = quantile(count, 0.25, na.rm = TRUE),
    Q3           = quantile(count, 0.75, na.rm = TRUE)
  ), by = .(gene, condition)]

# View summary
print(head(statsum_dt8))

# Reshape to wide format (one row per gene, columns for each condition)
wide_dt8 <- dcast(statsum_dt8, gene ~ condition, value.var = "mean_count")

# Filter by keeping only genes where treated mean ≥ 2 × control mean
filtered_dt8 <- wide_dt8[
  !is.na(.SD$treated) &
    !is.na(.SD$control) &
    .SD$treated >= 2 * .SD$control
]

# Return results
return(list(
  statsum_dt8  = statsum_dt8,
  wide_dt8     = wide_dt8,
  filtered_dt8 = filtered_dt8
  ))
}

#' TASK 8 - DATA.FRAME VERSION
#'
#' @description
#' This function performs grouped summary statistics using base R and tidyverse tools:
#' it merges counts and metadata, calculates mean, median, Q1, and Q3 per gene and condition,
#' converts to wide format, and selects genes where the treated mean is at least
#' two times higher than the control mean.
#'
#' @param counts_file Path to CSV with columns containing gene, sample_id, count
#' @param meta_file Path to CSV with columns containing sample_id, condition, batch, patient_id, and timepoint
#' @return list(statsum_df8, wide_df8, filtered_df8)
#' @export

task8_df <- function(counts_file = "bulk_counts_long.csv",
                     meta_file   = "sample_metadata.csv") {

# Load libraries
library(dplyr)
library(tidyr)
library(readr)

# Read CSV files
counts_df8 <- read_csv(counts_file)
meta_df8   <- read_csv(meta_file)

# Merge counts with metadata by sample_id
merged_df8 <- merge(counts_df8, meta_df8, by = "sample_id", sort = FALSE)

# Compute per-condition summary stats for each gene
statsum_df8 <- merged_df8 %>%
group_by(gene, condition) %>%
    summarise(
      mean_count   = mean(count, na.rm = TRUE),
      median_count = median(count, na.rm = TRUE),
      Q1           = quantile(count, 0.25, na.rm = TRUE),
      Q3           = quantile(count, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

# View summary
print(head(statsum_df8))

# Reshape from long → wide so we can compare control vs treated
wide_df8 <- statsum_df8 %>%
  select(gene, condition, mean_count) %>%
  pivot_wider(names_from = condition, values_from = mean_count)

# Filter genes where treated mean ≥ 2 × control mean
filtered_df8 <- wide_df8 %>%
filter(!is.na(treated) & !is.na(control) & treated >= 2 * control)

# Return results
return(list(
statsum_df8  = statsum_df8,
  wide_df8     = wide_df8,
   filtered_df8 = filtered_df8
  ))
}





