
#' TASK 1 DATA.TABLE VERSION
#'
#' @description
#' This function joins count data and metadata by `sample_id`, filters treated samples
#' and genes starting with "GENE_00", and summarizes counts per gene and condition.
#'
#'
#' @param counts_file Path to CSV with columns containing gene, sample_id, count
#' @param meta_file Path to CSV with columns containing sample_id, condition, batch, patient_id, and timepoint
#' @return list(merged, filtered, gene_summary, per_condition_mean)
#' @import data.table
#' @export

task1_dt <- function(counts_file = "bulk_counts_long.csv",
                     meta_file   = "sample_metadata.csv") {
# Load libraries
  library(data.table)

# Read th CSV files
counts <- fread(counts_file)
meta   <- fread(meta_file)

# Ensure sample_id is character
counts[, sample_id := as.character(sample_id)]
meta[,   sample_id := as.character(sample_id)]

# Join on sample_id
merged <- merge(counts, meta, by = "sample_id", sort = FALSE)

# Keep only treated samples and genes starting with "GENE_00"
filtered <- merged[
  condition == "treated" & grepl("^GENE_00", gene)
]


# Mean and median count by gene
gene_summary <- filtered[, .(
  mean_count   = mean(count, na.rm = TRUE),
  median_count = median(count, na.rm = TRUE)
), by = gene][order(gene)]


# Per-condition mean counts by gene
per_condition_mean <- merged[, .(
  mean_count = mean(count, na.rm = TRUE)
), by = .(gene, condition)][order(gene, condition)]


list(
  merged = merged,
  filtered = filtered,
  gene_summary = gene_summary,
  per_condition_mean = per_condition_mean
)
}

#' TASK DATA.FRAME VERSION
#'
#' @description
#' This function reproduces the same operations as `task1_dt()` but uses only base R.
#' @param counts_file Path to CSV with columns containing gene, sample_id, count
#' @param meta_file Path to CSV with columns containing sample_id, condition, batch, patient_id, and timepoint
#' @return list(merged_df, filtered_df, gene_summary_df, per_condition_mean_df)
#' @export

task1_df <- function(counts_file = "bulk_counts_long.csv",
                     meta_file   = "sample_metadata.csv") {

# Read CSV files
counts_df <- read.csv(counts_file, stringsAsFactors = FALSE)
meta_df   <- read.csv(meta_file,  stringsAsFactors = FALSE)

# Ensure sample_id is character
counts_df$sample_id <- as.character(counts_df$sample_id)
meta_df$sample_id   <- as.character(meta_df$sample_id)

# Join on sample_id
merged_df <- merge(counts_df, meta_df, by = "sample_id", sort = FALSE)

# Filter treated + genes starting with "GENE_00"
filtered_df <- merged_df[
  merged_df$condition == "treated" & grepl("^GENE_00", merged_df$gene),
]

# Mean and median by gene (on filtered)
mean_by_gene   <- aggregate(count ~ gene, data = filtered_df, FUN = function(x) mean(x, na.rm = TRUE))
median_by_gene <- aggregate(count ~ gene, data = filtered_df, FUN = function(x) median(x, na.rm = TRUE))

# Merge summaries and align column names/order to data.table version

gene_summary_df <- merge(mean_by_gene, median_by_gene, by = "gene")

names(gene_summary_df) <- c("gene", "mean_count", "median_count")

gene_summary_df <- gene_summary_df[order(gene_summary_df$gene), ]

# Per-condition mean counts by gene (on all merged)
per_condition_mean_df <- aggregate(
  count ~ gene + condition,
  data = merged_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)
per_condition_mean_df <- per_condition_mean_df[order(per_condition_mean_df$gene, per_condition_mean_df$condition), ]

list(
  merged_df = merged_df,
  filtered_df = filtered_df,
  gene_summary_df = gene_summary_df,
  per_condition_mean_df = per_condition_mean_df
)
}






