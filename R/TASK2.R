#TASK 2 : Add QC-style derived columns without copying.

#'TASK 2 DATA.TABLE VERSION
#'
#' @description
#' This function adds QC-style columns to bulk count data:
#' it computes `log2_count = log2(count + 1)` and a logical flag `high`
#' marking rows where the count is above the gene-wise median.
#'
#' @param counts_file path to CSV with columns containing gene, sample_id, count
#' @return data.table with added columns: log2_count and high
#' @import data.table
#' @export

task2_dt <- function(
    counts_file = "bulk_counts_long.csv") {

# Load the library
library(data.table)

# Read the CSV file
counts_dt <- fread(counts_file)

# Adding log2 counts coulumn which is applied to all rows
counts_dt[, log2_count := log2(count + 1)]

# Adding a binary flag high if count >100
counts_dt[, high := count > 100]

# Overwrite 'high' if the count is higher thant the median value for a certain gene
counts_dt[, high := count > median(count), by = gene]

# Return result
return(counts_dt)
}

#' TASK 2 DATA.FRAME VERSION
#'
#' @description
#' This function adds QC-style columns to bulk count data:
#' it computes `log2_count = log2(count + 1)` and a logical flag `high`
#' marking rows where the count is above the gene-wise median.
#' @param counts_file Path to CSV with columns containing gene, sample_id, count
#' @return data.frame with added columns: log2_count and high
#' @export

task2_df <- function(counts_file = "bulk_counts_long.csv") {


#  Read CSV file
counts_df <- read.csv(path, stringsAsFactors = FALSE)

# Add a new column: log2_count
counts_df$log2_count <- log2(counts_df$count + 1)

# Add a column: high = TRUE if count > 100
counts_df$high <- counts_df$count > 100

# Compute median count per gene
med_by_gene <- tapply(counts_df$count, counts_df$gene, median)

# Overwrite 'high' using gene-wise median
counts_df$high <- counts_df$count > med_by_gene[counts_df$gene]

return(counts_df)
}

