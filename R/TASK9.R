# TASK 9: Go wide‚ to long‚ to wide for downstream plotting.

#' TASK 9 - DATA.TABLE VERSION
#'
#' @description
#' This function reshapes count data between wide and long formats using `data.table`:
#' it converts a wide gene × sample matrix to long format, adds per-sample totals,
#' assigns condition labels, computes mean counts per gene and condition,
#' and returns a wide summary table for downstream plotting.
#'
#' @param wide_file Path to CSV with columns containing gene and sample count columns
#' @return list(counts_long_dt, counts_mean_dt, counts_summary_dt)
#' @import data.table
#' @export

task9_dt <- function(wide_file = "bulk_counts_wide.csv") {

# Load library
library(data.table)

# Read CSV file
counts_wide_dt <- fread(wide_file)   # columns: gene, S01, S02, …

str(counts_wide_dt)

# Identify numeric columns (sample count columns)
numeric_columns <- setdiff(names(counts_wide_dt), "gene")
numeric_columns  # should list S01, S02, …

# Convert from wide to long format
counts_long_dt <- melt(
  counts_wide_dt,
  id.vars = "gene",
  measure.vars = numeric_columns,
  variable.name = "sample_id",
  value.name = "count"
  )

setDT(counts_long_dt)

# Compute per-sample totals
counts_long_dt[, sample_total := sum(count, na.rm = TRUE), by = sample_id]

# Assign condition labels (example rule)
counts_long_dt[, condition := ifelse(
  grepl("treat", sample_id, ignore.case = TRUE),
    "treated", "control"
  )]

# Compute mean count per gene × condition
counts_mean_dt <- counts_long_dt[
    , .(mean_count = mean(count, na.rm = TRUE)),
    by = .(gene, condition)
  ][order(gene, condition)]

# Convert back to wide format (gene × condition)
counts_summary_dt <- dcast(
  counts_mean_dt,
  gene ~ condition,
  value.var = "mean_count"
  )

# Set outputs as data.table
setDT(counts_long_dt)
setDT(counts_mean_dt)
setDT(counts_summary_dt)

 # Return all key tables
list(
  counts_long_dt    = counts_long_dt,
  counts_mean_dt    = counts_mean_dt,
  counts_summary_dt = counts_summary_dt
  )
}

#' TASK 9 - DATA.FRAME VERSION
#'
#' @description
#' This function reshapes count data using base R and `reshape2`:
#' it transforms a wide gene × sample matrix to long format, calculates per-sample totals,
#' assigns conditions, computes mean counts by gene and condition, and converts
#' the result back to wide format for visualization or analysis.
#'
#' @param wide_file Path to CSV with columns containing gene and sample count columns
#' @return list(counts_long_df, counts_mean_df, counts_summary_df)
#' @export

task9_df <- function(wide_file = "bulk_counts_wide.csv") {

# Load library
library(reshape2)

# Read CSV file
counts_wide_df <- read.csv(wide_file, stringsAsFactors = FALSE)
head(counts_wide_df)
str(counts_wide_df)

# Identify numeric columns (sample count columns)
numeric_columns <- setdiff(names(counts_wide_df), "gene")
numeric_columns

# Convert from wide to long format
counts_long_df <- melt(
counts_wide_df,
id.vars = "gene",
measure.vars = numeric_columns,
variable.name = "sample_id",
value.name = "count"
  )

# Compute per-sample totals
counts_long_df$sample_total <- ave(
counts_long_df$count,
counts_long_df$sample_id,
FUN = function(x) sum(x, na.rm = TRUE)
  )


# Assign condition labels (example rule)
counts_long_df$condition <- ifelse(
  grepl("treat", counts_long_df$sample_id, ignore.case = TRUE),
  "treated", "control"
  )


# Compute mean count per gene × condition
counts_mean_df <- aggregate(
  count ~ gene + condition,
  data = counts_long_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)

names(counts_mean_df)[3] <- "mean_count"

counts_mean_df <- counts_mean_df[order(counts_mean_df$gene,
                                       counts_mean_df$condition), ]

# Convert back to wide format (gene × condition)
counts_summary_df <- reshape(
  counts_mean_df,
  timevar = "condition",
  idvar   = "gene",
  direction = "wide"
)

#Modify names
names(counts_summary_df) <- gsub("^mean_count\\.", "",
                                 names(counts_summary_df))

#Sort
col_order <- c("gene", setdiff(sort(names(counts_summary_df)), "gene"))
counts_summary_df <- counts_summary_df[, col_order]

# Return all key tables
list(
  counts_long_df    = counts_long_df,
  counts_mean_df    = counts_mean_df,
  counts_summary_df = counts_summary_df
  )
}

