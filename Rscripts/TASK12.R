# TASK 12 - Combine cohorts safely

#DATA.TABLE VERSION 
#' @param cohortA_file Path to cohortA_samples.csv
#' @param cohortB_file Path to cohortB_samples.csv
#' @param counts_file Path to bulk_counts_long.csv (columns: gene, sample_id, count)
#' @return list(combined_samples, merged_counts, top_genes_var, summary_counts, column_check)
#' @export

task12_dt <- function(cohortA_file = "cohortA_samples.csv",
                      cohortB_file = "cohortB_samples.csv",
                      counts_file  = "bulk_counts_long.csv") {
  
# Load library
  library(data.table)

# Read cohort files
A <- fread(cohortA_file)
B <- fread(cohortB_file)
  
# Combine cohorts
combined <- rbindlist(list(A, B), use.names = TRUE, fill = TRUE)
  
# Order by cohort, condition, sample_id
order_cols <- intersect(c("cohort", "condition", "sample_id"), names(combined))
if (length(order_cols) > 0) setorderv(combined, order_cols)
  
# Read counts and join by sample_id
counts <- fread(counts_file)
merged <- merge(counts, combined, by = "sample_id", allow.cartesian = TRUE)
  
# Compute variance per gene (sorted descending)
var_dt <- merged[, .(variance = var(count, na.rm = TRUE)), by = gene]
var_dt <- var_dt[order(-variance, gene)]
top_genes <- var_dt$gene[1:100]
  
# Compute per-cohort/per-condition mean counts
summary_dt <- merged[gene %in% top_genes,
                       .(mean_count = mean(count, na.rm = TRUE)),
                       by = .(gene, cohort, condition)]
setorderv(summary_dt, c("cohort", "condition", "gene"))
  
# Column alignment summary (silent)
alignment_summary <- data.table(
  Status = c("In both", "Only in A", "Only in B"),
  Count  = c(length(intersect(names(A), names(B))),
               length(setdiff(names(A), names(B))),
               length(setdiff(names(B), names(A))))
  )
  
# Return all results
return(list(
    combined_samples = combined,
    merged_counts    = merged,
    top_genes_var    = var_dt,
    summary_counts   = summary_dt,
    column_check     = alignment_summary
  ))
}

# SHOW RESULTS
res_dt12 <- task12_dt("cohortA_samples.csv", "cohortB_samples.csv", "bulk_counts_long.csv")

head(res_dt12$combined_samples, 5)
head(res_dt12$top_genes_var, 5)
head(res_dt12$summary_counts, 5)

# DATA.FRAME VERSION
#' @param cohortA_file Path to cohortA_samples.csv
#' @param cohortB_file Path to cohortB_samples.csv
#' @param counts_file Path to bulk_counts_long.csv (columns: gene, sample_id, count)
#' @return list(combined_samples, merged_counts, top_genes_var, summary_counts_df, column_check)
#' @export

task12_df <- function(cohortA_file = "cohortA_samples.csv",
                      cohortB_file = "cohortB_samples.csv",
                      counts_file  = "bulk_counts_long.csv") {
  
# Read and align cohorts
A <- read.csv(cohortA_file, stringsAsFactors = FALSE)
B <- read.csv(cohortB_file, stringsAsFactors = FALSE)
  
all_cols <- union(names(A), names(B))
for (col in setdiff(all_cols, names(A))) A[[col]] <- NA
for (col in setdiff(all_cols, names(B))) B[[col]] <- NA
  
combined_df <- rbind(A[all_cols], B[all_cols])
  
# Order by cohort, condition, sample_id
  order_cols <- intersect(c("cohort", "condition", "sample_id"), names(combined_df))
  if (length(order_cols) > 0)
    combined_df <- combined_df[do.call(order, combined_df[order_cols]), ]
  
# Read counts and join by sample_id
  counts <- read.csv(counts_file, stringsAsFactors = FALSE)
  merged <- merge(counts, combined_df, by = "sample_id")
  
# Compute variance per gene and select top 100
var_df <- aggregate(count ~ gene, data = merged, FUN = var, na.rm = TRUE)
names(var_df)[2] <- "variance"
var_df <- var_df[order(-var_df$variance, var_df$gene), ]
top_genes <- var_df$gene[1:100]
  
# Compute per-cohort/per-condition mean counts
filtered <- merged[merged$gene %in% top_genes, ]
summary_counts_df <- aggregate(count ~ gene + cohort + condition,
                                 data = filtered,
                                 FUN = function(x) mean(x, na.rm = TRUE))
names(summary_counts_df)[4] <- "mean_count"

summary_counts_df <- summary_counts_df[order(summary_counts_df$cohort,
                                               summary_counts_df$condition,
                                               summary_counts_df$gene), ]
  
# Column alignment summary
  alignment_summary <- data.frame(
    Status = c("In both", "Only in A", "Only in B"),
    Count  = c(length(intersect(names(A), names(B))),
               length(setdiff(names(A), names(B))),
               length(setdiff(names(B), names(A))))
  )
  
return(list(
    combined_samples  = combined_df,
    merged_counts     = merged,
    top_genes_var     = var_df,
    summary_counts_df = summary_counts_df,
    column_check      = alignment_summary
  ))
}

# SHOW RESULTS 
res_df12 <- task12_df("cohortA_samples.csv",
                      "cohortB_samples.csv", "bulk_counts_long.csv")
head(res_df12$combined_samples, 5)
head(res_df12$top_genes_var, 5)
head(res_df12$summary_counts_df, 5)

