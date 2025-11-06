# TASK 4 :Annotate counts with sample and patient info.


# DATA TABLE VERSION 

#' @param counts_file Path to CSV with columns: gene, sample_id, count
#' @param sample_metadata Path to CSV with columns: sample_id, condition, batch, patient_id, timepoint
#' @return list(counts_annotated, counts_patient_total, counts_top10_by_condition)
#' @export

task4_dt <- function(counts_file = "bulk_counts_long.csv",
                     sample_metadata = "sample_metadata.csv") {
# Load library
library(data.table)
  
# Read both files
counts_dt4 <- fread(counts_file)
meta_dt4   <- fread(sample_metadata)
  
# Annotate counts with metadata (join by sample_id)
counts_annotated <- merge(counts_dt4, meta_dt4, by = "sample_id", all.x = TRUE)
  
# Compute per-patient total counts (then sort by patient_id)
counts_patient_total <- counts_annotated[, .(total_count = sum(count)), by = patient_id]
setorder(counts_patient_total, patient_id)   
  
# Find top 10 genes by average count within each condition
counts_top_genes <- counts_annotated[, .(avg_count = mean(count)), by = .(condition, gene)]
  
# Order descending and pick top 10 per condition
counts_top10_by_condition <- counts_top_genes[order(condition, -avg_count), .SD[1:10], by = condition]

# Return as list
  return(list(
    counts_annotated = counts_annotated,
    counts_patient_total = counts_patient_total,
    counts_top10_by_condition = counts_top10_by_condition
  ))
}

# SHOW RESULTS
res_dt4 <- task4_dt("bulk_counts_long.csv", "sample_metadata.csv")
head(res_dt4$counts_patient_total, 5)
head(res_dt4$counts_top10_by_condition, 5)

# DATA.FRAME VERSION 
#' @param counts_file Path to CSV with columns: gene, sample_id, count
#' @param sample_metadata Path to CSV with columns: sample_id, condition, batch,
#'  patient_id, timepoint
#' @return list(counts_annotated, counts_patient_totals, counts_top10_genes)
#' @export

task4_df <- function(counts_file = "bulk_counts_long.csv",
                     sample_metadata = "sample_metadata.csv") {
  
# Read data
counts_df4 <- read.csv(counts_file, stringsAsFactors = FALSE)
meta_df4   <- read.csv(sample_metadata, stringsAsFactors = FALSE)
  
# Annotate counts with metadata (prevent automatic sorting)
counts_annotated <- merge(counts_df4, meta_df4, by = "sample_id", all.x = TRUE, sort = FALSE)
  
# Compute total counts per patient
counts_patient_totals <- aggregate(count ~ patient_id, data = counts_annotated, FUN = sum)
colnames(counts_patient_totals)[2] <- "total_count"
  
# Sort the data
counts_patient_totals <- counts_patient_totals[order(counts_patient_totals$patient_id), ]
  
# Find top 10 genes by average count within each condition
counts_avg_by_condition <- aggregate(count ~ condition + gene, data = counts_annotated, FUN = mean)
colnames(counts_avg_by_condition)[3] <- "avg_count"
  
# Order and extract top 10 per condition
counts_top10_genes <- do.call(rbind, lapply(split(counts_avg_by_condition, counts_avg_by_condition$condition), function(x) {
    x <- x[order(-x$avg_count), ]
    head(x, 10)
  }))
  
rownames(counts_top10_genes) <- NULL   # remove control.* rownames
  
  return(list(
    counts_annotated = counts_annotated,
    counts_patient_totals = counts_patient_totals,
    counts_top10_genes = counts_top10_genes
  ))
}

# SHOW RESULTS
res_df4 <- task4_df("bulk_counts_long.csv", "sample_metadata.csv")
head(res_df4$counts_patient_totals, 5)
head(res_df4$counts_top10_genes, 5)

