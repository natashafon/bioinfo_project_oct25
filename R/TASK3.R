
# TASK 3 — Speed up joins/lookups

#' TASK 3 - DATA.TABLE VERSION
#'
#' @description
#' This function speeds up joins and lookups using `data.table`:
#' it sets a key on `sample_metadata` by `sample_id`, performs an
#' equi-join with the long counts table, adds a secondary index on
#' `(gene, sample_id)`, and benchmarks a subset query before and after
#' indexing to measure speed improvement.
#'
#' @param counts_file Path to CSV with columns: gene, sample_id, count
#' @param sample_metadata Path to CSV with columns: sample_id, condition, batch, patient_id, timepoint
#' @return list(merged_dt3, benchmark3)
#' @export

task3_dt <- function(counts_file = "bulk_counts_long.csv",
                     sample_metadata = "sample_metadata.csv") {

# Load needed libraries
library(data.table)
library(microbenchmark)

# Read the files
counts_dt3 <- fread(counts_file)
meta_dt3   <- fread(sample_metadata)

# Set key and join by sample_id
setkey(meta_dt3, sample_id)
merged_dt3 <- counts_dt3[meta_dt3, on = "sample_id", nomatch = 0L]


# Benchmark before and after index
gene_query   <- merged_dt3$gene[1]
sample_query <- merged_dt3$sample_id[1]

before_time <- microbenchmark(
no_index = merged_dt3[gene == gene_query & sample_id == sample_query],
times = 20L
)

# Add secondary index
setindexv(merged_dt3, c("gene", "sample_id"))
after_time <- microbenchmark(
    with_index = merged_dt3[gene == gene_query & sample_id == sample_query],
    times = 20L
)

# Combine benchmark results
benchmark3 <- rbind(
    data.table(expr = "no_index",  median_ms = median(before_time$time) / 1e6),
    data.table(expr = "with_index", median_ms = median(after_time$time) / 1e6)
  )
  benchmark3[, speedup := round(benchmark3$median_ms[1] / benchmark3$median_ms, 2)]

return(list(
merged_dt3 = merged_dt3,
benchmark3 = benchmark3
  ))
}

# DATA.FRAME VERSION
#' TASK 3 - DATA.FRAME VERSION
#'
#' @description
#' This function performs the same join and lookup benchmark using
#' base R data.frames: it merges count and metadata tables by `sample_id`,
#' runs subset queries on selected `gene` and `sample_id`, and measures
#' execution time for comparison with the `data.table` approach.
#'
#' @param counts_file Path to CSV with columns: gene, sample_id, count
#' @param sample_metadata Path to CSV with columns: sample_id, condition, batch, patient_id, timepoint
#' @return list(merged_df3, benchmark_df3, subset_df3)
#' @export

task3_df <- function(counts_file = "bulk_counts_long.csv",
                     sample_metadata = "sample_metadata.csv") {

# Read data
counts_df3 <- read.csv(counts_file, stringsAsFactors = FALSE)
meta_df3   <- read.csv(sample_metadata, stringsAsFactors = FALSE)


# Merge using base R
merged_df3 <- merge(counts_df3, meta_df3, by = "sample_id", all = FALSE)

# Define query
gene_query   <- merged_df3$gene[1]
sample_query <- merged_df3$sample_id[1]


# Benchmark
benchmark_df3 <- microbenchmark(
  subset_lookup = merged_df3[merged_df3$gene == gene_query &
                               merged_df3$sample_id == sample_query, ],
  times = 100L
)

# Subset
subset_df3 <- merged_df3[merged_df3$gene == gene_query & merged_df3$sample_id == sample_query, ]

return(list(
  merged_df3 = merged_df3,
  benchmark_df3 = benchmark_df3,
  subset_df3 = subset_df3))
}



