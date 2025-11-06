# TASK 7 - Slice genomics windows efficiently

#' TASK 7 - DATA.TABLE VERSION
#'
#' @description
#' This function slices genomic peak data efficiently using `data.table`:
#' it extracts peaks located on chromosome 2 with start positions between
#' 2 Mb and 4 Mb, orders them by descending score, and returns the top 50
#' peaks within that region.
#'
#' @param peaks_file Path to CSV file containing ATAC-seq peaks
#'   (expected columns: chr, start, end, score, etc.)
#' @return list(chr2_window, top50_peaks)
#' @import data.table
#' @export

task7_dt <- function(peaks_file = "atac_peaks.bed.csv") {
# Load library
library(data.table)

# Read the peaks file
peaks_dt <- fread(peaks_file)

# Subset peaks on chr2 in 2–4 Mb region
chr2_window <- peaks_dt[
  chr == "chr2" & start >= 2e6 & start <= 4e6]

# Order by descending score
setorder(chr2_window, -score)

# Select top 50 peaks (safe even if <50)
top50_peaks <- head(chr2_window, 50)

# Return both objects
  return(list(
    chr2_window = chr2_window,
    top50_peaks = top50_peaks
  ))
}

#' TASK 7 - DATA.FRAME VERSION
#'
#' @description
#' This function processes genomic peak data using base R data.frames:
#' it selects peaks on chromosome 2 within the 2–4 Mb window, sorts them
#' by score in descending order, and extracts the top 50 peaks for analysis.
#'
#' @param peaks_file Path to CSV file containing ATAC-seq peaks
#'   (expected columns: chr, start, end, score, etc.)
#' @return list(chr2_window_df, top50_peaks_df)
#' @export

task7_df <- function(peaks_file = "atac_peaks.bed.csv") {
# Load library
library(readr)

# Read the peaks file
peaks_df <- read_csv(peaks_file)

# Show first few lines
head(peaks_df)

# Subset peaks on chr2 in 2–4 Mb region
chr2_window_df <- subset(
    peaks_df,
    chr == "chr2" & start >= 2e6 & start <= 4e6
  )

# Order by descending score
chr2_window_df <- chr2_window_df[order(-chr2_window_df$score), ]

# Select top 50 peaks
top50_peaks_df <- head(chr2_window_df, 50)

return(list(
    chr2_window_df = chr2_window_df,
    top50_peaks_df = top50_peaks_df
  ))
}

