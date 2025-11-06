
# TASK 5 - Classify values against reference intervals

#' TASK 5 - DATA.TABLE VERSION
#'
#' @description
#' This function classifies clinical lab values against reference intervals:
#' it merges lab results with reference ranges, labels each value as
#' `"normal"` or `"out_of_range"`, and summarizes abnormal rates by
#' patient and by lab using `data.table` operations.
#'
#' @param labs_file Path to CSV with columns: patient_id, time_iso, lab, value
#' @param ref_file  Path to CSV with columns: lab, lower, upper
#' @return list(classified_dt5, abnormal_by_patient_dt5, abnormal_by_lab_dt5)
#' @export

task5_dt <- function(labs_file = "clinical_labs.csv",
                     ref_file  = "lab_reference_ranges.csv"
                     ) {

#Load library
  library(data.table)

# Read data
  labs_dt <- fread(labs_file)
  ref_dt  <- fread(ref_file)

# Remove duplicate M/F reference ranges
  ref_dt <- unique(ref_dt[, .(lab, lower, upper)])

# Merge by lab and classify using vectorised conditions
  merged <- merge(labs_dt, ref_dt, by = "lab", all.x = TRUE)

  merged[, status := fifelse(value >= lower & value <= upper,
                             "normal", "out_of_range")]
#Summaries
  abnormal_by_patient_dt5 <- merged[, .(
    total_labs = .N,
    abnormal_n = sum(status == "out_of_range"),
    abnormal_rate = mean(status == "out_of_range")
  ), by = patient_id]

  abnormal_by_lab_dt5 <- merged[, .(
    total_patients = .N,
    abnormal_n = sum(status == "out_of_range"),
    abnormal_rate = mean(status == "out_of_range")
  ), by = lab]

# Order
setorder(merged, patient_id, time_iso, lab)
setorder(abnormal_by_patient_dt5, -abnormal_rate, patient_id)
setorder(abnormal_by_lab_dt5, -abnormal_rate, lab)


return(list(
    classified_dt5 = merged[, .(patient_id, time_iso, lab, value, status)],
    abnormal_by_patient_dt5 = abnormal_by_patient_dt5,
    abnormal_by_lab_dt5 = abnormal_by_lab_dt5
  ))
}

#' TASK 5 - DATA.FRAME VERSION
#'
#' @description
#' This function classifies clinical lab values against reference intervals:
#' it merges lab results and reference ranges, flags each measurement as
#' `"normal"` or `"out_of_range"`, and computes abnormal rates per patient
#' and per lab using base R data.frame functions.
#'
#' @param labs_file Path to CSV with columns: patient_id, time_iso, lab, value
#' @param ref_file  Path to CSV with columns: lab, lower, upper
#' @return list(classified_dt5, abnormal_by_patient_dt5, abnormal_by_lab_dt5)
#' @export

task5_df <- function(labs_file = "clinical_labs.csv",
                     ref_file  = "lab_reference_ranges.csv",
                     meta_file = "sample_metadata.csv") {

# Read data
  labs_df <- read.csv(labs_file, stringsAsFactors = FALSE)
  ref_df  <- read.csv(ref_file,  stringsAsFactors = FALSE)

# Remove duplicate M/F reference ranges
  ref_df <- unique(ref_df[, c("lab", "lower", "upper")])

# Merge and classify vectorised
  merged <- merge(labs_df, ref_df, by = "lab", all.x = TRUE)
  merged$status <- ifelse(merged$value >= merged$lower &
                            merged$value <= merged$upper,
                          "normal", "out_of_range")

# Summaries
abnormal_by_patient_df5 <- aggregate(status ~ patient_id, data = merged,
                                       FUN = function(x) mean(x == "out_of_range"))
 names(abnormal_by_patient_df5)[2] <- "abnormal_rate"

abnormal_by_patient_df5$total_labs <- as.numeric(table(merged$patient_id))
abnormal_by_patient_df5$abnormal_n <- round(
  abnormal_by_patient_df5$abnormal_rate * abnormal_by_patient_df5$total_labs, 0)

abnormal_by_lab_df5 <- aggregate(status ~ lab, data = merged,
                                   FUN = function(x) mean(x == "out_of_range"))
names(abnormal_by_lab_df5)[2] <- "abnormal_rate"
abnormal_by_lab_df5$total_patients <- as.numeric(table(merged$lab))
abnormal_by_lab_df5$abnormal_n <- round(
abnormal_by_lab_df5$abnormal_rate * abnormal_by_lab_df5$total_patients, 0)

# Order
  merged <- merged[order(merged$patient_id, merged$time_iso, merged$lab), ]

abnormal_by_patient_df5 <- abnormal_by_patient_df5[
    order(-abnormal_by_patient_df5$abnormal_rate, abnormal_by_patient_df5$patient_id), ]

abnormal_by_lab_df5 <- abnormal_by_lab_df5[
    order(-abnormal_by_lab_df5$abnormal_rate, abnormal_by_lab_df5$lab), ]


# Order the columns
abnormal_by_patient_df5 <- abnormal_by_patient_df5[, c("patient_id", "total_labs", "abnormal_n", "abnormal_rate")]
abnormal_by_lab_df5     <- abnormal_by_lab_df5[, c("lab", "total_patients", "abnormal_n", "abnormal_rate")]

return(list(
    classified_df5 = merged[, c("patient_id", "time_iso", "lab", "value", "status")],
    abnormal_by_patient_df5 = abnormal_by_patient_df5,
    abnormal_by_lab_df5 = abnormal_by_lab_df5
  ))
}
