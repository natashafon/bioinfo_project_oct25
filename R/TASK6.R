# TASK 6 - Nearest-time Matching of Vitals to Lab Draws


#' TASK 6 - DATA.TABLE VERSION
#'
#' @description
#' This function performs nearest-time matching of vitals to lab draws:
#' for each lab record, it finds the closest heart rate (HR) and systolic
#' blood pressure (SBP) measurement in time for the same patient using
#' a rolling join, computes the time lag in minutes, and summarizes the
#' correlation between CRP and HR/SBP by patient.
#'
#' @param labs_file Path to CSV with columns: patient_id, time_iso, lab, value
#' @param vitals_file Path to CSV with columns: patient_id, time_iso, vital (HR/SBP), value
#' @return list(merged_dt6, corr_summary_dt6)
#' @export

task6_dt <- function(labs_file   = "clinical_labs.csv",
                     vitals_file = "vitals_time_series.csv") {

# Load library
library(data.table)

# Read files
labs   <- fread(labs_file)
vitals <- fread(vitals_file)

# Ensure time columns exist and convert
if (!"time_iso" %in% names(labs))   stop("Column 'time_iso' missing in labs file")
if (!"time_iso" %in% names(vitals)) stop("Column 'time_iso' missing in vitals file")

labs[,   lab_time   := as.POSIXct(time_iso)]
vitals[, vital_time := as.POSIXct(time_iso)]

# Reshape vitals to wide format (HR + SBP)
vitals_w <- dcast(vitals, patient_id + vital_time ~ vital, value.var = "value")

# For each lab, find nearest vital record manually
results_list <- lapply(1:nrow(labs), function(i) {
pid <- labs$patient_id[i]
t0  <- labs$lab_time[i]
subset_v <- vitals_w[vitals_w$patient_id == pid]

if (nrow(subset_v) == 0) {
    return(data.table(
        patient_id   = pid,
        lab          = labs$lab[i],
        value        = labs$value[i],
        lab_time     = t0,
        vital_time   = NA,
        HR           = NA_real_,
        SBP          = NA_real_,
        time_lag_min = NA_real_
      ))
    }

subset_v[, diff_min := abs(as.numeric(difftime(vital_time, t0, units = "mins")))]
nearest <- subset_v[which.min(diff_min)]
    nearest[, `:=`(
      patient_id   = pid,
      lab          = labs$lab[i],
      value        = labs$value[i],
      lab_time     = t0,
      time_lag_min = diff_min
    )]

nearest[, .(patient_id, lab, value, lab_time, vital_time, HR, SBP, time_lag_min)]
  })

# Combine all results
  merged_dt6 <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

# Correlation summary (CRP only)
crp_dt <- merged_dt6[lab == "CRP"]
corr_summary_dt6 <- crp_dt[, .(
cor_CRP_HR  = cor(value, HR,  use = "complete.obs"),
cor_CRP_SBP = cor(value, SBP, use = "complete.obs")
  ), by = patient_id]

  return(list(
    merged_dt6 = merged_dt6,
    corr_summary_dt6 = corr_summary_dt6
  ))
}

#' TASK 6 - DATA.FRAME VERSION
#'
#' @description
#' This function performs nearest-time matching of vitals to lab draws:
#' it merges each lab with the temporally nearest HR and SBP values for
#' the same patient, calculates time differences in minutes, and computes
#' per-patient correlations between CRP and vital signs using base R data.frames.
#'
#' @param labs_file Path to CSV with columns: patient_id, time_iso, lab, value
#' @param vitals_file Path to CSV with columns: patient_id, time_iso, vital, value
#' @return list(merged_df6, corr_summary_df6)
#' @export

task6_df <- function(labs_file   = "clinical_labs.csv",
                     vitals_file = "vitals_time_series.csv") {

# Read CSVs
labs   <- read.csv(labs_file, stringsAsFactors = FALSE)
vitals <- read.csv(vitals_file, stringsAsFactors = FALSE)

# Convert time columns
labs$lab_time     <- as.POSIXct(labs$time_iso)
vitals$vital_time <- as.POSIXct(vitals$time_iso)

# Reshape vitals to wide (HR + SBP)
vitals_w <- reshape(
vitals[, c("patient_id", "vital_time", "vital", "value")],
    timevar = "vital",
    idvar   = c("patient_id", "vital_time"),
    direction = "wide"
  )
  names(vitals_w) <- gsub("value.", "", names(vitals_w))

# For each lab, find nearest vital record
  nearest_row <- function(pid, t0) {
    d <- subset(vitals_w, patient_id == pid)
    if (nrow(d) == 0)
      return(data.frame(time_lab = NA, HR = NA, SBP = NA, time_lag_min = NA))
    diffs <- abs(as.numeric(difftime(d$vital_time, t0, units = "mins")))
    k <- which.min(diffs)
    data.frame(time_lab = d$vital_time[k],
               HR = d$HR[k],
               SBP = d$SBP[k],
               time_lag_min = diffs[k])
  }

merged_list <- mapply(nearest_row, labs$patient_id, labs$lab_time, SIMPLIFY = FALSE)
merged_df6  <- cbind(labs[, c("patient_id", "lab", "value", "lab_time")],
                       do.call(rbind, merged_list))

# Correlation summary (CRP only)
crp_df <- subset(merged_df6, lab == "CRP")
corr_summary_df6 <- do.call(rbind, lapply(split(crp_df, crp_df$patient_id), function(sub) {
    data.frame(
      patient_id  = unique(sub$patient_id),
      cor_CRP_HR  = cor(sub$value, as.numeric(sub$HR),  use = "complete.obs"),
      cor_CRP_SBP = cor(sub$value, as.numeric(sub$SBP), use = "complete.obs")
    )
  }))


return(list(
    merged_df6 = merged_df6,
    corr_summary_df6 = corr_summary_df6
  ))
}
