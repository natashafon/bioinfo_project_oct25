
# TASK FINAL - Associate cell type to integration clusters (Normal vs Tumor)

#' TASK FINAL - DATA.TABLE VERSION
#'
#' @description
#' Associates cell types to Seurat integration clusters and tissue status (N/T):
#' it merges integration and annotation tables by `cell`, outputs the combined
#' table, counts cell types per cluster, builds a summary with normalized
#' percentages by cluster and sample type (N vs T), and produces bar plots
#' of raw counts and normalized % distributions.
#'
#' @param integration_file Path to annotated integration CSV
#'   (columns: cell, integration_cluster)
#' @param annotation_file Path to annotation CSV
#'   (columns: cell, cell_type, sample_type)
#' @return list(merged_dt13, count_dt13, summary_dt13, plot_counts_dt13, plot_percent_dt13)
#' @import data.table
#' @export

task13_dt <- function(integration_file = "annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv",
                      annotation_file  = "nt_combined_clustering.output.csv") {

# Load libraries
library(data.table)
library(ggplot2)

# Prepare output directory
if (!dir.exists("output")) dir.create("output")

# Read input files
integ_dt <- fread(integration_file)
annot_dt <- fread(annotation_file)

# Clean IDs
integ_dt[, cell := gsub("_X_", "", cell)]

# Merge
merged_dt13 <- merge(integ_dt, annot_dt, by = "cell", sort = FALSE)

# Count per cluster and cell type
count_dt13 <- merged_dt13[, .(count = .N), by = .(integration_cluster, cell_type)]
setorder(count_dt13, integration_cluster, cell_type)

# Summary table: add tissue type and normalized percentages
summary_dt13 <- merged_dt13[, .(count = .N),
                              by = .(integration_cluster, cell_type, sample_type)]
summary_dt13[, total := sum(count), by = .(integration_cluster, sample_type)]
summary_dt13[, percent := round((count / total) * 100, 2)]

setorder(summary_dt13, integration_cluster, sample_type, cell_type)

# Order
summary_dt13[, integration_cluster := factor(integration_cluster,
                                             levels = sort(unique(integration_cluster)))]
summary_dt13[, sample_type := factor(sample_type, levels = c("N", "T"))]

# Plot 1: Distribution plot
plot_counts_dt13 <- ggplot(summary_dt13,
                             aes(x = factor(integration_cluster),
                                 y = count,
                                 fill = cell_type)) +
geom_bar(stat = "identity") +
facet_wrap(~sample_type) +
theme_minimal() +
labs(title = "Cell Type Distribution by Cluster and Tissue Type",
         x = "Integration Cluster",
         y = "Number of Cells",
         fill = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Distribution plot with normalized percentages
plot_percent_dt13 <- ggplot(summary_dt13,
                              aes(x = factor(integration_cluster),
                                  y = percent,
                                  fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~sample_type) +
    theme_minimal() +
    labs(title = "Cell Distribution by Cluster and Tissue Type (Normalized %)",
         x = "Integration Cluster",
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save outputs
fwrite(merged_dt13, "output/merged_dt13.csv")
fwrite(count_dt13, "output/count_dt13.csv")
fwrite(summary_dt13, "output/summary_dt13.csv")
ggsave("output/plot_counts_dt13.pdf", plot_counts_dt13, width = 8, height = 5)
ggsave("output/plot_percent_dt13.pdf", plot_percent_dt13, width = 8, height = 5)


  return(list(
    merged_dt13       = merged_dt13,
    count_dt13        = count_dt13,
    summary_dt13      = summary_dt13,
    plot_counts_dt13  = plot_counts_dt13,
    plot_percent_dt13 = plot_percent_dt13
  ))
}

#' TASK FINAL - DATA.FRAME VERSION
#'
#' @description
#' Associates cell types to Seurat integration clusters and tissue status (N/T):
#' it merges integration and annotation data by `cell`, returns the merged table,
#' counts cell types per cluster, creates a summary with normalized percentages
#' by cluster and sample type (N vs T), and generates bar plots for counts and
#' normalized % distributions.
#'
#' @param integration_file Path to annotated integration CSV
#'   (columns: cell, integration_cluster)
#' @param annotation_file Path to annotation CSV
#'   (columns: cell, cell_type, sample_type)
#' @return list(merged_df13, count_df13, summary_df13, plot_counts_df13, plot_percent_df13)
#' @export

task13_df <- function(integration_file = "annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv",
                      annotation_file  = "nt_combined_clustering.output.csv") {

# Prepare output directory
if (!dir.exists("output")) dir.create("output")

# Read input files
integ_df <- read.csv(integration_file, stringsAsFactors = FALSE)
annot_df <- read.csv(annotation_file, stringsAsFactors = FALSE)

# Clean IDs
integ_df$cell <- gsub("_X_", "", integ_df$cell)

# Merge
merged_df13 <- merge(integ_df, annot_df, by = "cell")

# Count per cluster
count_df13 <- as.data.frame(table(
    integration_cluster = merged_df13$integration_cluster,
    cell_type = merged_df13$cell_type
  ))
names(count_df13)[3] <- "count"
count_df13$count <- as.integer(count_df13$count)
count_df13 <- count_df13[count_df13$count > 0, ]
count_df13 <- count_df13[order(count_df13$integration_cluster, count_df13$cell_type), ]

# Summary with sample_type and normalized %
summary_df13 <- as.data.frame(table(
    integration_cluster = merged_df13$integration_cluster,
    cell_type = merged_df13$cell_type,
    sample_type = merged_df13$sample_type
  ))

names(summary_df13)[4] <- "count"
summary_df13$count <- as.integer(summary_df13$count)
summary_df13 <- summary_df13[summary_df13$count > 0, ]
summary_df13$total <- ave(summary_df13$count,
                            interaction(summary_df13$integration_cluster, summary_df13$sample_type),
                            FUN = sum)
summary_df13$percent <- round((summary_df13$count / summary_df13$total) * 100, 2)
summary_df13 <- summary_df13[order(summary_df13$integration_cluster,
                                     summary_df13$sample_type,
                                     summary_df13$cell_type), ]
# Order
summary_df13$integration_cluster <- factor(summary_df13$integration_cluster,
                                           levels = sort(unique(summary_df13$integration_cluster)))
summary_df13$sample_type <- factor(summary_df13$sample_type,
                                   levels = c("N", "T"))

# Plot 1: Distribution plot
  plot_counts_df13 <- ggplot(summary_df13,
                             aes(x = factor(integration_cluster),
                                 y = count,
                                 fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~sample_type) +
    theme_minimal() +
    labs(title = "Cell Type Distribution by Cluster and Tissue Type",
         x = "Integration Cluster",
         y = "Number of Cells",
         fill = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Distribution plot with normalized percentages
  plot_percent_df13 <- ggplot(summary_df13,
                              aes(x = factor(integration_cluster),
                                  y = percent,
                                  fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~sample_type) +
    theme_minimal() +
    labs(title = "Cell Distribution by Cluster and Tissue Type (Normalized %)",
         x = "Integration Cluster",
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save outputs
write.csv(merged_df13, "output/merged_df13.csv", row.names = FALSE)
write.csv(count_df13, "output/count_df13.csv", row.names = FALSE)
write.csv(summary_df13, "output/summary_df13.csv", row.names = FALSE)
ggsave("output/plot_counts_df13.pdf", plot_counts_df13, width = 8, height = 5)
ggsave("output/plot_percent_df13.pdf", plot_percent_df13, width = 8, height = 5)


return(list(
    merged_df13       = merged_df13,
    count_df13        = count_df13,
    summary_df13      = summary_df13,
    plot_counts_df13  = plot_counts_df13,
    plot_percent_df13 = plot_percent_df13
  ))
}

