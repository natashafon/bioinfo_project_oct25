# TASK 11 - Map SNPs to Genes

#' TASK 11 - DATA.TABLE VERSION
#'
#' @description
#' This function maps SNP variants to genes using genomic interval overlap:
#' it converts variant positions to 1-bp ranges, identifies overlaps with gene
#' coordinates, counts HIGH-impact variants per gene and sample, and lists genes
#' containing HIGH-impact variants across all samples.
#'
#' @param variants_file Path to CSV file of variants (columns: chr, pos, sample_id, impact)
#' @param genes_file Path to CSV file of gene annotation (columns: chr, start, end, gene)
#' @return list(overlaps_dt11, counts_per_gene_dt, genes_high_all_dt)
#' @export

task11_dt <- function(variants_file = "variants.csv",
                      genes_file    = "gene_annotation.bed.csv") {

#Load library
  library(data.table)

# Read data
variants <- fread(variants_file, colClasses = "character")
genes    <- fread(genes_file,    colClasses = "character")

# Convert numeric columns
variants[, pos := as.integer(pos)]
genes[, c("start","end") := .(as.integer(start), as.integer(end))]

#Create overlaps
overlaps_list <- lapply(1:nrow(variants), function(i) {
    v <- variants[i]
    subset <- genes[chr == v$chr & start <= v$pos & end >= v$pos]
    if (nrow(subset) > 0) {
      data.table(chr = v$chr,
                 gene = subset$gene,
                 sample_id = v$sample_id,
                 impact = v$impact)
    } else NULL
  })

overlaps <- rbindlist(overlaps_list, use.names = TRUE, fill = TRUE)

# Count HIGH-impact variants
counts <- overlaps[impact == "HIGH", .(n_high = .N), by = .(gene, sample_id)]
setorder(counts, gene, sample_id)

# Genes with HIGH-impact variants in all samples
n_samples <- uniqueN(variants$sample_id)
genes_all <- counts[, .N, by = gene][N == n_samples, gene]


return(list(
    overlaps_dt11      = overlaps,
    counts_per_gene_dt = counts,
    genes_high_all_dt  = genes_all
  ))
}

#' TASK 11 - DATA.FRAME VERSION
#'
#' @description
#' This function maps SNPs to genes using base R:
#' it compares variant positions against gene coordinates, counts HIGH-impact
#' variants per gene and sample, and reports genes that show HIGH-impact
#' variants in every sample.
#'
#' @param variants_file Path to CSV file of variants (columns: chr, pos, sample_id, impact)
#' @param genes_file Path to CSV file of gene annotation (columns: chr, start, end, gene)
#' @return list(overlaps_df11, counts_per_gene_df, genes_high_all_df)
#' @export

task11_df <- function(variants_file = "variants.csv",
                      genes_file    = "gene_annotation.bed.csv") {

#Read data
variants <- read.csv(variants_file, stringsAsFactors = FALSE)
genes    <- read.csv(genes_file, stringsAsFactors = FALSE)

variants$chr <- as.character(variants$chr)
variants$pos <- as.integer(variants$pos)
genes$chr    <- as.character(genes$chr)
genes$start  <- as.integer(genes$start)
genes$end    <- as.integer(genes$end)
genes$gene   <- as.character(genes$gene)

# Create overlaps
overlaps <- data.frame(chr=character(), gene=character(),
                         sample_id=character(), impact=character(),
                         stringsAsFactors = FALSE)

for (i in 1:nrow(variants)) {
    v <- variants[i, ]
    matches <- subset(genes, chr == v$chr &
                        start <= v$pos &
                        end >= v$pos)
if (nrow(matches) > 0) {
      new_rows <- data.frame(chr = v$chr,
                             gene = matches$gene,
                             sample_id = v$sample_id,
                             impact = v$impact,
                             stringsAsFactors = FALSE)
      overlaps <- rbind(overlaps, new_rows)
    }
  }

# Count HIGH-impact variants
  high <- subset(overlaps, impact == "HIGH")
  if (nrow(high) == 0) {
    counts <- data.frame(gene=character(), sample_id=character(), n_high=integer())
  } else {
    counts <- aggregate(
      list(n_high = rep(1, nrow(high))),
      by = list(gene = high$gene, sample_id = high$sample_id),
      FUN = sum
    )
    counts <- counts[order(counts$gene, counts$sample_id), ]
  }

# Genes with HIGH-impact variants in all samples
  n_samples <- length(unique(variants$sample_id))
  if (nrow(counts) == 0) {
    genes_all <- character(0)
  } else {
    gene_counts <- aggregate(sample_id ~ gene, data = counts,
                             FUN = function(x) length(unique(x)))
    genes_all <- sort(gene_counts$gene[gene_counts$sample_id == n_samples])
  }


return(list(
    overlaps_df11      = overlaps,
    counts_per_gene_df = counts,
    genes_high_all_df  = genes_all
  ))
}
