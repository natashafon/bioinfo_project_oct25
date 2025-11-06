# TASK 10 — ATAC-to-Gene Mapping

# DATA.TABLE VERSION
#' @param peaks_file Path to CSV file of ATAC-seq peaks (columns: chr, start, end, score, etc.)
#' @param genes_file Path to CSV file of gene annotation (columns: chr, gene_start, gene_end, gene)
#' @return list(overlaps_dt10, peaks_per_gene_dt10, top20_dt10)
#' @export

task10_dt <- function(peaks_file = "atac_peaks.bed.csv",
                      genes_file = "gene_annotation.bed.csv") {
  
library(data.table)
  
#Read data
peaks <- fread(peaks_file)
genes <- fread(genes_file)
  
#Standardize column names
setnames(peaks, 1:3, c("chr", "start", "end"))
setnames(genes, 1:4, c("chr", "gene_start", "gene_end", "gene"))
  
# Set keys
setkey(peaks, chr, start, end)
setkey(genes, chr, gene_start, gene_end)
  
#Intersect peaks with gene bodies
overlaps <- foverlaps(peaks, genes,
                        by.x = c("chr", "start", "end"),
                        type = "any", nomatch = 0L)
  
#Compute overlap length (bp)
overlaps[, overlap_bp := pmin(end, gene_end) - pmax(start, gene_start)]
overlaps <- overlaps[overlap_bp > 0]
  
#Count peaks and sum overlap per gene
peaks_per_gene <- overlaps[, .(
    n_peaks = .N,
    total_overlap_bp = sum(overlap_bp)
  ), by = gene]
  
#Top 20 genes by total overlap
top20 <- peaks_per_gene[order(-total_overlap_bp)][1:20]
  
return(list(
    overlaps_dt10 = overlaps,
    peaks_per_gene_dt10 = peaks_per_gene,
    top20_dt10 = top20
  ))
}

# SHOW RESULTS 
res_dt10 <- task10_dt("atac_peaks.bed.csv", "gene_annotation.bed.csv")
head(res_dt10$overlaps_dt10)
head(res_dt10$peaks_per_gene_dt10)
res_dt10$top20_dt10

#DATA.FRAME VERSION 
#' @param peaks_file Path to CSV file of ATAC-seq peaks (columns: chr, start, end, score, etc.)
#' @param genes_file Path to CSV file of gene annotation (columns: chr, gene_start, gene_end, gene)
#' @return list(overlaps_df10, peaks_per_gene_df10, top20_df10)
#' @export

task10_df <- function(peaks_file = "atac_peaks.bed.csv",
                      genes_file = "gene_annotation.bed.csv") {
  
# Read data
peaks <- read.csv(peaks_file)
genes <- read.csv(genes_file)

names(peaks)[1:3] <- c("chr", "start", "end")
names(genes)[1:4] <- c("chr", "gene_start", "gene_end", "gene")
  
#Container for overlaps
overlap_list <- list()
  
#Loop per chromosome
  for (ch in intersect(unique(peaks$chr), unique(genes$chr))) {
    p_sub <- subset(peaks, chr == ch)
    g_sub <- subset(genes, chr == ch)
    
#For each gene, find overlapping peaks
for (i in seq_len(nrow(g_sub))) {
      g <- g_sub[i, ]
      hits <- which(p_sub$end > g$gene_start & p_sub$start < g$gene_end)
      
if (length(hits) > 0) {
        tmp <- p_sub[hits, ]
        tmp$gene <- g$gene
        tmp$overlap_bp <- pmin(tmp$end, g$gene_end) - pmax(tmp$start, g$gene_start)
        tmp <- tmp[tmp$overlap_bp > 0, ]
        overlap_list[[length(overlap_list) + 1]] <- tmp
      }
    }
  }
  
# Combine overlaps
overlaps <- do.call(rbind, overlap_list)

# Order 
overlaps <- overlaps[order(overlaps$chr, overlaps$start, overlaps$end), ]

  
#Summarize per gene
peaks_per_gene <- do.call(rbind,
                            lapply(split(overlaps, overlaps$gene), function(x)
                              data.frame(
                                gene = unique(x$gene),
                                n_peaks = nrow(x),
                                total_overlap_bp = sum(x$overlap_bp)
                              ))
  )
  
# Top 20 genes
top20 <- head(peaks_per_gene[order(-peaks_per_gene$total_overlap_bp), ], 20)
  
return(list(
    overlaps_df10 = overlaps,
    peaks_per_gene_df10 = peaks_per_gene,
    top20_df10 = top20
  ))
}


# SHOW RESULTS 
res_df10 <- task10_df("atac_peaks.bed.csv", "gene_annotation.bed.csv")
head(res_df10$overlaps_df10)
head(res_df10$peaks_per_gene_df10)
res_df10$top20_df10



