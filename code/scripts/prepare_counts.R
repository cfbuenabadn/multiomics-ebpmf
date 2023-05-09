library(dplyr)
library(tidyverse)
library(data.table)
library(bedr)
library(readr)

get_counts <- function(file_name, chrom, start, end){
    coord_str <- paste0(chrom, ':', as.character(start), '-', as.character(end))
    print(coord_str)
    print(file_name)
    gene_df <- tabix(coord_str, file_name, check.chr=FALSE, verbose=FALSE)
    colnames(gene_df) <- c('chrom', 'start', 'end', 'counts')
    counts <- get_region_counts(gene_df, chrom, start, end)
    return(counts)
}

get_region_counts <- function(gene_df, chrom, start, end) {
  coord_list <- c()
  sample_counts <- c()
  
  for (i in 1:nrow(gene_df)) {
    start_ <- as.integer(gene_df$start[i])
    end_ <- as.integer(gene_df$end[i])
    if (i == 1) {
      start_ <- start
    }
    if (i == nrow(gene_df)) {
      end_ <- end
    }
    
    segment_length <- end_ - start_
    
    segment_counts <- as.integer(gene_df$counts[i])
    
    segment_to_add <- rep(segment_counts, segment_length)
    
    sample_counts <- c(sample_counts, segment_to_add)
    coord_list <- c(coord_list, seq(start_, end_-1, by=1))
  }
    
  if (!startsWith(chrom, "chr")){
      chrom <- paste0('chr', chrom)
  }
  
  coord_list <- paste(chrom, coord_list, sep=":")
  
  return(list(sample_counts = sample_counts, coord_list = coord_list))
}

process_and_write <- function(assay, chrom, start, end, region_name){
    
    assay_files <- list.files(paste0('coverage/bed/', assay, '/'))
    sample_list <- unique(gsub("\\.bed.*", "", assay_files))
    
    #df <- data.frame()
    count_list <- list()
    for (sample in sample_list){
        file_name <- paste0('coverage/bed/', assay, '/', sample, '.bed.bgz')
        sample_counts <- get_counts(file_name, chrom, start, end)
        sample_counts_to_add <- sample_counts$sample_counts %>% as.numeric() #%>% rbind(df, .)
        print(paste('processed', file_name))
        count_list[[sample]] = sample_counts_to_add
        #df %>% dim() %>% print()
    }
    
    df <- data.frame(count_list) %>% t()
    
    df %>% dim() %>% print()

    colnames(df) <- sample_counts$coord_list
    rownames(df) <- paste0(sample_list, '.', assay)
                           
    out_name <- paste0("coverage/counts/", assay, "/", region_name, ".csv.gz")

    df %>% as.data.frame() %>% rownames_to_column(., var='Sample_ID') %>% write_csv(., out_name, col_names = TRUE)
}

args = commandArgs(trailingOnly=TRUE)
chrom = args[1]
start = as.integer(args[2])
end = as.integer(args[3])
region = args[4]

chrom_rnaseq <- sub("^chr", "", chrom)

process_and_write('rnaseq', chrom_rnaseq, start, end, region)
process_and_write('atacseq', chrom, start, end, region)
process_and_write('K4me3', chrom, start, end, region)
process_and_write('K4me1', chrom, start, end, region)
process_and_write('K27ac', chrom, start, end, region)
process_and_write('K27me3', chrom, start, end, region)