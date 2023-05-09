library(dplyr)
library(tidyverse)
library(data.table)
library(bedr)
library(readr)
#library(ebpmf)

library(NNLM)
#library(Matrix)
library(CountClust)
library(fastTopics)
library(FNN)

library(robustbase)

library(smashr)
library(ebpmf)
library(Matrix)


load_coverage_counts <- function(assay, gene){
    file_name <- paste0("coverage/counts/", assay, "/", gene, ".csv.gz")
    df <- read_csv(file_name) %>%
    column_to_rownames(var = "Sample_ID") 
    
    rownames(df) <- sub(paste0("\\.", assay), "", rownames(df))
    
    colnames(df) <- paste0(colnames(df), '.', assay)
    
    return (df)
    
}

group_matrix <- function(df, k=10){
    n_groups <- ncol(df) %/% k

    # Split dataframe into groups of 10 columns and sum them
    output <- lapply(seq_len(n_groups), function(i) {
      rowSums(df[, (k*(i-1)+1):(k*i)], na.rm = TRUE)
    })

    # Combine the output into a single dataframe
    output_df <- do.call(cbind, output)
    
    return(output_df)
}

plot_structure <- function(fit, gene_name, kfactors, annotation_){
    EL = fit$EL
    indis = rownames(EL)
    tissue_label  <- c()
    
    colnames(EL) <- sprintf("factor%d",seq(1:kfactors))
    colores = RColorBrewer::brewer.pal(kfactors+1,  "Paired")
    
    if (kfactors==2){
        colores <- c('#1f77b4','#ff7f0e')
    }
    
   StructureGGplot(EL, annotation = annotation_,
                      palette = colores, figure_title = gene_name, yaxis_label='samples',
                      axis_tick = list(axis_ticks_length = 0.1, 
                                       axis_ticks_lwd_y = 1, 
                                       axis_ticks_lwd_x = 1, 
                                       axis_label_size = 12, axis_label_face = "bold"),
                         legend_title_size = 12, legend_key_size = 1, legend_text_size = 12)
}

args = commandArgs(trailingOnly=TRUE)
region_name = args[1]
output = args[3]
output_tab = args[4]
K = as.integer(args[2])

set.seed(123)

rnaseq <- load_coverage_counts('rnaseq', region_name) %>% group_matrix(20)
K4me3 <- load_coverage_counts('K4me3', region_name) %>% group_matrix(20)
K4me1 <- load_coverage_counts('K4me1', region_name) %>% group_matrix(20)
K27ac <- load_coverage_counts('K27ac', region_name) %>% group_matrix(20)
atacseq <- load_coverage_counts('atacseq', region_name) %>% group_matrix(20)
K27me3 <- load_coverage_counts('K27me3', region_name) %>% group_matrix(20)

shared_idx <- rownames(rnaseq)[rownames(rnaseq) %in% rownames(K4me3)]

nsamples <- length(shared_idx)
counts <- cbind(rpois(nsamples, 1), rnaseq[shared_idx,], 0, K4me3[shared_idx,], 0, K27ac[shared_idx,], 0, atacseq[shared_idx,], rpois(nsamples, 1))

fit_ebpmf = ebpmf_identity(counts,K=K,tol = 1e-3,maxiter = 100,init = 'kmeans')


# remove the XX from the strings
sample_annot <- gsub("\\d+", "", rownames(counts))

# create a second vector with only "Flu" and "NI"
sample_superannot <- gsub(".*_(.*)", "\\1", rownames(counts))


indis <- counts %>% rownames() 
    
annotation_ = data.frame(sample_id = indis,tissue_label = factor(sample_annot))

colnames(fit_ebpmf$EL) <- paste0('factor', 1:K)
rownames(fit_ebpmf$EL) <- indis


fit_ebpmf$EL %>% head()
annotation_ %>% head()

pdf(file=paste0('plots/', region_name, '.K', as.character(K), '.structure.pdf'))
plot_structure(fit_ebpmf, region_name, K, annotation_)
dev.off()

saveRDS(list(region=region_name,
             annotation = annotation_,
             fit_ebpmf = fit_ebpmf,
             coords = colnames(rnaseq),
             samples = rownames(counts)
            ),
        file=output
       )

EF <- as.data.frame(fit_ebpmf$EF)

colnames(EF) <- paste0('factor', 1:K)

write_tsv(EF, output_tab)
