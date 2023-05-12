library(dplyr)
library(tidyverse)
library(edgeR)
library(bedr)
library(readr)
library(tidyverse)
library(magrittr)
library(edgeR)
library(RNOmni)

rename_cols <- function(MyString){
    return(
           gsub('.*/(\\w+)_(\\w+).*', '\\1_\\2', MyString)
    )
}

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}

featureCounts_FileIn <- 'featureCounts/rnaseq/Counts.txt'

counts <- read_tsv(featureCounts_FileIn, comment = "#") %>%
    rename_with(rename_cols, starts_with("/project2/"))

flu_cols <- counts %>% colnames() %>% grep("_Flu$", .)
ni_cols <- counts %>% colnames() %>% grep("_NI$", .)

gtf <- read_tsv('Annotations/gencode.v19.annotation.gtf', comment="#", n_max=Inf, 
                col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )) %>%
    filter(feature=="gene") %>%
    select(attribute)
gtf$gene_id <- unlist(lapply(gtf$attribute, extract_attributes, "gene_id"))
gtf$gene_type <- unlist(lapply(gtf$attribute, extract_attributes, "gene_type"))
ProteinCodingGenes <- filter(gtf, gene_type=="protein_coding") %>% pull(gene_id)

counts['Chr'] <- counts$Chr %>% sapply(., function(x) unique(strsplit(x, ';')[[1]])) %>% 
                                       unlist() %>% as.vector() 

counts['#Chr'] <- counts$Chr %>% sapply(., function(x) unique(strsplit(x, ';')[[1]])) %>% 
                                        unlist() %>% as.vector() %>% gsub('chr', '', .)

vec <- counts$Start
split_vec <- strsplit(vec, ';')
num_vec <- lapply(split_vec, as.numeric)
min_vec <- sapply(num_vec, min)
counts['start'] <- min_vec
                                        
vec <- counts$End
split_vec <- strsplit(vec, ';')
num_vec <- lapply(split_vec, as.numeric)
max_vec <- sapply(num_vec, max)
counts['end'] <- max_vec
                                        
counts['strand'] <- counts$Strand %>% sapply(., function(x) unique(strsplit(x, ';')[[1]])) %>% unlist() %>% as.vector()

counts['pid'] <- counts$Geneid
counts['gid'] <- counts$Geneid

dat.cpm <- counts %>%
    filter(Geneid %in% ProteinCodingGenes  ) %>%
    filter(Chr %in% paste0("chr", 1:22)) %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length", "#Chr", 'start', 'end', 'strand', 'pid', 'gid')) %>%
    cpm(log=T, prior.count=0.1)
                                             
MedCpm <- sort(apply(dat.cpm, 1, median), decreasing=T)
NGenesToInclude <- 14000
GenesToInclude <- MedCpm[1:NGenesToInclude] %>% names()
dat.cpm.filtered <- dat.cpm[GenesToInclude,] %>% as.data.frame()
                                             
flu_cols <- dat.cpm.filtered %>% colnames() %>% grep("_Flu$", .)
ni_cols <- dat.cpm.filtered %>% colnames() %>% grep("_NI$", .)

samples <- c("AF04", "AF06", "AF08", "AF10", "AF12", "AF14", "AF16", "AF18", "AF20", "AF22", "AF24", 
             "AF26", "AF28", "AF30", "AF34", "AF36", "AF38", "EU03", "EU05", "EU07", "EU09", "EU13", 
             "EU15", "EU17", "EU19", "EU21", "EU25", "EU27", "EU29", "EU33", "EU37", "EU39", "EU41", 
             "EU43", "EU47")

transform_string <- function(x) {
  if (x == "AF14") {
    return("Epi_AF14_Mock_WGS_1_Epi_AF14_Mock_WGS_1")
  } else if (x == "EU09") {
    return("Epi_EU09_Mock_WGS_1_Epi_EU09_Mock_WGS_1")
  } else {
    return(paste0("Epi_", x, "_NI_WGS_1_Epi_", x, "_NI_WGS_1"))
  }
}

dat.cpm.flu <- dat.cpm.filtered %>% select(all_of(flu_cols))
colnames(dat.cpm.flu) <- dat.cpm.flu %>% colnames() %>% gsub('_Flu$', '', .) 

dat.cpm.flu <- dat.cpm.flu %>% select(all_of(samples))
colnames(dat.cpm.flu) <- colnames(dat.cpm.flu) %>% sapply(., transform_string) %>% as.vector()

dat.cpm.ni <- dat.cpm.filtered %>% select(all_of(ni_cols))
colnames(dat.cpm.ni) <- dat.cpm.ni %>% colnames() %>% gsub('_NI$', '', .) 

dat.cpm.ni <- dat.cpm.ni %>% select(all_of(samples))
colnames(dat.cpm.ni) <- colnames(dat.cpm.ni) %>% sapply(., transform_string) %>% as.vector()

dat.standardized.flu <- dat.cpm.flu %>% t() %>% scale() %>% t()
dat.qqnormed.flu <- apply(dat.standardized.flu, 2, RankNorm)

dat.standardized.ni <- dat.cpm.ni %>% t() %>% scale() %>% t()
dat.qqnormed.ni <- apply(dat.standardized.ni, 2, RankNorm)

bed <- counts %>% 
filter(Geneid %in% rownames(dat.cpm.filtered)) %>%
select('#Chr', 'start', 'end', 'pid', 'gid', 'strand')

bed %>%
    inner_join(
               (dat.qqnormed.flu %>% as.data.frame() %>% rownames_to_column("pid")),
               by = "pid"
    ) %>%
write_tsv('QTLs/rnaseq.flu.unsort.bed.gz')


bed %>%
    inner_join(
               (dat.qqnormed.ni %>% as.data.frame() %>% rownames_to_column("pid")),
               by = "pid"
    ) %>%
write_tsv('QTLs/rnaseq.ni.unsort.bed.gz')

                                             
dat.rpkm <- counts %>%
    filter(Geneid %in% ProteinCodingGenes  ) %>%
    filter(Chr %in% paste0("chr", 1:22)) %>%
    column_to_rownames("Geneid") %>%
    select(everything(), -c("Chr", "Start", "End", "Strand", "Length", "#Chr", 'start', 'end', 'strand', 'pid', 'gid')) %>%
    rpkm(prior.count=0.1, log=T, gene.length=counts$Length)
                                             
dat.rpkm[GenesToInclude,] %>% as.data.frame() %>% rownames_to_column(., var = "Geneid") %>% write_tsv('QTLs/rnaseq.logRPKM.tab.gz')