import pandas as pd
import numpy as np
import os
import glob

samples = pd.read_csv(config["samples"],sep='\t')

Phenotype_list = ['rnaseq', 'K4me1', 'K4me3', 'K27ac', 'K27me3', 'atacseq']

def GetBamForPhenotype(wildcards):
    bam_files = list(samples.loc[(samples.Phenotype == wildcards.Phenotype)].BAM_file)
    return bam_files
    
def GetBamForSample(wildcards):
    if wildcards.Phenotype == 'rnaseq':
        return "/project2/mstephens/multiomics/{Phenotype}/{IndID}.sorted.mdup.bam"
    elif wildcards.Phenotype == 'atacseq':
        return "/project2/mstephens/multiomics/{Phenotype}/{IndID}.sorted.dup.bam"
    else:
        return "/project2/mstephens/multiomics/{Phenotype}/{IndID}_H3{Phenotype}.sorted.dup.bam"