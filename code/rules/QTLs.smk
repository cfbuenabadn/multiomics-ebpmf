rule PrepareBedForQTLTools:
    input:
        'featureCounts/rnaseq/Counts.txt',
        'Annotations/gencode.v19.annotation.gtf'
    output:
        'QTLs/rnaseq.flu.unsort.bed.gz',
        'QTLs/rnaseq.ni.unsort.bed.gz',
        'QTLs/rnaseq.logRPKM.tab.gz'
    log:
        'logs/prepareQTLsPhenotypes.log'
    shell:
        """
        Rscript scripts/PrepareQTLPhenotypes.R &> {log}
        """
        
rule SortQTLtoolsPhenotypeTable:
    input:
        "QTLs/rnaseq.{phenotype}.unsort.bed.gz",
    output:
        bed = "QTLs/rnaseq.{phenotype}.bed.gz",
        tbi = "QTLs/rnaseq.{phenotype}.bed.gz.tbi",
    log:
        "logs/SortQTLtoolsPhenotypeTable/{phenotype}.log"
    resources:
        mem_mb = 24000
    shell:
        """
        (bedtools sort -header -i {input} | bgzip /dev/stdin -c > {output.bed}) &> {log}
        (tabix -p bed {output.bed}) &>> {log}
        """


rule PhenotypePCs:
    """
    QTLtools format expression PCs as covariates
    including the number of PCs that explain more
    variance then when the phenotype table is
    permuted
    """
    input:
        "QTLs/rnaseq.{phenotype}.bed.gz",
    output:
        "QTLs/rnaseq.{phenotype}.bed.pca",
    log:
        "logs/PhenotypePCs.{phenotype}.log"
    resources:
        mem_mb = 24000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PermuteAndPCA.R {input} {output} &> {log}
        """

def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 4000
    else:
        return 52000
        
def GetQTLtoolsPassFlags(wildcards):
    if wildcards.Pass == "PermutationPass":
        return "--permute 1000"
    elif wildcards.Pass == "NominalPass":
        return "--nominal 1"
        
N_PermutationChunks = 50
ChunkNumbers = range(0, 1+N_PermutationChunks) 

rule QTLtools_generalized:
    input:
        vcf = 'Annotations/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.v4.vcf.bgz',
        tbi = 'Annotations/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.v4.vcf.bgz.tbi',
        bed = "QTLs/rnaseq.{phenotype}.bed.gz",
        bed_tbi = "QTLs/rnaseq.{phenotype}.bed.gz.tbi",
        cov = "QTLs/rnaseq.{phenotype}.bed.pca"
    output:
        temp("QTLs/Chunks/{Pass}/{QTLTools_chunk_n}.{phenotype}.txt")
    log:
        "logs/QTLtools_cis_permutation_pass/{phenotype}.{Pass}.{QTLTools_chunk_n}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    envmodules:
        "gsl/2.5"
    params:
        WindowFlag = "--window 100000",
        OtherFlags = "",
        PassFlags = GetQTLtoolsPassFlags,
        ExcFlag = "",
    wildcard_constraints:
        phenotype = 'ni|flu',
        n = "|".join(str(i) for i in ChunkNumbers),
        Pass = "PermutationPass|NominalPass"
    shell:
        """
        {config[QTLtools]} cis --std-err --chunk {wildcards.QTLTools_chunk_n} {N_PermutationChunks} --vcf {input.vcf} --bed {input.bed} --cov {input.cov} --out {output} {params.OtherFlags} {params.WindowFlag} {params.PassFlags} {params.ExcFlag} &> {log}
        if [ ! -f {output} ]
        then
            touch {output}
        fi
        """

rule Gather_QTLtools_cis_pass:
    input:
        expand("QTLs/Chunks/{{Pass}}/{QTLTools_chunk_n}.{{phenotype}}.txt", QTLTools_chunk_n=ChunkNumbers )
    output:
        "QTLs/{Pass}.{phenotype}.txt.gz"
    log:
        "logs/QTLs/{Pass}.{phenotype}.log"
    wildcard_constraints:
        phenotype = 'ni|flu',
        Pass = "PermutationPass|NominalPass"
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """
        