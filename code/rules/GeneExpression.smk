rule DownloadHg19Ref:
    output:
        primary_gtf = "Annotations/gencode.v19.annotation.gtf"
    shell:
        """
        wget -O- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz | zcat > {output}
        """
        
rule featureCounts:
    input:
        bam = GetBamForPhenotype,
        annotations = "Annotations/gencode.v19.annotation.gtf"
    output:
        "featureCounts/{Phenotype}/Counts.txt"
    threads:
        8
    wildcard_constraints:
        Phenotype = "|".join(Phenotype_list)
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{Phenotype}.log"
    shell:
        """
        featureCounts -p -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """

rule GetBigWigs:
    input:
        GetBamForSample
    output:
        "coverage/bigwigs/{Phenotype}/{IndID}.bw"
    resources:
        mem_mb = 42000,
    log:
        "logs/bigwigs/{Phenotype}.{IndID}.log"
    shell:
        """
        (bamCoverage -b {input} -o {output}) &> {log}
        """
        
rule GetBeds:
    input:
        GetBamForSample
    output:
        "coverage/bed/{Phenotype}/{IndID}.bed.bgz",
    resources:
        mem_mb = 42000,
    log:
        "logs/bed_and_tabix/{Phenotype}.{IndID}.log"
    shell:
        """
        (bedtools genomecov -5 -bga -ibam {input} | bgzip > {output}) &> {log};
        """
        
rule GetTabix:
    input:
        "coverage/bed/{Phenotype}/{IndID}.bed.bgz",
    output:
        "coverage/bed/{Phenotype}/{IndID}.bed.bgz.tbi",
    resources:
        mem_mb = 12000,
    log:
        "logs/bed_and_tabix/{Phenotype}.{IndID}.log"
    shell:
        """
        (tabix -p bed {input}) &>> {log}
        """
