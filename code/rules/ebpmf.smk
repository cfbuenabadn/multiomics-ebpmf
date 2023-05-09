def GetRegionCoords(wildcards):
    if wildcards.Region == 'IL27':
        return 'chr16 28474000 28542000'
    elif wildcards.Region == 'IDO1':
        return 'chr8 39750000 39890000'
    elif wildcards.Region == 'TAP1':
        return 'chr6 32785000 32924000'
    elif wildcards.Region == 'CMC2':
        return 'chr16 80987000 81113000'
    elif wildcards.Region == 'AIF1L':
        return 'chr9 133954000 134078000'
    elif wildcards.Region == 'RETN':
        return 'chr19 7729000 7754000'
    elif wildcards.Region == 'GTF3C6':
        return 'chr6 111263000 111364000'
    elif wildcards.Region == 'IRF2':
        return 'chr4 185109000 185596000'
    elif wildcards.Region == 'CMC2_extended':
        return 'chr16 80814000 81241000'
    elif wildcards.Region == 'FAM220A_extended':
        return 'chr7 6169000 6588600'
    elif wildcards.Region == 'GTF3C6_extended':
        return 'chr6 111079900 111489100'
                
    elif wildcards.Region == 'HDDC2_extended':
        return 'chr6 125418600 125823300'
        
    elif wildcards.Region == 'HLA-DRB5_extended':
        return 'chr6 32320700 32727800'
        
    elif wildcards.Region == 'IRF2_extended':
        return 'chr4 185108800 185595800'
        
    elif wildcards.Region == 'MRPL18_extended':
        return 'chr6 160011400 160419500'
        
    elif wildcards.Region == 'RP11-661A12.5_extended':
        return 'chr8 144424100 144831900'
        
    elif wildcards.Region == 'YBEY_extended':
        return 'chr21 47506200 47917700'
        
    elif wildcards.Region == 'ZFAND2A_extended':
        return 'chr7 991700 1399900'
        
    elif wildcards.Region == 'RETN_extended':
        return 'chr19 7680200 7776000'
        
        

rule CreateRegionCounts:
    input:
        expand("coverage/bed/{Phenotype}/{IndID}.bed.bgz", zip, Phenotype=samples.Phenotype, IndID=samples.Sample),
        expand("coverage/bed/{Phenotype}/{IndID}.bed.bgz.tbi", zip, Phenotype=samples.Phenotype, IndID=samples.Sample),
    output:
        expand("coverage/counts/{Phenotype}/{{Region}}.csv.gz", Phenotype=Phenotype_list)
    resources:
        mem_mb = 24000,
    log:
        "logs/counts/{Region}.log"
    wildcard_constraints:
        Region = '|'.join(['IL27', 'IDO1', 'TAP1', 'CMC2', 'AIF1L', 'RETN', 'GTF3C6', 'IRF2', 'CMC2_extended',
                           'FAM220A_extended', 'GTF3C6_extended', 'HDDC2_extended', 'HLA-DRB5_extended', 'IRF2_extended',
                           'MRPL18_extended', 'RP11-661A12.5_extended', 'YBEY_extended', 'ZFAND2A_extended', 'RETN_extended'])
    params:
        GetRegionCoords
    shell:
        """
        python scripts/prepare_counts.py {params} {wildcards.Region} &> {log}
        """
        
rule run_ebpmf:
    input:
        expand("coverage/counts/{Phenotype}/{{Region}}.csv.gz", Phenotype=Phenotype_list)
    output:
        "ebpmf_models/RDS/{Region}.K{K}.rds",
        "ebpmf_models/Factors/{Region}.K{K}.tab.gz",
    resources:
        mem_mb = 24000,
    log:
        "logs/ebpmf_models/{Region}.K{K}.log"
    wildcard_constraints:
        Region = '|'.join(['IL27', 'IDO1', 'TAP1', 'CMC2', 'AIF1L', 'RETN', 'GTF3C6', 'IRF2', 'CMC2_extended',
                           'FAM220A_extended', 'GTF3C6_extended', 'HDDC2_extended', 'HLA-DRB5_extended', 'IRF2_extended',
                           'MRPL18_extended', 'RP11-661A12.5_extended', 'YBEY_extended', 'ZFAND2A_extended', 'RETN_extended']),
        K = '2|3|4|5|8|10'
    shell:
        """
        {config[Rscript]} scripts/run_ebpmf.R {wildcards.Region} {wildcards.K} {output} &> {log}
        """
        
        
        
#{config[Rscript]} scripts/prepare_counts.R {params} {wildcards.Region} &> {log}