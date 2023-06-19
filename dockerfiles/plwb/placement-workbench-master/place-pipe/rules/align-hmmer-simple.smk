# =================================================================================================
#     Setup Dependencies
# =================================================================================================

include: "align-hmmer-common.smk"

# =================================================================================================
#     Alignment with hmmer
# =================================================================================================

rule hmmer_search:
    input:
        msa      = config["data"]["reference-alignment"],
        hmmfile  = "{outdir}/hmmer/profile.hmm",
        seqfile  = "{outdir}/{clusterer}/samples/{sample}/queries.fa"
    output:
        "{outdir}/{clusterer}/filtered/{sample}/filtered.fa"
    params:
        noali = True
    log:
        "{outdir}/{clusterer}/filtered/{sample}/hmmsearch.log"
    threads:
        get_threads( "hmmer" )
    conda:
        "../envs/hmmer.yaml"
    script:
        "../scripts/hmmsearch.py"

rule hmmer_align:
    input:
        msa      = config["data"]["reference-alignment"],
        hmmfile  = "{outdir}/hmmer/profile.hmm",
        seqfile  = rules.hmmer_search.output[0]
    output:
        "{outdir}/{clusterer}/aligned/{sample}/queries.afa"
    params:
        outformat   = "afa",
        states      = hmmer_datatype_string
    log:
        "{outdir}/{clusterer}/aligned/{sample}/hmmalign.log"
    threads:
        get_threads( "hmmer" )
    conda:
        "../envs/hmmer.yaml"
    script:
        "../scripts/hmmalign.py"