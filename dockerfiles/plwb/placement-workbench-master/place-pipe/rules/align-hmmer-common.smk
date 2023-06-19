# =================================================================================================
#     Auxilliary Rules and Functions
# =================================================================================================

rule hmmer_build:
    input:
        config["data"]["reference-alignment"]
    output:
        "{outdir}/hmmer/profile.hmm"
    params:
        extra   = config["params"]["hmmer"]["hmmbuild"]["extra"],
        states  = hmmer_datatype_string
    log:
        "{outdir}/hmmer/build.log"
    threads:
        get_threads( "hmmer" )
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmbuild --cpu {threads} --{params.states} {params.extra} "
        '{output} "{input}" > {log} 2>&1'
