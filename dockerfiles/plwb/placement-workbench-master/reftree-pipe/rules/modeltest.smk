# =================================================================================================
#     Optional step to determine model parameters for treesearch
# =================================================================================================

rule modeltest:
    input:
       "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/trimmed.afa"
    params:
        datatype = config["settings"]["datatype"]
    threads:
        get_threads( "modeltest-ng" )
    output:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/modeltest-ng/model.file"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/modeltest-ng/modeltest.log"
    conda:
        "../envs/modeltest-ng.yaml"
    script:
        "../scripts/modeltest-ng.py"
