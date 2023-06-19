# =================================================================================================
#     Setup Dependencies
# =================================================================================================

include: "place-epa-ng-common.smk"

# =================================================================================================
#     Placement with epa-ng
# =================================================================================================

rule epa_ng_place:
    group: "placement"
    input:
        tree    = config["data"]["reference-tree"],
        msa     = config["data"]["reference-alignment"],
        query   = "{outdir}/{clusterer}/aligned/{sample}/queries.afa",
        # add the model file for DAG resolution only
        model   = rules.raxml_ng_model_eval.output if use_evaluate else []
    output:
        jplace  = "{outdir}/{clusterer}/placed/{sample}.jplace"
    params:
        model   = model_params,
        redo    = True
    log:
        "{outdir}/{clusterer}/place/{sample}.log"
    threads:
        get_threads( "epa-ng" )
    conda:
        "../envs/epa-ng.yaml"
    script:
        "../scripts/epa-ng.py"
