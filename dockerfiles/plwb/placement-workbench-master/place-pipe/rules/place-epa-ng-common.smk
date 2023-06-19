# =================================================================================================
#     Auxilliary Rules and Functions
# =================================================================================================

# Rule to infer the best ML model parameters given the MSA and tree that we are going to use
# for placement.
rule raxml_ng_model_eval:
    input:
        tree    = config["data"]["reference-tree"],
        msa     = config["data"]["reference-alignment"]
    output:
        "{outdir}/model/model_eval.raxml.bestModel"
    params:
        model       = config["params"]["epa-ng"]["model"],
        model_dir   = lambda wildcards, input, output: os.path.dirname(output[0])
    log:
        "{outdir}/model/model_eval.log"
    threads:
        get_threads( "raxml-ng" )
    conda:
        "../envs/raxml-ng.yaml"
    shell:
        "raxml-ng --evaluate"
        ' --msa "{input.msa}"'
        ' --tree "{input.tree}"'
        " --model {params.model}"
        " --prefix {params.model_dir}/model_eval"
        " --threads {threads}"
        " > {log} 2>&1"

def model_params( wildcards, input ):
    if use_evaluate:
        # Use the raxml best model file
        return input.model
    else:
        # Use the provided file / model string
        return config["params"]["epa-ng"]["model-params"]
