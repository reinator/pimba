# =================================================================================================
#     Helper Functions
# =================================================================================================
#
# These functions are meant to be extended via wildcards if and when more tree inference tools are added
#
def bootstrap_params( wildcards ):
    num_trees   = get_highest_override( ['raxml-ng', 'treesearch'], "bs-trees" )
    auto_bs     = get_highest_override( ['raxml-ng', 'treesearch'], "auto-bootstrap" )

    if num_trees:
        if auto_bs:
            return f"autoMRE{{{{{num_trees}}}}}"
        else:
            return f"{num_trees}"
    return ""

def model_params( wildcards, input ):
    datatype    = config["settings"]["datatype"]
    model       = get_highest_override( ['raxml-ng', 'treesearch'], "model" )

    if use_auto_model:
        return input.model_file
    elif model:
        return model
    elif datatype and datatype in ['nt','aa']:
        return "GTR+G" if datatype == 'nt' else  "LG+G"
    else:
        util.fail("'datatype' field is required, and must be either 'nt' or 'aa'.")

def datatype( wildcards ):
    dt = config["settings"]["datatype"]
    if dt == "nt":
        return "DNA"
    elif dt == "aa":
        return "AA"
    else:
        util.fail("'datatype' field is required, and must be either 'nt' or 'aa'.")

# =================================================================================================
#     Tree Search with RAxML-ng
# =================================================================================================

rule get_constraint:
    input:
        tree    = config['data']['constraint_tree'],
        fasta   = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/trimmed.afa"
    output:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/constraint.newick"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/get_constraint.log"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/constraintify.py"

rule treesearch_raxmlng:
    input:
        msa             = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/trimmed.afa",
        model_file      = rules.modeltest.output[0] if use_auto_model else [],
        tree_constraint = rules.get_constraint.output[0] if config['data']['constraint_tree'] else []
    params:
        model       = model_params,
        pars_trees  = get_highest_override( ['raxml-ng', 'treesearch'], "parsimony-starting-trees"),
        rand_trees  = get_highest_override( ['raxml-ng', 'treesearch'], "random-starting-trees"),
        bs_metric   = get_highest_override( ['raxml-ng', 'treesearch'], "bootstrap-metric" ),
        bs_trees    = bootstrap_params,
        data_type   = datatype,
        redo        = True
    threads:
        get_threads( ['raxml-ng', 'treesearch'] )
    output:
        best_tree       = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.newick",
        best_model      = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.model",
        support_tree    = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/bootstrap.newick",
        ml_trees        = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees.newick",
        bs_trees        = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/bs_trees.newick"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/search.log"
    conda:
        "../envs/raxml-ng.yaml"
    script:
        "../scripts/raxml-ng-search.py"

# =================================================================================================
#     Consensus Tree with RAxML-ng
# =================================================================================================

rule treesearch_consensus:
    input:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees.newick"
    output:
        mr      = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/consensusTreeMR.newick",
        mre     = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/consensusTreeMRE.newick"
    params:
        prefix  = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees"
    threads:
        get_threads( ['raxml-ng', 'treesearch'] )
    log:
        mr      = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/mr.log",
        mre     = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/mre.log"
    conda:
        "../envs/raxml-ng.yaml"
    shell:
        "raxml-ng --consense MR  --tree {input} --prefix {params.prefix} --threads {threads} > {log.mr}  2>&1 && "
        "mv {params.prefix}.raxml.consensusTreeMR {output.mr} && "
        "raxml-ng --consense MRE --tree {input} --prefix {params.prefix} --threads {threads} > {log.mre} 2>&1 && "
        "mv {params.prefix}.raxml.consensusTreeMRE {output.mre}"

rule determine_best_run:
    """Compare all best trees of a sample and create a symlink to that folder"""
    input:
        expand(
            "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.newick",
            autoref=autoref_list,
            aligner=aligner_list,
            trimmer=trimmer_list,
            allow_missing=True
            ),
        expand(
            "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible.consensusTreeMR.newick",
            autoref=autoref_list,
            aligner=aligner_list,
            trimmer=trimmer_list,
            allow_missing=True
            ),
        expand(  
            "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible.consensusTreeMRE.newick",
            autoref=autoref_list,
            aligner=aligner_list,
            trimmer=trimmer_list,
            allow_missing=True
            )
    output:
        "{outdir}/result/{sample}/best_result/raxml-ng/tree/best.newick",
        "{outdir}/result/{sample}/best_result/raxml-ng/post/plausible.consensusTreeMR.newick",
        "{outdir}/result/{sample}/best_result/raxml-ng/post/plausible.consensusTreeMRE.newick"
    script:
        "../scripts/symlink-best-result.py"

localrules: determine_best_run, get_constraint
