# =================================================================================================
#     Statistical tests
# =================================================================================================

def modelstring_params( wildcards, input, output ):
    import re
    modelstring=""
    with open(input.best_model) as modelfile:
        lines = modelfile.readlines()
        if len(lines) != 1:
            raise InputError("Modelfile has weird number of lines")
        modelstring=re.sub( r",.*", "", re.sub(r"{.*}", "", lines[0]) ).rstrip()
        # fix for nonexistent "G4m"-like
        modelstring=re.sub( r"m$", "", modelstring )
    return modelstring

# IQTREE fixes non-unicode chars in the MSA, but not in the tree... so here we go
rule sanitize_labels:
    group: "postsearch"
    input:
        best_tree   = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.newick",
        ml_trees    = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees.newick",
        msa         = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/trimmed.afa"
    output:
        best_tree   = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/sanitized_best.newick",
        ml_trees    = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/sanitized_ml_trees.newick",
        msa         = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/sanitized_msa.afa",
    params:
        basedir = workflow.current_basedir
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/sanitizer.log"
    conda:
        "../envs/biopython.yaml"
    shell:
        "{params.basedir}/../scripts/sanitize_labels.py {input.best_tree} {output.best_tree} newick > {log}"
        " && {params.basedir}/../scripts/sanitize_labels.py {input.ml_trees} {output.ml_trees} newick > {log}"
        " && {params.basedir}/../scripts/sanitize_labels.py {input.msa} {output.msa} fasta > {log}"

rule iqtree_stats_test:
    group: "postsearch"
    input:
        msa         = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/sanitized_msa.afa",
        best_tree   = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/sanitized_best.newick",
        best_model  = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.model",
        ml_trees    = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/sanitized_ml_trees.newick"
    output:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/stats.iqtree"
    threads:
        get_threads( "iqtree" )
    params:
        workdir     = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post",
        modelstring = modelstring_params
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/iqtree.log"
    conda:
        "../envs/iqtree.yaml"
    shell:
        "iqtree -s {input.msa} -te {input.best_tree} -z {input.ml_trees}"
        " -m {params.modelstring}"
        " -pre {params.workdir}/stats -T {threads}"
        " -n 0 -zb 1000 -zw -au -treediff -redo > {log}"

rule summarize_iqtree_stats_test:
    group: "postsearch"
    input:
        iqtree_stats    = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/stats.iqtree",
        ml_trees        = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees.newick"
    output:
        summary         = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/significance.txt",
        plausible_trees = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible_trees.newick"
    script:
        "../scripts/iqtree_test_summarize.py"

rule plausible_consensus:
    group: "postsearch"
    input:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible_trees.newick"
    output:
        mr  = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible.consensusTreeMR.newick",
        mre = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible.consensusTreeMRE.newick"
    params:
        prefix  = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible"
    log:
        mr  = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/mr.log",
        mre = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/mre.log"
    conda:
        "../envs/raxml-ng.yaml"
    shell:
        "raxml-ng --consense MR  --tree {input} --prefix {params.prefix} --redo > {log.mr}  2>&1 && "
        "mv {params.prefix}.raxml.consensusTreeMR {output.mr} && "
        "raxml-ng --consense MRE --tree {input} --prefix {params.prefix} --redo > {log.mre} 2>&1 && "
        "mv {params.prefix}.raxml.consensusTreeMRE {output.mre}"

localrules: summarize_iqtree_stats_test, plausible_consensus, sanitize_labels
