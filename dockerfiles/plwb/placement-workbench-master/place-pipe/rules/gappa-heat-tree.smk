# =================================================================================================
#     gappa heat-tree
# =================================================================================================

# Rule to create a heat tree visualization for all samples combined
rule gappa_heat_tree_all:
    group: "postplacement"
    input:
        expand( "{outdir}/{clusterer}/placed/{sample}.jplace",
                allow_missing=True,
                sample=sample_names
                )
    output:
        expand( "{outdir}/{clusterer}/heat-tree.{ext}",
                allow_missing=True,
                ext=config["params"]["gappa"]["heat-tree"]["formats"]
                )
    params:
        allow_file_overwriting = True
    log:
        "{outdir}/{clusterer}/heat_tree_all.log"
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-heat-tree.py"

# Rule to create individual heat tree visualizations per sample
rule gappa_heat_tree:
    group: "postplacement"
    input:
        rules.epa_ng_place.output
    output:
        expand( "{outdir}/{clusterer}/heat-trees/{sample}/heat-tree.{ext}",
                allow_missing=True,
                ext=config["params"]["gappa"]["heat-tree"]["formats"]
                )
    params:
        allow_file_overwriting = True
    log:
        "{outdir}/{clusterer}/heat_trees/{sample}.log"
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-heat-tree.py"
