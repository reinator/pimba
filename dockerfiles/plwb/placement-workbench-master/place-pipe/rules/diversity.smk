# =================================================================================================
#     Calculate Diversity Metrics
# =================================================================================================

rule diversity_guppy:
    group: "postplacement"
    input:
        "{outdir}/{clusterer}/placed/{sample}.jplace"
    output:
        "{outdir}/{clusterer}/diversity/guppy_fpd/{sample}.csv"
    params:
        csv = True
    log:
        "{outdir}/{clusterer}/diversity/guppy_fpd/{sample}.log"
    conda:
        "../envs/pplacer.yaml"
    script:
        "../scripts/guppy-fpd.py"

# # No need to execute this on the cluster computed nodes.
# localrules: diversity_guppy
