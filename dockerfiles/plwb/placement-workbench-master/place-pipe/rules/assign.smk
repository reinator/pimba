# =================================================================================================
#     Perform Taxonomic Assignment
# =================================================================================================

rule assign:
    group: "postplacement"
    input:
        jplace      = "{outdir}/{clusterer}/placed/{sample}.jplace",
        taxon_file  = config["data"]["taxonomy-file"]
    output:
        "{outdir}/{clusterer}/taxonomic_assignment/gappa_assign/{sample}/profile.tsv",
        "{outdir}/{clusterer}/taxonomic_assignment/gappa_assign/{sample}/krona.profile"
        if make_krona_plots else []
    params:
        allow_file_overwriting = True
    log:
        "{outdir}/{clusterer}/taxonomic_assignment/gappa_assign/{sample}/gappa-assign.log"
    threads:
        get_threads( "gappa" )
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-assign.py"

rule krona_plot:
    group: "postplacement"
    input:
        "{outdir}/{clusterer}/taxonomic_assignment/gappa_assign/{sample}/krona.profile"
    output:
        "{outdir}/{clusterer}/taxonomic_assignment/gappa_assign/{sample}/krona.html"
    log:
        "{outdir}/{clusterer}/taxonomic_assignment/gappa_assign/{sample}/krona.log"
    conda:
        "../envs/krona.yaml"
    shell:
        "ktImportText -o {output} {input} > {log} 2>&1"

# No need to execute this on the cluster computed nodes.
localrules: assign
