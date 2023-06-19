# =================================================================================================
#     Chunkify and Unchunkify
# =================================================================================================

# We want the sample names of the output files to fit the naming given by the user.
# In case our input for the samples is a table with two columns, the sample names (first column)
# can be different from the file names. However, the gappa chunkify command uses file names as
# sample names. For now, it is easiest to just symlink to the files to get a list of properly
# named files. In the future, we might add an option to rename samples in gappa chunkify,
# to make this step a bit less convoluted...
rule chunkify_sample_prep:
    input:
        fasta = "{outdir}/{clusterer}/samples/{sample}/queries.fa"
    output:
        "{outdir}/{clusterer}/chunkify/samples/{sample}.fasta"
    log:
        "{outdir}/{clusterer}/chunkify/samples/{sample}.log"
    script:
        "../../common/symlink.py"
# No need to execute this on the cluster computed nodes.
localrules: chunkify_sample_prep

# The rule to chunkify input fasta samples (query sequences) into chunks of equal size without
# duplicate sequences.
# We need a snakemake checkpoint here, because we cannot predict the number of chunks being produced.
checkpoint chunkify:
    input:
        # Request renamed samples, using the rule above, to get chunkify to use proper sample names.
        expand( "{outdir}/{clusterer}/chunkify/samples/{sample}.fasta",
                outdir=outdir,
                sample=sample_names,
                allow_missing=True
                )
    output:
        abundances  = expand(   "{outdir}/{clusterer}/chunkify/abundances/abundances_{sample}.json",
                                outdir=outdir,
                                sample=sample_names,
                                allow_missing=True
                                ),
        chunks_dir  = directory(expand( "{outdir}/{clusterer}/chunkify/chunks",
                                        outdir=outdir,
                                        allow_missing=True
                                        )
                                )
    params:
        chunks_dir      = expand(   "{outdir}/{clusterer}/chunkify/chunks",
                                    outdir=outdir,
                                    allow_missing=True
                                    ),
        abundances_dir  = expand(   "{outdir}/{clusterer}/chunkify/abundances",
                                    outdir=outdir,
                                    allow_missing=True
                                    ),
        hashfunction    = config["params"]["chunkify"]["hash-function"],
        minabun         = config["params"]["chunkify"]["min-abundance"],
        chunksize       = config["params"]["chunkify"]["chunk-size"]
    log:
        expand( "{outdir}/{clusterer}/chunkify/chunkify.log",
                outdir=outdir,
                allow_missing=True
                )
    threads:
        get_threads( "gappa" )
    conda:
        "../envs/gappa.yaml"
    shell:
        "mkdir -p {params.chunks_dir} ;"
        " mkdir -p {params.abundances_dir} ;"
        " gappa prepare chunkify"
        " --fasta-path {input}"
        " --chunks-out-dir {params.chunks_dir}"
        " --abundances-out-dir {params.abundances_dir}"
        " --hash-function {params.hashfunction}"
        " --min-abundance {params.minabun}"
        " --chunk-size {params.chunksize}"
        " --threads {threads}"
        " > {log} 2>&1"

# Following the documentation tutorial here:
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
def aggregate_chunkify_chunks(wildcards):
    # Here we get the path to the chunks dir from the output of the chunkify rule.
    # Since we have wildcards in play, we have to pass those that are relevant in the 
    # directory path when we use the get function
    chunks     = checkpoints.chunkify.get(**wildcards).output["chunks_dir"][0]

    # Next we return the paths to all chunks
    return expand(  "{outdir}/{clusterer}/chunkify/placed/{chunk}.jplace",
                    chunk = glob_wildcards( os.path.join(chunks, "{chunk}.fasta")).chunk,
                    outdir=outdir,
                    allow_missing=True
                    )

rule unchunkify:
    input:
        aggregate_chunkify_chunks,
        expand( "{outdir}/{clusterer}/chunkify/abundances/abundances_{sample}.json",
                outdir=outdir,
                sample=sample_names,
                allow_missing=True
                )
    output:
        protected(  expand( "{outdir}/{clusterer}/placed/{sample}.jplace",
                            outdir=outdir,
                            sample=sample_names,
                            allow_missing=True
                            )
                    )
    params:
        hash_function = config["params"]["chunkify"]["hash-function"],
        base_dir        = expand(   "{outdir}/{clusterer}",
                                    outdir=outdir,
                                    allow_missing=True
                                    ),
        # clusterer = expand("{clusterer}", clusterer=clusterer_list)[0]
    log:
        expand( "{outdir}/{clusterer}/chunkify/unchunkify.log",
                outdir=outdir,
                allow_missing=True
                )
    threads:
        get_threads( "gappa" )
    conda:
        "../envs/gappa.yaml"
    shell:
        "gappa prepare unchunkify"
        " --abundances-path {params.base_dir}/chunkify/abundances"
        " --chunk-file-expression {params.base_dir}/chunkify/placed/chunk_@.jplace"
        " --hash-function {params.hash_function}"
        " --out-dir {params.base_dir}/placed"
        " --threads {threads}"
        " > {log} 2>&1"
