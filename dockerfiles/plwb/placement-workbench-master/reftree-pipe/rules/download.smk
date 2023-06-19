# =================================================================================================
#     Obtaining the input sequences from genbank
# =================================================================================================

rule download_sequences:
    input:
        get_accessions
    # params:
    #     # tell the download script whether it should expect to find tsv files alongside the csv files, 
    #     # specifying a PhAT-conformant taxonomy per input taxon
    #     copy_tsvs_over = use_phat
    output:
    	"{outdir}/result/{sample}/download/seqs.fa"
    log:
        "{outdir}/result/{sample}/download/dl.log"
    conda:
        "../envs/download.yaml"
    script:
        "../scripts/download_fasta.py"
localrules: download_sequences
