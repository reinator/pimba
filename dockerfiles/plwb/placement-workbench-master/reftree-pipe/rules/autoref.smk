rule no_phat:
    input:
        get_fasta
    output:
        "{outdir}/result/{sample}/no_autoref/ref_candidates.fa"
    log:
        "{outdir}/result/{sample}/no_autoref/log.txt"
    script:
        "../../common/symlink.py"
localrules: no_phat

# Rule to run gappa phat algorithm from a raw database of potential reference sequences
rule gappa_phat_raw:
    group: "alignment"
    input:
        # the taxonomy file is expected to be passed via the config, for now.
        taxonomy_file=get_taxonomy_file,
        # special call to get_fasta, as we are calling it from the optional PhAT step
        sequence_file=get_fasta
    output:
        "{outdir}/result/{sample}/phat/ref_candidates.fa"
    log:
        "{outdir}/result/{sample}/phat/log.txt"
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-phat.py"
