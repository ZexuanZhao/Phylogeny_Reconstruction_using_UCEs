rule nhmmer:
    conda:
        os.path.join(workflow.basedir,"envs/hmmer.yaml")
    input:
        unpack(get_asm),
        hmm=os.path.join(config["outdir"], "hmm_merged", "uce-all.hmm"),
    output:
        os.path.join(config["outdir"], "hmm_hits", "{sample}.hits.txt")
    log:
        os.path.join(config["outdir"], "logs", "nhmmer", "{sample}_nhmmer.log")
    threads:
        config["cpus_per_nhmmer"]
    shell:
        """
            nhmmer \
                --cpu {threads} \
                --tblout {output} \
                {input.hmm} {input.asm} \
                >{log} 2>{log}
        """

rule filter_hmm_hits:
    conda:
        os.path.join(workflow.basedir,"envs/R.yaml")
    input:
        unpack(get_asm),
        hits = os.path.join(config["outdir"], "hmm_hits", "{sample}.hits.txt")
    output:
        os.path.join(config["outdir"], "hmm_hits_filtered", "{sample}.done.txt")
    params:
        script="scripts/filter_hmm_hits.R",
        sample="{sample}",
        outdir=os.path.join(config["outdir"], "hmm_hits_filtered"),
        figuredir=os.path.join(config["outdir"], "logs", "hmm_filter",),
        window_size=config["window_size"],
        ratio_cutoff=config["ratio_cutoff"],
        min_score=config["min_score"],
        min_size=config["min_size"],
    log:
        os.path.join(config["outdir"], "logs", "hmm_filter", "{sample}_hmm_filter.log")
    threads:
        1
    shell:
        """
            Rscript {params.script} \
                {input.hits} {input.asm} {params.sample} {params.outdir} {params.figuredir} \
                {params.window_size} {params.ratio_cutoff} {params.min_score} {params.min_size} \
                >{log} 2>{log}
            touch {output}
        """