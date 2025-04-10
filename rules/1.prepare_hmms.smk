rule nexus_to_fasta:
    output:
        os.path.join(config["outdir"], "fasta", "done.txt")
    params:
        nexus_dir=config["uce_dir"],
        outdir=os.path.join(config["outdir"], "fasta"),
        script="scripts/nexus_to_fasta.sh"
    log:
        os.path.join(config["outdir"],"logs", "nexus_to_fasta.log")
    threads:
        config["cpus"]
    shell:
        """
            find {params.nexus_dir} -type f -name '*.nexus' -print0 | \
                xargs -0 -P {threads} \
                    -I {{}} \
                    bash -c \
                        'bash {params.script} "{{}}" {params.nexus_dir} {params.outdir}' \
                         >{log} 2>{log}
            touch {output}
        """

rule build_hmm:
    conda:
        os.path.join(workflow.basedir,"envs/hmmer.yaml")
    input:
        os.path.join(config["outdir"], "fasta", "done.txt")
    output:
        os.path.join(config["outdir"], "hmm", "done.txt")
    params:
        fasta_dir=os.path.join(config["outdir"], "fasta"),
        outdir=os.path.join(config["outdir"], "hmm"),
        script="scripts/build_hmm.sh"
    log:
        os.path.join(config["outdir"],"logs", "build_hmm.log")
    threads:
        config["cpus"]
    shell:
        """
            find {params.fasta_dir} -type f -name '*.fasta' -print0 | \
                xargs -0 -P {threads} \
                    -I {{}} \
                    bash -c \
                        'bash {params.script} "{{}}" {params.fasta_dir} {params.outdir}' \
                         >{log} 2>{log}
            touch {output}
        """

rule merge_hmm:
    input:
        os.path.join(config["outdir"], "hmm", "done.txt")
    output:
        os.path.join(config["outdir"], "hmm_merged", "uce-all.hmm")
    params:
        hmm_dir=os.path.join(config["outdir"], "hmm")
    log:
        os.path.join(config["outdir"],"logs", "merge_hmm.log")
    threads:
        1
    shell:
        """
            cat {params.hmm_dir}/*.hmm > {output} 2> {log}
        """