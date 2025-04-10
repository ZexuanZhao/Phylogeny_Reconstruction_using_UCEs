rule astral:
    input:
        os.path.join(config["outdir"],"gene_trees", "merged_gene_trees.txt")
    output:
        os.path.join(config["outdir"],"species_tree", "speciesTree.txt")
    log:
        os.path.join(config["outdir"],"logs","astral.log")
    params:
        astral = config["astral"]
    threads:
        config["cpus"]
    shell:
        """
            bin/{params.astral} \
                -t {threads} \
                -o {output} \
                {input} \
                >{log} 2>{log}
        """