#!/usr/bin/env python
import os.path
import pandas as pd

configfile: "configuration/config.yaml"
sample_sheet = pd.read_csv(config["sample_sheet"],
    dtype=str,
    names = ["sample", "asm"]).set_index("sample")

wildcard_constraints:
    sample = "|".join(sample_sheet.index)

include: "rules/utils.smk"
include: "rules/1.prepare_hmms.smk"
include: "rules/2.extract_UCEs_and_filter.smk"
include: "rules/3.organize_sequences.smk"
include: "rules/4.reconstruct_UCE_trees.smk"
include: "rules/5.reconstruct_species_tree.smk"

rule all:
    input:
        os.path.join(config["outdir"],"species_tree", "speciesTree.txt")
    shell:
        """
        echo "Job done!"
        """