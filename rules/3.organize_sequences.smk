rule merge_sequences:
    input:
        expand(os.path.join(config["outdir"], "hmm_hits_filtered", "{sample}.done.txt"), sample = sample_sheet.index)
    output:
        os.path.join(config["outdir"], "hmm_hits_filtered_merged", "done.txt")
    params:
        fasta_dir=os.path.join(config["outdir"], "fasta"),
        hits_dir=os.path.join(config["outdir"], "hmm_hits_filtered"),
        merged_dir=os.path.join(config["outdir"], "hmm_hits_filtered_merged"),
        samples=" ".join(sample_sheet.index)
    log:
        os.path.join(config["outdir"], "logs", "merge_sequences.log")
    threads:
        1
    shell:
        """
            IFS=' ' read -r -a s_array <<< "{params.samples}"
            for file in $(ls {params.fasta_dir}/*.fasta); do
              file=$(basename "${{file}}")
              cat {params.fasta_dir}/${{file}} > {params.merged_dir}/${{file%.fasta}}_merged.fasta
              for s in "${{s_array[@]}}"; do
                if [ -f {params.hits_dir}/${{s}}/${{file}} ]; then
                  cat {params.hits_dir}/${{s}}/${{file}} >> {params.merged_dir}/${{file%.fasta}}_merged.fasta
                else
                  echo ${{file}} not found in ${{s}} >{log}
                fi
              done
            done
            touch {output}
        """

