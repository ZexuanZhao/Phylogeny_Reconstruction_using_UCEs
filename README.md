# Phylogeny_Reconstruction_using_UCEs

## Description:
 - A snakemake pipeline to reconstruct phylogeny based on genome assemblies and UCE alignment files in nexus format.
 - Use `hmmer` to extract UCE loci in genome assemblies, and build UCE trees and a species tree.

## Files to prepare:
 - A sample sheet - sample_sheet.csv: a comma delimited file with 2 columns (no column name):
   - `sample` as the tip labels in the final species tree, `path to the genome assembly`
 - Modify configuration file - `configuration/config.yaml`:
   - `project`: a name for your project as the prefix of the final gene tree file.
     
   - `sample_sheet`: path to the sample sheet prepared above
   - `outdir`: path to the output directory
   - `uce_dir`: path to the directory that has the UCE alignments in nexues format
     
   - `cpus`: number of CPUs to use
   - `cpus_per_nhmmer`: number of CPUs to use for `hmmer` on each genome assembly
   - `n_parallel_gblocks`: number of jobs to run `gblocks` in parallel
   - `n_parallele_mafft`: number of jobs to run `mafft` in parallel
   - `n_parallele_iqtree`: number of jobs to run `iqtree` in parallel
   - `astral`: name of the astral excutable in the `\bin` folder
  
   - `window_size`: window size to visualize the hmm hits
   - `ratio_cutoff`: if the score of second hmmer hit is more than `ratio_cutoff`% of that of the first hit, drop both hits.
   - `min_score`: minimum bitscore to include a hmmer hit
   - `min_size`: minimum size of a hmmer hit to be used for building trees.
  
   - `p_Ns`: a sequence should have no more than `p_Ns` columns as Ns for building trees.

## Environment:
 - Make sure snakemake is installed in current environment.
   
## Usage:
`snakemake --use-conda --cores [ncpu]`
