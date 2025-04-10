# Install Biostrings if you haven't already
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("Biostrings")
}

# Load the Biostrings package
library(Biostrings)

# Define the function to filter sequences based on N content
filter_by_n_proportion <- function(input_fasta, max_n_proportion, output_dir) {
  # Check if the input file exists
  if (!file.exists(input_fasta)) {
    stop(paste("Input FASTA file not found:", input_fasta))
  }

  # Read the FASTA file into a DNAStringSet object
  sequences <- readDNAStringSet(input_fasta)

  # Calculate the proportion of Ns in each sequence
  n_proportions <- vcountPattern("N", sequences) / width(sequences)

  # Identify sequences to keep (where N proportion is less than or equal to the threshold)
  keep_indices <- which(n_proportions <= max_n_proportion)
  filtered_sequences <- sequences[keep_indices]

  # Construct the output file name
  base_name <- sub(".fsa$", "", basename(input_fasta))
  output_fasta <- file.path(output_dir, paste0(base_name, ".filterN.fasta"))

  # Write the filtered sequences to the output FASTA file
  writeXStringSet(filtered_sequences, file = output_fasta)

  cat(paste("Filtered", length(sequences) - length(filtered_sequences),
            "sequences with more than", max_n_proportion,
            "proportion of Ns from", basename(input_fasta), "\n"))
  cat(paste("Filtered sequences written to:", output_fasta, "\n"))
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 3) {
  stop("Usage: Rscript filter_sequences_with_Ns.R <max_n_proportion> <input_fasta_file> <output_directory>")
}

# Assign arguments to variables
max_n_proportion_input <- as.numeric(args[1])
input_fasta_file <- args[2]
output_directory <- args[3]

# Basic input validation for the N proportion
if (is.na(max_n_proportion_input) || max_n_proportion_input < 0 || max_n_proportion_input > 1) {
  stop("The maximum N proportion must be a decimal number between 0 and 1.")
}

# Check if the output directory exists and create it if it doesn't
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
  cat(paste("Output directory created:", output_directory, "\n"))
}

# Call the filtering function
filter_by_n_proportion(input_fasta_file, max_n_proportion_input, output_directory)