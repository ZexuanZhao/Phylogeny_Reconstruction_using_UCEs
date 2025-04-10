#!/bin/bash

# Check if input parameter is provided
if [ -z "$1" ]; then
    echo "Error: No input file specified" >&2
    echo "Usage: $0 <nexus_file> <input_dir> <output_dir>" >&2
    exit 1
fi

# Check if input parameter is provided
if [ -z "$2" ]; then
    echo "Error: No input dir specified" >&2
    echo "Usage: $0 <nexus_file> <input_dir> <output_dir>" >&2
    exit 1
fi

# Check if input parameter is provided
if [ -z "$3" ]; then
    echo "Error: No output dir specified" >&2
    echo "Usage: $0 <nexus_file> <input_dir> <output_dir>" >&2
    exit 1
fi

# Set input and output dir
nexus_dir=$2
outdir=$3

# Extract locus name from input file
locus=$(basename "$1")
locus="${locus%.nexus}"

# Set input and output paths
input="${nexus_dir}/${locus}.nexus"
output="${outdir}/${locus}.fasta"

# Check if input file exists
if [ ! -f "$input" ]; then
    echo "Error: Input file $input not found" >&2
    exit 1
fi

# Process the file
awk 'FNR >5' "$input" | \
  tr -d "'" | \
  tr " " "\n" | \
  sed -e 's/^FW/>FW/g' \
      -e 's/^W/>W/g' \
      -e 's/?/N/g' \
      -e 's/^;$//g' \
      -e 's/^end;//g' \
      -e '/^$/d' > "$output"

# Check if processing was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to process $input" >&2
    exit 1
fi

echo "Successfully converted $input to $output"