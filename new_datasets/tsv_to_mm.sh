#!/bin/bash

# Exit immediately if a command fails
set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file.tsv> <output_file.mtx>" >&2
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found." >&2
    exit 1
fi

if file "$INPUT_FILE" | grep -q 'gzip compressed data'; then
    echo "$INPUT_FILE is a gzip file. Extracting..."
    gunzip "$INPUT_FILE"
fi


echo "Processing matrix..."

#Get matrix dimensions using GNU utilities
TOTAL_ROWS=$(wc -l < "$INPUT_FILE")
TOTAL_COLS=$(head -n 1 "$INPUT_FILE" | awk -F'\t' '{print NF}')

#GNU-optimized non-zero count
TOTAL_NZ=$(tr -s '\t ' '\n' < "$INPUT_FILE" | grep -v -x '0' | grep -v -x '' | wc -l)

#Create/overwrite the output file with the MatrixMarket header
echo "%%MatrixMarket matrix coordinate integer general" > "$OUTPUT_FILE"
echo "$TOTAL_ROWS $TOTAL_COLS $TOTAL_NZ" >> "$OUTPUT_FILE"

#Stream data coordinates directly into the output file
awk -F'\t' -v OFS=" " '
{
  for (c=1; c<=NF; c++) {
    if ($c != 0 && $c != "") {
      print NR, c, $c
    }
  }
}' "$INPUT_FILE" >> "$OUTPUT_FILE"

echo "Done! Matrix successfully written to '$OUTPUT_FILE'"

