#!/bin/bash

usage() {
	echo "Usage: $0 -r <reference_fasta_path> -s <simulated_seqs_fasta_path>"
	echo "	<reference_fasta_path> Path to reference fasta"
	echo "	<simulated_seqs_fasta_path> Path to fasta output of simulation"
	exit 1
}

ref_path=''
seq_path=''

while getopts 'r:s:' flag; do
	case "${flag}" in
		r) ref_path=$OPTARG ;;
		s) seq_path=$OPTARG ;;
		*) usage ;;
	esac
done

# extract the "core_names" of ref_path and seq_path (i.e. strip extensions for downstream savenames)
ref_path_base="${ref_path%.fasta}"
seq_path_base="${seq_path%.fasta}"

# checks if paths are empty strings
if [ -z "${ref_path}" ] || [ -z "${seq_path}" ]; then
		usage

fi
if [ ! -f "$ref_path" ] || [ ! -f "$seq_path" ]; then
		echo "Both input paths must point to existing files." >&2
		exit 1
fi

# helper function that converts from fasta to fastq
fasta_to_fastq() {
  local input_fasta_path=$1
  local output_fastq_path=$2

	  awk '
	    function emit_record(quality) {
	      if (header == "") return
	      quality = sequence
	      gsub(/./, "I", quality)
	      print "@" header
	      print sequence
	      print "+"
	      print quality
	    }
	    /^>/ {
	      emit_record()
	      header = substr($0, 2)
	      sequence = ""
	      next
	    }
	    {
	      gsub(/[[:space:]]/, "", $0)
	      sequence = sequence $0
	    }
	    END {
	      emit_record()
	    }
	  ' "$input_fasta_path" > "$output_fastq_path"

}

# Call the function with input and output file paths
# convert_fasta_to_fastq "${ref_path}" "updated_output.fastq"


echo "Prepping IGV input files for ${ref_path}..."

# index the reference:
samtools faidx "$ref_path"

# convert simulation fasta to fastq using null quality scores
fastq_path="${seq_path_base}.fastq"
sam_path="${seq_path_base}.sam"
bam_path="${seq_path_base}.bam"
igv_input_bam_path="${seq_path_base}_IGV_input.bam"

# convert simulation output fasta to fastq
fasta_to_fastq "$seq_path" "$fastq_path"

# align
bwa mem "${ref_path}" "$fastq_path" > "$sam_path"

# sort
samtools view -h "$sam_path" | samtools view -Sb - > "$bam_path"
samtools sort "$bam_path" > "$igv_input_bam_path"
samtools index "$igv_input_bam_path"


