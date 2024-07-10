#! /usr/bin/env bash

# dumb approach go through all

RST_DIR="./lammps_docs_cleaned"
MD_DIR="./lammps_docs_md"


# Has a flat structure so easy to convert

# TODO: Use GNU-parrallel 
for file in "$RST_DIR"/*.rst ; do
	echo "$file"
	# Just the file name, not full path
	filename=$(basename "$file")
	# bash expansion
	file_prefix="${filename%.*}"
	out="${MD_DIR}/${file_prefix}.md"

	pandoc -o "$out"  "$file"
done




