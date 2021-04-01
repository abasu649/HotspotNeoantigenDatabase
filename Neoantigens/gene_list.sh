#!/bin/bash

for filename in /Users/anamikabasu/Code/Database/Neoantigens/file_col_with_mut_name_1-5/*.tsv; do
	awk 'NR==2{print $1}' $filename >> gene_list.tsv
done
