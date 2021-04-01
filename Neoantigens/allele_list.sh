#!/bin/bash

for filename in /Users/anamikabasu/Code/Database/Neoantigens/file_col_with_mut_name_1-5/*.tsv; do
	cut -f8 $filename | uniq >> allele_list.tsv
done
