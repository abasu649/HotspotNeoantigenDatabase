#!/bin/bash

count=0
while read -r LINE 
do
	if [ $count -eq 0 ]
	then
	cat ./file_allele_freq_col_added/$LINE > /Users/anamikabasu/Code/Database/Neoantigens/merged_files/merged_file.tsv
	else 
	tail +2 ./file_allele_freq_col_added/$LINE >> /Users/anamikabasu/Code/Database/Neoantigens/merged_files/merged_file.tsv
	fi
	count=$((count+1))
done

