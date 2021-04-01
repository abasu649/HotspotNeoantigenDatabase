import os
import sys
import re
import pandas as pd

input_dir = "/Users/anamikabasu/Code/Database/Neoantigens/variant_tsvs_1-5"
output_dir = "/Users/anamikabasu/Code/Database/Neoantigens/file_col_with_mut_name_1-5"

open_files = [f for f in sorted(os.listdir(
    input_dir)) if f.endswith("_TUMOR.all_epitopes.tsv")]

for file1 in open_files:
    df_input = pd.read_csv(f'{input_dir}/{file1}', sep='\t')
    df_input["Hotspot Mutation"] = (file1.split(
        '_TUMOR.all_epitopes.tsv')[0]).split('filtered_')[1]
    df_input = df_input.drop(
        ['Chromosome', 'Start', 'Stop', 'Reference', 'Variant'], axis=1)
    df_input = df_input.rename({'Mutation': 'Protein Change'}, axis=1)
    df_input = df_input[['Gene Name', 'Hotspot Mutation', 'Protein Change', 'Protein Position', 'Variant Type', "MT Epitope Seq", "WT Epitope Seq", "HLA Allele",
                         "Median MT Score", "Median WT Score", "Median Fold Change", "Peptide Length", "Sub-peptide Position", "Mutation Position"]]
    output_file = f'{output_dir}/{file1}'
    df_input.to_csv(output_file, sep='\t', index=False)
