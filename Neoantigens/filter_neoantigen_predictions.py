"""
Usage: python3 filter_neoantigen_predictions.py <input_dir> <output_dir>
Parses hotspot pVACtools results files for further analysis
"""

import os
import sys
import re
import pandas as pd


def identify_neoantigens(df):
    """Identify class I and II neoantigens from pVACtools total predictions df
    Args:
        df (dataframe): pVACtools total predictions df
    Returns:
        df: neoantigens df
    """
    df = df.drop_duplicates()
    df_neo = df.loc[df["Median MT Score"] <= 500]
    for row_index, row in df_neo.iterrows():
        if "HLA-A" in row['HLA Allele'] or "HLA-B" in row['HLA Allele'] or "HLA-C" in row['HLA Allele']:
            wt_score, fold_change, mutation_position, seq_length = (
                row["Median WT Score"], row["Median Fold Change"], row["Mutation Position"], row["Peptide Length"])
            # class I - first and last 3 sequence positions are considered anchor positions
            anchor_positions = list(range(
                1, seq_length+1))[0:3] + list(range(1, seq_length+1))[(seq_length-3):seq_length]
            if fold_change <= 1:
                # if mutation in an anchor position, drop neo
                if mutation_position in anchor_positions:
                    df_neo = df_neo.drop(row_index)
            # mt binds better than wt
            if fold_change > 1:
                # if mutation is in anchor position
                if mutation_position in anchor_positions:
                    # but wt score is a good binder, drop neo
                    if wt_score < 500:
                        df_neo = df_neo.drop(row_index)
    return df_neo


input_dir = sys.argv[1]
output_dir = sys.argv[2]

open_files = [f for f in sorted(os.listdir(
    input_dir)) if f.endswith("_TUMOR.all_epitopes.tsv")]
for file in open_files:
    print(file)
    df_input = pd.read_csv(f'{input_dir}/{file}', sep='\t')
    df_input = df_input.drop(
        df_input.columns[[5, 6, 7, 12, 13, 20, 21, 22, 23]], axis=1)
    df_filtered = identify_neoantigens(df_input)
    output_file = f'{output_dir}/filtered_{file}'
    df_filtered.to_csv(output_file, index=False, sep='\t')
