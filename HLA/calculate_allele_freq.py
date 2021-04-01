import pandas as pd
import sys
import os
import numpy as np


# gene_allele_dict: {gene : {allele: count}, {allele: count} ..., gene: ...}
# gene_counts_dict: {gene: count, gene: count, ...}

gene_allele_dict = {}
gene_counts_dict = {}

input_dir = '/Users/anamikabasu/Code/Database/HLA/database_files'

open_files = [f for f in sorted(os.listdir(input_dir)) if f.endswith(".csv")]

# for csv file in /Users/anamikabasu/Code/Database/HLA/database_files
for file in open_files:
    # if csv file is not empty
    if os.path.getsize(f'{input_dir}/{file}') == 2:
        print(f'{file} is empty')
    else:
        df = pd.read_csv(f'{input_dir}/{file}', index_col=0, header=None)
        # for row in csv file
        for row in range(0, len(df.index)):
            # for col in csv file
            for col in range(0, len(df.columns)):
                # if cell has at least 1 identifier
                if ":" in str(df.iloc[row, col]):
                    # split after the second identifier (won't affect alleles with only 1 identifier)
                    allele_with_id = df.iloc[row, col]
                    allele_formatted = ":".join(
                        allele_with_id.split(":")[:2])
                    gene = allele_formatted.split("*")[0]
                    # update gene_count dictionary
                    if gene not in list(gene_allele_dict.keys()):
                        gene_allele_dict.update({gene: {}})
                        gene_counts_dict.update({gene: 0})
                    # if the allele is an inner key in the gene_allele_dict
                    gene_counts_dict.update({gene: gene_counts_dict[gene] + 1})
                    if allele_formatted in gene_allele_dict[gene].keys():
                        # update the count of the allele in gene_allele_dict
                        gene_allele_dict[gene].update(
                            {allele_formatted: gene_allele_dict[gene][allele_formatted] + 1})
                    else:
                        # add the allele to gene_allele_dict and then update the count
                        gene_allele_dict[gene].update({allele_formatted: 1})

print(gene_counts_dict)

# alleles = []
# freqs = []

# # for every gene in gene_allele_dict
# for gene in gene_allele_dict.keys():
#     # for every allele in gene
#     for allele in gene_allele_dict[gene].keys():
#         # calculate allele counts/total gene counts
#         freq = (gene_allele_dict[gene][allele]/gene_counts_dict[gene]) * 100
#         # add allele to list of alleles
#         alleles.append(allele)
#         # add frequency to list of frequencies
#         freqs.append(freq)


# # create a dataframe with list of alleles and list of frequencies
# output_dict = {'allele': alleles, 'frequency': freqs}
# df_output = pd.DataFrame(data=output_dict)

# df_output.to_csv(
#     '/Users/anamikabasu/Code/Database/HLA/global_hla_freqs.tsv', index=False, sep='\t')
# # convert dataframe to excel file
