import pandas as pd
import os

# read in file_allele_freqs.tsv into dictionary called allele_freqs_dict
allele_freqs_df = pd.read_csv(
    '/Users/anamikabasu/Code/Database/Neoantigens/file_allele_freqs.tsv', header=0, sep='\t')
allele_freqs_dict = dict(
    zip(allele_freqs_df.allele, allele_freqs_df.frequency))

input_dir = '/Users/anamikabasu/Code/Database/Neoantigens/file_col_with_mut_name_1-5'

open_files = [f for f in sorted(os.listdir(
    input_dir)) if f.endswith(".tsv")]

# for all files in the /Users/anamikabasu/Code/Database/Neoantigens/file_col_with_mut_name_1-5
for file in open_files:
    # upload file to dataframe called df_input
    df_input = pd.read_csv(f'{input_dir}/{file}', sep='\t')
    # empty list of freqs called file_freqs_list
    file_freqs_list = []
    # for all rows in the file
    for row in range(0, len(df_input.index)):
        allele = df_input.iloc[row, 7]
        file_freqs_list.append(allele_freqs_dict[allele])

    # add a new column to the dataframe with file_freqs_list
    df_input.insert(loc=8, column="HLA Allele Frequency",
                    value=file_freqs_list)

    # convert the file_df to tsv file stored in a new folder
    output_name = f'/Users/anamikabasu/Code/Database/Neoantigens/file_allele_freq_col_added/{file}'
    df_input.to_csv(output_name, index=False, sep='\t')
