import csv
import pandas as pd
import math

# helper methods

# gets dict key from value


def get_key(val):
    # sample name in gene_expression_matrix.tsv: TCGA-OR-A5J1-01A-11R-A29S-07
    # sample name in nontumor_mapping.tsv: TCGA-OR-A5J1-01
    sampleName = "-".join(val.split("-", 4)[:4])
    valueToCheck = sampleName[0:len(sampleName) - 1]
    for key, value in sample_dict.items():
        if valueToCheck in value:
            return key
    return "key doesn't exist"


def avg_list(lst):
    return sum(lst) / len(lst)


def read_sample_file(filename: str) -> dict:
    """Group together samples of same cancer type in dictionary
    Args:
        string: filename of file containing samples and corresponding cancer types
    Returns:
        dict: dictionary {cancer: [sample, sample,...], cancer: [sample, sample,...]...}
    """
    sample_dict = dict()
    with open(filename) as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')
        for row in reader:
            # Check if acronym already exists as a key
            if row[1] in sample_dict.keys():
                # Add sample as value to acronym key
                sample_dict[row[1]].append(row[0])
            else:
                # If acronym doesn't already exist as a key, add acronym and sample as key-value pair
                sample_dict[row[1]] = [row[0]]
    return sample_dict


sample_dict = read_sample_file(
    '/Users/anamikabasu/Code/Database/TCGA/nontumor_mapping.tsv')
cancer_list = list(sample_dict.keys())
gene_df = pd.read_csv(
    "/Users/anamikabasu/Code/Database/Neoantigens/gene_list.tsv", header=None, sep='\t')

gene_list = list(gene_df.iloc[:, 0])
gene_list = list(set(gene_list))


gene_expression_df = pd.read_csv(
    "/Users/anamikabasu/Code/Database/TCGA/gene_expression_matrix.tsv", sep='\t', header=0)

print("necessary files in dfs")


####gene_expression_df data cleaning: replacing cells with log2(cell+1), reformatting genes from "A2M|2" to A2M####

# reformatting genes from ?| 100130426 to 100130426 and A2M | 2 to A2M
for r in range(0, len(gene_expression_df.index)):
    if "?" in gene_expression_df.iloc[r, 0]:
        gene_expression_df.iloc[r, 0] = (
            gene_expression_df.iloc[r, 0]).split('|')[1]
    else:
        gene_expression_df.iloc[r, 0] = (
            gene_expression_df.iloc[r, 0]).split('|')[0]

print("reformatted gene names")


# # gene_expression_df data cutting
gene_expression_df = gene_expression_df[gene_expression_df['gene_id'].isin(
    gene_list)]

gene_expression_df = gene_expression_df.set_index("gene_id")

print("gene_expression_df cut")

# replacing cells with log2(cell+1)
for r in range(0, len(gene_expression_df.index)):
    for c in range(0, len(gene_expression_df.columns)):
        gene_expression_df.iat[r, c] = math.log(
            float(gene_expression_df.iloc[r, c]) + 1, 2)


print("gene_expression_df log transformed")

####populating avg_expression_df####

# creating empty avg_expression_df
avg_expression_df = pd.DataFrame(columns=cancer_list, index=gene_list)

avg_expression_df.to_csv(
    '/Users/anamikabasu/Code/Database/TCGA/average.tsv', index=True, sep='\t')


# for each cell in gene_expression_df, added it to correct key in exp_dict, which is averaged and placed in avg_expression_df
for gene in gene_expression_df.index:
    # exp_dict = {cancer1: None, cancer2: None, ... }
    exp_dict = dict.fromkeys(cancer_list)
    for sample in gene_expression_df.columns:
        if exp_dict[get_key(sample)] == None:
            exp_dict[get_key(sample)] = [
                gene_expression_df.loc[[gene], [sample]].iloc[0, 0]]
        else:
            exp_dict[get_key(sample)].append(
                gene_expression_df.loc[[gene], [sample]].iloc[0, 0])
    # exp_dict = {cancer1: [exp1, exp2, ...], cancer2: [exp1, exp2, ...], ... }
    for cancer_key in exp_dict.keys():
        if exp_dict[cancer_key] != None:
            avg_exp = avg_list(exp_dict[cancer_key])
            avg_expression_df.loc[[gene], [cancer_key]] = avg_exp


print("populating avg_expression_df done")

####converting avg_expression_df to tsv file####
avg_expression_df.to_csv(
    '/Users/anamikabasu/Code/Database/TCGA/average_gene_expression_matrix_chr1-5.tsv', index=True, sep='\t')


#####Don't need anymore####
# stripping sample names of quotes (Don't need this)
# stripped_header = [gene_expression_df.iloc[0, c].strip(
#     '"') for c in range(0, len(gene_expression_df.columns))]
# gene_expression_df.iloc[0, :] = stripped_header

# only reformatting genes from A2M|2 to A2M (Don't need this)
# stripped_column = [(gene_expression_df.iloc[r, 0]).split('|')[
#     0] for r in range(0, len(gene_expression_df.index))]
# gene_expression_df.iloc[:, 0] = stripped_column

# getting gene_expression_df ready for indexing
# set the first row to be the header
# gene_expression_df.columns = gene_expression_df.iloc[0, :]

# remove the first row after the header
# gene_expression_df = gene_expression_df.iloc[1:]
# set the first column to be index
# gene_expression_df.set_index('gene_id', inplace=True)
