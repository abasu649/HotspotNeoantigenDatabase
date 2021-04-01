import pandas as pd

# read in unique_allele_list.tsv into file_alleles_df
file_alleles_df = pd.read_csv(
    '/Users/anamikabasu/Code/Database/Neoantigens/unique_allele_list.tsv', header=0, sep='\t')

# read in global_hla_freqs.tsv into database_alleles_dict
database_alleles_df = pd.read_csv(
    '/Users/anamikabasu/Code/Database/HLA/global_hla_freqs.tsv', header=0, sep='\t')
database_alleles_dict = dict(
    zip(database_alleles_df.allele, database_alleles_df.frequency))

# empty list of freqs to add into file_alleles_df
freqs_to_be_added = []
count = 0
# for row in file_alleles_df
for row in range(0, len(file_alleles_df.index)):
    allele = file_alleles_df.iloc[row, 0]
    if "HLA-" in file_alleles_df.iloc[row, 0]:
        allele = allele.split("HLA-")[1]
    if "-" in allele:
        allele1 = allele.split("-")[0]
        allele2 = allele.split("-")[1]
        if allele1 not in database_alleles_dict.keys() or allele2 not in database_alleles_dict.keys():
            count += 1
            freqs_to_be_added.append("NA")
        else:
            allele1_freq = database_alleles_dict[allele1]
            allele2_freq = database_alleles_dict[allele2]
            paired_freq = (allele1_freq * allele2_freq) / 100
            freqs_to_be_added.append(str(paired_freq))
    elif "\\" in allele:
        count += 1
        freqs_to_be_added.append("NA")
    else:
        if allele in database_alleles_dict.keys():
            freqs_to_be_added.append(str(database_alleles_dict[allele]))
        else:
            count += 1
            freqs_to_be_added.append("NA")


file_alleles_df.insert(loc=1, column="frequency", value=freqs_to_be_added)

file_alleles_df.to_csv(
    '/Users/anamikabasu/Code/Database/Neoantigens/file_allele_freqs.tsv', index=False, sep='\t')
