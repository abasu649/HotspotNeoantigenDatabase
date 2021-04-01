import csv
import json
from typing import Tuple


def read_codes_file(filename: str, abbrevs: dict) -> Tuple[dict, set]:
    result_dict = dict()
    result_set = set()
    with open(filename) as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')
        for row in reader:
            result_dict[row[0]] = abbrevs[row[2]]
            result_set.add(row[2])
    return result_dict, result_set


def read_abbrevs_file(filename: str) -> set:
    result_dict = dict()
    with open(filename) as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')
        for row in reader:
            result_dict[row[1]] = row[0]
    return result_dict


def create_tumor_mapping(inputfilename, outputfilename):
    with open(inputfilename) as inputfile:
        with open(outputfilename, 'w') as outputfile:
            reader = csv.reader(inputfile, delimiter='\t')
            for row in reader:
                sample = row[0].split('.')[0]
                sample_code = sample.split('-')[1]
                tumor_type = code_dict[f'{sample_code}']
                outputfile.write(f'{sample}\t{tumor_type}\n')


study_dict = read_abbrevs_file(
    '/Users/anamikabasu/Code/Database/TCGA/StudyAbbrevs.tsv')
code_dict, code_set = read_codes_file(
    '/Users/anamikabasu/Code/Database/TCGA/TSS_codes.tsv', study_dict)
create_tumor_mapping('/Users/anamikabasu/Code/Database/TCGA/nonprimary_tumors.tsv',
                     '/Users/anamikabasu/Code/Database/TCGA/nontumor_mapping.tsv')
print(code_dict['EY'])
