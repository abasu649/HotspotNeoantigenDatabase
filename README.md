# HotspotNeoantigenDatabase

**Anchor Position Filtering** 

filter_neoantigen_prediction.py

**AAdding filenames as a column for hotspot mutation**A 

add_file_name_col.py

**ATCGA averaged gene expression table**A 

ID_to_cancertype.py // got cancer code for each sample from pancan datafile 

Gene_list.sh //gets list of genes from neoantigen files, duplication

average_gene_expression_matrix.py //creates 33)(cancers) x 64 (genes in 1-5) averaged gene 

**AAdding HLA Frequencies to Columns**A 

curl http://www.allelefrequencies.net/BrowseGenotype.aspx -Ls | 
  pup '.table01 tbody tr td:nth-child(1) text{}' | xargs -I {} wget http://www.allelefrequencies.net/tools/getrawdata.asp'?'pop_id={}'&resolved=true' -O {}.csv

gene_list.sh //gets list of genes from neoantigen files, duplication

allele_list.sh //gets list of alleles from neoantigen files, duplication 

calculate_allele_freq.py //calculates frequencies for alleles in database 

add_freqs_to_unique_alleles_list.py //gets frequencies for all alleles that in the neoantigen files

add_hla_freqs.py //adds column of hla allele frequencies to all files 

**AMerging Files**

merge_files.sh
