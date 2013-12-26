import string
import sys
import os
import glob
import re
import csv, logging

PROTEOMIC_DATA_FILE = '../data/hsa_nagaraj/hsa_actin_tubulin_abundance.csv'
KO_TO_UNIPROT_FILE = '../data/KO_gene_hierarchy/KO_gene_hierarchy_organism_mapping/hsa_mapping.csv'
OUTPUT_PATH = '../res/hsa_nagaraj/'
if not os.path.exists(OUTPUT_PATH):
    os.mkdir(OUTPUT_PATH)


uniprot_known = set()
for row in csv.reader(open(KO_TO_UNIPROT_FILE ,'r'), delimiter='\t'):
    uniprot_known.add(row[0])
    
data_csv = csv.reader(open(PROTEOMIC_DATA_FILE, 'r'), delimiter='\t')
output_csv = csv.writer(open(OUTPUT_PATH + 'hsa_actin_tubulin_abundance.csv', 'w'), delimiter='\t')

not_mapped = 0
for row in data_csv:
    uniprot_semicolon, abundance = row
    abundance = float(abundance)
    
    uniprot_unique = None
    for uniprot in uniprot_semicolon.split(';'):
        if uniprot in uniprot_known:
            uniprot_unique = uniprot
            break

    if uniprot_unique is None:
        not_mapped += abundance
    else:
        output_csv.writerow([uniprot_unique, abundance])

output_csv.writerow(['NotMapped', abundance])
   
