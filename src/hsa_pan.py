import string
import sys
import os
import glob
import re
import csv, logging

ENSEMBLE_TO_UNIPROT_FILE = '../data/hsa_pan/ensemble_to_uniprot.csv'
KO_TO_UNIPROT_FILE = '../data/KO_gene_hierarchy/KO_gene_hierarchy_organism_mapping/hsa_mapping.csv'
PRIMATE_PROTEOME_FILE = '../data/hsa_pan/hsa_pan_abundance.csv'
OUTPUT_PATH = '../res/hsa_pan/'

uniprot_known = set()
for row in csv.reader(open(KO_TO_UNIPROT_FILE ,'r'), delimiter='\t'):
    uniprot_known.add(row[0])

ensg2uniprot_list = {}
for row in csv.reader(open(ENSEMBLE_TO_UNIPROT_FILE, 'r'), delimiter='\t'):
    ensg, uniprot = row
    ensg2uniprot_list.setdefault(ensg, []).append(uniprot)

ensg2uniprot_semicolon = {}
ensg2uniprot_unique = {}
for ensg, uniprot_list in ensg2uniprot_list.iteritems():
    ensg2uniprot_semicolon[ensg] = ';'.join(uniprot_list)
    for uniprot in uniprot_list:
        if uniprot in uniprot_known:
            ensg2uniprot_unique[ensg] = uniprot
            break


data_csv = csv.reader(open(PRIMATE_PROTEOME_FILE, 'r'), delimiter='\t')
titles = data_csv.next()

output_csv = csv.writer(open(OUTPUT_PATH + 'all_data_with_uniprot.csv', 'w'), delimiter='\t')
output_csv.writerow(['UniProt', 'full UniProt list'] + titles)

dataset_names = titles[1:]
uniprot2abundance = dict([(n, {}) for n in titles])

for row in data_csv:
    ensg = row.pop(0)
    if ensg not in ensg2uniprot_unique:
        logging.warning('The EnsembleID [%s] is not in the mapping table', ensg)

    uniprot_unique = ensg2uniprot_unique.get(ensg, 'NotMapped')
    uniprot_semicolon = ensg2uniprot_semicolon.get(ensg, '')
    output_csv.writerow([uniprot_unique, uniprot_semicolon, ensg] + row)

    for i, dataset_name in enumerate(dataset_names):
        uniprot2abundance[dataset_name].setdefault(uniprot_unique, 0)
        uniprot2abundance[dataset_name][uniprot_unique] += float(row[i])
        
for dataset_name in dataset_names:
    output_csv = csv.writer(open(OUTPUT_PATH + dataset_name + '.csv', 'w'), delimiter='\t')
    output_csv.writerows(sorted(uniprot2abundance[dataset_name].iteritems()))
    
