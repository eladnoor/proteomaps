import sys
import os
import subprocess
import re

PYTHON_PATH = '../python/'
GENERAL_HIERARCHY_FNAME = '../genomic_data/KO_gene_hierarchy/KO_gene_hierarchy_general.tms'
CHANGES_FNAME = '../genomic_data/KO_gene_hierarchy/KO_gene_hierarchy_changes_v1.0.csv'
MODIFIED_HIERARCHY_FNAME = '../res/hierarchy_modified.tms'
MAPPINGS_PATH = '../genomic_data/KO_gene_hierarchy/KO_gene_hierarchy_organism_mapping/'

p = subprocess.Popen(['python',
                      PYTHON_PATH + 'modify_ko_hierarchy.py',
                      GENERAL_HIERARCHY_FNAME,
                      CHANGES_FNAME,
                      MODIFIED_HIERARCHY_FNAME])
p.wait()

for f in os.listdir(MAPPINGS_PATH):
    for org in re.findall('(\w+)_mapping\.csv', f):
        subprocess.Popen(['python',
                          PYTHON_PATH + 'extend_hierarchy_with_genes.py',
                          MODIFIED_HIERARCHY_FNAME,
                          MAPPINGS_PATH + org + '_mapping.csv',
                          '--root_name', org],
                          stdout=open('../res/' + org + '_hierarchy.tms', 'w'))
    

