import sys
import modify_ko_hierarchy, extend_hierarchy_with_genes
import re
import os

try:
    os.mkdir('../res')
except OSError:
    pass

PYTHON_PATH = './'
GENERAL_HIERARCHY_FNAME = '../../../Dropbox/proteomaps_data/genomic_data/KO_gene_hierarchy/KO_gene_hierarchy_general.tms'
CHANGES_FNAME = '../data/KO_gene_hierarchy/KO_gene_hierarchy_changes.csv'
MAPPINGS_PATH = '../data/KO_gene_hierarchy/KO_gene_hierarchy_organism_mapping/'
MODIFIED_HIERARCHY_FNAME = '../res/hierarchy_modified.tms'

if not os.path.exists(GENERAL_HIERARCHY_FNAME):
    sys.stderr.write('ERROR: cannot find the general hierarchy file')
    
modify_ko_hierarchy.main(open(GENERAL_HIERARCHY_FNAME, 'r'),
                         open(CHANGES_FNAME, 'r'),
                         open(MODIFIED_HIERARCHY_FNAME, 'w'))

for f in os.listdir(MAPPINGS_PATH):
    for org in re.findall('(\w+)_mapping\.csv', f):
        extend_hierarchy_with_genes.extend(hierarchy_file=open(MODIFIED_HIERARCHY_FNAME, 'r'),
                                           mapping_file=open(MAPPINGS_PATH + org + '_mapping.csv', 'r'),
                                           output_file=open('../res/' + org + '_hierarchy.tms', 'w'),
                                           root_name=org)
    

