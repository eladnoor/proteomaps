import sys
import re
import os
import hierarchy

try:
    os.mkdir('../res')
except OSError:
    pass

PYTHON_PATH = './'
GENERAL_HIERARCHY_FNAME = '../data/KO_gene_hierarchy/KO_gene_hierarchy_general.tms'
CHANGES_FNAME = '../data/KO_gene_hierarchy/KO_gene_hierarchy_changes.csv'
MAPPINGS_PATH = '../data/KO_gene_hierarchy/KO_gene_hierarchy_organism_mapping/'
MODIFIED_HIERARCHY_FNAME = '../res/hierarchy_modified.tms'

if not os.path.exists(GENERAL_HIERARCHY_FNAME):
    sys.stderr.write('ERROR: cannot find the general hierarchy file')
    
for f in os.listdir(MAPPINGS_PATH):

    try:
        org = re.findall('(\w+)_mapping\.csv', f)[0]
    except IndexError:
        continue
    
    if org != 'hsa':
        continue

    sys.stderr.write('Mapping orgnaism: %s...\n' % org)
    KO_tree = hierarchy.Hierarchy.load(open(GENERAL_HIERARCHY_FNAME, 'r'))
    KO_tree.apply_updates(open(CHANGES_FNAME, 'r'), org)
    KO_tree.extend(open(MAPPINGS_PATH + org + '_mapping.csv', 'r'), org)
    KO_tree.save(open('../res/' + org + '_hierarchy.tms', 'w'))


