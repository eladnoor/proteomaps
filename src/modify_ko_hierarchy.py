# ------------------------------------------------------
# This is a simplified version of filter_ko_hierarchy.py
# ------------------------------------------------------

# USAGE modify_ko_hierarchy_keep_all_occurences.py [ORIGINAL_HIERARCHY.tms] [REANNOTATION_FILE.csv] 

# 1. read KO hierarchy from tms file
# 2. add all entries requested in reannotation
# 3. add categories for non-mapped genes
# 4. Write new new hierachy file

import os
import glob
import re
import sys
import argparse
import csv
import logging
import hierarchy

def main(hierarchy_file, modification_file, output_file, tree_depth=4):
    # read the original hierarchy file into a tree structure
    
    ko_tree = hierarchy.Hierarchy.load(hierarchy_file)
    ko_tree.create_node('Not mapped', (1, 'Not mapped'), parent=ko_tree.root)
    for level in xrange(2, 5):
        ko_tree.create_node('Not mapped', (level, 'Not mapped'), parent=(level-1, 'Not mapped'))

    # -----------------------------------------------------------
    ## read file KO_gene_hierarchy_changes.csv containing genes to be added to the hierarchy
    ## and move the KO to the new pathway.
    ## Note: the last mapping in the modification file is the one that will be used

    ko_to_pathway_dict = {}
    for row in csv.reader(modification_file, delimiter='\t'):
        ko, pathway = row[3:5]
        try:
            ko_tree.move_node(ko, pathway)
        except KeyError:
            logging.debug('The KO [%s] is not in the original hierarchy' % ko)

    ko_tree.save(output_file)
    output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='a simplified version of filter_ko_hierarchy.py')
    parser.add_argument('hierarchy_file', type=file, help='path to hierarchy file')
    parser.add_argument('modification_file', type=file, help='path to modification file')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='path to modification file')
    args = parser.parse_args()
    main(args.hierarchy_file, args.modification_file, args.output_file)

