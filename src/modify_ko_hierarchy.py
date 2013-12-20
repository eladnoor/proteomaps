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

parser = argparse.ArgumentParser(description='a simplified version of filter_ko_hierarchy.py')
parser.add_argument('hierarchy_fname', help='path to hierarchy file')
parser.add_argument('modification_fname', help='path to modification file')
parser.add_argument('output_fname', help='path to modification file')
args = parser.parse_args()

# -----------------------------------------------------------
## read file modify_KO_hierarchy.tsv containing genes to be added to the hierarchy
## (overriding possibly existing entries)

added_ko_2_pathway = set()
modified_ko_dictionary = {} 
for row in csv.reader(open(args.modification_fname, 'r'), delimiter='\t'):
    ko, pathway = row[3:5]

    # if there are several pathways for the same KO, we use only the last one
    added_ko_2_pathway.add(ko)
    modified_ko_dictionary.setdefault(pathway, set()).add(ko)

# -----------------------------------------------------------

## go through hierarchy and write all lines that are either no KO numbers or 
## KO numbers that appear in the relevant list; also exclude all KO Numbers that
## already appeared further up

ko_used = set()
outlines = []
for row in csv.reader(open(args.hierarchy_fname, 'r'), delimiter='\t'):
    # in a hierarchy file, only one of the cells in each row is not empty
    # we store the index of that cell in 'col'
    level = 0
    while not row[level]:
        level += 1
    
    if level == 3:
        # if the KO is in one of the changed pathways, skip this line since
        # we are moving it to a different location in the tree
        if row[level] in added_ko_2_pathway:
           continue

        # if the KO has already been used somewhere in the tree before this line
        # remove the line as well, since we cannot have duplicate entries
        if row[level] in ko_used:
            continue
        ko_used.add(row[level])
    
    outlines.append(row)
    
    # if this is a level-2 line, and the label matches one of the
    # dictionary entries, drop all the relevant KOs in this location
    if level == 2:
        for ko in modified_ko_dictionary.get(row[level], []):
            outlines.append(['', '', '', ko])
            ko_used.add(ko)

for i in xrange(4):
    outlines.append([''] * i + ['Not mapped'])
    
fo = csv.writer(open(args.output_fname, 'w'), delimiter='\t')
fo.writerows(outlines)


