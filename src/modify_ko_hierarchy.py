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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='a simplified version of filter_ko_hierarchy.py')
    parser.add_argument('hierarchy_file', type=file, help='path to hierarchy file')
    parser.add_argument('modification_file', type=file, help='path to modification file')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='path to modification file')
    args = parser.parse_args()
    KO_tree = hierarchy.Hierarchy.load(args.hierarchy_file)
    KO_tree.apply_updates(args.modification_file)
    KO_tree.save(args.output_file)
    output_file.close()

