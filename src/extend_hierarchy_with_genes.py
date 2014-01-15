import logging
import sys
import csv
import argparse
import hierarchy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extend a tree using node mapping.')
    parser.add_argument('hierarchy_file', type=file, help='path to hierarchy file')
    parser.add_argument('mapping_file', type=file, help='path to gene mapping file')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='path to modification file')
    parser.add_argument('--organism', metavar='o', default=None, help='name of organism, to be used as the root node')
    
    args = parser.parse_args()
    KO_tree = hierarchy.Hierarchy.load(args.hierarchy_file)
    KO_tree.extend(args.mapping_file, organism_name=args.organism_name)
    KO_tree.save(args.output_file)
    output_file.close()
