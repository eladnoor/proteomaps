import logging
import sys
import csv
import argparse
import hierarchy

def syntax():
    print "input error"
    sys.exit(-1)

def extend(hierarchy_file, mapping_file, output_file, organism_name=None):
    
    ko_tree = hierarchy.Hierarchy.load(hierarchy_file)
    ko_tree.get_node((0, 'KO')).tag = organism_name
    mapping = {'Not mapped': 'Not mapped:NotMapped'}
    
    for row in csv.reader(mapping_file, delimiter='\t'):
        systematic, gene, KO = row
        if KO in mapping:
            # this is only a warning but the mapping is still changed to the
            # new gene. therefore, the mapping will always be to the last
            # entry.
            logging.debug('KO %s is mapped to two different genes' % KO)
        mapping[KO] = '%s:%s' % (gene, systematic)

    # get all nodes that are leaves in the tree (i.e. genes)
    for child in ko_tree.all_leaves():
        if child.tag in mapping.keys():
            child.tag = mapping[child.tag]
        else:
            ko_tree.remove_node(child.identifier)

    ko_tree.save(output_file)
    output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extend a tree using node mapping.')
    parser.add_argument('hierarchy_file', type=file, help='path to hierarchy file')
    parser.add_argument('mapping_file', type=file, help='path to gene mapping file')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='path to modification file')
    parser.add_argument('--organism', metavar='o', default=None, help='name of organism, to be used as the root node')
    
    args = parser.parse_args()
    extend(args.hierarchy_file, args.mapping_file, args.output_file, args.organism_name)
