import logging
import sys
import csv
import argparse

def syntax():
    print "input error"
    sys.exit(-1)

def extend(hierarchy_fname, mapping_fname, root_name):
    
    mapping = {'Not mapped': 'Not Mapped:NotMapped'}
    for row in csv.reader(open(mapping_fname, 'r'), delimiter='\t'):
        systematic, gene, KO = row
        if KO in mapping:
            logging.debug('KO %s is mapped to two different genes' % KO)
        mapping[KO] = ':'.join([gene, systematic])

    output = csv.writer(sys.stdout, delimiter='\t')
    if root_name is not None:
        output.writerow([root_name])
        KO_level = 4
    else:
        KO_level = 5

    for row in csv.reader(open(hierarchy_fname, 'r'), delimiter='\t'):
        if len(row) == KO_level:
            KO = row[-1]
            if KO not in mapping:
                logging.debug('KO %s appears in the hierarchy but is not mapped' % KO)
                continue
            row = row[:-1] + [mapping[KO]]
        if root_name is not None:
            row = [''] + row
        output.writerow(row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extend a tree using node mapping.')
    parser.add_argument('hierarchy_fname', help='path to hierarchy file')
    parser.add_argument('mapping_fname', help='path to gene mapping file')
    parser.add_argument('--root_name', metavar='r', default=None, help='name of root node to be added')
    
    args = parser.parse_args()
    extend(args.hierarchy_fname, args.mapping_fname, args.root_name)
