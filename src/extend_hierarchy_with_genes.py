import logging
import sys
import csv
import argparse

def syntax():
    print "input error"
    sys.exit(-1)

def extend(hierarchy_file, mapping_file, output_file, root_name=None):
    
    mapping = {'Not mapped': 'Not Mapped:NotMapped'}
    for row in csv.reader(mapping_file, delimiter='\t'):
        systematic, gene, KO = row
        if KO in mapping:
            logging.debug('KO %s is mapped to two different genes' % KO)
        mapping[KO] = ':'.join([gene, systematic])

    output = csv.writer(output_file, delimiter='\t')
    if root_name is not None:
        output.writerow([root_name])
        KO_level = 4
    else:
        KO_level = 5

    for row in csv.reader(hierarchy_file, delimiter='\t'):
        if len(row) == KO_level:
            KO = row[-1]
            if KO not in mapping:
                logging.debug('KO %s appears in the hierarchy but is not mapped' % KO)
                continue
            row = row[:-1] + [mapping[KO]]
        if root_name is not None:
            row = [''] + row
        output.writerow(row)
    output_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extend a tree using node mapping.')
    parser.add_argument('hierarchy_file', type=file, help='path to hierarchy file')
    parser.add_argument('mapping_file', type=file, help='path to gene mapping file')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='path to modification file')
    parser.add_argument('--root_name', metavar='r', default=None, help='name of root node to be added')
    
    args = parser.parse_args()
    extend(args.hierarchy_file, args.mapping_file, args.output_file, args.root_name)
