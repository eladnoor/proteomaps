#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 10:40:32 2014

@author: eladnoor
"""

import urllib
from src.obo_parser import GODag
import logging
import os
import gzip
import csv

if not os.path.exists('res'):
    os.mkdir('res')

# file locations and URLs
GO_OBO_URL = 'http://geneontology.org/ontology/go.obo'
GO_OBO_FNAME = 'res/go.obo'

GENE_ASSOCIATION_URL = 'http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.sgd.gz?rev=HEAD'
GENE_ASSOCIATION_FNAME = 'res/gene_association.sgd.gz'

GO_TMS_FNAME = 'res/go.tms'
GO_MAPPED_TMS_FNAME = 'res/go_sce.tms'

GO_ROOT = 'GO:0008150' # biological_process
#GO_ROOT = 'GO:0003674' # molecular_function

GO_PRIORIY_FILE = 'data/GO/priority.csv'

###############################################################################

if not os.path.exists(GO_OBO_FNAME):
    logging.info('Reading OBO file from "%s"' % GO_OBO_URL)
    f = open(GO_OBO_FNAME, 'w')
    f.write(urllib.urlopen(GO_OBO_URL).read())
    f.close()

if not os.path.exists(GENE_ASSOCIATION_FNAME):
    logging.info('Reading Gene Association file from "%s"' 
                 % GENE_ASSOCIATION_URL)
    f = open(GENE_ASSOCIATION_FNAME, 'w')
    f.write(urllib.urlopen(GENE_ASSOCIATION_URL).read())
    f.close()

logging.info('Parsing OBO file "%s"' % GO_OBO_FNAME)
dag = GODag(GO_OBO_FNAME)

logging.info('Loading GO priorities from "%s"' % GO_PRIORIY_FILE)
priority_list = [row[1] for row in csv.reader(open(GO_PRIORIY_FILE, 'r'), delimiter='\t')]

logging.info('Converting OBO directed acyclic graph to a tree')
GOtree = dag.to_tree(GO_ROOT, priority_list=priority_list)

logging.info('Writing GO hierarchical tree to "%s"' % GO_TMS_FNAME)
GOtree.save(open(GO_TMS_FNAME, 'w'))

logging.info('Reading gene association mapping from "%s"'
             % GENE_ASSOCIATION_FNAME)

# save the set of all GOs that every gene can be mapped to
gene_to_go = {}
gene_to_display_name = {}
for i, row in enumerate(csv.reader(gzip.open(GENE_ASSOCIATION_FNAME, 'r'),
                                   delimiter='\t')):
    if len(row) <= 10 or row[0][0] == '!':
        continue
    
    display_name = row[2]
    go = row[4]
    systematic_name = row[10].split('|')[0]
    
    if not systematic_name:
        logging.warning('There is no systematic name in line #%d' % i)
        continue
    if not systematic_name:
        short_name = systematic_name
    node = GOtree.get_node(go)
    
    # check that this GO is a leaf in the tree (i.e. has no children)
    if node is not None and len(node.fpointer) == 0:
        gene_to_go.setdefault(systematic_name, []).append(go)
        if systematic_name not in gene_to_display_name:
            gene_to_display_name[systematic_name] = display_name

# go through all the genes in the dictionary and order the list of GO IDs
# according to the priority_list (if relevant)
for systematic_name in gene_to_go.keys():
    go_ids = set(gene_to_go[systematic_name])
    
    gene_to_go[systematic_name] = []
    for go_id in priority_list:
        if go_id in go_ids:
            gene_to_go[systematic_name].append(go_id)
            go_ids.remove(go_id)
    gene_to_go[systematic_name] += list(go_ids)
    
# Now, we need to select only one GO for each gene, and add that gene to
# the hierarchy. Currently, we simply choose the first GO that appears in
# the association file.
logging.info('Mapping specific genes into the GO hierarchy as leaves')
for systematic_name, go_list in gene_to_go.iteritems():
    display_name = gene_to_display_name[systematic_name]
    gene = display_name + ":" + systematic_name
    go = go_list[0]
    GOtree.create_node(gene, gene, parent=go)

logging.info('Writing mapped GO hierarchy to "%s"' % GO_MAPPED_TMS_FNAME)
GOtree.save(open(GO_MAPPED_TMS_FNAME, 'w'))
