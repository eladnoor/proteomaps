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
import sys

if True:
    GO_OBO_URL = 'http://geneontology.org/ontology/go.obo'
    GO_OBO_FNAME = 'res/go.obo'
    GO_TMS_FNAME = 'res/go.tms'
    
    GENE_ASSOCIATION_URL = 'http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.sgd.gz?rev=HEAD'
    GENE_ASSOCIATION_FNAME = 'res/gene_association.sgd.gz'
    
    GO_ROOT = 'GO:0008150' # biological_process
    #GO_ROOT = 'GO:0003674' # molecular_function
else:
    GO_OBO_URL = 'http://www.geneontology.org/GO_slims/goslim_yeast.obo'
    GO_OBO_FNAME = 'res/goslim_yeast.obo'
    GO_TMS_FNAME = 'res/goslim_yeast.tms'

    GENE_ASSOCIATION_URL = 'http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.sgd.gz?rev=HEAD'
    GENE_ASSOCIATION_FNAME = '../res/gene_association.sgd.gz'

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

logging.info('Converting OBO directed acyclic graph to a tree')
GOtree = dag.to_tree(GO_ROOT)

logging.info('Reading gene association mapping from "%s"'
             % GENE_ASSOCIATION_FNAME)

# save the set of all GOs that every gene can be mapped to
gene_to_go = {}
for row in csv.reader(gzip.open(GENE_ASSOCIATION_FNAME, 'r'), delimiter='\t'):
    if len(row) <= 10 or row[0][0] == '!':
        continue
    
    go = row[4]
    gene = row[10].split('|')[0]
    node = GOtree.get_node(go)
    
    # check that this GO is a leaf in the tree (i.e. has no children)
    if node is not None and node.is_leaf():
        gene_to_go.setdefault(gene, set([])).add(go)

# now, we need to select only one GO for each gene, and add that gene to
# the hierarchy
logging.info('Mapping specific genes into the GO hierarchy as leaves')
for gene, go_list in gene_to_go.iteritems():
    go = list(go_list)[0]
    GOtree.create_node(gene, gene, parent=go)

logging.info('Writing mapped GO hierarchy to "%s"'
             % GO_TMS_FNAME)
GOtree.save(open(GO_TMS_FNAME, 'w'))
