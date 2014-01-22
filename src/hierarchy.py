from treelib import Node, Tree
import csv
import logging

class Hierarchy(Tree):

    def __init__(self):
        Tree.__init__(self)

    def save(self, f, nid=None, level=Tree.ROOT, idhidden=True, filter=None, cmp=None, key=None, reverse=False):
        leading = ''
        lasting = ''
        nid = self.root if (nid is None) else Node.sanitize_id(nid)
        label = ("{0}".format(self[nid].tag)) if idhidden else ("{0}[{1}]".format(self[nid].tag, self[nid].identifier))
        filter = (self._real_true) if (filter is None) else filter

        if level == self.ROOT:
            f.write(label + '\n')
        else:
            leading += '\t' * level
            f.write("{0}{1}{2}\n".format(leading, lasting, label))

        if filter(self[nid]) and self[nid].expanded:
            queue = [self[i] for i in self[nid].fpointer if filter(self[i])]
            key = (lambda x: x) if (key is None) else key
            queue.sort(cmp=cmp, key=key, reverse=reverse)
            level += 1
            for element in queue:
                self.save(f, element.identifier, level, idhidden, filter, cmp, key, reverse)

    @staticmethod
    def load(f):
        h = Hierarchy()
        csv_reader = csv.reader(f, delimiter='\t')
    
        branch = {}
        for row in csv_reader:
            level = 0
            while row[level] == '': level += 1

            branch[level] = row[level]
            
            tag = branch[level]
            identifier = (level, branch[level])
            
            if h.get_node(identifier) is not None:
                logging.debug('The label [%s] appears twice in the hierarchy' 
                              % branch[level])
            elif level == 0:
                h.create_node(tag, identifier, parent=None)
            else:
                h.create_node(tag, identifier, parent=(level-1, branch[level-1]))
        return h

    def all_leaves(self):
        """
            get all nodes that are leaves in the tree (i.e. have no children)
        """
        for node in self.all_nodes():
            if len(node.fpointer) == 0:
                yield node
        
    def apply_updates(self, update_file, organism_name=None, KO_level=4):
        # read the original hierarchy file into a tree structure
        
        self.create_node('Not mapped', (1, 'Not mapped'), parent=self.root)
        for level in xrange(2, 5):
            self.create_node('Not mapped', (level, 'Not mapped'), parent=(level-1, 'Not mapped'))
        self.create_node('Not mapped:NotMapped', (5, 'Not mapped'), parent=(4, 'Not mapped'))

        # -----------------------------------------------------------
        ## read file KO_gene_hierarchy_changes.csv containing genes to be added to the hierarchy
        ## and move the KO to the new pathway.
        ## Note: the last mapping in the modification file is the one that will be used

        for row in csv.reader(update_file, delimiter='\t'):
            if row[0] != organism_name:
                continue
            KO, pathway = row[3:5]

            pathway_node = self.get_node((KO_level-1, pathway))
            if pathway_node is None:
                logging.debug('pathway %s is not in the general hierarchy file' % pathway)
                continue
            
            KO_node = self.get_node((KO_level, KO))
            if KO_node is None:
                logging.debug('KO %s is not in the general hierarchy file, creating new one' % KO)
                self.create_node(KO, (KO_level, KO), parent=pathway_node.identifier)
            else:
                self.move_node(KO_node.identifier, pathway_node.identifier)
                
    def extend(self, mapping_file, organism_name=None, KO_level=4):
        self.get_node((0, 'KO')).tag = organism_name
        used_systematic_names = set()    
        for row in csv.reader(mapping_file, delimiter='\t'):
            systematic, gene, KO = row
            if systematic in used_systematic_names:
                logging.debug('The gene %s appears more than once, ignoring all '
                              'duplicate occurances' % gene)
                continue
            used_systematic_names.add(systematic)
            
            KO_node = self.get_node((KO_level, KO))
            if KO_node is None:
                logging.debug('The KO %s is not in the general hierarchy file' % KO)
                continue

            leaf_tag = '%s:%s' % (gene, systematic)
            self.create_node(leaf_tag, leaf_tag, parent=KO_node.identifier)
