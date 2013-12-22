from treelib import Node, Tree, tree
import csv
import logging
import sys

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
        
