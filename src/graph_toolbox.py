# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 22:52:31 2014

@author: eladn
"""

from treelib import Node, Tree

def dag_to_tree(dag, root):
    """
        Use a BFS from the root to decide the level of each node in the graph.
        For each child, keep only one of the edges leading to it so it will
        have a single parent, from the previous level. Choose these edges 
        in a way which will keep the tree balanced (probably not optimally 
        though). This will assure that only one of the shortest paths from 
        the root to every node is kept.
    """
    
    tree = Tree()
    tree.create_node(root, root, parent=None)
    
    parents = set([root])
    
    while parents:
        children = set()
        
        child_to_parent_list = {}
        for parent in parents:
            for child in dag[parent]:
                if tree.get_node(child) is not None:
                    continue
                children.add(child)
                child_to_parent_list.setdefault(child, []).append(parent)

        # the following list will contain pairs, counting the number of children
        # that each parent has.
        child_count_parent_pairs_list = [[0, p] for p in parents]
        for child, parent_list in child_to_parent_list.iteritems():
            # give this child to the parent with the least children
            for i, (counter, parent) in enumerate(child_count_parent_pairs_list):
                if parent in parent_list:
                    break
            tree.create_node(child, child, parent)
            child_count_parent_pairs_list[i][0] += 1
            child_count_parent_pairs_list.sort()
        parents = children
        
    tree.show()
    
if __name__ == "__main__":
    dag = {'a' : ['b', 'c'],
           'b' : ['d', 'f'],
           'c' : ['d', 'e', 'f', 'g'],
           'd' : ['e'],
           'e' : [],
           'f' : [],
           'g' : []}

    dag_to_tree(dag, 'a')