import dendropy
import numpy as np

def construct_precedes_map(taxon_order): # [x, y] = True if x precedes y TODO: Do something that isn't quadratic
    precedes = {}
    for i in range(len(taxon_order)):
        for j in range(i):
            precedes[taxon_order[i], taxon_order[j]] = False
        for j in range(i + 1, len(taxon_order)):
            precedes[taxon_order[i], taxon_order[j]] = True
    return precedes

def construct_leaf_descendant_map(tree): # [x] = y only if y is a leaf descendant of x
    leaf_descendants = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            leaf_descendants[node] = node
        else:
            leaf_descendants[node] = leaf_descendants[node.child_nodes()[0]]
    return leaf_descendants

def postorder_node_iter(tree, taxon_order):
    if taxon_order is None:
        return tree.postorder_node_iter()
    else:
        precedes = construct_precedes_map(taxon_order)
        leaf_desc_map = construct_leaf_descendant_map(tree)
    stack = [tree.seed_node]
    reverse = []
    while len(stack) > 0:
        reverse.append(stack.pop())
        if not reverse[-1].is_leaf():
            left, right = reverse[-1].child_nodes()
            if precedes[leaf_desc_map[left].taxon.label, leaf_desc_map[right].taxon.label]:
                stack.append(left)
                stack.append(right)
            else:
                stack.append(right)
                stack.append(left)
    
    return reverse[::-1]

def parse_beast_trees(tree_file, format="nexus", metadata_keys=[], taxon_order=None):
    trees = dendropy.TreeList.get(path=tree_file, schema=format, rooting="default-rooted", preserve_underscores=True)
    taxon_count = len(trees.taxon_namespace)
    branch_lengths = np.zeros((len(trees), 2 * taxon_count - 2))
    metadata = { key: np.zeros((len(trees), 2 * taxon_count - 2)) for key in metadata_keys }

    label_length = max([len(taxon.label) for taxon in trees.taxon_namespace])
    label_dtype = ">U{0}".format(label_length)
    taxon_names = np.empty((len(trees), taxon_count), dtype=label_dtype)

    for i, tree in enumerate(trees):
        leaf_index = 0
        node_index = taxon_count
        for node in list(postorder_node_iter(tree, taxon_order))[:-1]:
            if node.is_leaf():
                branch_lengths[i, leaf_index] = node.edge.length
                taxon_names[i, leaf_index] = node.taxon.label
                for key in metadata_keys:
                    metadata[key][i, leaf_index] = float(node.annotations[key].value)
                leaf_index += 1
            else:
                branch_lengths[i, node_index] = node.edge.length
                for key in metadata_keys:
                    metadata[key][i, node_index] = float(node.annotations[key].value)
                node_index += 1

    return dict(branch_lengths=branch_lengths, taxon_names=taxon_names, metadata=metadata)