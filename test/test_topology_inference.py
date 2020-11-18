import io
import pytest
import Bio.Phylo
import treeflow_pipeline.topology_inference
import numpy as np
from numpy.testing import assert_allclose

test_newick = "(((A_2.5:1.0,B_2.0:0.4999):0.0,C_2.25:0.75):1.5,D_3.0:3.0)"
@pytest.fixture
def test_tree():
    with io.StringIO(test_newick) as f:
        tree = next(Bio.Phylo.parse(f, format="newick"))
    return tree

def leaf_heights(tree):
    distances = np.array([tree.distance(x) for x in tree.get_terminals()])
    return max(distances) - distances

def leaf_height_dict(tree):
    heights = treeflow_pipeline.topology_inference.get_node_heights(tree)
    return { clade.name: heights[clade] for clade in tree.get_terminals() }

def test_adjust_zero_branches(test_tree):
    init_heights = leaf_heights(test_tree)
    treeflow_pipeline.topology_inference.adjust_zero_branches(test_tree.clade)
    final_heights = leaf_heights(test_tree)
    assert_allclose(final_heights, init_heights)
    for node in test_tree.find_clades():
        if node is not test_tree.clade:
            assert node.branch_length > 0.0

def test_fix_leaf_dates(test_tree):
    name_height_dict = treeflow_pipeline.topology_inference.parse_leaf_heights([leaf.name for leaf in test_tree.get_terminals()])
    treeflow_pipeline.topology_inference.fix_leaf_dates(test_tree)
    tree_height_dict = leaf_height_dict(test_tree)
    for name in name_height_dict:
        assert_allclose(tree_height_dict[name], name_height_dict[name])
    
    