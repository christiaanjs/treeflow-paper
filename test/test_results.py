import treeflow_pipeline.results
import treeflow.tree_processing
import treeflow.sequences
from numpy.testing import assert_allclose


reference_file = "test/data/test-tree.newick"
reordered_file = "test/data/test-tree-shuffled.newick"

def test_beast_branch_length_parsing_same_file():
    reference_tree, taxon_names = treeflow.tree_processing.parse_newick(reference_file)
    parsed_tree = treeflow_pipeline.results.parse_beast_trees(reference_file, format="newick")
    reference_blens = treeflow.sequences.get_branch_lengths(reference_tree)
    assert_allclose(reference_blens.numpy(), parsed_tree['branch_lengths'][0])


def test_beast_parsing_different_order():
    reference_tree, taxon_names = treeflow.tree_processing.parse_newick(reference_file)
    parsed_tree = treeflow_pipeline.results.parse_beast_trees(reordered_file, taxon_order=taxon_names, format="newick")
    reference_blens = treeflow.sequences.get_branch_lengths(reference_tree)
    assert_allclose(reference_blens.numpy(), parsed_tree['branch_lengths'][0])