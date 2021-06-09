import pytest
import pathlib


@pytest.fixture
def test_data_dir():
    return pathlib.Path("test") / "data"


@pytest.fixture
def test_newick_file(test_data_dir):
    return test_data_dir / "test-tree.newick"


@pytest.fixture
def test_fasta_file(test_data_dir):
    return test_data_dir / "test-sequences.fasta"
