import pytest
import pathlib


@pytest.fixture
def test_newick_file():
    return pathlib.Path("test") / "data" / "test-tree.newick"
