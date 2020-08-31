import pathlib
from treeflow_pipeline.util import build_shell_command

OUT_PATH = pathlib.Path("out")

rule dengue:
    input:
        "out/dengue/RAxML_info.raxml"

rule ml_topology:
    input:
        "data/{dataset}.fasta"
    output:
        "out/{dataset}/RAxML_info.raxml"
        "out/{dataset}/RAxML_bestTree.raxml"
    shell:
        build_shell_command(
            "raxmlHPC",
            ["HKY"],
            dict(
                p=123, # seed
                n="raxml",
                s="{input}",
                w=(OUT_PATH / "{wildcards.dataset}").resolve(), # output file,
                m="GTRGAMMAX",  # GTR Gamma substitution, estimate frequencies
            )
        )