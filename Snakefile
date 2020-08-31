import pathlib
import treeflow_pipeline.topology_inference as top
from treeflow_pipeline.util import build_shell_command

OUT_PATH = pathlib.Path("out")

rule dengue:
    input:
        "out/dengue/lsd-tree.nexus"

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
                n="raxml", # ID
                s="{input}", # input file
                w=(OUT_PATH / "{wildcards.dataset}").resolve(), # output file
                m="GTRGAMMAX",  # GTR Gamma substitution, estimate frequencies
            )
        )

rule lsd_dates:
    input:
        "data/{dataset}.fasta"
    output:
        "out/{dataset}/lsd.dates"
    run:
        top.build_lsd_date_file(input[0], 'fasta', output[0])

rule root_topology:
    input:
        dates = "out/{dataset}/lsd.dates",
        tree = "out/{dataset}/RAxML_bestTree.raxml"
    output:
        "out/{dataset}/lsd-tree.date.nexus",
        "out/{dataset}/lsd-tree.nexus",
        log="out/{dataset}/lsd-tree"
    shell:
        build_shell_command(
            "lsd",
            ["c"],
            dict(
                r="a",
                i="{input.tree}",
                d="{input.dates}",
                o="{output.log}"
            ),
            arg_prefix="-"
        )
