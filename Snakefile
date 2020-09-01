import pathlib
import treeflow_pipeline.topology_inference as top
import treeflow_pipeline.templating as tem
from treeflow_pipeline.util import build_shell_command, yaml_input, yaml_output, sequence_input, text_input, text_output
 
OUT_PATH = pathlib.Path("out")

rule dengue:
    input:
        "out/dengue/beast-strict-estimate.log"

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
        top.build_lsd_date_file(sequence_input(input[0]), output[0])

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

rule starting_values:
    input:
        date_tree = "out/{dataset}/lsd-tree.date.nexus",
        distance_tree = "out/{dataset}/lsd-tree.nexus",
        ml_info = "out/{dataset}/RAxML_info.raxml"
    output:
        "out/{dataset}/starting-values.yaml"
    run:
        yaml_output(top.get_starting_values(
            input.date_tree,
            input.distance_tree,
            input.ml_info
        ), output[0])

rule prepare_tree:
    input:
        "out/{dataset}/lsd-tree.date.nexus"
    output:
        "out/{dataset}/lsd-tree.date.newick"
    run:
        top.convert_tree(input[0], 'nexus', 'newick', strip_data=True, allow_zero_branches=False)

rule beast_xml:
    input:
        fasta = "data/{dataset}.fasta",
        tree = "out/{dataset}/lsd-tree.date.newick",
        starting_values = "out/{dataset}/starting-values.yaml",
        prior_params = "config/prior-params.yaml",
        beast_config = "config/beast-config.yaml"
    output:
        "out/{dataset}/beast-{clock}-{tree_type}.xml"
    run:
        text_output(tem.build_beast_analysis(
            wildcards.clock,
            wildcards.tree_type,
            sequence_input(input.fasta),
            text_input(input.tree),
            yaml_input(input.starting_values),
            yaml_input(input.prior_params),
            yaml_input(input.beast_config),
            output[0]
        ), output[0])

rule beast_run:
    input:
        "out/{dataset}/beast-{model}.xml"
    output:
        "out/{dataset}/beast-{model}.trees",
        "out/{dataset}/beast-{model}.log"
    shell:
        "beast {input}" 
