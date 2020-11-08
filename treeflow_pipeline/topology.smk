import pathlib
import treeflow_pipeline.topology_inference as top
from treeflow_pipeline.util import yaml_input, yaml_output, sequence_input, text_input, text_output, pickle_output

wd = pathlib.Path(config["working_directory"])

unrooted_tree_filename = dict(
    raxml="RAxML_bestTree.raxml"
)[config["tree_method"]]

rooted_tree_filename, rooted_tree_format = dict(
    lsd=("lsd-tree.date.nexus", "nexus")
)[config["rooting_method"]]

rule raxml_topology:
    input:
        config["alignment"]
    output:
        wd / "RAxML_info.raxml",
        wd / "RAxML_bestTree.raxml"
    shell:
        top.build_raxml_command("{input}", wd.resolve(), config["seed"], config["site_model"], config["subst_model"])

rule lsd_dates:
    input:
        config["alignment"]
    output:
        wd / "lsd.dates"
    run:
        top.build_lsd_date_file(sequence_input(input[0]), output[0])

rule root_topology_lsd:
    input:
        dates = wd / "lsd.dates",
        tree = wd / unrooted_tree_filename
    output:
        wd / "lsd-tree.date.nexus",
        wd / "lsd-tree.nexus",
        log= (wd / "lsd-tree")
    shell:
        "lsd -c -r a -i {input.tree} -d {input.dates} -o {output.log}"

rule starting_values_raxml:
    input:
        wd / "RAxML_info.raxml"
    output:
        wd / "starting-values-raxml.yaml"
    run:
        yaml_output(top.get_starting_values_raxml(input[0], config["subst_model"]), output[0])

rule starting_values_lsd:
    input:
        date_tree = wd / "lsd-tree.date.nexus",
        distance_tree = wd / "lsd-tree.nexus"
    output:
        wd / "starting-values-lsd.yaml"
    run:
        yaml_output(top.get_starting_values_lsd(
            input.date_tree,
            input.distance_tree
        ), output[0])

rule starting_values:
    input:
        tree_starting_values = wd / "starting-values-{0}.yaml".format(config["tree_method"]),
        rooting_starting_values = wd / "starting-values-{0}.yaml".format(config["rooting_method"])
    output:
        wd / "starting-values.yaml"
    run:
        yaml_output({
            **yaml_input(input.tree_starting_values),
            **yaml_input(input.rooting_starting_values)
        }, output[0])

rule tree:
    input:
        wd / rooted_tree_filename
    output:
        config["output"]
    run:
        top.convert_tree(input[0], rooted_tree_format, output[0], 'newick', strip_data=True, allow_zero_branches=False)
