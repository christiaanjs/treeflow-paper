import pathlib
import treeflow_pipeline.topology_inference as top
from treeflow_pipeline.util import yaml_input, yaml_output, sequence_input, text_input, text_output, pickle_output
from Bio import SeqIO

wd = pathlib.Path(config["working_directory"])

# Compute sequence length for LSD2
_alignment_format = config.get("alignment_format", "fasta")
_first_seq = next(SeqIO.parse(config["alignment"], _alignment_format))
seq_length = len(_first_seq.seq)

unrooted_tree_filename = dict(
    raxml="RAxML_bestTree.raxml"
)[config["tree_method"]]

lsd_output_format = config["lsd_output_format"]

# LSD2 v2.4+ outputs .nwk for distance tree and .date.nexus for date tree
rooted_tree_filename, rooted_tree_format, fix_dates = {
    "lsd": ("lsd-tree.nwk", "newick", False),
    "lsd-dates": ("lsd-dates-tree.date.nexus", "nexus", True)
}[config["rooting_method"]]

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
        tree = wd / unrooted_tree_filename
    output:
        date_tree = wd / "lsd-tree.date.nexus",
        distance_tree = wd / "lsd-tree.nwk",
        log = (wd / "lsd-tree")
    shell:
        f"lsd -r a -s {seq_length} -i {{input.tree}} -o {{output.log}}"

rule root_topology_lsd_dates:
    input:
        dates = wd / "lsd.dates",
        tree = wd / unrooted_tree_filename
    output:
        date_tree = wd / "lsd-dates-tree.date.nexus",
        distance_tree = wd / "lsd-dates-tree.nwk",
        log= (wd / "lsd-dates-tree")
    shell:
        f"lsd -r a -s {seq_length} -i {{input.tree}} -d {{input.dates}} -o {{output.log}}"

rule starting_values_raxml:
    input:
        wd / "RAxML_info.raxml"
    output:
        wd / "starting-values-raxml.yaml"
    run:
        yaml_output(top.get_starting_values_raxml(input[0], config["subst_model"]), output[0])

rule starting_values_lsd:
    input:
        date_tree = rules.root_topology_lsd.output.date_tree,
        distance_tree = rules.root_topology_lsd.output.distance_tree
    output:
        wd / "starting-values-lsd.yaml"
    run:
        yaml_output(top.get_starting_values_lsd(
            input.date_tree,
            input.distance_tree,
            "nexus",
            config["tree_model"],
            distance_tree_format="newick",
        ), output[0])

rule starting_values_lsd_dates:
    input:
        date_tree = rules.root_topology_lsd_dates.output.date_tree,
        distance_tree = rules.root_topology_lsd_dates.output.distance_tree
    output:
        wd / "starting-values-lsd-dates.yaml"
    run:
        yaml_output(top.get_starting_values_lsd(
            input.date_tree,
            input.distance_tree,
            "nexus",
            config["tree_model"],
            distance_tree_format="newick",
        ), output[0])

rule starting_values:
    input:
        tree_starting_values = wd / "starting-values-{0}.yaml".format(config["tree_method"]),
        rooting_starting_values = wd / "starting-values-{0}.yaml".format(config["rooting_method"])
    output:
        wd / "starting-values.yaml"
    run:
        yaml_output(top.finalise_starting_values(
            yaml_input(input.tree_starting_values),
            yaml_input(input.rooting_starting_values),
            config["clock_model"]
        ), output[0])

rule tree:
    input:
        wd / rooted_tree_filename
    output:
        config["output"]
    run:
        top.convert_tree(input[0], rooted_tree_format, output[0], 'newick', strip_data=True, allow_zero_branches=False, fix_dates=fix_dates)
