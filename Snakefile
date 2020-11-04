include: "workflow/SimSnakefile"

import pathlib
import treeflow_pipeline.topology_inference as top
import treeflow_pipeline.templating as tem
import treeflow_pipeline.model as mod
import treeflow_pipeline.results as res
from treeflow_pipeline.util import build_shell_command, yaml_input, yaml_output, sequence_input, text_input, text_output, pickle_output
 
OUT_PATH = pathlib.Path("out")
sequence_lengths =  [100, 1000, 5000, 20000] 
rule sim:
    input:
        ["out/sim/sequence_length{sequence_length}/beast-relaxed-fixed.pickle".format(sequence_length=sequence_length) for sequence_length in sequence_lengths],
        ["out/sim/sequence_length{sequence_length}/variational-relaxed-{approx}.pickle".format(sequence_length=sequence_length, approx=approx) for sequence_length in sequence_lengths for approx in ["mean_field", "scaled", "tuneable"]]

rule ml_topology:
    input:
        "data/{dataset}.fasta"
    output:
        "out/{dataset}/RAxML_info.raxml",
        "out/{dataset}/RAxML_bestTree.raxml"
    shell:
        build_shell_command(
            "raxmlHPC",
            ["--HKY85", "-V"],
            dict(
                p=123, # seed
                n="raxml", # ID
                s="{input}", # input file
                w=(OUT_PATH / "{wildcards.dataset}").resolve(), # output file
                m="GTRCATX",  # GTR substitution, estimate frequencies
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
        "lsd -c -r a -i {input.tree} -d {input.dates} -o {output.log}"

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

rule beast_xml: # TODO: Should template be input?
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

rule variational_fit:
    input:
        fasta = "data/{dataset}.fasta",
        tree = "out/{dataset}/lsd-tree.date.newick",
        starting_values = "out/{dataset}/starting-values.yaml",
        prior_params = "config/prior-params.yaml",
        vi_config = "config/vi-config.yaml"
    output:
        "out/{dataset}/variational-{clock}-{approx}.pickle"
    run:
        pickle_output(mod.get_variational_fit(
            input.tree,
            input.fasta,
            yaml_input(input.starting_values),
            yaml_input(input.prior_params),
            yaml_input(input.vi_config),
            wildcards.clock,
            wildcards.approx
        ), output[0])
        
rule beast_results:
    input:
        topology = "out/{dataset}/lsd-tree.date.newick",
        trees = "out/{dataset}/beast-{clock}-fixed.trees",
        trace = "out/{dataset}/beast-{clock}-fixed.log",
        beast_config = "config/beast-config.yaml"
    output:
        "out/{dataset}/beast-{clock}-fixed.pickle"
    run:
        pickle_output(
            res.process_beast_results(
                input.trees,
                input.trace,
                input.topology,
                yaml_input(input.beast_config),
                wildcards.clock
            ),
            output[0]
        )
# TODO: Split out tree parsing
# TODO: Same trees for different sequence length

rule relaxed_plot:
    input:
        topology = "out/{dataset}/lsd-tree.date.newick",
        trees = "out/{dataset}/beast-relaxed-fixed.trees",
        trace = "out/{dataset}/beast-relaxed-fixed.log",
        beast_config = "config/beast-config.yaml",
        prior_params = "config/prior-params.yaml",
        variational_fit = "out/{dataset}/variational-relaxed-{approx}.pickle",
        notebook = "notebook/plot-posterior-relaxed.ipynb"
    output:
        notebook = "out/{dataset}/plot-posterior-relaxed-{approx}.ipynb",
        plot = "out/{dataset}/rate-correlations-{approx}.png"
    shell:
        """
        papermill {input.notebook} {output.notebook} \
            -p trace_file {input.trace} \
            -p tree_file {input.trees} \
            -p topology_file {input.topology} \
            -p beast_config_file {input.beast_config} \
            -p prior_params_file {input.prior_params} \
            -p variational_fit_file {input.variational_fit} \
            -p clock_approx {wildcards.approx} \
            -p plot_out_file {output.plot}
        """

rule relaxed_report:
    input:
        "out/{dataset}/plot-posterior-relaxed-{approx}.ipynb"
    output:
        "out/{dataset}/plot-posterior-relaxed-{approx}.html"
    shell:
        "jupyter nbconvert --to html {input}"