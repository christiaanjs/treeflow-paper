include: "workflow/sim.smk"

import pathlib
import treeflow_pipeline.templating as tem
import treeflow_pipeline.model as mod
import treeflow_pipeline.results as res
from treeflow_pipeline.util import build_shell_command, yaml_input, yaml_output, sequence_input, text_input, text_output, pickle_output
 
OUT_PATH = pathlib.Path("out")
sequence_lengths =  [100, 1000, 5000, 20000] 
rule sim:
    input:
        
        #["out/sim/sequence_length{sequence_length}/plot-posterior-relaxed-{approx}.html".format(sequence_length=sequence_length, approx=approx) for sequence_length in sequence_lengths for approx in ["mean_field", "scaled", "tuneable"]]
        #["out/sim/sequence_length{sequence_length}/beast-relaxed-fixed.pickle".format(sequence_length=sequence_length) for sequence_length in sequence_lengths],
        #["out/sim/sequence_length{sequence_length}/variational-relaxed-{approx}.pickle".format(sequence_length=sequence_length, approx=approx) for sequence_length in sequence_lengths for approx in ["mean_field", "scaled", "tuneable"]]

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
        beast_result = "out/{dataset}/beast-relaxed-fixed.pickle",
        beast_config = "config/beast-config.yaml",
        prior_params = "config/prior-params.yaml",
        variational_fit = "out/{dataset}/variational-relaxed-{approx}.pickle",
        notebook = "notebook/plot-posterior-relaxed.ipynb"
    output:
        notebook = "out/{dataset}/plot-posterior-relaxed-{approx}.ipynb",
        correlation_plot = "out/{dataset}/rate-correlations-{approx}.png",
        marginal_plot = "out/{dataset}/marginals-{approx}.png",
        rate_marginal_plot = "out/{dataset}/rate-marginals-{approx}.png"
    shell:
        """
        papermill {input.notebook} {output.notebook} \
            -p beast_result_file {input.beast_result} \
            -p topology_file {input.topology} \
            -p beast_config_file {input.beast_config} \
            -p prior_params_file {input.prior_params} \
            -p variational_fit_file {input.variational_fit} \
            -p clock_approx {wildcards.approx} \
            -p correlation_plot_out_file {output.correlation_plot} \
            -p marginal_plot_out_file {output.marginal_plot} \
            -p rate_marginal_plot_out_file {output.rate_marginal_plot}
        """

rule relaxed_report:
    input:
        "out/{dataset}/plot-posterior-relaxed-{approx}.ipynb"
    output:
        "out/{dataset}/plot-posterior-relaxed-{approx}.html"
    shell:
        "jupyter nbconvert --to html {input}"