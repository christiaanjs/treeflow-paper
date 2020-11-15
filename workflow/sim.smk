from treeflow_pipeline.util import yaml_input, yaml_output, text_input, text_output, beast_log_input, sequence_input, pickle_output
import treeflow_pipeline.simulation as sim
import treeflow_pipeline.templating as tem
import treeflow_pipeline.model as mod
import treeflow_pipeline.topology_inference as top
import treeflow_pipeline.results as res
import pathlib

configfile: "config/sim-config.yaml"

wd = pathlib.Path(config["working_directory"])
model = mod.Model(yaml_input(config["model_file"]))

SEQUENCE_LENGTHS = [1000]#[100, 1000, 10000]
APPROXES = ["scaled"]#["mean_field", "scaled"]

rule test_sim:
    input:
        expand("out/sim/10taxa/2seed/{sequence_length}sites/plot-posterior-relaxed-{clock_approx}.html", sequence_length=SEQUENCE_LENGTHS, clock_approx=APPROXES)

taxon_dir = "{taxon_count}taxa"

rule sampling_times:
    output:
        wd / taxon_dir / "sampling-times.yaml"
    run:
        yaml_output(sim.get_sampling_times(config, int(wildcards.taxon_count)), output[0])

seed_dir = "{seed}seed"

rule sample_prior:
    input:
        sampling_times = wd / taxon_dir / "sampling-times.yaml"
    output:
        wd / taxon_dir / seed_dir / "prior-sample.yaml"
    run:
        yaml_output(sim.sample_prior(yaml_input(input.sampling_times), model, int(wildcards.seed)), output[0])
    
rule tree_sim_xml:
    input:
        sampling_times = wd / taxon_dir / "sampling-times.yaml",
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml"
    output:
        wd / taxon_dir / seed_dir / "tree-sim.xml"
    run:
        text_output(
            tem.build_tree_sim(
                config,
                yaml_input(input.sampling_times),
                yaml_input(input.prior_sample),
                output[0]
            ),
            output[0]
        )

rule tree_sim:
    input:
        wd / taxon_dir / seed_dir / "tree-sim.xml"
    output:
        wd / taxon_dir / seed_dir / "tree-sim.trees"
    shell:
        "beast -seed {config[seed]} {input}"

rule tree_sim_newick:
    input:
        "{dir}/tree-sim.trees"
    output:
        "{dir}/tree-sim.newick"
    run:
        top.convert_tree(input[0], 'nexus', output[0], 'newick')

rule branch_rate_sim_xml:
    input:
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml",
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick"
    output:
        wd / taxon_dir / seed_dir / "branch-rate-sim.xml"
    run:
        text_output(
            tem.build_branch_rate_sim(
                text_input(input.tree),
                yaml_input(input.prior_sample),
                output[0]
            ),
            output[0]
        )

rule branch_rate_sim:
    input:
        wd / taxon_dir / seed_dir / "branch-rate-sim.xml"
    output:
        wd / taxon_dir / seed_dir / "branch-rate-sim.log",
        wd / taxon_dir / seed_dir / "branch-rate-sim.trees"
    shell:
        "beast -seed {config[seed]} {input}"

rule convert_branch_rates:
    input:
        wd / taxon_dir / seed_dir / "branch-rate-sim.log"
    output:
        wd / taxon_dir / seed_dir / "branch-rates.yaml"
    run:
        yaml_output(
            sim.parse_branch_rates(beast_log_input(input[0])),
            output[0]
        )

sequence_dir = "{sequence_length}sites"

rule sequence_sim_xml:
    input:
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick",
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml",
        branch_rates = wd / taxon_dir / seed_dir / "branch-rates.yaml",
        sampling_times = wd / taxon_dir / "sampling-times.yaml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "sim-seq.xml"
    run:
        text_output(
            tem.build_sequence_sim(
                config,
                text_input(input.tree),
                yaml_input(input.prior_sample),
                yaml_input(input.branch_rates),
                yaml_input(input.sampling_times),
                int(wildcards.sequence_length),
                output[0]
            ),
            output[0]
        )

rule sequence_sim:
    input:
        wd / taxon_dir / seed_dir / sequence_dir / "sim-seq.xml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "sequences.xml"
    shell:
        "beast -seed {config[seed]} {input}"

rule fasta_sim:
    input:
        "{wd}/sequences.xml"
    output:
        "{wd}/sequences.fasta"
    run:
        sim.convert_simulated_sequences(input[0], output[0], 'fasta')

rule sim_starting_values:
    input:
        prior_sample = wd / taxon_dir / seed_dir / "prior-sample.yaml"
    output:
        wd / taxon_dir / seed_dir / "starting-values.yaml"
    run:
        yaml_output(dict(
            frequencies=config["frequencies"],
            kappa=config["kappa"],
            **yaml_input(input.prior_sample)
        ), output[0])

rule beast_xml:
    input:
        fasta = wd / taxon_dir / seed_dir / sequence_dir / "sequences.fasta",
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick",
        starting_values = wd / taxon_dir / seed_dir / "starting-values.yaml",
        beast_config = "config/beast-config.yaml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "beast.xml"
    run:
        text_output(tem.build_beast_analysis(
            sequence_input(input.fasta),
            text_input(input.tree),
            yaml_input(input.starting_values),
            model,
            yaml_input(input.beast_config),
            output[0]
        ), output[0])

rule beast_run:
    input:
        wd / taxon_dir / seed_dir / sequence_dir / "beast.xml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "beast.log",
        wd / taxon_dir / seed_dir / sequence_dir / "beast.trees"
    shell:
        "beast -seed {config[seed]} {input}"

rule beast_results:
    input:
        topology = wd / taxon_dir / seed_dir / "tree-sim.newick",
        trees = wd / taxon_dir / seed_dir / sequence_dir / "beast.trees",
        trace = wd / taxon_dir / seed_dir / sequence_dir / "beast.log",
        beast_config = "config/beast-config.yaml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "beast.pickle"
    run:
        pickle_output(
            res.process_beast_results(
                input.trees,
                input.trace,
                input.topology,
                yaml_input(input.beast_config),
                model
            ),
            output[0]
        )

rule variational_fit:
    input:
        fasta = wd / taxon_dir / seed_dir / sequence_dir / "sequences.fasta",
        tree = wd / taxon_dir / seed_dir / "tree-sim.newick",
        starting_values = wd / taxon_dir / seed_dir / "starting-values.yaml"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "variational-fit-{clock_approx}.pickle"
    shell:
        """
        treeflow_pipeline -s {config[seed]} \
            {input.fasta} {config[model_file]} {output} \
            variational-fit \
            -t {input.tree} \
            -s {input.starting_values} \
            -c {wildcards.clock_approx} \
            --config {config[vi_config]}
        """

rule relaxed_plot:
    input:
        topology = wd / taxon_dir / seed_dir / "tree-sim.newick",
        beast_result = wd / taxon_dir / seed_dir / sequence_dir / "beast.pickle",
        variational_fit = wd / taxon_dir / seed_dir / sequence_dir / "variational-fit-{clock_approx}.pickle",
        notebook = "notebook/plot-posterior-relaxed.ipynb"
    output:
        notebook = wd / taxon_dir / seed_dir / sequence_dir / "plot-posterior-relaxed-{clock_approx}.ipynb",
        correlation_plot = wd / taxon_dir / seed_dir / sequence_dir / "rate-correlations-{clock_approx}.png",
        marginal_plot = wd / taxon_dir / seed_dir / sequence_dir / "marginals-{clock_approx}.png",
        rate_marginal_plot = wd / taxon_dir / seed_dir / sequence_dir / "rate-marginals-{clock_approx}.png"
    shell:
        """
        papermill {input.notebook} {output.notebook} \
            -p beast_result_file {input.beast_result} \
            -p topology_file {input.topology} \
            -p model_file {config[model_file]} \
            -p variational_fit_file {input.variational_fit} \
            -p clock_approx {wildcards.clock_approx} \
            -p correlation_plot_out_file {output.correlation_plot} \
            -p marginal_plot_out_file {output.marginal_plot} \
            -p rate_marginal_plot_out_file {output.rate_marginal_plot}
        """

rule relaxed_report:
    input:
        wd / taxon_dir / seed_dir / sequence_dir / "plot-posterior-relaxed-{clock_approx}.ipynb"
    output:
        wd / taxon_dir / seed_dir / sequence_dir / "plot-posterior-relaxed-{clock_approx}.html"
    shell:
        "jupyter nbconvert --to html {input}"